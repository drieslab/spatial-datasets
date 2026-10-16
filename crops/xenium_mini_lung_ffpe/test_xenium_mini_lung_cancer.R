#!/usr/bin/env Rscript
# =============================================================================
# Acceptance tests for the xenium_mini_lung_cancer fixture.
# All six must pass for the fixture to be usable by GiottoAgent / Giotto Suite.
# =============================================================================
suppressMessages({
    library(reticulate); library(arrow); library(data.table)
    library(jsonlite); library(terra)
})

OUT_DIR <- "/Users/rubendries/Documents/Datasets/Cropped_datasets/xenium_mini_lung_cancer"
PNG_OUT <- "/Users/rubendries/Documents/Datasets/Cropped_datasets/xenium_mini_lung_cancer_alignment.png"

# Giotto's OME-TIFF reader needs these importable from the reticulate python
py_require(c("tifffile", "imagecodecs", "leidenalg", "python-igraph", "numpy"))
invisible(import("tifffile")); invisible(import("imagecodecs"))

suppressMessages(library(Giotto))

RESULTS <- list()
ok <- function(name, pass, detail = "") {
    RESULTS[[length(RESULTS) + 1L]] <<- list(name = name, pass = isTRUE(pass), detail = detail)
    cat(sprintf("%-4s %-42s %s\n", if (isTRUE(pass)) "PASS" else "FAIL", name, detail))
}

# ---------------------------------------------------------------- test 1 ----
cat("\n== 1. importer, default arguments ==\n")
g <- createGiottoXeniumObject(xenium_dir = OUT_DIR)
print(g)
sl <- list_spatial_locations(g); pl <- list_spatial_info(g)
fi <- list_feature_info(g);      im <- list_images(g)
ex <- list_expression(g)
ok("1. createGiottoXeniumObject(defaults)",
   !is.null(g) && !is.null(ex) && !is.null(fi) && !is.null(im) &&
   nrow(pl) >= 2 && nrow(im) >= 4,
   sprintf("%d poly layers, %d feat layers, %d images, %d expr", nrow(pl), nrow(fi), nrow(im), nrow(ex)))

g0 <- g   # the un-filtered object, used for the registration checks
# the importer writes a tif_exports/ cache into the bundle; keep the fixture clean
te <- file.path(OUT_DIR, "morphology_focus", "tif_exports")

# ---------------------------------------------------------------- test 4a ---
# extent agreement, done on the gobject (Giotto flips y, so compare |y|)
cells_pq <- as.data.table(read_parquet(file.path(OUT_DIR, "cells.parquet")))
tif <- file.path(OUT_DIR, "morphology_focus", "morphology_focus_0000.ome.tif")
tx <- as.data.table(read_parquet(file.path(OUT_DIR, "transcripts.parquet"),
                                 col_select = c("x_location", "y_location", "cell_id")))
img <- getGiottoImage(g, name = list_images(g)$name[1])
ie <- ext(img)
sp <- getSpatialLocations(g, output = "data.table")
exp_meta <- fromJSON(file.path(OUT_DIR, "experiment.xenium"))
mani <- fromJSON(file.path(OUT_DIR, "crop_manifest.json"))
W <- mani$crop$fixture_extent_um
tol <- exp_meta$pixel_size
img_ok <- abs(ie[1] - 0) < tol && abs(ie[2] - W$x_max) < tol &&
          abs(ie[3] + W$y_max) < tol && abs(ie[4] - 0) < tol
ok("4a. image extent == window (<1 px)", img_ok,
   sprintf("image ext [%.4f,%.4f]x[%.4f,%.4f] vs window %.3f x %.3f um",
           ie[1], ie[2], ie[3], ie[4], W$x_max, W$y_max))

# ---------------------------------------------------------------- test 2 ----
cat("\n== 2. downstream pipeline ==\n")
g <- filterGiotto(g, expression_threshold = 1,
                  feat_det_in_min_cells = 3, min_det_feats_per_cell = 5)
g <- normalizeGiotto(g)
g <- addStatistics(g)
g <- calculateHVF(g, show_plot = FALSE, save_plot = FALSE, return_plot = FALSE)
g <- runPCA(g, feats_to_use = NULL, scale_unit = TRUE, center = TRUE, ncp = 30)
g <- createNearestNetwork(g, dimensions_to_use = 1:15, k = 15)
g <- tryCatch(doLeidenCluster(g, resolution = 0.5, n_iterations = 200),
              error = function(e) { message("  leidenalg unavailable (", conditionMessage(e),
                                            "); using doLeidenClusterIgraph()")
                                    doLeidenClusterIgraph(g, resolution = 0.5) })
cm <- pDataDT(g)
clcol <- grep("leiden", colnames(cm), value = TRUE)[1]
nclust <- uniqueN(cm[[clcol]])
ok("2. filter->norm->HVF->PCA->kNN->Leiden", nclust >= 3,
   sprintf("%d cells after filtering, %d Leiden clusters (%s)", nrow(cm), nclust, clcol))

# ---------------------------------------------------------------- test 3 ----
hits <- grep("xenium", c(basename(OUT_DIR),
                         list.files(OUT_DIR, recursive = TRUE, include.dirs = TRUE)),
             ignore.case = TRUE, value = TRUE)
ok("3. detect_dataset_type: >=2 'xenium' hits", length(hits) >= 2,
   paste(hits, collapse = ", "))

# ---------------------------------------------------------------- test 4b ---
# Registration, in the coordinate convention a Giotto user actually sees:
# tables are flipped (sdimy = -y_um) and the image carries negative_y = TRUE.
cat("\n== 4. registration + alignment plot ==\n")
sl <- getSpatialLocations(g0, spat_unit = "cell", output = "data.table")
rimg <- getGiottoImage(g0, name = list_images(g0)$name[1])@raster_object
set.seed(1)
v_cell <- terra::extract(rimg, cbind(sl$sdimx, sl$sdimy))[, 1]
v_rand <- terra::extract(rimg, cbind(runif(nrow(sl), 0, W$x_max),
                                     runif(nrow(sl), -W$y_max, 0)))[, 1]
v_mis  <- terra::extract(rimg, cbind(sl$sdimx, -W$y_max - sl$sdimy))[, 1]
ratio  <- mean(v_cell, na.rm = TRUE) / mean(v_rand, na.rm = TRUE)
ok("4b. DAPI at centroids >> background", ratio > 2 && mean(v_cell) > mean(v_mis),
   sprintf("mean DAPI: centroids %.0f, random %.0f, mis-flipped %.0f (ratio %.2fx)",
           mean(v_cell, na.rm = TRUE), mean(v_rand, na.rm = TRUE),
           mean(v_mis, na.rm = TRUE), ratio))

# ---------------------------------------------------------------- test 4c ---
# Pixel-exact check against the untouched source image at the original coords.
PYM <- import_main()
PYM$cx <- cells_pq$x_centroid; PYM$cy <- cells_pq$y_centroid
PYM$src_tif <- file.path(mani$source$dir, "morphology_focus", "morphology_focus_0000.ome.tif")
PYM$dst_tif <- tif
PYM$ox <- mani$crop$window_px_original$x_min; PYM$oy <- mani$crop$window_px_original$y_min
PYM$pxs <- exp_meta$pixel_size
py_require("zarr")
py_run_string('
import numpy as np, tifffile, zarr
cx = np.asarray(cx); cy = np.asarray(cy)
crop = tifffile.imread(dst_tif)
crop = crop[0] if crop.ndim == 3 else crop
rr = np.round(cy / pxs).astype(int); cc = np.round(cx / pxs).astype(int)
v_crop = crop[rr, cc]
with tifffile.TiffFile(src_tif) as f:
    st = f.pages[0].aszarr(); z = zarr.open(st, mode="r")
    tile = np.asarray(z[int(oy):int(oy) + crop.shape[0], int(ox):int(ox) + crop.shape[1]])
    st.close()
v_src = tile[rr, cc]
identical_px = bool(np.array_equal(crop, tile))
identical_at_cells = bool(np.array_equal(v_crop, v_src))
mean_at_cells = float(v_crop.mean())
')
ok("4c. cropped image == source window (pixel-exact)",
   isTRUE(PYM$identical_px) && isTRUE(PYM$identical_at_cells),
   sprintf("all %d x %d pixels identical to source[%d:%d, %d:%d]; mean DAPI at centroids %.1f",
           2000L, 2000L, mani$crop$window_px_original$y_min, mani$crop$window_px_original$y_max,
           mani$crop$window_px_original$x_min, mani$crop$window_px_original$x_max,
           PYM$mean_at_cells))

# ---------------------------------------------------------------- test 4d ---
png(PNG_OUT, width = 1900, height = 980, res = 130)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
stretch <- function(x) {   # quantile stretch so faint tissue stays visible
    q <- quantile(values(x), c(0.01, 0.995), na.rm = TRUE)
    clamp(x, q[1], q[2], values = TRUE)
}
rs <- stretch(rimg)
plot(rs, col = grey.colors(256, start = 0, end = 1), legend = FALSE,
     main = "DAPI + cell centroids", mar = c(4, 4, 3, 1))
points(sl$sdimx, sl$sdimy, pch = 16, cex = .20, col = "#ff3b30cc")
plot(rs, col = grey.colors(256, start = 0, end = 1), legend = FALSE,
     main = "DAPI + transcripts (2% sample); blue = UNASSIGNED", mar = c(4, 4, 3, 1))
sm <- tx[sample(.N, round(.N * 0.02))]
points(sm$x_location, -sm$y_location, pch = 16, cex = .14,
       col = ifelse(sm$cell_id == "UNASSIGNED", "#00b0ffaa", "#ffd60aaa"))
dev.off()
r <- rast(tif)
rext <- c(0, ncol(r) * exp_meta$pixel_size, 0, nrow(r) * exp_meta$pixel_size)
tab_ok <- abs(rext[2] - W$x_max) < tol && abs(rext[4] - W$y_max) < tol &&
          max(cells_pq$x_centroid) < W$x_max && max(cells_pq$y_centroid) < W$y_max
ok("4d. raster size == table extent", tab_ok,
   sprintf("raster %.4f x %.4f um; centroid max %.3f, %.3f", rext[2], rext[4],
           max(cells_pq$x_centroid), max(cells_pq$y_centroid)))
cat("   plot:", PNG_OUT, "\n")

# ---------------------------------------------------------------- test 5 ----
if (dir.exists(te)) unlink(te, recursive = TRUE)   # drop the importer's cache
sz <- sum(file.size(list.files(OUT_DIR, recursive = TRUE, full.names = TRUE)))
ok("5. fixture size <= 50 MB", sz <= 50 * 1024^2, sprintf("%.2f MB", sz / 1024^2))

# ---------------------------------------------------------------- test 6 ----
cat("\n== 6. consistency audit ==\n")
cb <- as.datatable <- as.data.table(read_parquet(file.path(OUT_DIR, "cell_boundaries.parquet")))
nb <- as.data.table(read_parquet(file.path(OUT_DIR, "nucleus_boundaries.parquet")))
h5bc <- rhdf5::h5read(file.path(OUT_DIR, "cell_feature_matrix.h5"), "/matrix/barcodes")
rhdf5::h5closeAll()
keep <- cells_pq$cell_id
a1 <- all(cb$cell_id %in% keep) && all(nb$cell_id %in% keep)
a2 <- identical(as.character(h5bc), as.character(keep))
txc <- unique(tx$cell_id)
a3 <- all(txc %in% c(keep, "UNASSIGNED"))
a4 <- min(cells_pq$x_centroid) >= 0 && min(cells_pq$y_centroid) >= 0 &&
      max(cells_pq$x_centroid) < W$x_max && max(cells_pq$y_centroid) < W$y_max &&
      min(tx$x_location) >= -tol && min(tx$y_location) >= -tol &&
      max(tx$x_location) <= W$x_max + tol && max(tx$y_location) <= W$y_max + tol
bnd <- 20   # cells legitimately straddle the window edge
a5 <- min(cb$vertex_x) > -bnd && max(cb$vertex_x) < W$x_max + bnd &&
      min(cb$vertex_y) > -bnd && max(cb$vertex_y) < W$y_max + bnd &&
      min(nb$vertex_x) > -bnd && max(nb$vertex_x) < W$x_max + bnd &&
      min(nb$vertex_y) > -bnd && max(nb$vertex_y) < W$y_max + bnd
a6 <- exp_meta$num_cells == nrow(cells_pq)
declared <- unlist(exp_meta$images)
a7 <- all(file.exists(file.path(OUT_DIR, declared))) && is.null(exp_meta$xenium_explorer_files)
h5shape <- rhdf5::h5read(file.path(OUT_DIR, "cell_feature_matrix.h5"), "/matrix/shape")
rhdf5::h5closeAll()
a8 <- h5shape[2] == nrow(cells_pq)
for (x in list(
    c("6a. boundary cell_ids subset of cells.parquet", a1),
    c("6b. h5 barcodes == cells.parquet (same order)", a2),
    c("6c. tx cell_ids in kept cells or UNASSIGNED",   a3),
    c("6d. centroids/transcripts inside window",       a4),
    c("6e. boundary vertices within 20um of window",   a5),
    c("6f. experiment.xenium num_cells matches",       a6),
    c("6g. no manifest entry points at a missing file", a7),
    c("6h. h5 shape[2] == n cells",                    a8)))
    ok(x[[1]], as.logical(x[[2]]))

# ------------------------------------------------------------------ summary -
cat("\n================ SUMMARY ================\n")
for (r in RESULTS) cat(sprintf("%-4s %s\n", if (r$pass) "PASS" else "FAIL", r$name))
np <- sum(!vapply(RESULTS, function(r) r$pass, logical(1)))
cat(sprintf("\n%d/%d passed\n", length(RESULTS) - np, length(RESULTS)))
if (np > 0) quit(status = 1)
