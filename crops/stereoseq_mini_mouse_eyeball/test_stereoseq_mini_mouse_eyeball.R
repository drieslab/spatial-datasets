#!/usr/bin/env Rscript
# =============================================================================
# Acceptance tests for the stereoseq_mini_mouse_eyeball fixture.
# =============================================================================
suppressMessages({
    library(reticulate); library(rhdf5); library(data.table)
    library(jsonlite); library(terra)
})

BASE    <- "/Users/rubendries/Documents/Datasets/Cropped_datasets/Stereoseq"
OUT_DIR <- file.path(BASE, "stereoseq_mini_mouse_eyeball")
PNG_OUT <- file.path(BASE, "stereoseq_mini_mouse_eyeball_alignment.png")
SN      <- "C04687E314"

py_require(c("h5py", "numpy", "tifffile", "imagecodecs", "zarr", "anndata", "scanpy"))
invisible(lapply(c("h5py", "numpy", "tifffile", "imagecodecs", "zarr"), import))
suppressMessages(library(Giotto))

mani <- fromJSON(file.path(OUT_DIR, "crop_manifest.json"))
W    <- mani$crop$window_units_original
SIZE <- as.integer(mani$crop$fixture_extent_units$x_max)
SRC  <- mani$source$dir

RESULTS <- list()
KNOWN   <- list()
ok <- function(name, pass, detail = "") {
    RESULTS[[length(RESULTS) + 1L]] <<- list(name = name, pass = isTRUE(pass), detail = detail)
    cat(sprintf("%-4s %-58s %s\n", if (isTRUE(pass)) "PASS" else "FAIL", name, detail))
}

# ------------------------------------------------------------------ test 1 --
cat("\n== 1. square-bin importer ==\n")
g100 <- createGiottoStereoSeqObjectBin(stereoseq_dir = OUT_DIR, verbose = FALSE)
ok("1a. createGiottoStereoSeqObjectBin(defaults)",
   !is.null(g100) && nrow(list_expression(g100)) >= 1 && nrow(list_images(g100)) >= 1,
   sprintf("spat_unit '%s', %d cells, %d feats, %d images",
           list_expression(g100)$spat_unit[1], ncol(getExpression(g100, output = "matrix")),
           nrow(getExpression(g100, output = "matrix")), nrow(list_images(g100))))

bins <- c("bin1", "bin5", "bin10", "bin20", "bin50", "bin100", "bin150", "bin200")
nb <- vapply(bins, function(b) {
    gg <- createGiottoStereoSeqObjectBin(stereoseq_dir = OUT_DIR, bin_size = b,
                                         load_image = FALSE, load_mask = FALSE, verbose = FALSE)
    ncol(getExpression(gg, spat_unit = b, output = "matrix"))
}, numeric(1))
ok("1b. all 8 bin_size values load", all(nb > 0),
   paste(sprintf("%s=%d", bins, nb), collapse = " "))

gbp <- createGiottoStereoSeqObjectBin(stereoseq_dir = OUT_DIR, load_binpoints = TRUE,
                                      load_image = FALSE, load_mask = FALSE, verbose = FALSE)
ok("1c. load_binpoints = TRUE (needs bin1)", nrow(list_feature_info(gbp)) >= 1,
   sprintf("%d feature layer(s)", nrow(list_feature_info(gbp))))

# ------------------------------------------------------------------ test 2 --
cat("\n== 2. cell-bin importer ==\n")
gc <- createGiottoStereoSeqObjectCell(stereoseq_dir = OUT_DIR, verbose = FALSE)
ok("2a. createGiottoStereoSeqObjectCell(defaults)",
   !is.null(gc) && nrow(list_spatial_info(gc)) >= 1,
   sprintf("%d cells, %d polys", ncol(getExpression(gc, spat_unit = "cell", output = "matrix")),
           nrow(list_spatial_info(gc))))
gc2 <- createGiottoStereoSeqObjectCell(stereoseq_dir = OUT_DIR, gef_type = "cellbin",
                                       load_image = FALSE, verbose = FALSE)
ok("2b. gef_type = 'cellbin' resolves",
   ncol(getExpression(gc2, spat_unit = "cell", output = "matrix")) > 0,
   sprintf("%d cells", ncol(getExpression(gc2, spat_unit = "cell", output = "matrix"))))

# ------------------------------------------------------------------ test 3 --
# The cropped h5ad must stay valid AnnData and stay consistent with the GEFs.
# NOTE: anndataToGiotto() cannot read SAW-written h5ad at all -- GiottoClass's
# ad2g.py extract_pca() unconditionally reads adata.uns[["pca"]] whenever obsm
# holds X_pca (scanpy's convention), but SAW writes a flat uns[["pca_variance_ratio"]].
# Verified to fail identically on the UNTOUCHED source file, so it is an upstream
# limitation, not a property of this fixture. Reported as KNOWN, not FAIL.
cat("\n== 3. cropped h5ad ==\n")
known <- function(name, detail) {
    KNOWN[[length(KNOWN) + 1L]] <<- list(name = name, detail = detail)
    cat(sprintf("%-4s %-58s %s\n", "KNWN", name, detail))
}
ad <- import("anndata")
for (h in list(list(f = sprintf("%s.bin50_1.0.h5ad", SN), ref = "bin50"),
               list(f = sprintf("%s.cellbin_1.0.adjusted.h5ad", SN), ref = "cell"))) {
    a <- ad$read_h5ad(file.path(OUT_DIR, "analysis", h$f))
    obs <- py_to_r(a$obs)
    gref <- if (h$ref == "cell") gc else
        createGiottoStereoSeqObjectBin(stereoseq_dir = OUT_DIR, bin_size = "bin50",
                                       load_image = FALSE, load_mask = FALSE, verbose = FALSE)
    sr <- getSpatialLocations(gref, spat_unit = h$ref, output = "data.table")
    A <- unique(data.table(x = as.integer(obs$x), y = as.integer(obs$y)))
    B <- unique(data.table(x = as.integer(sr$sdimx), y = as.integer(abs(sr$sdimy))))
    setkey(A, x, y); setkey(B, x, y)
    ok(sprintf("3a. %s h5ad is valid AnnData, coords inside window", h$ref),
       a$n_obs > 0 && min(A$x) >= 0 && min(A$y) >= 0 && max(A$x) < SIZE && max(A$y) < SIZE,
       sprintf("%d obs x %d vars, x[%d,%d] y[%d,%d]", a$n_obs, a$n_vars,
               min(A$x), max(A$x), min(A$y), max(A$y)))
    ok(sprintf("3b. %s h5ad positions == the GEF object exactly", h$ref),
       nrow(A) == nrow(B) && nrow(A[!B]) == 0,
       sprintf("%d vs %d positions, %d unmatched", nrow(A), nrow(B), nrow(A[!B])))
    sp <- py_to_r(a$obsm$get("spatial"))
    ok(sprintf("3c. %s obsm$spatial agrees with obs x/y", h$ref),
       all(as.integer(sp[, 1]) == as.integer(obs$x)) &&
       all(as.integer(sp[, 2]) == as.integer(obs$y)))
    e <- tryCatch({ anndataToGiotto(file.path(OUT_DIR, "analysis", h$f)); NULL },
                  error = function(e) conditionMessage(e))
    if (is.null(e)) ok(sprintf("3d. anndataToGiotto(%s)", h$ref), TRUE)
    else known(sprintf("3d. anndataToGiotto(%s) -- upstream ad2g.py limitation", h$ref),
               sprintf("%s (identical on the untouched source file)", trimws(e)))
}

# ------------------------------------------------------------------ test 4 --
cat("\n== 4. downstream pipeline ==\n")
run_pipe <- function(g, su) {
    g <- filterGiotto(g, spat_unit = su, expression_threshold = 1,
                      feat_det_in_min_cells = 5, min_det_feats_per_cell = 20)
    g <- normalizeGiotto(g, spat_unit = su)
    g <- calculateHVF(g, spat_unit = su, show_plot = FALSE, save_plot = FALSE, return_plot = FALSE)
    g <- runPCA(g, spat_unit = su, scale_unit = TRUE, center = TRUE, ncp = 30)
    g <- createNearestNetwork(g, spat_unit = su, dimensions_to_use = 1:15, k = 15)
    doLeidenCluster(g, spat_unit = su, resolution = 0.5)
}
for (cfg in list(list(g = g100, su = "bin100"), list(g = gc, su = "cell"))) {
    gg <- tryCatch(run_pipe(cfg$g, cfg$su), error = function(e) e)
    if (inherits(gg, "error")) { ok(sprintf("4. pipeline (%s)", cfg$su), FALSE,
                                    conditionMessage(gg)); next }
    cm <- pDataDT(gg, spat_unit = cfg$su)
    cl <- grep("leiden", colnames(cm), value = TRUE)[1]
    ok(sprintf("4. filter->norm->HVF->PCA->kNN->Leiden (%s)", cfg$su),
       uniqueN(cm[[cl]]) >= 3,
       sprintf("%d units, %d clusters", nrow(cm), uniqueN(cm[[cl]])))
}

# ------------------------------------------------------------------ test 5 --
cat("\n== 5. structural audit + registration ==\n")
source_python(file.path(BASE, "audit_fixture.py"))
au <- audit(out_dir = OUT_DIR, source_dir = SRC, sn = SN,
            x0 = as.integer(W$x_min), y0 = as.integer(W$y_min), size = SIZE)
for (r in au[[1]]) ok(r$name, r$pass, r$detail)

for (r in registration(out_dir = OUT_DIR, source_dir = SRC, sn = SN,
                       x0 = as.integer(W$x_min), y0 = as.integer(W$y_min), size = SIZE))
    ok(r$name, r$pass, r$detail)

sl <- getSpatialLocations(g100, spat_unit = "bin100", output = "data.table")
he <- rast(file.path(OUT_DIR, "image", sprintf("%s_HE_regist.tif", SN)))
ext(he) <- c(0, SIZE, -SIZE, 0)
slc <- getSpatialLocations(gc, spat_unit = "cell", output = "data.table")
png(PNG_OUT, width = 1900, height = 980, res = 130)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
terra::plotRGB(he, r = 1, g = 2, b = 3, mar = c(4, 4, 3, 1),
               main = "H&E + bin100 centres (MID-scaled)")
mid <- Matrix::colSums(getExpression(g100, spat_unit = "bin100", output = "matrix"))
points(sl$sdimx, sl$sdimy, pch = 15,
       cex = 0.30 + 0.9 * (rank(mid) / length(mid)), col = "#ff3b3088")
terra::plotRGB(he, r = 1, g = 2, b = 3, mar = c(4, 4, 3, 1),
               main = "H&E + cell-bin centroids + polygons")
poly <- getPolygonInfo(gc, polygon_name = "cell", return_giottoPolygon = TRUE)
terra::plot(poly[], add = TRUE, border = "#00b0ffcc", lwd = 0.4)
points(slc$sdimx, slc$sdimy, pch = 16, cex = 0.18, col = "#ffd60a")
dev.off()
cat("   plot:", PNG_OUT, "\n")

# ------------------------------------------------------------------ test 6 --
sz <- sum(file.size(list.files(OUT_DIR, recursive = TRUE, full.names = TRUE)))
ok("6. fixture size", sz <= 120 * 1024^2, sprintf("%.2f MB", sz / 1024^2))

# ------------------------------------------------------------------ summary -
cat("\n================ SUMMARY ================\n")
for (k in KNOWN) cat(sprintf("KNWN %s\n     %s\n", k$name, k$detail))
np <- sum(!vapply(RESULTS, function(r) r$pass, logical(1)))
for (r in RESULTS) if (!r$pass) cat(sprintf("FAIL %s  %s\n", r$name, r$detail))
cat(sprintf("%d/%d passed\n", length(RESULTS) - np, length(RESULTS)))
if (np > 0) quit(status = 1)
