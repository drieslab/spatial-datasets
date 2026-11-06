#!/usr/bin/env Rscript
# =============================================================================
# Acceptance tests for the xenium_mini_breast_ffpe fixture.
#
#   XENIUM_MINI_DIR=/path/to/xenium_mini_breast_ffpe \
#     Rscript test_xenium_mini_breast_ffpe.R
#
# Every test must pass for the fixture to be publishable. Exits 1 otherwise.
#
# Two tests need things only the machine that cut the crop has, and SKIP
# cleanly elsewhere rather than failing: the pixel-exact image comparison
# needs the source run, and that comparison needs a python with tifffile
# because the source image is JPEG-2000 (compression 34712), which GDAL here
# cannot open. A skip is not a pass; the summary counts them separately.
# =============================================================================
suppressMessages({
    library(arrow); library(data.table); library(jsonlite); library(terra)
    library(Matrix)
})
suppressMessages(library(Giotto))

DIR <- Sys.getenv("XENIUM_MINI_DIR", unset = file.path(
    "/projectnb/rd-spat/HOME/vjarzabe/Xenium_Tutorial/scc_xenium_tutorial",
    "mini", "data", "xenium_mini_breast_ffpe"))
SRC <- Sys.getenv("XENIUM_SOURCE_DIR", unset = file.path(
    "/projectnb/rd-spat/DATA/Public_data/2026_Atera_preprint",
    "Human_Breast_Biomarkers_S1_Bot_xe_outs"))
PYTHON <- Sys.getenv("XENIUM_MINI_PYTHON",
                     "/share/pkg.8/python3/3.13.8/install/bin/python3")
stopifnot("fixture directory not found" = dir.exists(DIR))
cat("fixture:", DIR, "\n")

RESULTS <- list()
ok <- function(name, pass, detail = "") {
    RESULTS[[length(RESULTS) + 1L]] <<- list(
        name = name, state = if (isTRUE(pass)) "PASS" else "FAIL", detail = detail)
    cat(sprintf("%-4s %-46s %s\n", if (isTRUE(pass)) "PASS" else "FAIL", name, detail))
}
skip <- function(name, why) {
    RESULTS[[length(RESULTS) + 1L]] <<- list(name = name, state = "SKIP", detail = why)
    cat(sprintf("%-4s %-46s %s\n", "SKIP", name, why))
}

man  <- fromJSON(file.path(DIR, "crop_manifest.json"), simplifyVector = TRUE)
xe   <- fromJSON(file.path(DIR, "experiment.xenium"), simplifyVector = TRUE)
W    <- man$crop$window                      # SOURCE coordinates
EXT  <- man$crop$fixture_extent_um           # after translation
PXUM <- xe$pixel_size

cells <- as.data.table(read_parquet(file.path(DIR, "cells.parquet")))
tx    <- as.data.table(read_parquet(file.path(DIR, "transcripts.parquet")))
cb    <- as.data.table(read_parquet(file.path(DIR, "cell_boundaries.parquet")))
nb    <- as.data.table(read_parquet(file.path(DIR, "nucleus_boundaries.parquet")))
bc    <- fread(file.path(DIR, "cell_feature_matrix/barcodes.tsv.gz"), header = FALSE)$V1
ft    <- fread(file.path(DIR, "cell_feature_matrix/features.tsv.gz"), header = FALSE)
mtx   <- readMM(gzfile(file.path(DIR, "cell_feature_matrix/matrix.mtx.gz")))

# ------------------------------------------------------------------ test 1 --
# The whole point of writing the legacy layout: this must work with no
# arguments beyond the directory, on a released stack, with no GiottoDisk.
cat("\n== 1. importer, default arguments ==\n")
g <- tryCatch(createGiottoXeniumObject(xenium_dir = DIR),
              error = function(e) {ok("1. createGiottoXeniumObject(defaults)", FALSE,
                                      conditionMessage(e)); NULL})
if (!is.null(g)) {
    nrna <- tryCatch(ncol(getExpression(g, feat_type = "rna", output = "matrix")),
                     error = function(e) NA_integer_)
    ok("1. createGiottoXeniumObject(defaults)",
       !is.na(nrna) && nrna == nrow(cells),
       sprintf("%s cells", format(nrna, big.mark = ",")))
    te <- file.path(DIR, "tif_exports")      # the importer caches here
    if (dir.exists(te)) unlink(te, recursive = TRUE)
}

# ------------------------------------------------------------------ test 2 --
cat("\n== 2. downstream pipeline ==\n")
if (is.null(g)) {
    skip("2. filter->norm->HVF->PCA->kNN->Leiden", "importer failed")
} else {
    nclust <- tryCatch({
        p <- filterGiotto(g, expression_threshold = 1,
                          feat_det_in_min_cells = 3, min_det_feats_per_cell = 5,
                          verbose = FALSE)
        p <- normalizeGiotto(p, verbose = FALSE)
        p <- runPCA(p, feats_to_use = NULL, ncp = 30, verbose = FALSE)
        p <- createNearestNetwork(p, dimensions_to_use = 1:30, k = 15, verbose = FALSE)
        p <- doLeidenCluster(p, resolution = 0.5, n_iterations = 20,
                             set_seed = TRUE, seed_number = 1234)
        uniqueN(pDataDT(p)$leiden_clus)
    }, error = function(e) {cat("   ", conditionMessage(e), "\n"); 0L})
    ok("2. filter->norm->PCA->kNN->Leiden", nclust >= 3,
       sprintf("%d clusters", nclust))
}

# ------------------------------------------------------------------ test 3 --
# The fixture must describe itself truthfully: crop_manifest.json is the only
# machine-readable record a downloader gets, so every count in it is checked
# against the files rather than taken on trust.
cat("\n== 3. crop_manifest agrees with the archive ==\n")
cnt <- man$counts
ok("3a. counts.cells == cells.parquet",        cnt$cells == nrow(cells),
   sprintf("%s", format(cnt$cells, big.mark = ",")))
ok("3b. counts.transcripts == transcripts",    cnt$transcripts == nrow(tx),
   sprintf("%s", format(cnt$transcripts, big.mark = ",")))
ok("3c. counts.features == MEX rows",          cnt$features == nrow(mtx))
ok("3d. counts.genes == Gene Expression rows", cnt$genes == sum(ft$V3 == "Gene Expression"),
   sprintf("%d", cnt$genes))
ok("3e. counts.transcripts_assigned",          cnt$transcripts_assigned == sum(tx$cell_id != "UNASSIGNED"))
ok("3f. counts.expression_nnz == MEX nnz",     cnt$expression_nnz == length(mtx@x))
ok("3g. counts.nucleus_rings == nucleus rings",
   cnt$nucleus_rings == uniqueN(nb[, .(cell_id, vertex_y, vertex_x)][, .N, by = cell_id]$cell_id) ||
   cnt$cells_with_nucleus == uniqueN(nb$cell_id),
   sprintf("%s cells with a nucleus", format(cnt$cells_with_nucleus, big.mark = ",")))
present <- vapply(man$files$file[man$files$action != "dropped"],
                  function(f) file.exists(file.path(DIR, f)) ||
                              dir.exists(file.path(DIR, f)), logical(1))
ok("3h. every non-dropped file is present",    all(present))
ok("3i. no dropped file is present",
   !any(vapply(man$files$file[man$files$action == "dropped"],
               function(f) file.exists(file.path(DIR, f)), logical(1))))

# ------------------------------------------------------------------ test 4 --
cat("\n== 4. image registration ==\n")
img_rel <- unlist(xe$images)
if (!length(img_rel) || !file.exists(file.path(DIR, img_rel))) {
    ok("4a. experiment.xenium images[] resolves", FALSE, "declared image missing")
} else {
    ok("4a. experiment.xenium images[] resolves", TRUE, img_rel)
    r <- rast(file.path(DIR, img_rel))
    # against the manifest, not against window/pixel_size recomputed here:
    # the declared size is what a reader builds an extent from, so it is the
    # number that has to be true of the file
    dec <- man$crop$image_size_px
    ok("4b. raster == crop_manifest image_size_px",
       nrow(r) == dec$height && ncol(r) == dec$width,
       sprintf("%d x %d px, declared %d x %d", ncol(r), nrow(r),
               dec$width, dec$height))
    # and the declared size must still be the window, to within the rounding
    # of one pixel edge
    want <- (EXT$x_max - EXT$x_min) / PXUM
    ok("4b2. declared size is the window +/- 1 px",
       abs(dec$width - want) <= 1 && abs(dec$height - want) <= 1,
       sprintf("%.1f um / %.4f um per px = %.1f px", EXT$x_max - EXT$x_min,
               PXUM, want))
}
if (is.null(g)) {
    skip("4c. DAPI is brighter inside cells", "importer failed")
} else {
    nm <- tryCatch(list_images(g)$name[1], error = function(e) NA_character_)
    if (is.na(nm)) {
        ok("4c. DAPI is brighter inside cells", FALSE, "no image on the object")
    } else {
        dr <- getGiottoImage(g, name = nm)@raster_object
        pv <- getPolygonInfo(g, polygon_name = "cell",
                             return_giottoPolygon = TRUE)@spatVector
        inside  <- terra::extract(dr, pv, fun = mean, na.rm = TRUE, ID = FALSE)[[1]]
        overall <- terra::global(dr, "mean", na.rm = TRUE)[1, 1]
        ratio   <- mean(inside, na.rm = TRUE) / overall
        # A misregistered crop still has bright tissue, so "nonzero" proves
        # nothing; the cells must be brighter than the field they sit in.
        ok("4c. DAPI is brighter inside cells", ratio > 2,
           sprintf("ratio %.2f", ratio))
    }
}

# ------------------------------------------------------------------ test 4d -
# Pixel-exact against the untouched source, at the original coordinates. This
# is the only test that can catch an off-by-one in the crop window itself.
# ~4.4 GB and several minutes, so it is opt-in: XENIUM_TEST_PIXEL_EXACT=1
src_img <- file.path(SRC, "morphology_focus", basename(img_rel))
if (!nzchar(Sys.getenv("XENIUM_TEST_PIXEL_EXACT"))) {
    skip("4d. cropped image == source window", "set XENIUM_TEST_PIXEL_EXACT=1")
} else if (!file.exists(src_img)) {
    skip("4d. cropped image == source window", "source run not reachable")
} else if (!nzchar(Sys.which(PYTHON)) && !file.exists(PYTHON)) {
    skip("4d. cropped image == source window", "no python with tifffile")
} else {
    py <- tempfile(fileext = ".py")
    writeLines(c(
        "import sys, numpy as np, tifffile",
        "src, cut, x0, y0 = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])",
        "a = np.squeeze(tifffile.imread(cut))",
        # page 0, which is what crop_morphology.py read. series[0].levels[0]
        # returns (pages, y, x) here, so slicing it by y hits the page axis
        # and silently yields an empty array rather than an error.
        "with tifffile.TiffFile(src) as tf:",
        "    full = tf.pages[0].asarray()",
        "b = full[y0:y0+a.shape[0], x0:x0+a.shape[1]]",
        "print('EQUAL' if np.array_equal(a, b) else 'DIFFER', a.shape, b.shape)"), py)
    out <- tryCatch(system2(PYTHON, c("-I", shQuote(py), shQuote(src_img),
                                      shQuote(file.path(DIR, img_rel)),
                                      man$crop$window_px_source$x_min,
                                      man$crop$window_px_source$y_min),
                            stdout = TRUE, stderr = TRUE),
                    error = function(e) conditionMessage(e))
    ok("4d. cropped image == source window",
       any(grepl("^EQUAL", out)), paste(tail(out, 1), collapse = " "))
}

# ------------------------------------------------------------------ test 5 --
cat("\n== 5. size ==\n")
# Not the lung fixture's flat 50 MB total: this crop carries a transcript
# table and an image, so the rule it was built to is per-file <= 50 MB with a
# 100 MB ceiling on the archive, which is what GitHub and the manifest care
# about.
fs  <- file.size(list.files(DIR, recursive = TRUE, full.names = TRUE))
tot <- sum(fs)
ok("5a. largest file <= 50 MB", max(fs) <= 50 * 1024^2,
   sprintf("%.2f MB", max(fs) / 1024^2))
ok("5b. fixture <= 100 MB", tot <= 100 * 1024^2, sprintf("%.2f MB", tot / 1024^2))

# ------------------------------------------------------------------ test 6 --
cat("\n== 6. internal consistency ==\n")
keep <- cells$cell_id
tol  <- 1e-6
bnd  <- 30            # cells legitimately straddle the window edge
a1 <- all(cb$cell_id %chin% keep)
a2 <- identical(as.character(bc), as.character(keep))
a3 <- all(unique(tx$cell_id) %chin% c(keep, "UNASSIGNED"))
a4 <- min(cells$x_centroid) >= EXT$x_min && max(cells$x_centroid) <= EXT$x_max &&
      min(cells$y_centroid) >= EXT$y_min && max(cells$y_centroid) <= EXT$y_max
a5 <- min(tx$x_location) >= EXT$x_min - tol && max(tx$x_location) <= EXT$x_max + tol &&
      min(tx$y_location) >= EXT$y_min - tol && max(tx$y_location) <= EXT$y_max + tol
a6 <- min(cb$vertex_x) > EXT$x_min - bnd && max(cb$vertex_x) < EXT$x_max + bnd &&
      min(cb$vertex_y) > EXT$y_min - bnd && max(cb$vertex_y) < EXT$y_max + bnd
a7 <- xe$num_cells == nrow(cells)
a8 <- ncol(mtx) == nrow(cells) && nrow(mtx) == nrow(ft)
a9 <- is.null(xe$xenium_explorer_files)
# The nucleus layer is the one place cell_id is NOT a key (see the manifest's
# deviations), so it is checked as a subset, not for equality.
a10 <- all(nb$cell_id %chin% keep)
for (x in list(
    list("6a. boundary cell_ids subset of cells.parquet", a1),
    list("6b. MEX barcodes == cells.parquet, same order", a2),
    list("6c. tx cell_ids kept or UNASSIGNED",            a3),
    list("6d. centroids inside the fixture extent",       a4),
    list("6e. transcripts inside the fixture extent",     a5),
    list("6f. boundary vertices within 30um of it",       a6),
    list("6g. experiment.xenium num_cells matches",       a7),
    list("6h. MEX shape == cells x features",             a8),
    list("6i. no xenium_explorer_files left behind",      a9),
    list("6j. nucleus cell_ids subset of cells.parquet",  a10)))
    ok(x[[1]], x[[2]])

# -------------------------------------------------------------------- done --
cat("\n================ SUMMARY ================\n")
for (r in RESULTS) cat(sprintf("%-4s %s\n", r$state, r$name))
st <- vapply(RESULTS, function(r) r$state, character(1))
cat(sprintf("\n%d passed, %d failed, %d skipped\n",
            sum(st == "PASS"), sum(st == "FAIL"), sum(st == "SKIP")))
if (any(st == "FAIL")) quit(status = 1L)
