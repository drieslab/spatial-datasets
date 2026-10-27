#!/usr/bin/env Rscript
# =============================================================================
# make_stereoseq_mini_mouse_eyeball.R
#
# Build a small, format-faithful Stereo-seq test fixture by SPATIALLY CROPPING
# the full SAW 8.2.0 run for chip C04687E314 (mouse eyeball).
#
#   SOURCE_DIR (a SAW `outs/`) -> OUT_DIR (an `outs/`-equivalent) in one run.
#
# This driver picks the window, provisions the python dependencies through
# reticulate/uv, calls stereoseq_crop.py to do the HDF5 / OME-TIFF / h5ad
# surgery, and writes crop_manifest.json.
#
# Requires: R  -- rhdf5, data.table, jsonlite, digest, reticulate
#           py -- h5py, numpy, tifffile, imagecodecs, zarr, anndata (via uv)
# =============================================================================

suppressMessages({
    library(rhdf5)
    library(data.table)
    library(jsonlite)
    library(reticulate)
})

# ---------------------------------------------------------------- config ----
SOURCE_DIR <- "/Users/rubendries/Documents/Datasets/stereoseq/C04687E314_backup/outs"
BASE_DIR   <- "/Users/rubendries/Documents/Datasets/Cropped_datasets/Stereoseq"
OUT_DIR    <- file.path(BASE_DIR, "stereoseq_mini_mouse_eyeball")
SN         <- "C04687E314"

BIN_GEF      <- sprintf("%s.tissue.gef", SN)
CELLBIN_GEFS <- sprintf(c("%s.cellbin.gef", "%s.adjusted.cellbin.gef"), SN)
H5ADS        <- sprintf(c("%s.bin50_1.0.h5ad", "%s.cellbin_1.0.adjusted.h5ad"), SN)

# lcm(1,5,10,20,50,100,150,200): the window origin and size must be multiples of
# this so that, after translation, every bin level's origins stay multiples of
# its own bin size.
LATTICE <- 600L
SIZE    <- 2400L                                  # 1.20 mm at 500 nm/unit
PIN     <- list(x0 = 4200L, y0 = 5400L)           # NULL to re-run the search

SCRIPT_PATH <- normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), mustWork = FALSE)

msg <- function(...) cat(sprintf(...), "\n", sep = "")

# ------------------------------------------------------- python bootstrap ----
py_require(c("h5py", "numpy", "tifffile", "imagecodecs", "zarr", "anndata"))
for (m in c("h5py", "numpy", "tifffile", "imagecodecs", "zarr", "anndata")) {
    invisible(import(m))
    stopifnot(py_module_available(m))
}
msg("[0] python ok: %s", py_config()$python)
source_python(file.path(BASE_DIR, "stereoseq_crop.py"))

# ------------------------------------------------------------ window ---------
h5closeAll()
gef <- file.path(SOURCE_DIR, "feature_expression", BIN_GEF)
b100 <- as.data.table(h5read(gef, "/geneExp/bin100/expression"))
h5closeAll()
occ <- b100[, .(mid = sum(count)), by = .(x, y)]

if (is.null(PIN)) {
    # densest SIZE x SIZE window on the LATTICE grid
    best <- NULL
    for (x0 in seq(2400L, max(occ$x) - SIZE, by = LATTICE)) {
        for (y0 in seq(4800L, max(occ$y) - SIZE, by = LATTICE)) {
            n <- occ[x >= x0 & x < x0 + SIZE & y >= y0 & y < y0 + SIZE, .N]
            if (is.null(best) || n > best$n) best <- list(x0 = x0, y0 = y0, n = n)
        }
    }
    W <- list(x0 = best$x0, y0 = best$y0)
} else {
    W <- PIN
}
stopifnot(W$x0 %% LATTICE == 0, W$y0 %% LATTICE == 0, SIZE %% LATTICE == 0)
n_occ <- occ[x >= W$x0 & x < W$x0 + SIZE & y >= W$y0 & y < W$y0 + SIZE, .N]
msg("[1] window chip units x:[%d, %d) y:[%d, %d)  (%d x %d units = %.2f x %.2f mm)",
    W$x0, W$x0 + SIZE, W$y0, W$y0 + SIZE, SIZE, SIZE, SIZE * 5e-4, SIZE * 5e-4)
msg("[1] occupied bin100 bins in window: %d of %d (%.0f%% fill)",
    n_occ, nrow(occ), 100 * n_occ / (SIZE / 100)^2)

# ------------------------------------------------------------- crop ----------
msg("[2] cropping ...")
t0 <- Sys.time()
rep <- main(source_dir = SOURCE_DIR, out_dir = OUT_DIR,
            x0 = W$x0, y0 = W$y0, size = SIZE,
            bin_gef = BIN_GEF, cellbin_gefs = as.list(CELLBIN_GEFS),
            h5ads = as.list(H5ADS))
msg("[2] done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))

bg <- rep$bin_gef[[BIN_GEF]]
for (b in names(bg$bins)) {
    x <- bg$bins[[b]]
    msg("    %-7s rows %8d -> %8d  occupied %6d  genes>0 %5d  MID %9d",
        b, x$rows_source, x$rows_kept, x$occupied_bins, x$genes_nonzero, x$mid_total)
}
for (f in names(rep$cellbin_gef)) {
    x <- rep$cellbin_gef[[f]]
    msg("    %-34s cells %5d -> %5d  nnz %8d -> %8d  blocks %dx%d",
        f, x$cells_source, x$cells_kept, x$nnz_source, x$nnz_kept,
        x$blocks[[1]], x$blocks[[2]])
}
for (f in names(rep$images)) {
    x <- rep$images[[f]]
    msg("    %-34s %s %s  %.2f MB", f, paste(unlist(x$shape), collapse = "x"),
        x$dtype, x$bytes / 1024^2)
}
for (f in names(rep$h5ad)) {
    x <- rep$h5ad[[f]]
    msg("    %-34s obs %6d -> %5d  vars %5d  %.2f MB (dropped uns: %s)",
        f, x$obs_source, x$obs_kept, x$vars, x$bytes / 1024^2,
        paste(unlist(x$uns_dropped), collapse = ", "))
}

# ---------------------------------------------------------- manifest ---------
fsize <- function(p) if (!file.exists(p)) 0 else if (dir.exists(p))
    sum(file.size(list.files(p, recursive = TRUE, full.names = TRUE))) else file.size(p)

act <- function(file, action, reason = NULL)
    c(list(file = file, action = action),
      if (!is.null(reason)) list(reason = reason) else NULL,
      list(bytes_source = unname(fsize(file.path(SOURCE_DIR, file))),
           bytes_fixture = unname(fsize(file.path(OUT_DIR, file)))))

actions <- list(
    act(file.path("feature_expression", BIN_GEF), "cropped"),
    act(file.path("feature_expression", CELLBIN_GEFS[1]), "cropped"),
    act(file.path("feature_expression", CELLBIN_GEFS[2]), "cropped"),
    act(file.path("image", sprintf("%s_HE_regist.tif", SN)), "cropped"),
    act(file.path("image", sprintf("%s_HE_mask.tif", SN)), "cropped"),
    act(file.path("analysis", H5ADS[1]), "cropped"),
    act(file.path("analysis", H5ADS[2]), "cropped"),
    act(file.path("feature_expression", sprintf("%s.gef", SN)), "dropped",
        "full-chip bin GEF; tissue.gef is the default gef_type and covers the bin workflow"),
    act(file.path("feature_expression", sprintf("%s.raw.gef", SN)), "dropped",
        "pre-tissue-mask bin GEF, not needed by the fixture"),
    act(file.path("feature_expression", sprintf("%s.merge.barcodeReadsCount.txt", SN)), "dropped",
        "read-level QC text, never read by Giotto"),
    act(file.path("feature_expression", sprintf("%s_raw_barcode_gene_exp.txt", SN)), "dropped",
        "raw barcode/gene text dump, duplicates the GEF"),
    act(file.path("analysis", sprintf("%s.bin20_1.0.h5ad", SN)), "dropped",
        "geometry already covered by bin50; largest of the four h5ad"),
    act(file.path("analysis", sprintf("%s.cellbin_1.0.h5ad", SN)), "dropped",
        "unadjusted cell-bin h5ad; the adjusted one is the importer default"),
    act(file.path("analysis", sprintf("%s.bin20_1.0.marker_features.csv", SN)), "dropped",
        "run-level differential expression that would be stale after cropping"),
    act(file.path("analysis", sprintf("%s.bin50_1.0.marker_features.csv", SN)), "dropped",
        "run-level differential expression that would be stale after cropping"),
    act(file.path("analysis", sprintf("%s.cellbin_1.0.marker_features.csv", SN)), "dropped",
        "run-level differential expression that would be stale after cropping"),
    act(file.path("analysis", sprintf("%s.cellbin_1.0.adjusted.marker_features.csv", SN)), "dropped",
        "run-level differential expression that would be stale after cropping"),
    act(file.path("image", sprintf("%s_HE_tissue_cut.tif", SN)), "dropped",
        "never matched by .stereoseq_image or .stereoseq_find_mask"),
    act(file.path("image", sprintf("%s_HE_mask_edm_dis_10.tif", SN)), "dropped",
        "distance map, never read by Giotto"),
    act("bam", "dropped", "not needed"),
    act(sprintf("%s.report.html", SN), "dropped", "run-level report that would be stale"),
    act("visualization.tar.gz", "dropped", "Stereo-seq Explorer bundle, never read by Giotto")
)

cfg <- tryCatch(readLines(file.path(dirname(SOURCE_DIR), "pipeline-logs", "config.yaml"),
                          warn = FALSE), error = function(e) character())
grab <- function(k) {
    h <- grep(sprintf("^\\s*%s:", k), cfg, value = TRUE)[1]
    if (is.na(h)) NA_character_ else trimws(sub("^[^:]*:", "", h))
}

manifest <- list(
    fixture_name = basename(OUT_DIR),
    created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    source = list(
        dir = SOURCE_DIR, sn = SN,
        saw_version = grab("SAW_VERSION"), kit_version = grab("kit_version"),
        organism = grab("organism"), tissue = grab("tissue"),
        resolution_nm = 500L, chip_size_units = c(11760L, 14700L)
    ),
    crop = list(
        method = paste("single rectangular window in chip units;",
                       "bins kept iff their ORIGIN is inside, cells kept iff their CENTROID is inside"),
        lattice = LATTICE,
        lattice_reason = paste("lcm(1,5,10,20,50,100,150,200) -- keeps every bin level's",
                               "origins on its own grid after translation"),
        window_units_original = list(x_min = W$x0, x_max = W$x0 + SIZE,
                                     y_min = W$y0, y_max = W$y0 + SIZE),
        translation_units = list(x = -W$x0, y = -W$y0),
        fixture_extent_units = list(x_min = 0L, x_max = SIZE, y_min = 0L, y_max = SIZE),
        fixture_extent_mm = SIZE * 5e-4,
        image_size_px = list(width = SIZE, height = SIZE)
    ),
    counts = rep[c("bin_gef", "cellbin_gef", "h5ad")],
    rebuilt = list(
        "geneExp/<bin>/gene offset+count rebuilt as an exclusive prefix sum; all genes retained, zero-count genes kept with count = 0",
        "wholeExp/<bin> dense grids recomputed from the cropped expression rows (validated to reproduce the source grids bit-for-bit before cropping)",
        "cellBin/cell offset rebuilt; cell ids renumbered to 0..n-1 because cellBin/geneExp$cellID indexes cells by row",
        "cellBin/geneExp rebuilt gene-major from the cell-major triplets (the two are the same triplet set, verified on the source)",
        "cellBin/geneExon recomputed -- note SAW's own value overflows a uint16 accumulator for 5 of 27820 genes",
        "cellBin/blockIndex + blockSize recomputed for the new 256-unit block grid, cells re-sorted into the new row-major block order",
        "every dataset attribute recording an extent, min/max/median/average or a count"
    ),
    deviations = list(
        "HDF5 datasets are written with gzip level 6; SAW writes them uncompressed. Transparent to rhdf5/h5py/gefpy.",
        "stat/ and contour/ are dropped from the bin GEF: stat/gene$E10 is a SAW-internal statistic whose formula is not recoverable, and tissueContour describes the whole section. Giotto reads neither.",
        "codedCellBlock/ is dropped from the cell-bin GEFs: an opaque serialised spatial index that cannot be regenerated and that Giotto never reads.",
        "root gef_area is scaled by the fraction of occupied bin1 positions retained -- it is an area measure and areas are additive, but SAW's exact definition is not recoverable from the file.",
        "cellBorder vertices are cell-relative with a 32767 sentinel, so they are NOT translated; only the dataset's absolute bbox attributes are recomputed.",
        "h5ad uns/rank_genes_groups, uns/gene_exp_leiden and uns/pca_variance_ratio are dropped as run-level results; var/n_cells and var/n_counts are recomputed exactly from the cropped counts."
    ),
    files = actions,
    script = list(
        driver = SCRIPT_PATH,
        module = file.path(BASE_DIR, "stereoseq_crop.py"),
        sha256 = list(
            driver = if (file.exists(SCRIPT_PATH)) digest::digest(file = SCRIPT_PATH, algo = "sha256") else NA_character_,
            module = digest::digest(file = file.path(BASE_DIR, "stereoseq_crop.py"), algo = "sha256")
        )
    )
)
write(toJSON(manifest, auto_unbox = TRUE, pretty = 4, digits = NA),
      file.path(OUT_DIR, "crop_manifest.json"))

total <- sum(file.size(list.files(OUT_DIR, recursive = TRUE, full.names = TRUE)))
msg("[3] crop_manifest.json written")
msg("[done] %s  total %.2f MB", OUT_DIR, total / 1024^2)
