#!/usr/bin/env Rscript
# =============================================================================
# Acceptance tests for the visiumhd_mini_crc fixture.
#
# 1. internal consistency: h5 == mtx, positions == barcodes, 8 / 16 um and the
#    cell matrix are exact aggregates of the 2 um bins, geojson == cells
# 2. fidelity: every kept column equals the same column in the source
# 3. readers: Giotto in memory and GiottoDisk on disk, bins and cells
# 4. alignment: spatlocs and polygons sit on the cropped image
#
# Env: VHD_OUT_DIR, the fixture to check (default: the local copy it was made in),
# and VHD_SOURCE_DIR, the extracted 10x `outs/` it was cut from.
# =============================================================================
suppressMessages({
    library(data.table); library(jsonlite); library(hdf5r); library(Matrix)
})
options(giotto.has_conda = FALSE, giotto.no_python_warn = TRUE,
        giottodisk.dgc_max_rows = Inf, giottodisk.dgc_max_cols = Inf)

BASE    <- "/Users/rubendries/Documents/Datasets/Cropped_datasets/VisiumHD"
OUT_DIR <- Sys.getenv("VHD_OUT_DIR", file.path(BASE, "visiumhd_mini_crc"))
PNG_OUT <- sub("/?$", "_alignment.png", OUT_DIR)
mani <- fromJSON(file.path(OUT_DIR, "crop_manifest.json"), simplifyVector = FALSE)
SRC  <- Sys.getenv("VHD_SOURCE_DIR",
    "/Users/rubendries/Documents/Datasets/visiumHD/Human_Colon_Cancer_v4/outs")

RESULTS <- list()
ok <- function(name, pass, detail = "") {
    RESULTS[[length(RESULTS) + 1L]] <<- list(name = name, pass = isTRUE(pass), detail = detail)
    cat(sprintf("%-4s %-60s %s\n", if (isTRUE(pass)) "PASS" else "FAIL", name, detail))
}
sq <- function(um) sprintf("square_%03dum", um)
bdir <- function(root, um) file.path(root, "binned_outputs", sq(um))

read_h5 <- function(f, cols = NULL) {
    h <- H5File$new(f, "r"); on.exit(h$close_all())
    m <- sparseMatrix(i = h[["matrix/indices"]][] + 1, p = as.numeric(h[["matrix/indptr"]][]),
                      x = as.numeric(h[["matrix/data"]][]), dims = h[["matrix/shape"]][],
                      dimnames = list(h[["matrix/features/id"]][], h[["matrix/barcodes"]][]))
    if (!is.null(cols)) m <- m[, cols, drop = FALSE]
    m
}
read_mex <- function(d) {
    m <- as(readMM(gzfile(file.path(d, "matrix.mtx.gz"))), "CsparseMatrix")
    dimnames(m) <- list(fread(file.path(d, "features.tsv.gz"), header = FALSE)$V1,
                        readLines(gzfile(file.path(d, "barcodes.tsv.gz"))))
    m
}
# group columns of m by `groups` (one per column; NA dropped)
agg <- function(m, groups) {
    keep <- !is.na(groups)
    lv <- unique(groups[keep])
    g <- sparseMatrix(i = which(keep), j = match(groups[keep], lv), x = 1,
                      dims = c(ncol(m), length(lv)), dimnames = list(NULL, lv))
    m %*% g
}
same <- function(a, b) identical(dim(a), dim(b)) && identical(dimnames(a), dimnames(b)) &&
    nnzero(drop0(a - b)) == 0L

# ------------------------------------------------------------------ test 1 --
cat("\n== 1. internal consistency ==\n")
m <- list(); tp <- list()
for (um in c(2L, 8L, 16L)) {
    d <- bdir(OUT_DIR, um)
    m[[sq(um)]] <- read_h5(file.path(d, "raw_feature_bc_matrix.h5"))
    tp[[sq(um)]] <- as.data.table(arrow::read_parquet(file.path(d, "spatial", "tissue_positions.parquet")))
    ok(sprintf("1a. %s h5 == mtx (raw)", sq(um)),
       same(m[[sq(um)]], read_mex(file.path(d, "raw_feature_bc_matrix"))),
       sprintf("%d x %d, nnz %d", nrow(m[[sq(um)]]), ncol(m[[sq(um)]]), length(m[[sq(um)]]@x)))
    fh <- read_h5(file.path(d, "filtered_feature_bc_matrix.h5"))
    # an FFPE run's filtered matrix also keeps only the filtered probe set
    ok(sprintf("1b. %s filtered == raw[features, in_tissue]", sq(um)),
       setequal(colnames(fh), tp[[sq(um)]][in_tissue == 1, barcode]) &&
           same(fh, m[[sq(um)]][rownames(fh), colnames(fh), drop = FALSE]) &&
           same(fh, read_mex(file.path(d, "filtered_feature_bc_matrix"))),
       sprintf("%d in tissue, %d features", ncol(fh), nrow(fh)))
    ok(sprintf("1c. %s positions == barcodes", sq(um)),
       setequal(tp[[sq(um)]]$barcode, colnames(m[[sq(um)]])))
}
bm <- as.data.table(arrow::read_parquet(file.path(OUT_DIR, "barcode_mappings.parquet")))
bm <- bm[match(colnames(m$square_002um), square_002um)]
for (um in c(8L, 16L)) {
    a <- agg(m$square_002um, bm[[sq(um)]])
    ok(sprintf("1d. %s == aggregate of 2 um bins", sq(um)),
       same(a[, colnames(m[[sq(um)]])], m[[sq(um)]]))
}
seg <- file.path(OUT_DIR, "segmented_outputs")
mc <- read_h5(file.path(seg, "raw_feature_cell_matrix.h5"))
ok("1e. cell h5 == mtx", same(mc, read_mex(file.path(seg, "raw_feature_cell_matrix"))),
   sprintf("%d cells", ncol(mc)))
ac <- agg(m$square_002um, ifelse(bm$in_cell, bm$cell_id, NA_character_))
ok("1f. cell matrix == aggregate of in-cell 2 um bins",
   setequal(colnames(ac), colnames(mc)) && same(ac[, colnames(mc)], mc))
for (gj in c("cell_segmentations", "nucleus_segmentations",
             "graphclust_annotated_cell_segmentations",
             "graphclust_annotated_nucleus_segmentations")) {
    g <- fromJSON(file.path(seg, paste0(gj, ".geojson")), simplifyVector = FALSE)
    ids <- sprintf("cellid_%09d-1", vapply(g$features, function(f) f$properties$cell_id, 1))
    ok(sprintf("1g. %s ids == cells", gj), setequal(ids, colnames(mc)))
}

# ------------------------------------------------------------------ test 2 --
cat("\n== 2. fidelity to the source ==\n")
for (um in c(2L, 8L, 16L)) {
    s <- read_h5(file.path(bdir(SRC, um), "raw_feature_bc_matrix.h5"), colnames(m[[sq(um)]]))
    ok(sprintf("2a. %s columns == source columns", sq(um)), same(s, m[[sq(um)]]))
}
s <- read_h5(file.path(SRC, "segmented_outputs", "raw_feature_cell_matrix.h5"), colnames(mc))
ok("2b. cell columns == source columns", same(s, mc))
stp <- as.data.table(arrow::read_parquet(file.path(bdir(SRC, 8L), "spatial", "tissue_positions.parquet")))
stp <- stp[match(tp$square_008um$barcode, barcode)]
tr <- mani$crop$translation
ok("2c. positions == source positions translated",
   isTRUE(all.equal(tp$square_008um$pxl_col_in_fullres, stp$pxl_col_in_fullres + tr$x)) &&
       isTRUE(all.equal(tp$square_008um$pxl_row_in_fullres, stp$pxl_row_in_fullres + tr$y)) &&
       identical(tp$square_008um$array_row, stp$array_row))

# ------------------------------------------------------------------ test 3 --
cat("\n== 3. readers ==\n")
suppressMessages({ library(Giotto); library(GiottoDisk) })
q <- function(x) suppressWarnings(suppressMessages(x))
bo <- file.path(OUT_DIR, "binned_outputs")
for (um in c(2L, 8L, 16L)) {
    g <- q(createGiottoVisiumHDObjectBin(bo, bin = um, verbose = FALSE))
    su <- sprintf("bin%03d", um)
    ok(sprintf("3a. in-memory bin %d um (with hires image)", um),
       ncol(getExpression(g, spat_unit = su, output = "matrix")) > 0 &&
           nrow(list_images(g)) == 1L,
       sprintf("%d bins", ncol(getExpression(g, spat_unit = su, output = "matrix"))))
}
n_tis <- sum(tp$square_008um$in_tissue)
gt <- q(createGiottoVisiumHDObjectBin(bo, bin = 8, tissue_only = TRUE,
                                      load_image = FALSE, verbose = FALSE))
gdt <- q(createGiottoVisiumHDObjectBin(bo, bin = 8, tissue_only = TRUE,
    load_image = FALSE, verbose = FALSE, backend = tempfile("vhd_tis_")))
ok("3a. tissue_only keeps exactly the in-tissue 8 um bins (memory, disk)",
   n_tis < nrow(tp$square_008um) &&
       setequal(spatIDs(gt, "bin008"), tp$square_008um[in_tissue == 1, barcode]) &&
       setequal(spatIDs(gdt, "bin008"), spatIDs(gt, "bin008")),
   sprintf("%d of %d", n_tis, nrow(tp$square_008um)))
gc <- q(createGiottoVisiumHDObjectCell(seg, load_transcripts = TRUE, verbose = FALSE))
ok("3b. in-memory cells + 2 um points",
   identical(sort(spatIDs(gc, "cell")), sort(colnames(mc))) &&
       inherits(getFeatureInfo(gc, return_giottoPoints = TRUE), "giottoBinPoints"))

proj <- tempfile("vhd_mini_")
gd <- q(importVisiumHDDisk(bo, backend = proj, bin = 8)$create_gobject(
    bin = c(2, 8, 16), load_bin_mapping = TRUE, load_transcripts = TRUE, verbose = FALSE))
gd <- q(importVisiumHDDisk(seg, backend = proj)$create_gobject(gobject = gd, verbose = FALSE))
ok("3c. disk object: bins + cells + parents",
   setequal(names(gd@expression), c("bin002", "bin008", "bin016", "cell")) &&
       all(c("bin008", "bin016", "cell") %in% names(pDataDT(gd, spat_unit = "bin002"))),
   toString(names(pDataDT(gd, spat_unit = "bin002"))))
md <- storeRead(getExpression(gd, spat_unit = "cell", output = "exprObj")[], output = "dgcmatrix")
mm <- getExpression(gc, spat_unit = "cell", output = "matrix")
ok("3d. disk cells == in-memory cells",
   identical(dim(md), dim(mm)) && max(abs(md[rownames(mm), colnames(mm)] - mm)) == 0)

# ------------------------------------------------------------------ test 4 --
cat("\n== 4. alignment ==\n")
g8 <- q(createGiottoVisiumHDObjectBin(bo, bin = 8, verbose = FALSE))
e <- terra::ext(getGiottoImage(g8, name = list_images(g8)$name[1]))
sl <- getSpatialLocations(g8, output = "data.table")
inside <- mean(sl$sdimx >= e$xmin & sl$sdimx <= e$xmax & sl$sdimy >= e$ymin & sl$sdimy <= e$ymax)
ok("4a. 8 um spatlocs inside the hires image extent", inside == 1, sprintf("%.3f", inside))
pe <- terra::ext(getPolygonInfo(gc, "cell", return_giottoPolygon = TRUE)[])
ok("4b. cell polygons inside the image extent",
   pe$xmin >= e$xmin && pe$xmax <= e$xmax && pe$ymin >= e$ymin && pe$ymax <= e$ymax)
g16 <- q(createGiottoVisiumHDObjectBin(bo, bin = 16, verbose = FALSE))
p1 <- spatPlot2D(g16, show_image = TRUE, point_size = 2.5, point_alpha = 0.5,
                 cell_color = "yellow", point_shape = "no_border", return_plot = TRUE)
p2 <- spatInSituPlotPoints(gc, show_image = TRUE, feats = NULL, show_polygon = TRUE,
                           polygon_feat_type = "cell", polygon_color = "yellow",
                           polygon_alpha = 0, polygon_line_size = 0.4, return_plot = TRUE)
grDevices::png(PNG_OUT, width = 1600, height = 800)
print(cowplot::plot_grid(p1, p2, ncol = 2))
grDevices::dev.off()
ok("4c. alignment plot written", file.exists(PNG_OUT), PNG_OUT)

# ------------------------------------------------------------------ summary -
total <- sum(file.size(list.files(OUT_DIR, recursive = TRUE, full.names = TRUE)))
ok("5. fixture size < 50 MB", total < 50 * 1024^2, sprintf("%.1f MB", total / 1024^2))
n_fail <- sum(!vapply(RESULTS, `[[`, TRUE, "pass"))
cat(sprintf("\n%d checks, %d failed\n", length(RESULTS), n_fail))
quit(status = as.integer(n_fail > 0))
