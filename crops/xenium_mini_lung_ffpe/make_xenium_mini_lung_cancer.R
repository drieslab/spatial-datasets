#!/usr/bin/env Rscript
# =============================================================================
# make_xenium_mini_lung_cancer.R
#
# Build a small, format-faithful Xenium test fixture by SPATIALLY CROPPING the
# full 10x Xenium human lung cancer FFPE run.
#
#   SOURCE_DIR -> OUT_DIR in one run.  No manual steps.
#
# The fixture keeps exactly one file format per information type (parquet for
# tables, h5 for expression), retains original filenames / schemas / dtypes /
# compression, translates every micron coordinate so the window origin is
# (0, 0), crops the morphology_focus OME-TIFFs to the matching pixel window,
# and rewrites experiment.xenium so no field or manifest entry is stale.
#
# Requires: R packages arrow, data.table, rhdf5, jsonlite, digest, reticulate
#           python (via reticulate/uv): tifffile, imagecodecs, zarr, h5py, numpy
# =============================================================================

suppressMessages({
    library(arrow)
    library(data.table)
    library(jsonlite)
    library(reticulate)
})

# ---------------------------------------------------------------- config ----
SOURCE_DIR <- "/Users/rubendries/Documents/Datasets/Xenium/LungCancer_10X/Xenium_V1_humanLung_Cancer_FFPE_outs"
OUT_DIR    <- "/Users/rubendries/Documents/Datasets/Cropped_datasets/xenium_mini_lung_cancer"
SCRIPT_PATH <- normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), mustWork = FALSE)

TARGET_CELLS <- 1500
SIDE_UM      <- 425      # exactly 2000 px at 0.2125 um/px, and a multiple of STEP_UM
STEP_UM      <- 25       # sliding-window stride for the auto search

# Window found by the Phase-2 auto search (deterministic).  Set to NULL to
# re-run the search from scratch; it will reproduce these same pixel bounds.
PIN_WINDOW_PX <- list(x_min = 11882L, x_max = 13882L, y_min = 6588L, y_max = 8588L)

msg <- function(...) cat(sprintf(...), "\n", sep = "")

# ------------------------------------------------------- python bootstrap ----
py_require(c("tifffile", "imagecodecs", "zarr", "h5py", "numpy"))
invisible(import("tifffile")); invisible(import("imagecodecs")); invisible(import("h5py"))
stopifnot(py_module_available("tifffile"), py_module_available("imagecodecs"),
          py_module_available("zarr"), py_module_available("h5py"))
PY <- import_main()   # assign R values into python __main__
msg("[0] python ok: %s", py_config()$python)

# ------------------------------------------------ experiment.xenium truth ----
exp_src <- fromJSON(file.path(SOURCE_DIR, "experiment.xenium"), simplifyVector = TRUE)
PX <- exp_src$pixel_size
stopifnot(is.numeric(PX), PX > 0)
msg("[0] XOA %s  pixel_size %s um/px  num_cells %d",
    exp_src$analysis_sw_version, format(PX), exp_src$num_cells)

if (dir.exists(OUT_DIR)) unlink(OUT_DIR, recursive = TRUE)
dir.create(file.path(OUT_DIR, "morphology_focus"), recursive = TRUE)

# =========================================================== PHASE 2 ========
# Choose the crop window.
# ============================================================================
cells_all <- as.data.table(read_parquet(
    file.path(SOURCE_DIR, "cells.parquet"),
    col_select = c("cell_id", "x_centroid", "y_centroid")))

auto_window <- function() {
    # graphclust labels are used only to score heterogeneity; they are not shipped
    td <- file.path(tempdir(), "xen_analysis"); dir.create(td, showWarnings = FALSE)
    untar(file.path(SOURCE_DIR, "analysis.tar.gz"), exdir = td,
          files = "analysis/clustering/gene_expression_graphclust/clusters.csv")
    cl <- fread(file.path(td, "analysis/clustering/gene_expression_graphclust/clusters.csv"))
    setnames(cl, c("cell_id", "cluster"))
    d <- merge(cells_all, cl, by = "cell_id", all.x = TRUE)

    B <- STEP_UM
    d[, bx := as.integer(floor(x_centroid / B))][, by_ := as.integer(floor(y_centroid / B))]
    H <- matrix(0L, max(d$bx) + 1L, max(d$by_) + 1L)
    tb <- d[, .N, by = .(bx, by_)]
    H[cbind(tb$bx + 1L, tb$by_ + 1L)] <- tb$N
    CS <- rbind(0L, cbind(0L, t(apply(apply(H, 2, cumsum), 1, cumsum))))
    boxcount <- function(x0, y0) {
        i0 <- x0 / B; i1 <- (x0 + SIDE_UM) / B; j0 <- y0 / B; j1 <- (y0 + SIDE_UM) / B
        CS[cbind(i1 + 1, j1 + 1)] - CS[cbind(i0 + 1, j1 + 1)] -
            CS[cbind(i1 + 1, j0 + 1)] + CS[cbind(i0 + 1, j0 + 1)]
    }
    # stay >=100 um inside the tissue bounding box (no empty margin)
    xs <- seq(ceiling((min(d$x_centroid) + 100) / B) * B,
              floor((max(d$x_centroid) - 100 - SIDE_UM) / B) * B, by = B)
    ys <- seq(ceiling((min(d$y_centroid) + 100) / B) * B,
              floor((max(d$y_centroid) - 100 - SIDE_UM) / B) * B, by = B)
    grid <- CJ(x0 = xs, y0 = ys)
    grid[, n := boxcount(x0, y0)]
    cand <- grid[abs(n - TARGET_CELLS) <= 250][order(abs(n - TARGET_CELLS))]
    msg("[2] %d candidate windows, %d within +/-250 of target", nrow(grid), nrow(cand))

    best <- NULL
    for (i in seq_len(min(nrow(cand), 4000L))) {
        x0 <- cand$x0[i]; y0 <- cand$y0[i]
        w <- d[x_centroid >= x0 & x_centroid < x0 + SIDE_UM &
               y_centroid >= y0 & y_centroid < y0 + SIDE_UM]
        if (!nrow(w)) next
        ct <- sort(table(w$cluster), decreasing = TRUE) / nrow(w)
        # (c) no empty margin: every 100 um sub-tile must hold >=1 cell
        ntiles <- uniqueN(paste(floor((w$x_centroid - x0) / 100),
                                floor((w$y_centroid - y0) / 100)))
        if (ntiles < ceiling(SIDE_UM / 100)^2) next
        # (b) >=3 clusters at reasonable proportions, none dominant
        if (sum(ct >= 0.05) < 3 || ct[1] > 0.70) next
        # (a) closest to TARGET_CELLS, mildly rewarding cluster diversity
        sc <- abs(nrow(w) - TARGET_CELLS) / TARGET_CELLS - 0.05 * sum(ct >= 0.05)
        if (is.null(best) || sc < best$sc) best <- list(x0 = x0, y0 = y0, sc = sc)
    }
    stopifnot(!is.null(best))
    list(x_min = floor(best$x0 / PX), y_min = floor(best$y0 / PX))
}

if (is.null(PIN_WINDOW_PX)) {
    a <- auto_window()
    WPX <- list(x_min = a$x_min, x_max = a$x_min + round(SIDE_UM / PX),
                y_min = a$y_min, y_max = a$y_min + round(SIDE_UM / PX))
} else {
    WPX <- PIN_WINDOW_PX
}
W <- list(x_min = WPX$x_min * PX, x_max = WPX$x_max * PX,
          y_min = WPX$y_min * PX, y_max = WPX$y_max * PX)
W$width_um  <- W$x_max - W$x_min
W$height_um <- W$y_max - W$y_min
IMG_W <- as.integer(WPX$x_max - WPX$x_min)
IMG_H <- as.integer(WPX$y_max - WPX$y_min)
msg("[2] window um  x:[%.4f, %.4f]  y:[%.4f, %.4f]  (%.3f x %.3f um)",
    W$x_min, W$x_max, W$y_min, W$y_max, W$width_um, W$height_um)
msg("[2] window px  x:[%d, %d]  y:[%d, %d]  (%d x %d px)",
    WPX$x_min, WPX$x_max, WPX$y_min, WPX$y_max, IMG_W, IMG_H)

# =========================================================== PHASE 3 ========
# 1. the one cell set that drives everything cell-indexed
# ============================================================================
keep <- cells_all[x_centroid >= W$x_min & x_centroid < W$x_max &
                  y_centroid >= W$y_min & y_centroid < W$y_max]
setorder(keep, cell_id)
keep_cells <- keep$cell_id
N_KEEP <- length(keep_cells)
msg("[3.1] kept cells: %d", N_KEEP)
stopifnot(N_KEEP > 0)

# --------------------------------------------------- parquet write helper ----
# Preserve the source schema (names, dtypes) and its embedded pandas metadata,
# updating only the row count that metadata records.
retag_pandas <- function(meta, nrow_new) {
    if (is.null(meta$pandas)) return(meta)
    j <- fromJSON(meta$pandas, simplifyVector = FALSE)
    if (length(j$index_columns) && is.list(j$index_columns[[1]]) &&
        identical(j$index_columns[[1]]$kind, "range")) {
        j$index_columns[[1]]$stop <- nrow_new
    }
    meta$pandas <- as.character(toJSON(j, auto_unbox = TRUE, null = "null"))
    meta
}
# arrow's dataset filter can emit row groups out of source order, which would
# split a cell's boundary vertices into two non-adjacent runs.  Read and filter
# one row group at a time instead, so source row order is preserved exactly.
filter_rowgroups <- function(path, mask_fn) {
    pr <- ParquetFileReader$create(path)
    out <- list()
    for (i in seq_len(pr$num_row_groups)) {
        rg <- pr$ReadRowGroup(i - 1L)
        m <- mask_fn(rg)
        if (any(m)) out[[length(out) + 1L]] <- rg$Filter(Array$create(m))
    }
    stopifnot(length(out) > 0L)
    if (length(out) == 1L) out[[1L]] else do.call(concat_tables, out)
}

set_col <- function(tb, name, values, type) {
    i <- match(name, tb$ColumnNames()) - 1L
    tb$SetColumn(i, Field$create(name, type),
                 ChunkedArray$create(Array$create(values)$cast(type)))
}
write_faithful <- function(tb, src_meta, path) {
    tb$metadata <- retag_pandas(src_meta, tb$num_rows)
    write_parquet(tb, path, compression = "zstd")
    invisible(tb)
}

# --------------------------------------------------- 2. cells.parquet --------
cells_path <- file.path(SOURCE_DIR, "cells.parquet")
meta_cells <- open_dataset(cells_path)$schema$metadata
tb <- filter_rowgroups(cells_path, function(rg) {
    x <- as.numeric(rg$GetColumnByName("x_centroid"))
    y <- as.numeric(rg$GetColumnByName("y_centroid"))
    x >= W$x_min & x < W$x_max & y >= W$y_min & y < W$y_max
})
# the source is sorted by cell_id and row order is preserved, so this holds
stopifnot(identical(as.vector(tb$GetColumnByName("cell_id")), keep_cells))
tb <- set_col(tb, "x_centroid", as.numeric(tb$GetColumnByName("x_centroid")) - W$x_min, float64())
tb <- set_col(tb, "y_centroid", as.numeric(tb$GetColumnByName("y_centroid")) - W$y_min, float64())
write_faithful(tb, meta_cells, file.path(OUT_DIR, "cells.parquet"))
KEEP_TX_COUNTS <- as.numeric(tb$GetColumnByName("transcript_counts"))
KEEP_CELL_AREA <- sum(as.numeric(tb$GetColumnByName("cell_area")))
msg("[3.2] cells.parquet: %d rows", tb$num_rows)

# ------------------------------------------- 3. boundary parquets ------------
keep_set <- keep_cells
for (bf in c("cell_boundaries.parquet", "nucleus_boundaries.parquet")) {
    bpath <- file.path(SOURCE_DIR, bf)
    m <- open_dataset(bpath)$schema$metadata
    b <- filter_rowgroups(bpath, function(rg)
        as.vector(rg$GetColumnByName("cell_id")) %chin% keep_set)
    # every cell's vertices must stay in one contiguous, correctly ordered run
    bid <- as.vector(b$GetColumnByName("cell_id"))
    stopifnot(length(rle(bid)$lengths) == uniqueN(bid), !is.unsorted(bid))
    b <- set_col(b, "vertex_x", as.numeric(b$GetColumnByName("vertex_x")) - W$x_min, float32())
    b <- set_col(b, "vertex_y", as.numeric(b$GetColumnByName("vertex_y")) - W$y_min, float32())
    write_faithful(b, m, file.path(OUT_DIR, bf))
    msg("[3.3] %s: %d vertices, %d cells", bf, b$num_rows,
        uniqueN(as.vector(b$GetColumnByName("cell_id"))))
}

# ------------------------------------------- 4. transcripts.parquet ----------
# Filter by TRANSCRIPT COORDINATE, not by cell_id, so in-window UNASSIGNED
# transcripts are retained the way a real run has them.  Then drop the few
# transcripts assigned to a cell whose centroid fell outside the window, so
# every referenced cell_id is covered by cells.parquet and the h5 matrix.
tx_path <- file.path(SOURCE_DIR, "transcripts.parquet")
meta_tx <- open_dataset(tx_path)$schema$metadata
tx <- filter_rowgroups(tx_path, function(rg) {
    x <- as.numeric(rg$GetColumnByName("x_location"))
    y <- as.numeric(rg$GetColumnByName("y_location"))
    x >= W$x_min & x < W$x_max & y >= W$y_min & y < W$y_max
})
tx_cell <- as.vector(tx$GetColumnByName("cell_id"))
n_pre <- tx$num_rows
mask <- tx_cell %in% keep_cells | tx_cell == "UNASSIGNED"
msg("[3.4] transcripts in window: %d; dropping %d assigned to out-of-window cells",
    n_pre, sum(!mask))
tx <- tx$Filter(Array$create(mask))
tx <- set_col(tx, "x_location", as.numeric(tx$GetColumnByName("x_location")) - W$x_min, float32())
tx <- set_col(tx, "y_location", as.numeric(tx$GetColumnByName("y_location")) - W$y_min, float32())
write_faithful(tx, meta_tx, file.path(OUT_DIR, "transcripts.parquet"))
N_TX <- tx$num_rows
N_TX_UNASSIGNED <- sum(as.vector(tx$GetColumnByName("cell_id")) == "UNASSIGNED")
N_TX_Q20 <- sum(as.numeric(tx$GetColumnByName("qv")) >= 20)
msg("[3.4] transcripts.parquet: %d rows (%d UNASSIGNED)", N_TX, N_TX_UNASSIGNED)
rm(tx, tx_cell, mask); invisible(gc())

# ------------------------------------------- 5. cell_feature_matrix.h5 -------
# Column-subset the 10x CSC matrix, preserving group layout, feature order,
# dataset dtypes (data int32 / indices int64 / indptr int64 / shape int32) and
# every root attribute.  h5py is used for exact int64 control.
PY$h5_src   <- file.path(SOURCE_DIR, "cell_feature_matrix.h5")
PY$h5_dst   <- file.path(OUT_DIR, "cell_feature_matrix.h5")
PY$keep_bc  <- keep_cells
py_run_string('
import h5py, numpy as np

with h5py.File(h5_src, "r") as f:
    m = f["matrix"]
    barcodes = m["barcodes"][:]
    data     = m["data"][:]
    indices  = m["indices"][:]
    indptr   = m["indptr"][:]
    shape    = m["shape"][:]
    feats    = {k: m["features"][k][:] for k in m["features"].keys()}
    root_attrs = {k: v for k, v in f.attrs.items()}
    dt_bc    = m["barcodes"].dtype
    dt_feat  = {k: m["features"][k].dtype for k in m["features"].keys()}
    dt_data, dt_idx, dt_ptr, dt_shape = (m["data"].dtype, m["indices"].dtype,
                                         m["indptr"].dtype, m["shape"].dtype)

bc_str = np.array([b.decode() if isinstance(b, bytes) else b for b in barcodes])
pos = {b: i for i, b in enumerate(bc_str)}
sel = np.array([pos[b] for b in keep_bc], dtype=np.int64)   # keeps keep_bc order

starts, ends = indptr[sel], indptr[sel + 1]
counts = (ends - starts).astype(np.int64)
new_indptr = np.concatenate(([0], np.cumsum(counts))).astype(dt_ptr)
if counts.sum():
    take = np.concatenate([np.arange(s, e) for s, e in zip(starts, ends)]).astype(np.int64)
else:
    take = np.zeros(0, dtype=np.int64)
new_data    = data[take].astype(dt_data)
new_indices = indices[take].astype(dt_idx)
new_shape   = np.array([shape[0], len(sel)], dtype=dt_shape)
new_bc      = np.array([b.encode() for b in keep_bc], dtype=dt_bc)

with h5py.File(h5_dst, "w") as g:
    for k, v in root_attrs.items():
        g.attrs[k] = v
    m = g.create_group("matrix")
    m.create_dataset("barcodes", data=new_bc,      dtype=dt_bc)
    m.create_dataset("data",     data=new_data,    dtype=dt_data,  maxshape=(None,))
    m.create_dataset("indices",  data=new_indices, dtype=dt_idx,   maxshape=(None,))
    m.create_dataset("indptr",   data=new_indptr,  dtype=dt_ptr,   maxshape=(None,))
    m.create_dataset("shape",    data=new_shape,   dtype=dt_shape, maxshape=(None,))
    fg = m.create_group("features")
    for k, v in feats.items():
        fg.create_dataset(k, data=v, dtype=dt_feat[k])

h5_nnz   = int(new_data.shape[0])
h5_ncell = int(new_shape[1])
h5_nfeat = int(new_shape[0])
')
N_NNZ <- PY$h5_nnz
msg("[3.5] cell_feature_matrix.h5: %d features x %d cells, nnz %d",
    PY$h5_nfeat, PY$h5_ncell, N_NNZ)
stopifnot(PY$h5_ncell == N_KEEP)

# ------------------------------------------- 6. morphology_focus crop --------
# Multi-file OME series: 4 files, one channel each (FirstC=N -> _000N), each
# carrying the shared 4-channel OME-XML.  Crop page 0 of each file with
# tile-wise zarr reads and re-emit the original XML with updated SizeX/SizeY.
PY$src_focus <- file.path(SOURCE_DIR, "morphology_focus")
PY$dst_focus <- file.path(OUT_DIR, "morphology_focus")
PY$x0 <- WPX$x_min; PY$x1 <- WPX$x_max; PY$y0 <- WPX$y_min; PY$y1 <- WPX$y_max
py_run_string('
import os, re, glob, tifffile, zarr, numpy as np

files = sorted(glob.glob(os.path.join(src_focus, "morphology_focus_*.ome.tif")))
focus_report = []
for p in files:
    with tifffile.TiffFile(p) as f:
        xml = f.ome_metadata
        page = f.pages[0]
        store = page.aszarr()
        z = zarr.open(store, mode="r")
        sub = np.asarray(z[int(y0):int(y1), int(x0):int(x1)])
        store.close()
    new_xml = re.sub(r\'SizeX="\\d+"\', \'SizeX="%d"\' % sub.shape[1], xml)
    new_xml = re.sub(r\'SizeY="\\d+"\', \'SizeY="%d"\' % sub.shape[0], new_xml)
    out = os.path.join(dst_focus, os.path.basename(p))
    # the OME-XML contains "\u00b5m"; write it as UTF-8 bytes exactly as the source does
    tifffile.imwrite(out, sub, description=new_xml.encode("utf-8"), ome=False, bigtiff=True,
                     photometric="minisblack", tile=(512, 512),
                     compression="zlib", metadata=None)
    focus_report.append((os.path.basename(p), list(sub.shape), str(sub.dtype),
                         int(sub.min()), int(sub.max()), float(sub.mean()),
                         os.path.getsize(out)))
')
for (r in PY$focus_report)
    msg("[3.6] %s -> %dx%d %s  range [%d, %d] mean %.1f  %.2f MB",
        r[[1]], r[[2]][[1]], r[[2]][[2]], r[[3]], r[[4]], r[[5]], r[[6]], r[[7]] / 1024^2)

# ------------------------------------------- 7. gene_panel.json --------------
stopifnot(file.copy(file.path(SOURCE_DIR, "gene_panel.json"), file.path(OUT_DIR, "gene_panel.json")))
msg("[3.7] gene_panel.json copied verbatim")

# ------------------------------------------- 8. experiment.xenium rewrite ----
exp_new <- fromJSON(file.path(SOURCE_DIR, "experiment.xenium"), simplifyVector = FALSE)
exp_new$num_cells <- N_KEEP
# same definitions 10x uses: median transcripts per cell, and high-quality
# (qv >= 20) transcripts per 100 um^2 of segmented cell area
exp_new$transcripts_per_cell <- as.integer(round(median(KEEP_TX_COUNTS)))
exp_new$transcripts_per_100um <- N_TX_Q20 / KEEP_CELL_AREA * 100
# manifests: drop every entry that points at a file we did not keep
exp_new$images <- list(morphology_focus_filepath = "morphology_focus/morphology_focus_0000.ome.tif")
exp_new$xenium_explorer_files <- NULL
exp_txt <- toJSON(exp_new, auto_unbox = TRUE, pretty = 4, digits = NA)
# jsonlite renders whole-number doubles without a decimal point; restore the
# source rendering for fields we did not change
exp_txt <- sub('"z_step_size": 3,', '"z_step_size": 3.0,', exp_txt, fixed = TRUE)
write(exp_txt, file.path(OUT_DIR, "experiment.xenium"))
msg("[3.8] experiment.xenium: num_cells %d, transcripts_per_cell %d, transcripts_per_100um %.4f",
    exp_new$num_cells, exp_new$transcripts_per_cell, exp_new$transcripts_per_100um)

# ------------------------------------------- 9. crop_manifest.json -----------
fsize <- function(p) if (file.exists(p)) file.size(p) else
    sum(file.size(list.files(p, recursive = TRUE, full.names = TRUE)))
actions <- list(
    list(file = "transcripts.parquet",       action = "kept-filtered"),
    list(file = "cells.parquet",             action = "kept-filtered"),
    list(file = "cell_boundaries.parquet",   action = "kept-filtered"),
    list(file = "nucleus_boundaries.parquet",action = "kept-filtered"),
    list(file = "cell_feature_matrix.h5",    action = "kept-subset"),
    list(file = "morphology_focus/",         action = "cropped"),
    list(file = "gene_panel.json",           action = "kept-verbatim"),
    list(file = "experiment.xenium",         action = "rewritten"),
    list(file = "cell_feature_matrix/",      action = "dropped", reason = "MEX duplicates the h5"),
    list(file = "cell_feature_matrix.tar.gz",action = "dropped", reason = "duplicates the h5"),
    list(file = "cells.csv.gz",              action = "dropped", reason = "duplicates cells.parquet"),
    list(file = "transcripts.csv.gz",        action = "dropped", reason = "duplicates transcripts.parquet"),
    list(file = "cell_boundaries.csv.gz",    action = "dropped", reason = "duplicates cell_boundaries.parquet"),
    list(file = "nucleus_boundaries.csv.gz", action = "dropped", reason = "duplicates nucleus_boundaries.parquet"),
    list(file = "cells.zarr.zip",            action = "dropped", reason = "zarr not supported by Giotto"),
    list(file = "transcripts.zarr.zip",      action = "dropped", reason = "zarr not supported by Giotto"),
    list(file = "cell_feature_matrix.zarr.zip", action = "dropped", reason = "zarr not supported by Giotto"),
    list(file = "analysis.zarr.zip",         action = "dropped", reason = "zarr not supported by Giotto"),
    list(file = "analysis.tar.gz",           action = "dropped", reason = "run-level analysis; read only to score the window"),
    list(file = "analysis_summary.html",     action = "dropped", reason = "run-level report, would be stale"),
    list(file = "aux_outputs.tar.gz",        action = "dropped", reason = "per-cycle raw images, not needed"),
    list(file = "morphology.ome.tif",        action = "dropped", reason = "full 3D z-stack, never read by the importer"),
    list(file = "metrics_summary.csv",       action = "dropped", reason = "run-level QC that would be stale; dropped rather than falsified")
)
for (i in seq_along(actions)) {
    actions[[i]]$bytes_source <- unname(fsize(file.path(SOURCE_DIR, actions[[i]]$file)))
    actions[[i]]$bytes_fixture <- if (grepl("dropped", actions[[i]]$action)) 0L else
        unname(fsize(file.path(OUT_DIR, actions[[i]]$file)))
}
manifest <- list(
    fixture_name = basename(OUT_DIR),
    created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    source = list(
        dir = SOURCE_DIR,
        run_name = exp_src$run_name, region_name = exp_src$region_name,
        slide_id = exp_src$slide_id, panel_name = exp_src$panel_name,
        analysis_sw_version = exp_src$analysis_sw_version,
        pixel_size_um = PX, num_cells = exp_src$num_cells
    ),
    crop = list(
        method = "spatial rectangular window; cells kept iff centroid strictly inside",
        selection = sprintf("auto, TARGET_CELLS=%d, %g um square, %g um stride, scored on graphclust diversity + no empty 100um sub-tile",
                            TARGET_CELLS, SIDE_UM, STEP_UM),
        window_um_original = list(x_min = W$x_min, x_max = W$x_max, y_min = W$y_min, y_max = W$y_max),
        window_px_original = list(x_min = WPX$x_min, x_max = WPX$x_max, y_min = WPX$y_min, y_max = WPX$y_max),
        translation_um = list(x = -W$x_min, y = -W$y_min),
        fixture_extent_um = list(x_min = 0, x_max = W$width_um, y_min = 0, y_max = W$height_um),
        image_size_px = list(width = IMG_W, height = IMG_H)
    ),
    counts = list(cells = N_KEEP, transcripts = N_TX, transcripts_unassigned = N_TX_UNASSIGNED,
                  transcripts_qv20 = N_TX_Q20,
                  features = PY$h5_nfeat, expression_nnz = N_NNZ),
    deviations = list(
        "morphology_focus OME-TIFFs are re-encoded from JPEG-2000 (compression 34712) to zlib so they are readable by GDAL/terra as well as tifffile",
        "pyramid sub-resolutions are not regenerated; only the full-resolution level is written, which is the only level the Giotto importer reads",
        "the multi-file OME-XML is preserved verbatim apart from SizeX/SizeY"
    ),
    files = actions,
    script = list(path = SCRIPT_PATH,
                  sha256 = if (nzchar(SCRIPT_PATH) && file.exists(SCRIPT_PATH))
                      digest::digest(file = SCRIPT_PATH, algo = "sha256") else NA_character_)
)
write(toJSON(manifest, auto_unbox = TRUE, pretty = 4, digits = NA),
      file.path(OUT_DIR, "crop_manifest.json"))

total <- sum(file.size(list.files(OUT_DIR, recursive = TRUE, full.names = TRUE)))
msg("[3.9] crop_manifest.json written")
msg("[done] %s  total %.2f MB", OUT_DIR, total / 1024^2)
