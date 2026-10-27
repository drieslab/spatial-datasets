#!/usr/bin/env Rscript
# =============================================================================
# make_visiumhd_mini_crc.R
#
# Build a small, format-faithful Visium HD fixture by SPATIALLY CROPPING a
# Space Ranger 4 `outs/` directory (default: 10x's Visium HD human colorectal
# cancer FFPE run, Space Ranger 4.0.1).
#
#   SOURCE_DIR (extracted `outs/`) -> OUT_DIR (an `outs/`-equivalent).
#
# One square window on the 2 um array, placed and sized on the 16 um grid so
# every 8 and 16 um bin it covers is whole. Bins are kept iff their array
# position is inside; cells iff ALL their 2 um bins are inside, so the shipped
# cell matrix stays the exact aggregate of the kept 2 um bins. Fullres pixel
# coordinates are translated so the cropped images start at pixel (0, 0),
# which is where Giotto places them. Barcodes and array indices are unchanged.
#
# Env: VHD_SOURCE_DIR (the extracted 10x `outs/`), VHD_OUT_DIR, and for dry
# runs VHD_SIZE_16 and VHD_PIN ("row16,col16"). The defaults are the paths
# the published fixture was made with.
#
# Requires: R -- hdf5r, arrow, data.table, jsonlite, png, digest; python3
# (standard library only) for visiumhd_crop_geojson.py.
# =============================================================================

suppressMessages({
    library(hdf5r)
    library(data.table)
    library(jsonlite)
})

# ---------------------------------------------------------------- config ----
SCRIPT_PATH <- normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]),
    mustWork = FALSE)
BASE_DIR   <- dirname(SCRIPT_PATH)                      # this crop's folder
SOURCE_DIR <- Sys.getenv("VHD_SOURCE_DIR",
    "/Users/rubendries/Documents/Datasets/visiumHD/Human_Colon_Cancer_v4/outs")
OUT_DIR    <- Sys.getenv("VHD_OUT_DIR",
    "/Users/rubendries/Documents/Datasets/Cropped_datasets/VisiumHD/visiumhd_mini_crc")
SIZE_16    <- as.integer(Sys.getenv("VHD_SIZE_16", "25"))   # 25 x 16 um = 400 um
PIN        <- Sys.getenv("VHD_PIN", "")                     # "" -> search
BINS       <- c(2L, 8L, 16L)
GEOJSONS   <- c("cell_segmentations", "nucleus_segmentations",
                "graphclust_annotated_cell_segmentations",
                "graphclust_annotated_nucleus_segmentations")

msg <- function(...) cat(sprintf(...), "\n", sep = "")

stopifnot(dir.exists(SOURCE_DIR), !dir.exists(OUT_DIR))
S2 <- 8L * SIZE_16                                  # window size, 2 um units
sq <- function(um) sprintf("square_%03dum", um)
bin_dir <- function(root, um) file.path(root, "binned_outputs", sq(um))

# ------------------------------------------------------------- helpers ------

# array row / col from a Visium HD barcode, s_NNNum_RRRRR_CCCCC-1
bc_rc <- function(bc) {
    p <- tstrsplit(sub("-1$", "", bc), "_", fixed = TRUE)
    list(r = as.integer(p[[3]]), c = as.integer(p[[4]]))
}
in_window <- function(bc, um, r0, c0) {
    f <- um %/% 2L
    rc <- bc_rc(bc)
    rc$r >= r0 %/% f & rc$r < (r0 + S2) %/% f &
        rc$c >= c0 %/% f & rc$c < (c0 + S2) %/% f
}

# Column subset of a 10x feature-barcode h5, read in contiguous column runs.
# Datasets keep their source dtypes; features are copied unchanged.
h5_subset <- function(src, dst, keep) {
    s <- H5File$new(src, "r")
    on.exit(s$close_all(), add = TRUE)
    m <- s[["matrix"]]
    bc <- m[["barcodes"]][]
    j <- which(bc %in% keep)
    ip <- as.numeric(m[["indptr"]][])
    runs <- split(j, cumsum(c(1L, diff(j) != 1L)))
    parts <- lapply(runs, function(r) {
        lo <- ip[r[1L]] + 1; hi <- ip[r[length(r)] + 1L]
        if (hi < lo) return(list(i = numeric(), x = numeric()))
        list(i = as.numeric(m[["indices"]][lo:hi]),
             x = as.numeric(m[["data"]][lo:hi]))
    })
    i <- unlist(lapply(parts, `[[`, "i"), use.names = FALSE)
    x <- unlist(lapply(parts, `[[`, "x"), use.names = FALSE)
    p <- c(0, cumsum(ip[j + 1L] - ip[j]))
    nfeat <- m[["shape"]][][1L]

    d <- H5File$new(dst, "w")
    on.exit(d$close_all(), add = TRUE)
    for (a in h5attr_names(s)) {
        sa <- s$attr_open(a)
        d$create_attr(a, robj = sa$read(), dtype = sa$get_type(),
                      space = sa$get_space())
    }
    g <- d$create_group("matrix")
    put <- function(name, val) {
        sd <- m[[name]]
        n <- length(val)
        pl <- H5P_DATASET_CREATE$new()
        if (n > 0) {
            pl$set_chunk(min(65536L, n))
            pl$set_shuffle()
            pl$set_deflate(4L)
        }
        g$create_dataset(name, robj = val, dtype = sd$get_type(),
                         space = H5S$new(dims = n, maxdims = n),
                         dataset_create_pl = pl, chunk_dims = NULL)
    }
    put("barcodes", bc[j])
    put("data", x)
    put("indices", i)
    put("indptr", p)
    put("shape", as.integer(c(nfeat, length(j))))
    d$obj_copy_from(src_loc = s, src_name = "matrix/features",
                    dst_name = "matrix/features")
    list(barcodes = bc[j], i = i, p = p, x = x, nfeat = nfeat)
}

# Matrix-market directory from the subset triplets, with the source header
write_mex <- function(src_dir, dst_dir, sub) {
    dir.create(dst_dir, recursive = TRUE)
    hdr <- readLines(gzfile(file.path(src_dir, "matrix.mtx.gz")), n = 2L)
    stopifnot(grepl("^%%MatrixMarket", hdr[1L]))
    hdr <- hdr[grepl("^%", hdr)]
    n <- length(sub$x)
    tmp <- file.path(dst_dir, "matrix.mtx")
    writeLines(c(hdr, sprintf("%d %d %d", sub$nfeat, length(sub$barcodes), n)), tmp)
    if (n) {
        jj <- rep.int(seq_along(sub$barcodes), diff(sub$p))
        fwrite(data.table(sub$i + 1, jj, sub$x), tmp, append = TRUE, sep = " ",
               col.names = FALSE, scipen = 100)
    }
    writeLines(sub$barcodes, file.path(dst_dir, "barcodes.tsv"))
    system2("gzip", c("-n", "-6", shQuote(tmp), shQuote(file.path(dst_dir, "barcodes.tsv"))))
    file.copy(file.path(src_dir, "features.tsv.gz"), dst_dir)
}

# Both feature matrices (raw, filtered) of an output directory, as h5 + mtx
crop_matrices <- function(src, dst, stem, keep) {
    out <- list()
    for (k in c("raw", "filtered")) {
        nm <- sprintf("%s_%s", k, stem)
        if (!file.exists(file.path(src, paste0(nm, ".h5")))) next
        sub <- h5_subset(file.path(src, paste0(nm, ".h5")),
                         file.path(dst, paste0(nm, ".h5")), keep)
        if (dir.exists(file.path(src, nm))) {
            write_mex(file.path(src, nm), file.path(dst, nm), sub)
        }
        out[[k]] <- list(barcodes = length(sub$barcodes), nnz = length(sub$x),
                         umi = sum(sub$x))
    }
    out
}

crop_png <- function(src, dst, x0, y0, w, h) {
    im <- png::readPNG(src)
    stopifnot(y0 + h <= dim(im)[1L], x0 + w <= dim(im)[2L])
    png::writePNG(im[(y0 + 1L):(y0 + h), (x0 + 1L):(x0 + w), , drop = FALSE], dst)
    c(w = w, h = h)
}

# ------------------------------------------------------------ window ---------
src16 <- bin_dir(SOURCE_DIR, 16L)
tp16 <- as.data.table(arrow::read_parquet(file.path(src16, "spatial", "tissue_positions.parquet")))
if (nzchar(PIN)) {
    w <- as.integer(strsplit(PIN, ",")[[1]])
    W <- list(r16 = w[1L], c16 = w[2L])
} else {
    cl <- fread(file.path(src16, "analysis", "clustering",
                          "gene_expression_graphclust", "clusters.csv"))
    tp16[cl, cluster := i.Cluster, on = c(barcode = "Barcode")]
    best <- NULL
    step <- max(1L, SIZE_16 %/% 2L)
    for (r in seq(0L, max(tp16$array_row) - SIZE_16, by = step)) {
        for (c in seq(0L, max(tp16$array_col) - SIZE_16, by = step)) {
            w <- tp16[array_row >= r & array_row < r + SIZE_16 &
                      array_col >= c & array_col < c + SIZE_16]
            # across the tissue edge, so in- and out-of-tissue bins both occur
            f <- mean(w$in_tissue)
            if (f < 0.6 || f > 0.8) next
            if (!any(!is.na(w$cluster))) next
            sh <- prop.table(table(w$cluster))
            sc <- sum(sh >= 0.05) + (-sum(sh * log(sh)))  / 10
            if (sum(sh >= 0.05) < 3L) next
            if (is.null(best) || sc > best$score) {
                best <- list(r16 = r, c16 = c, score = sc, n_clusters = sum(sh >= 0.05))
            }
        }
    }
    stopifnot(!is.null(best))
    W <- best[c("r16", "c16")]
    msg("[1] search: %d clusters >= 5%% in the best window", best$n_clusters)
}
r0 <- 8L * W$r16
c0 <- 8L * W$c16
msg("[1] window 2um array rows [%d, %d) cols [%d, %d)  (%d um square)",
    r0, r0 + S2, c0, c0 + S2, 2L * S2)

# ------------------------------------------------ pixel offset / image box ---
sf <- fromJSON(file.path(bin_dir(SOURCE_DIR, 2L), "spatial", "scalefactors_json.json"))
hs <- sf$tissue_hires_scalef
ls <- sf$tissue_lowres_scalef
mpp <- sf$microns_per_pixel
tp2 <- as.data.table(arrow::read_parquet(file.path(bin_dir(SOURCE_DIR, 2L), "spatial", "tissue_positions.parquet")))
w2 <- tp2[array_row >= r0 & array_row < r0 + S2 & array_col >= c0 & array_col < c0 + S2]
stopifnot(nrow(w2) == S2^2)
half <- 1 / mpp                                     # 1 um (half a 2 um bin) in px
bb <- list(x0 = min(w2$pxl_col_in_fullres) - half, x1 = max(w2$pxl_col_in_fullres) + half,
           y0 = min(w2$pxl_row_in_fullres) - half, y1 = max(w2$pxl_row_in_fullres) + half)
# hires origin on a multiple of hires/lowres so the lowres crop is whole pixels
ratio <- hs / ls
step_h <- if (abs(ratio - round(ratio)) < 1e-6) round(ratio) else 1
hx0 <- floor(bb$x0 * hs / step_h) * step_h
hy0 <- floor(bb$y0 * hs / step_h) * step_h
hw <- ceiling(bb$x1 * hs) - hx0
hh <- ceiling(bb$y1 * hs) - hy0
dx <- hx0 / hs                                      # fullres translation
dy <- hy0 / hs
lx0 <- round(dx * ls); ly0 <- round(dy * ls)
lw <- ceiling(bb$x1 * ls) - lx0; lh <- ceiling(bb$y1 * ls) - ly0
lowres_exact <- abs(dx * ls - lx0) < 1e-6 && abs(dy * ls - ly0) < 1e-6
msg("[1] fullres translation (-%.4f, -%.4f) px; hires box %dx%d at (%d, %d); lowres %dx%d%s",
    dx, dy, hw, hh, hx0, hy0, lw, lh, if (lowres_exact) "" else " (sub-pixel shift)")

# --------------------------------------------------------------- crop --------
dir.create(OUT_DIR, recursive = TRUE)
counts <- list()

crop_spatial <- function(src, dst) {
    dir.create(dst, recursive = TRUE)
    file.copy(file.path(src, "scalefactors_json.json"), dst)
    crop_png(file.path(src, "tissue_hires_image.png"),
             file.path(dst, "tissue_hires_image.png"), hx0, hy0, hw, hh)
    crop_png(file.path(src, "tissue_lowres_image.png"),
             file.path(dst, "tissue_lowres_image.png"), lx0, ly0, lw, lh)
    invisible()
}

for (um in BINS) {
    src <- bin_dir(SOURCE_DIR, um)
    dst <- bin_dir(OUT_DIR, um)
    dir.create(dst, recursive = TRUE)
    tbl <- arrow::read_parquet(file.path(src, "spatial", "tissue_positions.parquet"),
                               as_data_frame = FALSE)
    tp <- as.data.table(tbl)
    tp <- tp[in_window(barcode, um, r0, c0)]
    stopifnot(nrow(tp) == (S2 %/% (um %/% 2L))^2)
    tp[, `:=`(pxl_row_in_fullres = pxl_row_in_fullres - dy,
              pxl_col_in_fullres = pxl_col_in_fullres - dx)]
    crop_spatial(file.path(src, "spatial"), file.path(dst, "spatial"))
    arrow::write_parquet(arrow::Table$create(tp, schema = tbl$schema),
                         file.path(dst, "spatial", "tissue_positions.parquet"))
    counts[[sq(um)]] <- c(list(bins = nrow(tp), in_tissue = sum(tp$in_tissue)),
                          crop_matrices(src, dst, "feature_bc_matrix", tp$barcode))
    msg("[2] %s: %d bins (%d in tissue), raw nnz %d, umi %d", sq(um), nrow(tp),
        sum(tp$in_tissue), counts[[sq(um)]]$raw$nnz, counts[[sq(um)]]$raw$umi)
}

# cells whose every 2 um bin is inside the window
bm_tbl <- arrow::read_parquet(file.path(SOURCE_DIR, "barcode_mappings.parquet"),
                              as_data_frame = FALSE)
bm <- as.data.table(bm_tbl)
bm[, inside := in_window(square_002um, 2L, r0, c0)]
cells_all <- bm[!is.na(cell_id), .(all_in = all(inside), any_in = any(inside)), by = cell_id]
keep_cells <- cells_all[all_in == TRUE, cell_id]
msg("[3] cells touching the window %d, wholly inside %d",
    cells_all[any_in == TRUE, .N], length(keep_cells))

bmw <- bm[inside == TRUE]
bmw[!cell_id %in% keep_cells, `:=`(cell_id = NA_character_, in_nucleus = FALSE, in_cell = FALSE)]
bmw[, inside := NULL]
arrow::write_parquet(arrow::Table$create(bmw, schema = bm_tbl$schema),
                     file.path(OUT_DIR, "barcode_mappings.parquet"))

seg_src <- file.path(SOURCE_DIR, "segmented_outputs")
seg_dst <- file.path(OUT_DIR, "segmented_outputs")
dir.create(seg_dst)
crop_spatial(file.path(seg_src, "spatial"), file.path(seg_dst, "spatial"))
counts$segmented_outputs <- crop_matrices(seg_src, seg_dst, "feature_cell_matrix", keep_cells)

ids_txt <- tempfile(fileext = ".txt")
writeLines(as.character(as.integer(sub("^cellid_0*([0-9]+)-1$", "\\1", keep_cells))), ids_txt)
for (gj in GEOJSONS) {
    f <- file.path(seg_src, paste0(gj, ".geojson"))
    if (!file.exists(f)) next
    r <- system2("python3", c(shQuote(file.path(BASE_DIR, "visiumhd_crop_geojson.py")),
        shQuote(f), shQuote(file.path(seg_dst, paste0(gj, ".geojson"))),
        shQuote(ids_txt), format(dx, digits = 15), format(dy, digits = 15)), stdout = TRUE)
    counts$geojson[[gj]] <- fromJSON(r)
}
msg("[3] segmented: %d cells, raw nnz %d; geojson kept %s", counts$segmented_outputs$raw$barcodes,
    counts$segmented_outputs$raw$nnz,
    paste(vapply(counts$geojson, function(x) x$features_kept, 1), collapse = "/"))

# ---------------------------------------------------------- manifest ---------
fsize <- function(p) if (!file.exists(p)) 0 else if (dir.exists(p))
    sum(file.size(list.files(p, recursive = TRUE, full.names = TRUE))) else file.size(p)
act <- function(file, action, reason = NULL)
    c(list(file = file, action = action),
      if (!is.null(reason)) list(reason = reason) else NULL,
      list(bytes_source = unname(fsize(file.path(SOURCE_DIR, file))),
           bytes_fixture = unname(fsize(file.path(OUT_DIR, file)))))

actions <- c(
    lapply(BINS, function(um) act(file.path("binned_outputs", sq(um)), "cropped")),
    list(
        act("segmented_outputs", "cropped"),
        act("barcode_mappings.parquet", "cropped"),
        act("binned_outputs/*/analysis, segmented_outputs/analysis", "dropped",
            "run-level PCA / clustering / DE that would be stale after cropping"),
        act("*/spatial/{cytassist_image.tiff, aligned_fiducials.jpg, aligned_tissue_image.jpg, detected_tissue_image.jpg}",
            "dropped", "QC images, never read by Giotto"),
        act("*/cloupe*.cloupe", "dropped", "Loupe browser files, never read by Giotto"),
        act("spatial", "dropped", "top-level QC images; readers use each output's own spatial/"),
        act("binned_outputs/square_002um/raw_probe_bc_matrix.h5", "dropped",
            "probe-level counts (FFPE probe set), never read by Giotto"),
        act("feature_slice.h5", "dropped", "never read by Giotto"),
        act("molecule_info.h5", "dropped", "never read by Giotto"),
        act("metrics_summary.csv", "dropped", "run-level QC that would be stale"),
        act("web_summary.html", "dropped", "run-level report that would be stale")
    )
)

src_metrics <- tryCatch(fread(list.files(dirname(SOURCE_DIR), "metrics_summary.csv$",
                                         full.names = TRUE, recursive = TRUE)[1]),
                        error = function(e) NULL)
h5a <- local({
    s <- H5File$new(file.path(bin_dir(SOURCE_DIR, 8L), "raw_feature_bc_matrix.h5"), "r")
    on.exit(s$close_all()); h5attributes(s)
})

manifest <- list(
    fixture_name = basename(OUT_DIR),
    created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    source = list(
        dir = SOURCE_DIR,
        sample_id = if (!is.null(src_metrics)) src_metrics[["Sample ID"]][1] else NA,
        software_version = h5a$software_version,
        chemistry = h5a$chemistry_description,
        microns_per_pixel = mpp,
        scalefactors = sf
    ),
    crop = list(
        method = paste("single square window on the 2 um array, placed and sized on the 16 um grid,",
                       "across the tissue edge (60-80% of 16 um bins in tissue);",
                       "bins kept iff their array position is inside, cells kept iff ALL their 2 um bins are inside"),
        window_array_2um = list(row_min = r0, row_max = r0 + S2, col_min = c0, col_max = c0 + S2),
        window_array_16um = list(row_min = W$r16, row_max = W$r16 + SIZE_16,
                                 col_min = W$c16, col_max = W$c16 + SIZE_16),
        size_um = 2L * S2,
        translation_fullres_px = list(x = -dx, y = -dy),
        hires_box_px = list(x0 = hx0, y0 = hy0, width = hw, height = hh),
        lowres_box_px = list(x0 = lx0, y0 = ly0, width = lw, height = lh,
                             exact = lowres_exact)
    ),
    counts = counts,
    rebuilt = list(
        "feature-barcode h5: barcodes / data / indices / indptr / shape rewritten for the kept columns in source order, source dtypes kept; features/ copied unchanged (all genes kept); root attributes copied",
        "matrix-market directories written from the same triplets with the source header lines; features.tsv.gz copied unchanged",
        "tissue_positions.parquet: kept rows in source order, pxl_row/col_in_fullres translated; same schema",
        "barcode_mappings.parquet: rows whose 2 um bin is inside; cell_id / in_cell / in_nucleus cleared for cells that are not wholly inside (they are not in the fixture)",
        "segmentation geojson: features of kept cells, vertices translated; properties unchanged",
        "tissue_hires/lowres_image.png cropped to the window; scalefactors_json.json unchanged"
    ),
    deviations = list(
        "HDF5 datasets are written with shuffle + deflate level 4 and 65536-element chunks; the source's exact filter settings are not reproduced. Transparent to readers.",
        "Fullres pixel coordinates are translated (see crop$translation_fullres_px); barcodes and array_row/array_col keep their source values, so the array grid does not start at 0.",
        "Geojson vertices are rounded to 6 decimals after translation."
    ),
    files = actions,
    # repo-relative: the scripts live in drieslab/spatial-datasets
    script = list(
        repo = "https://github.com/drieslab/spatial-datasets",
        driver = "crops/visiumhd_mini_crc/make_visiumhd_mini_crc.R",
        module = "crops/visiumhd_mini_crc/visiumhd_crop_geojson.py",
        sha256 = list(
            driver = if (file.exists(SCRIPT_PATH)) digest::digest(file = SCRIPT_PATH, algo = "sha256") else NA_character_,
            module = digest::digest(file = file.path(BASE_DIR, "visiumhd_crop_geojson.py"), algo = "sha256")
        )
    )
)
write(toJSON(manifest, auto_unbox = TRUE, pretty = 4, digits = NA, null = "null"),
      file.path(OUT_DIR, "crop_manifest.json"))

total <- sum(file.size(list.files(OUT_DIR, recursive = TRUE, full.names = TRUE)))
msg("[4] crop_manifest.json written")
msg("[done] %s  total %.2f MB", OUT_DIR, total / 1024^2)
