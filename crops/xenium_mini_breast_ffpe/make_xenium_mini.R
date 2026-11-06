#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Build the Xenium mini fixture: a square crop of the S1_Bot breast section,
# written in 10x's LEGACY parquet/MEX layout so it loads under stock
# Giotto::createGiottoXeniumObject() with no GiottoDisk and no zarr reader.
#
#   Rscript make_xenium_mini.R [window_um]
#
# Legacy, not zarr: a fixture is only useful if it loads on a released stack.
# Subsetting a zarr export would mean reproducing 10x's 7-level tile grid and
# reading it back would need GiottoDisk dev plus Rarr.
#
# ACCEPTED TRADE-OFF: the mini does NOT exercise the zarr reader. Say so when
# handing it to anyone; the full tutorial covers that path.
#
# The crop is TRANSLATED TO THE ORIGIN, spanning [0, window] in both axes.
# That is not cosmetic: readers anchor a morphology image at (0, 0), so a
# cropped image beside untranslated coordinates is misregistered by exactly
# the crop offset. crop_manifest.json records the shift under crop$translation.
#
# What it contains, and how each file was made, is recorded per-file in the
# crop_manifest.json this script writes. Short version: cells, cell and
# nucleus boundaries, transcripts, the vendor MEX matrix, one morphology
# channel (DAPI), gene_panel.json and a crop-scoped experiment.xenium.
# ---------------------------------------------------------------------------
suppressMessages({
  library(arrow); library(GiottoDisk); library(data.table)
  library(Matrix); library(jsonlite)
})

HERE   <- "/projectnb/rd-spat/HOME/vjarzabe/Xenium_Tutorial/scc_xenium_tutorial"
XD     <- Sys.getenv("XENIUM_ZARR_DIR", file.path(HERE, "xenium_zarr"))
OUTDIR <- Sys.getenv("XENIUM_MINI_OUTDIR",
                     file.path(HERE, "mini", "data", "xenium_mini_breast_ffpe"))
WORK   <- file.path(HERE, "mini", ".work")
SIZE_CAP_MB <- 50

# The morphology channel to ship, and the interpreter that can decode it.
# Only DAPI: all four channels re-encode to 89.1 MB, which would put the
# archive at 123 MB, over the 100 MB ceiling. DAPI alone is 23.6 MB.
# crop_morphology.py explains why this cannot be done in R.
MORPH_CH <- "ch0000_dapi.ome.tif"
PYTHON   <- Sys.getenv("XENIUM_MINI_PYTHON",
                       "/share/pkg.8/python3/3.13.8/install/bin/python3")

# The zarr reader fixes this script used to patch in locally are upstream as of
# GiottoDisk #86 (2026-10-02), so nothing is sourced here any more.

args <- commandArgs(trailingOnly = TRUE)
WINDOW_UM <- if (length(args)) as.numeric(args[[1]]) else 1200

dir.create(WORK, recursive = TRUE, showWarnings = FALSE)
# The shipped fixture is a published artifact whose sha256 is quoted in the
# GiottoData manifest, so refuse to overwrite one. Build elsewhere with
# XENIUM_MINI_OUTDIR, or set XENIUM_MINI_OVERWRITE=1 and re-record the hash.
if (dir.exists(OUTDIR) && length(list.files(OUTDIR)) &&
    !nzchar(Sys.getenv("XENIUM_MINI_OVERWRITE"))) {
  stop(sprintf(paste0(
    "refusing to overwrite an existing fixture at\n  %s\n",
    "set XENIUM_MINI_OUTDIR to build elsewhere, or XENIUM_MINI_OVERWRITE=1 ",
    "to replace it."), OUTDIR), call. = FALSE)
}

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
ns <- asNamespace("GiottoDisk")

step <- function(...) cat(sprintf("[mini] %s\n", paste0(...)))

# --- 1. full-section cells + boundaries as parquet -------------------------
# xeniumZarrToParquet() writes them under 10x's own legacy names, so these are
# the files that get cropped rather than anything bespoke.
if (!file.exists(file.path(WORK, "cells.parquet"))) {
  step("converting cells + cell_boundaries from zarr")
  xeniumZarrToParquet(xenium_dir = XD, output_dir = WORK,
                      what = c("cells", "cell_boundaries"),
                      overwrite = TRUE, verbose = TRUE)
} else step("reusing cached cells/cell_boundaries in .work/")

cells <- as.data.table(read_parquet(file.path(WORK, "cells.parquet")))
step(sprintf("full section: %s cells, x %.0f-%.0f um, y %.0f-%.0f um",
             format(nrow(cells), big.mark = ","),
             min(cells$x_centroid), max(cells$x_centroid),
             min(cells$y_centroid), max(cells$y_centroid)))

# --- 2. vendor cell x gene matrix, in memory -------------------------------
# 280 genes x 254,480 cells of int counts is small; the transcripts are the
# only thing here that needs streaming.
step("reading vendor cell_feature_matrix from zarr")
cfm_path <- ns$.zarr_find_product(XD, "cell_feature_matrix")
inp <- tenxZarrInput(cfm_path, feature_id_col = 2L)
itr <- storeRead(inp)
trip <- list(); repeat { b <- itr$next_batch(); if (is.null(b)) break
                         if (nrow(b)) trip[[length(trip) + 1L]] <- b }
itr$close()
trip <- rbindlist(trip)
cfm <- sparseMatrix(i = trip$col_id, j = trip$row_id, x = trip$value,
                    dims = c(inp@n_genes, inp@n_cells),
                    dimnames = list(inp@feat_ids, inp@cell_ids))
feat_types_all <- inp@feat_types
step(sprintf("vendor matrix: %d features x %s cells, %s nonzero",
             nrow(cfm), format(ncol(cfm), big.mark = ","),
             format(length(cfm@x), big.mark = ",")))

# --- 3. choose the window --------------------------------------------------
# Scored on cell count AND lineage diversity, so the fixture contains tumour,
# stroma and immune rather than a dense block of one thing. Each cell is given
# a crude argmax label over a few canonical panel genes -- good enough to pick
# a window, and never written out.
probe <- list(
  epithelial = c("KRT8", "KRT18", "EPCAM", "ELF3", "FOXA1"),
  tumor      = c("ERBB2", "MKI67", "TOP2A", "CCNE1"),
  stroma     = c("COL1A1", "COL1A2", "DCN", "LUM", "POSTN"),
  vascular   = c("PECAM1", "VWF", "CDH5", "RGS5", "ACTA2"),
  immune     = c("PTPRC", "CD68", "CD3D", "MS4A1", "IGKC", "LYZ"))
probe <- lapply(probe, function(g) intersect(g, rownames(cfm)))
probe <- probe[lengths(probe) > 0]

score_mat <- sapply(probe, function(g)
  Matrix::colSums(cfm[g, , drop = FALSE]))
lab <- colnames(score_mat)[max.col(score_mat, ties.method = "first")]
lab[rowSums(score_mat) == 0] <- NA_character_
cells[, probe_label := lab[match(cell_id, colnames(cfm))]]

cells[, `:=`(gx = floor(x_centroid / WINDOW_UM),
             gy = floor(y_centroid / WINDOW_UM))]
shannon <- function(x) { x <- x[!is.na(x)]; if (!length(x)) return(0)
                         p <- prop.table(table(x)); -sum(p * log(p)) }
grid <- cells[, .(n = .N, H = shannon(probe_label)), by = .(gx, gy)]
grid <- grid[n >= 500]
stopifnot("no candidate window has >= 500 cells" = nrow(grid) > 0)
# rank on diversity first, cell count as the tiebreak
grid[, rank_score := scale(H)[, 1] + 0.5 * scale(log10(n))[, 1]]
setorder(grid, -rank_score)
win <- grid[1]

XMIN <- win$gx * WINDOW_UM; XMAX <- XMIN + WINDOW_UM
YMIN <- win$gy * WINDOW_UM; YMAX <- YMIN + WINDOW_UM
step(sprintf("window %.0f um at x [%.0f, %.0f], y [%.0f, %.0f]: %d cells, H = %.2f",
             WINDOW_UM, XMIN, XMAX, YMIN, YMAX, win$n, win$H))
print(cells[gx == win$gx & gy == win$gy, .N, by = probe_label][order(-N)])

# --- 4. crop ---------------------------------------------------------------
keep_cells <- cells[x_centroid >= XMIN & x_centroid < XMAX &
                    y_centroid >= YMIN & y_centroid < YMAX]
keep_ids <- keep_cells$cell_id
step(sprintf("kept %s cells", format(length(keep_ids), big.mark = ",")))

out_cells <- copy(keep_cells)
out_cells[, c("probe_label", "gx", "gy") := NULL]
# Shift the crop to the origin. A reader anchors every image at (0, 0)
# (`rescale(img, micron, x0 = 0, y0 = 0)`) and has no way to learn an offset,
# so source-frame coordinates would put the image 2400 um left and 3600 um
# above the cells.
out_cells[, `:=`(x_centroid = x_centroid - XMIN,
                 y_centroid = y_centroid - YMIN,
                 nucleus_centroid_x = nucleus_centroid_x - XMIN,
                 nucleus_centroid_y = nucleus_centroid_y - YMIN)]
write_parquet(out_cells, file.path(OUTDIR, "cells.parquet"), compression = "zstd")

bounds <- as.data.table(read_parquet(file.path(WORK, "cell_boundaries.parquet")))
bounds <- bounds[cell_id %chin% keep_ids]
bounds[, `:=`(vertex_x = vertex_x - XMIN, vertex_y = vertex_y - YMIN)]
write_parquet(bounds, file.path(OUTDIR, "cell_boundaries.parquet"),
              compression = "zstd")
step(sprintf("boundaries: %s vertex rows over %s cells",
             format(nrow(bounds), big.mark = ","),
             format(uniqueN(bounds$cell_id), big.mark = ",")))

# Nucleus boundaries come from cells.zarr polygon_sets/0: the zarr->parquet
# conversion emits cell rings only. The published xenium_mini_lung_ffpe ships
# this file, so this one does.
# A cell may own several nucleus rings, and Giotto then falls back to
# synthetic poly_IDs, so the layer does not join to cell_id. Upstream
# behaviour; the lung fixture does the same. Quantified in the manifest.
nuc_stats_json <- file.path(WORK, "nucleus_stats.json")
rc <- system2(PYTHON, c("-I", shQuote(file.path(HERE, "mini", "extract_nucleus.py")),
                        "--cells-zarr", shQuote(file.path(XD, "cells.zarr.zip")),
                        "--keep-cells", shQuote(file.path(OUTDIR, "cells.parquet")),
                        "--xmin", XMIN, "--ymin", YMIN,
                        "--out", shQuote(file.path(OUTDIR, "nucleus_boundaries.parquet")),
                        "--stats-out", shQuote(nuc_stats_json)))
if (rc != 0L) stop("extract_nucleus.py failed (exit ", rc, ")", call. = FALSE)
nuc_stats <- fromJSON(nuc_stats_json, simplifyVector = TRUE)
step(sprintf("nucleus: %s rings over %s cells",
             format(nuc_stats$rings, big.mark = ","),
             format(nuc_stats$cells_with_nucleus, big.mark = ",")))

# Transcripts are cropped on COORDINATE, not on cell_id: the transcripts zarr
# stores no per-transcript assignment at all (its tile arrays are
# codeword_identity, gene_identity, gene_offset, id, location, quality_score,
# status, uuid, valid). The assignment is recovered from the segmentation
# masks in step 4b below -- see recover_assignment.py.
TXC <- file.path(HERE, "full", "vault", "artifacts", "zarr_cache")
tx_dir <- list.files(TXC, pattern = "^transcripts$", recursive = TRUE,
                     include.dirs = TRUE, full.names = TRUE)
stopifnot("no cached transcripts conversion found -- run verify_ingest.R first" =
            length(tx_dir) > 0)
tx <- as.data.table(dplyr::collect(dplyr::filter(
  arrow::open_dataset(tx_dir[[1]]),
  x_location >= XMIN, x_location < XMAX,
  y_location >= YMIN, y_location < YMAX)))
# transcript_id breaks ties so the row order -- and therefore the file's
# sha256 -- is reproducible. 4,683 detections share an exact (x, y) with
# another, and without a unique third key arrow's fragment read order decides
# which comes first, so two builds of identical data hashed differently.
tx[, `:=`(x_location = x_location - XMIN, y_location = y_location - YMIN)]
setorder(tx, y_location, x_location, transcript_id)
write_parquet(tx, file.path(OUTDIR, "transcripts.parquet"), compression = "zstd")
step(sprintf("transcripts: %s detections, %d distinct features",
             format(nrow(tx), big.mark = ","), uniqueN(tx$feature_name)))

# --- 4b. recover per-transcript cell assignment ----------------------------
# Without this every transcript ships UNASSIGNED, which is not what a Xenium
# export looks like (the vendor reports fraction_transcripts_assigned 0.759
# for this section). cells.zarr carries masks/1, the cell-segmentation label
# raster whose value at a pixel is the 1-based row index into the cells array;
# masks/0 gives overlaps_nucleus the same way. Point-in-polygon is 10x further
# off (0.962 against 0.996). Step 5a measures the result. nucleus_distance
# stays NA: absent from the zarr and not documented precisely enough.
tx_path <- file.path(OUTDIR, "transcripts.parquet")
# The temp sits beside the output, not in WORK: file.rename() cannot cross a
# filesystem, and it returns FALSE rather than throwing, so the build once
# carried on and shipped the unassigned file.
tx_tmp  <- file.path(OUTDIR, ".transcripts_assigned.parquet")
asg_json <- file.path(WORK, "assignment_stats.json")
rc <- system2(PYTHON, c("-I", shQuote(file.path(HERE, "mini", "recover_assignment.py")),
                        "--cells-zarr", shQuote(file.path(XD, "cells.zarr.zip")),
                        "--transcripts", shQuote(tx_path),
                        "--keep-cells",  shQuote(file.path(OUTDIR, "cells.parquet")),
                        "--xmin", XMIN, "--ymin", YMIN, "--window", WINDOW_UM,
                        "--out", shQuote(tx_tmp), "--stats-out", shQuote(asg_json)))
if (rc != 0L) stop("recover_assignment.py failed (exit ", rc, ")", call. = FALSE)
stopifnot("recover_assignment.py wrote no output" = file.exists(tx_tmp))
asg_stats <- fromJSON(asg_json, simplifyVector = TRUE)
if (!file.rename(tx_tmp, tx_path)) {
    stopifnot("could not move the recovered transcripts into place" =
                file.copy(tx_tmp, tx_path, overwrite = TRUE))
    unlink(tx_tmp)
}

# Re-read the two recovered columns so every count below is taken from the
# file that actually ships, not from the pre-assignment table in memory.
recovered <- as.data.table(read_parquet(tx_path,
                                        col_select = c("cell_id", "overlaps_nucleus")))
stopifnot("recovered transcript count changed" = nrow(recovered) == nrow(tx))
tx[, `:=`(cell_id = recovered$cell_id, overlaps_nucleus = recovered$overlaps_nucleus)]
frac_assigned <- mean(tx$cell_id != "UNASSIGNED")
step(sprintf("assignment: %s of %s transcripts in a cell (%.1f%%), %s nucleus-overlapping",
             format(sum(tx$cell_id != "UNASSIGNED"), big.mark = ","),
             format(nrow(tx), big.mark = ","), 100 * frac_assigned,
             format(sum(tx$overlaps_nucleus > 0), big.mark = ",")))
# Shipping an all-UNASSIGNED fixture is the exact defect this step exists to
# prevent, so refuse to build one. The vendor reports 0.759 for the section.
if (frac_assigned < 0.5)
    stop(sprintf("only %.1f%% of transcripts were assigned -- recovery did not land",
                 100 * frac_assigned), call. = FALSE)

# --- 5. vendor matrix for the kept cells, as a MEX triplet -----------------
sub <- cfm[, keep_ids, drop = FALSE]
mex <- file.path(OUTDIR, "cell_feature_matrix")
dir.create(mex, recursive = TRUE, showWarnings = FALSE)
mm <- file.path(mex, "matrix.mtx")
invisible(writeMM(as(sub, "dgCMatrix"), mm))
R.utils::gzip(mm, destname = paste0(mm, ".gz"), overwrite = TRUE, remove = TRUE)
fwrite(data.table(id = rownames(sub), name = rownames(sub),
                  type = ifelse(feat_types_all == "rna", "Gene Expression",
                                feat_types_all)),
       file.path(mex, "features.tsv.gz"), sep = "\t", col.names = FALSE)
fwrite(data.table(keep_ids), file.path(mex, "barcodes.tsv.gz"),
       col.names = FALSE)
step(sprintf("MEX: %d features x %s cells", nrow(sub),
             format(ncol(sub), big.mark = ",")))

# --- 5a. measure the recovery against the vendor's own matrix --------------
# Per-cell TOTAL detections at qv >= 20 over all features: the filter the
# vendor matrix applies. Counting every detection gives ratio 1.21, which
# measures the qv filter instead. Written into crop_manifest.json.
val <- local({
    a <- tx[qv >= 20 & cell_id != "UNASSIGNED", .N, by = cell_id]
    v <- Matrix::colSums(sub)[a$cell_id]
    list(definition = "per-cell total detections, qv >= 20, all features",
         cells_compared = nrow(a),
         overall_ratio = sum(a$N) / sum(v),
         frac_cells_exact = mean(a$N == v),
         median_cell_ratio = median(a$N / v),
         pearson_r = cor(a$N, v))
})
step(sprintf("vendor check: ratio %.5f, %.1f%% of cells exact, median %.5f, r %.6f",
             val$overall_ratio, 100 * val$frac_cells_exact,
             val$median_cell_ratio, val$pearson_r))

# --- 5b. morphology image --------------------------------------------------
# Cropped and re-encoded by crop_morphology.py, which needs tifffile and
# imagecodecs. If that interpreter is not available the build still produces a
# usable fixture, just without the image, so this warns rather than stops.
morph_src <- file.path(XD, "morphology_focus", MORPH_CH)
morph_out <- file.path(OUTDIR, "morphology_focus", MORPH_CH)
morph_px  <- fromJSON(file.path(XD, "experiment.xenium"))$pixel_size
have_morph <- FALSE
if (!file.exists(morph_src)) {
  warning("no morphology source at ", morph_src, "; fixture ships no image",
          call. = FALSE)
} else {
  step("cropping + re-encoding ", MORPH_CH)
  morph_geom_json <- file.path(WORK, "morph_geom.json")
  rc <- system2(PYTHON, c("-I", shQuote(file.path(HERE, "mini",
                                                  "crop_morphology.py")),
                          "--src", shQuote(morph_src),
                          "--out", shQuote(morph_out),
                          "--xmin", XMIN, "--xmax", XMAX,
                          "--ymin", YMIN, "--ymax", YMAX,
                          "--pixel-size", morph_px,
                          "--geom-out", shQuote(morph_geom_json)))
  have_morph <- rc == 0L && file.exists(morph_out)
  if (have_morph) morph_geom <- fromJSON(morph_geom_json, simplifyVector = TRUE)
  if (!have_morph) {
    warning("crop_morphology.py failed (status ", rc,
            "); fixture ships no image", call. = FALSE)
  }
}

# --- 6. sidecars -----------------------------------------------------------
invisible(file.copy(file.path(XD, "gene_panel.json"),
                    file.path(OUTDIR, "gene_panel.json"), overwrite = TRUE))
xe <- fromJSON(file.path(XD, "experiment.xenium"), simplifyVector = TRUE)

# experiment.xenium stays a VENDOR file: correct the fields that no longer
# describe the shipped data, drop the ones naming files the fixture does not
# contain, add nothing. A field that is now false is worse than one that is
# absent. Crop provenance goes in crop_manifest.json, not here. Applied by
# patch_experiment.py because jsonlite would re-round every float. The method
# and area mappings were verified by reproducing the vendor's own
# section-level values from the full cells.parquet first.
meth <- keep_cells$segmentation_method
frac_meth <- function(pattern) sum(grepl(pattern, meth, fixed = TRUE)) / length(meth)
nuc_meth <- nuc_stats$methods
nuc_named <- function(pattern) {
  hit <- names(nuc_meth)[grepl(pattern, names(nuc_meth), fixed = TRUE)]
  if (!length(hit)) 0L else sum(unlist(nuc_meth[hit]))
}
nuc_total <- max(1L, sum(unlist(nuc_meth)))

xe_patch <- list(
  set = list(
    num_cells = length(keep_ids),
    num_transcripts = nrow(tx),
    region_name = paste0(xe$region_name, "_mini"),
    fraction_transcripts_assigned = sum(tx$cell_id != "UNASSIGNED") / nrow(tx),
    total_cell_area = sum(keep_cells$cell_area),
    non_zero_matrix_entries = length(sub@x),
    segmented_cell_boundary_frac = frac_meth("boundary stain (ATP1A1"),
    segmented_cell_interior_frac = frac_meth("interior stain"),
    segmented_cell_nuc_expansion_frac = frac_meth("nucleus expansion"),
    segmented_cell_stain_frac = frac_meth("boundary stain (ATP1A1") +
                                frac_meth("interior stain"),
    segmented_nucleus_focus_count = nuc_named("multi-focus"),
    segmented_nucleus_focus_frac = nuc_named("multi-focus") / nuc_total,
    segmented_nucleus_zstack_count = nuc_named("z stack"),
    segmented_nucleus_zstack_frac = nuc_named("z stack") / nuc_total,
    images = if (have_morph) {
      # name only what the fixture ships: the 3D stack is not included
      list(morphology_focus_filepath = file.path("morphology_focus", MORPH_CH))
    } else NULL
  ),
  # I() keeps this an array: auto_unbox would collapse a one-element vector
  # to a bare string, and the patcher would then iterate its characters and
  # silently delete nothing.
  delete = I(c(
    "xenium_explorer_files",          # names files the fixture does not ship
    # Section-scoped and NOT derivable for a crop. region_area is a tissue-mask
    # area with no crop equivalent, and the two rates are counts per unit of
    # it. num_transcripts_high_quality and transcripts_per_cell have no
    # definition that reproduces the vendor's own numbers from the files, so
    # recomputing them would be a guess; thickness_... is a run-level optic.
    "region_area",
    "transcripts_per_100um",
    "nuclear_transcripts_per_100um",
    "num_transcripts_high_quality",
    "thickness_of_high_quality_decoded_transcripts",
    "transcripts_per_cell",
    if (!have_morph) "images"))
)
xe_patch$set <- Filter(Negate(is.null), xe_patch$set)
patch_json <- file.path(WORK, "experiment_patch.json")
write(toJSON(xe_patch, auto_unbox = TRUE, digits = I(17)), patch_json)
rc <- system2(PYTHON, c("-I", shQuote(file.path(HERE, "mini",
                                                "patch_experiment.py")),
                        "--src", shQuote(file.path(XD, "experiment.xenium")),
                        "--out", shQuote(file.path(OUTDIR, "experiment.xenium")),
                        "--patch", shQuote(patch_json)))
stopifnot("patch_experiment.py failed" = rc == 0L)
xe <- fromJSON(file.path(OUTDIR, "experiment.xenium"), simplifyVector = TRUE)

# crop_manifest.json, GiottoData crop manifest SCHEMA V1:
#   giotto-suite/GiottoData .claude/skills/add-spatial-dataset/crop_manifest.schema.json
# Check a built archive against it before publishing anything:
#   Rscript check_crop_manifest.R xenium_mini_breast_ffpe.zip xenium_mini_breast_ffpe
# The top level is CLOSED (additionalProperties false), so only the keys the
# schema names may appear; source, crop and counts accept platform extras.
# The schema deliberately does not want panel, slide or pixel size repeated
# here -- they are already in the experiment.xenium this archive carries.
xe_src  <- fromJSON(file.path(XD, "experiment.xenium"), simplifyVector = TRUE)
bfiles  <- list.files(OUTDIR, recursive = TRUE, full.names = TRUE)
relpath <- function(f) sub(paste0(OUTDIR, "/"), "", f, fixed = TRUE)
ssize   <- function(...) unname(file.size(file.path(XD, ...)))
fsize   <- function(rel) {
  f <- file.path(OUTDIR, rel)
  if (dir.exists(f)) sum(file.size(list.files(f, recursive = TRUE,
                                              full.names = TRUE)))
  else unname(file.size(f))
}
# bytes_source is null wherever no single source file corresponds to the
# shipped one. Three of these tables come out of the one cells.zarr.zip, and
# naming its size three times would read as three source files.
fentry <- function(file, action, bytes_source, reason = NULL) {
  c(list(file = file, action = action),
    if (!is.null(reason)) list(reason = reason),
    list(bytes_source = bytes_source, bytes_fixture = fsize(file)))
}
morph_other <- setdiff(list.files(file.path(XD, "morphology_focus")), MORPH_CH)

manifest <- list(
  schema_version = "1",
  fixture_name = basename(OUTDIR),
  # SOURCE_DATE_EPOCH pins this so two builds of the same data are
  # bit-identical; left unset it records the real build time.
  created = local({
    sde <- Sys.getenv("SOURCE_DATE_EPOCH")
    tm <- if (nzchar(sde)) as.POSIXct(as.numeric(sde), origin = "1970-01-01",
                                      tz = "UTC") else Sys.time()
    format(tm, "%Y-%m-%dT%H:%M:%S%z")
  }),
  source = list(
    technology = "xenium",
    organism = tolower(xe_src$panel_organism),
    tissue = "breast cancer",
    preservation = "FFPE",
    dataset = xe_src$run_name,
    run_id = xe_src$region_name,
    # the public 10x dataset page for the full run. The download is behind a
    # free 10x registration, which the schema treats as public: the point of
    # this field is that someone can get the run the crop was cut from.
    url = "https://www.10xgenomics.com/datasets/xenium-ffpe-human-breast-biomarkers",
    software = list(name = "xenium onboard analysis",
                    version = xe_src$analysis_sw_version),
    format = "zarr export (.zarr.zip), not a legacy _outs bundle",
    num_cells = nrow(cells)),
  crop = list(
    method = paste("single square window placed by hand on the section;",
                   "cells kept iff their centroid is strictly inside, so a",
                   "cell straddling the edge is either wholly in or wholly",
                   "out and its boundary is never clipped"),
    window = list(units = "um", x_min = XMIN, x_max = XMAX,
                  y_min = YMIN, y_max = YMAX),
    translation = list(units = "um", x = -XMIN, y = -YMIN),
    size_um = WINDOW_UM,
    fixture_extent_um = list(x_min = 0, x_max = WINDOW_UM,
                             y_min = 0, y_max = WINDOW_UM),
    # measured from the image that was written, not recomputed: the two
    # disagree by a pixel, and a declared size the file does not have is the
    # kind of thing nobody checks until it misregisters something
    window_px_source = if (have_morph) {
      morph_geom[c("x_min", "x_max", "y_min", "y_max")]
    } else NULL,
    image_size_px = if (have_morph) {
      list(width = morph_geom$width, height = morph_geom$height)
    } else NULL),
  counts = list(
    cells = length(keep_ids),
    # the whole feature axis, as the MEX ships it: panel genes plus every
    # control row. `genes` is the panel size a reader means by "280-gene".
    features = nrow(sub),
    # inp@feat_types already holds the vendor strings ("Gene Expression",
    # "Negative Control Probe", ...), which is why the ifelse above it is a
    # no-op. Accept both spellings so this does not silently read 0 again.
    genes = sum(feat_types_all %in% c("rna", "Gene Expression")),
    features_detected = uniqueN(tx$feature_name),
    transcripts = nrow(tx),
    transcripts_assigned = asg_stats$assigned,
    transcripts_unassigned = nrow(tx) - asg_stats$assigned,
    transcripts_in_cell_per_mask = asg_stats$in_cell_per_mask,
    transcripts_cell_outside_window = asg_stats$cell_outside_window,
    transcripts_overlapping_nucleus = asg_stats$overlaps_nucleus,
    transcripts_qv20 = sum(tx$qv >= 20),
    cells_with_transcripts = asg_stats$distinct_cells,
    nucleus_rings = nuc_stats$rings,
    cells_with_nucleus = nuc_stats$cells_with_nucleus,
    expression_nnz = length(sub@x)),
  rebuilt = c(
    paste("the whole archive is written in 10x's LEGACY parquet/MEX layout",
          "from a zarr export, so every table is a new structure rather than",
          "a sliced copy of a source file"),
    sprintf(paste("cells.parquet, cell_boundaries.parquet and",
                  "nucleus_boundaries.parquet all come out of the one",
                  "%s-byte cells.zarr.zip (the cells array, polygon_sets/1",
                  "and polygon_sets/0); bytes_source is null for each",
                  "because no single source file corresponds to them"),
            format(ssize("cells.zarr.zip"), scientific = FALSE)),
    paste("transcripts.parquet: every detection in the window, qv unfiltered,",
          "coordinates translated; cell_id and overlaps_nucleus recovered",
          "from the vendor masks (see deviations)"),
    paste("cell_feature_matrix/: the VENDOR matrix restricted to the kept",
          "cells, written as an MEX triplet. Counts are the vendor's own, not",
          "re-aggregated from the transcripts, so the fixture can run the",
          "vendor-vs-aggregated cross-check the full tutorial does"),
    paste("experiment.xenium: crop-scoped fields recomputed from the kept",
          "cells, section-scoped fields that no longer describe anything",
          "dropped. Nothing is invented; see files[] and deviations"),
    paste("morphology_focus/: cropped to the window and re-encoded from",
          "JPEG-2000 to deflate so GDAL and terra can read it without a",
          "codec; no pyramid levels written")),
  deviations = c(
    sprintf(paste("cell_id and overlaps_nucleus are RECOVERED from the vendor",
                  "segmentation masks (cells.zarr masks/1 and masks/0),",
                  "because the transcripts zarr stores no per-transcript",
                  "assignment. Measured against the vendor matrix on this",
                  "window (%s): ratio %.5f, %.1f%% of %s cells exact, median",
                  "per-cell ratio %.5f, r %.6f. The lookup is 2D while the",
                  "vendor's is z-aware, so a small residual remains."),
            val$definition, val$overall_ratio, 100 * val$frac_cells_exact,
            format(val$cells_compared, big.mark = ","),
            val$median_cell_ratio, val$pearson_r),
    sprintf(paste("transcripts whose mask label is a cell OUTSIDE the window",
                  "are left UNASSIGNED rather than pointing at a cell the",
                  "fixture does not ship. That is %s detections, %.2f%%."),
            format(asg_stats$cell_outside_window, big.mark = ","),
            100 * asg_stats$cell_outside_window / nrow(tx)),
    paste("nucleus_distance is NA for every transcript: it is absent from the",
          "zarr and 10x does not document the convention precisely enough to",
          "reproduce it"),
    paste("written in the legacy parquet and MEX layout, so the fixture loads",
          "on a released stack and does NOT exercise the zarr reader the full",
          "tutorial uses"),
    if (have_morph) {
      paste("morphology_focus carries DAPI only. It is therefore NOT",
            "byte-identical to the vendor image, the fixture does not",
            "exercise the JPEG-2000 read path, and no pyramid levels are",
            "written. The other three channels are omitted on size")
    } else {
      "no morphology images"
    },
    sprintf(paste("nucleus_boundaries.parquet IS shipped, but %s cells own",
                  "more than one nucleus ring, so Giotto names the nucleus",
                  "polygons nucleus_1, nucleus_2 rather than by cell_id. The",
                  "nucleus layer therefore does not join to cell_id and",
                  "subsetGiotto() drops it with a warning. The published",
                  "xenium_mini_lung_ffpe fixture behaves identically; this is",
                  "upstream behaviour, not a property of this crop"),
            format(nuc_stats$rings - nuc_stats$cells_with_nucleus,
                   big.mark = ",")),
    sprintf(paste("coordinates ARE translated to the origin: translated by",
                  "(%g, %g) um, i.e. the window origin is subtracted, so the",
                  "fixture spans [0, %g] um in both",
                  "axes and lines up with the cropped morphology image, which",
                  "every reader anchors at (0,0). Add the shift back to",
                  "compare against the source run"),
            -XMIN, -YMIN, WINDOW_UM)),
  # crop_manifest.json itself is not listed: it describes the archive rather
  # than coming from the source run, and the published minis omit it too.
  files = c(
    list(
      fentry("cells.parquet", "cropped", NA,
             "centroid, area and nucleus counts of the kept cells, translated; from cells.zarr"),
      fentry("cell_boundaries.parquet", "cropped", NA,
             "segmentation ring vertices of the kept cells, translated; from cells.zarr polygon_sets/1"),
      fentry("nucleus_boundaries.parquet", "cropped", NA,
             "nucleus ring vertices, translated; from cells.zarr polygon_sets/0, which the legacy conversion omits"),
      fentry("transcripts.parquet", "cropped", ssize("transcripts.zarr.zip"),
             "detections inside the window, qv unfiltered, cell_id and overlaps_nucleus recovered from the masks"),
      fentry("cell_feature_matrix", "cropped", ssize("cell_feature_matrix.zarr.zip"),
             "vendor matrix restricted to the kept cells, as an MEX triplet"),
      fentry("experiment.xenium", "rewritten", ssize("experiment.xenium"),
             "crop-scoped fields recomputed, section-scoped fields dropped"),
      fentry("gene_panel.json", "verbatim", ssize("gene_panel.json"))),
    if (have_morph) list(
      fentry(file.path("morphology_focus", MORPH_CH), "cropped",
             ssize("morphology_focus", MORPH_CH),
             "DAPI only, cropped to the window and re-encoded from JPEG-2000 to deflate")),
    lapply(morph_other, function(ch) list(
      file = file.path("morphology_focus", ch), action = "dropped",
      reason = "non-DAPI channel; all four would put the archive over the 100 MB ceiling",
      bytes_source = ssize("morphology_focus", ch), bytes_fixture = 0L)),
    list(
      list(file = "analysis.zarr.zip", action = "dropped",
           reason = "run-level clustering, PCA and UMAP that would be stale after cropping",
           bytes_source = ssize("analysis.zarr.zip"), bytes_fixture = 0L),
      list(file = "analysis_summary.html", action = "dropped",
           reason = "run-level report that would be stale after cropping",
           bytes_source = ssize("analysis_summary.html"), bytes_fixture = 0L))),
  # which scripts produced this, and exactly which version of each. The repo
  # is where they are PUBLISHED, beside the other crops -- the sha256s pin the
  # content whether or not a commit is named.
  script = list(
    repo = "https://github.com/drieslab/spatial-datasets",
    files = lapply(c("make_xenium_mini.R", "recover_assignment.py",
                     "extract_nucleus.py", "crop_morphology.py",
                     "patch_experiment.py", "package_mini.sh"),
                   function(f) list(
                     path = file.path("crops", basename(OUTDIR), f),
                     sha256 = paste0(openssl::sha256(
                       file(file.path(HERE, "mini", f))))))))

write(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE, digits = NA,
             na = "null"),
      file.path(OUTDIR, "crop_manifest.json"))

# --- 7. report and enforce the size cap ------------------------------------
files <- list.files(OUTDIR, recursive = TRUE, full.names = TRUE)
rep <- data.table(file = sub(paste0(OUTDIR, "/"), "", files, fixed = TRUE),
                  mb = round(file.size(files) / 1024^2, 2))
setorder(rep, -mb)
cat("\n=== bundle ===\n"); print(rep)
cat(sprintf("total: %.2f MB   largest file: %.2f MB (cap %d MB)\n",
            sum(rep$mb), max(rep$mb), SIZE_CAP_MB))
if (max(rep$mb) > SIZE_CAP_MB) {
  cat(sprintf("\nWARNING: largest file exceeds the %d MB cap.\n", SIZE_CAP_MB))
  cat(sprintf("Re-run with a smaller window, e.g. %.0f um.\n",
              WINDOW_UM * sqrt(SIZE_CAP_MB / max(rep$mb)) * 0.95))
}
saveRDS(list(window_um = WINDOW_UM, xmin = XMIN, xmax = XMAX,
             ymin = YMIN, ymax = YMAX, n_cells = length(keep_ids),
             n_tx = nrow(tx), report = rep),
        file.path(HERE, "mini", "mini_build_info.rds"))
cat("\nMINI BUILD OK ->", OUTDIR, "\n")
