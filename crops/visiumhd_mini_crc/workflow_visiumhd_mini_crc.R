#!/usr/bin/env Rscript
# =============================================================================
# workflow_visiumhd_mini_crc.R
#
# Every Visium HD command we support, run on the visiumhd_mini_crc fixture, in
# memory (Giotto) and on disk (GiottoDisk), with one report line per step.
# Where the two backends should agree, the agreement is checked as its own
# line. Nothing stops at the first failure: the report is the result.
#
#   A  reader loaders             (each @calls loader, both readers)
#   B  create_gobject options     (bins, filters, subsets, sources, tar input)
#   C  wrappers                   (createGiottoVisiumHDObjectBin / Cell)
#   D  one object, several units  (2 / 8 / 16 um + cells, parent units)
#   E  downstream analysis        (bin008 and cell: the website tutorial's
#                                  pipeline, then markers, enrichment, spatial
#                                  stats, subsets, save / load)
#   F  aggregation                (the tutorial's path: 2 um points into cells
#                                  and hexagons, then the pipeline on hexagons)
#
# Env: VHD_MINI_DIR, the fixture to run on (default: the local copy it was made in).
# Writes <fixture>_workflow_report.csv and <fixture>_workflow_plots.pdf.
# =============================================================================

options(giotto.has_conda = FALSE, giotto.no_python_warn = TRUE, giotto.use_conda = FALSE,
        giotto.warn_sequential = FALSE, tilework.warn_sequential = FALSE,
        giottodisk.dgc_max_rows = Inf, giottodisk.dgc_max_cols = Inf)
suppressMessages({
    library(data.table)
    library(Giotto)
    library(GiottoDisk)
})

BASE <- "/Users/rubendries/Documents/Datasets/Cropped_datasets/VisiumHD"
DIR  <- Sys.getenv("VHD_MINI_DIR", file.path(BASE, "visiumhd_mini_crc"))
BO   <- file.path(DIR, "binned_outputs")
SO   <- file.path(DIR, "segmented_outputs")
OUT  <- sub("/?$", "", DIR)
WORK <- tempfile("vhd_workflow_")
dir.create(WORK)
PROJ <- function(tag) file.path(WORK, paste0("proj_", tag))

# ------------------------------------------------------------- reporting ----
REPORT <- list()
.log <- function(section, name, backend, status, secs = NA_real_, detail = "") {
    REPORT[[length(REPORT) + 1L]] <<- data.table(section, step = name, backend,
        status, secs = round(secs, 2), detail = substr(detail, 1, 400))
    cat(sprintf("%-2s %-5s %-7s %-46s %6s  %s\n", section, status, backend,
                name, if (is.na(secs)) "" else sprintf("%.1fs", secs),
                substr(detail, 1, 110)))
}
# run `expr`; log ok / ERROR; return its value (NULL on error)
step <- function(section, name, backend, expr, detail = NULL) {
    t0 <- Sys.time()
    v <- tryCatch(suppressWarnings(suppressMessages(force(expr))),
                  error = identity)
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (inherits(v, "error")) {
        .log(section, name, backend, "ERROR", secs, conditionMessage(v))
        if (nzchar(Sys.getenv("VHD_WORKFLOW_DEBUG"))) print(v)
        return(NULL)
    }
    d <- if (is.function(detail)) tryCatch(detail(v), error = function(e) "") else ""
    .log(section, name, backend, "ok", secs, d)
    v
}
# a comparison between the backends (or against ground truth)
agree <- function(section, name, cond, detail = "") {
    ok <- tryCatch(isTRUE(cond), error = function(e) FALSE)
    .log(section, name, "both", if (ok) "AGREE" else "DIFFER", NA, detail)
}
skip_step <- function(section, name, backend, why) .log(section, name, backend, "skip", NA, why)

BACKENDS <- c("memory", "disk")
be <- function(b, tag) if (b == "disk") PROJ(tag) else NULL
mat <- function(g, su, values = "raw") {
    x <- GiottoClass::getExpression(g, spat_unit = su, values = values, output = "exprObj")[]
    if (inherits(x, "parquetExprBase")) x <- storeRead(x, output = "dgcmatrix")
    methods::as(x, "CsparseMatrix")
}
same_mat <- function(a, b) {
    identical(dim(a), dim(b)) && setequal(colnames(a), colnames(b)) &&
        setequal(rownames(a), rownames(b)) &&
        max(abs(a[rownames(b), colnames(b)] - b)) < 1e-8
}
nbins <- function(g, su) length(spatIDs(g, spat_unit = su))
# features of one unit's expression (featIDs(g, x) would read x as a feat_type)
feats <- function(g, su) featIDs(GiottoClass::getExpression(g, spat_unit = su, output = "exprObj"))
has_python <- function() isTRUE(tryCatch(reticulate::py_available(initialize = FALSE) &&
    reticulate::py_module_available("community"), error = function(e) FALSE))
truth8 <- as.data.table(arrow::read_parquet(file.path(BO, "square_008um", "spatial",
                                                      "tissue_positions.parquet")))

cat(sprintf("fixture: %s\nwork dir: %s\n\n", DIR, WORK))

# ================================================================ A ==========
cat("== A. reader loaders ==\n")
rd <- list(memory = importVisiumHD(BO, bin = 8),
           disk = importVisiumHDDisk(BO, backend = PROJ("A"), bin = 8))
rs <- list(memory = importVisiumHD(SO),
           disk = importVisiumHDDisk(SO, backend = PROJ("A")))
A <- list()
for (b in BACKENDS) {
    r <- rd[[b]]
    A[[b]]$sl <- step("A", "load_tissue_position", b, r$load_tissue_position(verbose = FALSE),
                      function(v) sprintf("%d bins", nrow(v[])))
    step("A", "load_scalefactors", b, r$load_scalefactors(verbose = FALSE))
    for (it in c("hires", "lowres")) {
        step("A", paste0("load_image(", it, ")"), b, r$load_image(image_type = it, verbose = FALSE))
    }
    A[[b]]$ex <- step("A", "load_expression(raw)", b, r$load_expression(verbose = FALSE),
                      function(v) sprintf("%d feat types", length(v)))
    A[[b]]$exf <- step("A", "load_expression(filtered)", b,
                       r$load_expression(expression_source = "filtered", verbose = FALSE),
                       function(v) sprintf("%d bins", length(spatIDs(v[[1]]))))
    A[[b]]$exe <- step("A", "load_expression(ensembl)", b,
                       r$load_expression(feature_id_type = "ensembl", verbose = FALSE))
    for (p in c("raw_feature_bc_matrix.h5", "raw_feature_bc_matrix")) {
        step("A", paste0("load_expression(path = ", p, ")"), b,
             r$load_expression(path = file.path(BO, "square_008um", p), verbose = FALSE))
    }
    step("A", "load_cellmeta", b, r$load_cellmeta(verbose = FALSE))
    A[[b]]$tx <- step("A", "load_transcripts (2 um bin points)", b,
                      r$load_transcripts(verbose = FALSE),
                      function(v) class(v[[1]])[1])
    step("A", "tessellate_polygon(hexagon 50)", b,
         r$tessellate_polygon(shape = "hexagon", shape_size = 50, verbose = FALSE))
    if (b == "disk") {
        step("A", "load_bin_mapping", b, r$load_bin_mapping(verbose = FALSE),
             function(v) toString(names(v)))
    }
    s <- rs[[b]]
    step("A", "segmented load_expression", b, s$load_expression(verbose = FALSE))
    for (pt in c("cell", "nucleus")) {
        for (ga in c(FALSE, TRUE)) {
            step("A", sprintf("load_polygon(%s%s)", pt, if (ga) ", graphclust" else ""), b,
                 s$load_polygon(type = pt, graphclust_annotated = ga, verbose = FALSE),
                 function(v) sprintf("%d polys", length(spatIDs(v))))
        }
    }
}
ex_m <- function(x) { y <- x[[1]][]; if (inherits(y, "parquetExprBase")) y <- storeRead(y, output = "dgcmatrix"); y }
if (!is.null(A$memory$ex) && !is.null(A$disk$ex)) {
    agree("A", "load_expression(raw): same matrix", same_mat(ex_m(A$disk$ex), ex_m(A$memory$ex)))
}
if (!is.null(A$memory$exf) && !is.null(A$disk$exf)) {
    agree("A", "load_expression(filtered): same matrix",
          same_mat(ex_m(A$disk$exf), ex_m(A$memory$exf)))
}
if (!is.null(A$memory$exe) && !is.null(A$disk$exe)) {
    agree("A", "ensembl ids: same features, not symbols",
          setequal(featIDs(A$disk$exe[[1]]), featIDs(A$memory$exe[[1]])) &&
              mean(grepl("ENSG", featIDs(A$memory$exe[[1]]))) > 0.99,
          "genes without a symbol keep their Ensembl ID under both")
}

# ================================================================ B ==========
cat("\n== B. create_gobject options ==\n")
variants <- list(
    "bin = 2"                 = list(bin = 2),
    "bin = 8"                 = list(bin = 8),
    "bin = 16"                = list(bin = 16),
    "tissue_only = TRUE"      = list(bin = 8, tissue_only = TRUE),
    "micron = TRUE"           = list(bin = 8, micron = TRUE),
    "expression_source = filtered" = list(bin = 8, expression_source = "filtered"),
    "feature_id_type = ensembl"    = list(bin = 8, feature_id_type = "ensembl"),
    "barcodes = first 100"    = list(bin = 8, barcodes = truth8$barcode[1:100]),
    "array_subset_row/col"    = list(bin = 8,
        array_subset_row = range(truth8$array_row)[1] + c(0, 19),
        array_subset_col = range(truth8$array_col)[1] + c(0, 19)),
    # rows are given negated, as documented: y is inverted before filtering
    "pxl_subset_row/col"      = list(bin = 8,
        pxl_subset_row = -rev(quantile(truth8$pxl_row_in_fullres, c(0, 0.5))),
        pxl_subset_col = quantile(truth8$pxl_col_in_fullres, c(0, 0.5))),
    "image_type = lowres"     = list(bin = 8, image_type = "lowres"),
    "create_tessellated_polys" = list(bin = 8, create_tessellated_polys = TRUE,
        tessellate_shape_size = 50),
    "load_transcripts = TRUE" = list(bin = 8, load_transcripts = TRUE)
)
for (v in names(variants)) {
    a <- variants[[v]]
    su <- sprintf("bin%03d", as.integer(a$bin))
    g <- list()
    for (b in BACKENDS) {
        r <- if (b == "memory") importVisiumHD(BO, bin = a$bin)
             else importVisiumHDDisk(BO, backend = PROJ(paste0("B", gsub("\\W", "", v))), bin = a$bin)
        g[[b]] <- step("B", v, b, do.call(r$create_gobject, c(a, verbose = FALSE)),
                       function(x) sprintf("%d bins", nbins(x, su)))
    }
    if (!is.null(g$memory) && !is.null(g$disk)) {
        agree("B", paste0(v, ": same bins and raw matrix"),
              setequal(spatIDs(g$disk, su), spatIDs(g$memory, su)) &&
                  same_mat(mat(g$disk, su), mat(g$memory, su)),
              sprintf("%d bins", nbins(g$memory, su)))
    }
    if (v == "tissue_only = TRUE" && !is.null(g$memory)) {
        agree("B", "tissue_only keeps exactly the in-tissue bins",
              setequal(spatIDs(g$memory, su), truth8[in_tissue == 1, barcode]),
              sprintf("%d of %d", truth8[, sum(in_tissue)], nrow(truth8)))
    }
}
# tar input: the 8 um output only as 008um_outputs.tar
tar_bo <- file.path(WORK, "tar_outs", "binned_outputs")
dir.create(tar_bo, recursive = TRUE)
withr::with_dir(BO, utils::tar(file.path(tar_bo, "008um_outputs.tar"), "square_008um"))
for (b in BACKENDS) {
    r <- if (b == "memory") importVisiumHD(tar_bo, bin = 8)
         else importVisiumHDDisk(tar_bo, backend = PROJ("Btar"), bin = 8)
    step("B", "NNNum_outputs.tar input", b, r$create_gobject(load_image = FALSE, verbose = FALSE),
         function(x) sprintf("%d bins", nbins(x, "bin008")))
}

# ================================================================ C ==========
cat("\n== C. wrappers ==\n")
for (b in BACKENDS) {
    step("C", "createGiottoVisiumHDObjectBin", b,
         createGiottoVisiumHDObjectBin(BO, bin = 8, backend = be(b, "C1"), verbose = FALSE))
    step("C", "createGiottoVisiumHDObjectCell(+ 2 um points)", b,
         createGiottoVisiumHDObjectCell(SO, load_transcripts = TRUE, backend = be(b, "C2"),
                                        verbose = FALSE))
}
step("C", "createGiottoVisiumHDObject (deprecated)", "memory",
     createGiottoVisiumHDObject(visiumhd_dir = BO, bin = 8, verbose = FALSE))

# ================================================================ D ==========
cat("\n== D. one object, several units ==\n")
gD <- step("D", "create_gobject(bin = c(2, 8, 16), load_bin_mapping)", "disk",
           importVisiumHDDisk(BO, backend = PROJ("D"), bin = 8)$create_gobject(
               bin = c(2, 8, 16), load_bin_mapping = TRUE, load_transcripts = TRUE,
               verbose = FALSE),
           function(x) toString(names(x@expression)))
if (!is.null(gD)) {
    gD <- step("D", "add segmented cells (gobject =)", "disk",
               importVisiumHDDisk(SO, backend = PROJ("D"))$create_gobject(gobject = gD,
                                                                        verbose = FALSE),
               function(x) toString(names(x@expression)))
}
if (!is.null(gD)) {
    c2 <- pDataDT(gD, spat_unit = "bin002")
    agree("D", "bin002 records bin008, bin016, cell, in_cell, in_nucleus",
          all(c("bin008", "bin016", "cell", "in_cell", "in_nucleus") %in% names(c2)))
    agree("D", "parent ids exist in their units",
          all(c2$bin008 %in% spatIDs(gD, "bin008")) &&
              all(stats::na.omit(c2$cell) %in% spatIDs(gD, "cell")))
    # parent units as ordinary metadata: 2 um counts summed by `cell`
    m2 <- mat(gD, "bin002")
    grp <- c2$cell[match(colnames(m2), c2$cell_ID)]
    keep <- !is.na(grp)
    G <- Matrix::sparseMatrix(i = which(keep), j = match(grp[keep], unique(grp[keep])), x = 1,
                              dims = c(ncol(m2), length(unique(grp[keep]))),
                              dimnames = list(NULL, unique(grp[keep])))
    s2 <- m2 %*% G
    mc <- mat(gD, "cell")
    r <- intersect(rownames(s2), rownames(mc))
    agree("D", "2 um bins summed by parent cell == cell matrix",
          setequal(colnames(s2), colnames(mc)) &&
              max(abs(s2[r, colnames(mc)] - mc[r, ])) == 0)
}

# ================================================================ E ==========
# The website tutorial's pipeline (VisiumHD_Human_Colorectal_Cancer.Rmd), same
# calls and parameters. One deviation on disk, where z-score scaling is
# refused (it densifies the store): scale_feats / scale_cells are off, so
# feature plots show normalized rather than scaled values.
cat("\n== E. downstream analysis (tutorial pipeline) ==\n")
pdf(paste0(OUT, "_workflow_plots.pdf"), width = 8, height = 7)

tutorial_pipeline <- function(g, su, b, section) {
    tag <- sprintf("[%s]", su)
    run <- function(name, f, detail = NULL) {
        g2 <- step(section, paste(tag, name), b, f(g), detail)
        if (!is.null(g2) && inherits(g2, "giotto")) g <<- g2
        g2
    }
    show <- function(name, f) {
        pl <- step(section, paste(tag, name), b, f())
        if (!is.null(pl)) try(print(pl), silent = TRUE)
    }
    scaled <- b == "memory"
    run("filterGiotto", function(g) filterGiotto(g, spat_unit = su,
        expression_threshold = 1, feat_det_in_min_cells = 5,
        min_det_feats_per_cell = 25, verbose = FALSE),
        function(x) sprintf("%d x %d", length(feats(x, su)), nbins(x, su)))
    run("normalizeGiotto", function(g) normalizeGiotto(g, spat_unit = su,
        scalefactor = 1000, scale_feats = scaled, scale_cells = scaled,
        verbose = FALSE))
    run("addStatistics", function(g) addStatistics(g, spat_unit = su, verbose = FALSE))
    show("spatPlot2D(nr_feats)", function() spatPlot2D(g, spat_unit = su,
        cell_color = "nr_feats", color_as_factor = FALSE, show_image = TRUE,
        return_plot = TRUE, show_plot = FALSE, save_plot = FALSE))
    run("calculateHVF", function(g) calculateHVF(g, spat_unit = su, zscore_threshold = 1,
        show_plot = FALSE, return_plot = FALSE, save_plot = FALSE))
    run("runPCA", function(g) runPCA(g, spat_unit = su, expression_values = "normalized",
        feats_to_use = "hvf", verbose = FALSE))
    show("screePlot", function() screePlot(g, spat_unit = su, ncp = 30,
        return_plot = TRUE, show_plot = FALSE, save_plot = FALSE))
    run("runUMAP", function(g) runUMAP(g, spat_unit = su, dimensions_to_use = 1:14))
    run("createNearestNetwork (sNN)", function(g) createNearestNetwork(g, spat_unit = su,
        dimensions_to_use = 1:14, k = 5))
    run("doLeidenClusterIgraph", function(g) doLeidenClusterIgraph(g, spat_unit = su,
        resolution = 0.5, n_iterations = 1000),
        function(x) sprintf("%d clusters", length(unique(pDataDT(x, spat_unit = su)$leiden_clus))))
    show("plotUMAP", function() plotUMAP(g, spat_unit = su, cell_color = "leiden_clus",
        point_size = 1, return_plot = TRUE, show_plot = FALSE, save_plot = FALSE))
    show("spatPlot2D(leiden_clus)", function() spatPlot2D(g, spat_unit = su,
        cell_color = "leiden_clus", show_image = TRUE, return_plot = TRUE,
        show_plot = FALSE, save_plot = FALSE))
    run("createSpatialNetwork (kNN, k = 8)", function(g) createSpatialNetwork(g,
        spat_unit = su, name = "kNN_network", method = "kNN", k = 8))
    hvf_genes <- tryCatch(fDataDT(g, spat_unit = su)[hvf == "yes"]$feat_ID,
                          error = function(e) character())
    ranktest <- step(section, paste(tag, "binSpect (rank)"), b, binSpect(g, spat_unit = su,
        subset_feats = hvf_genes, bin_method = "rank", calc_hub = FALSE,
        do_fisher_test = TRUE, spatial_network_name = "kNN_network"),
        function(x) sprintf("%d feats", nrow(x)))
    if (!is.null(ranktest)) {
        top <- head(ranktest$feats, 4)
        show("spatFeatPlot2D(top spatial genes)", function() spatFeatPlot2D(g,
            spat_unit = su, expression_values = if (scaled) "scaled" else "normalized",
            feats = top, gradient_style = "sequential", cow_n_col = 2, point_size = 1,
            return_plot = TRUE, show_plot = FALSE, save_plot = FALSE))
        spat_genes <- ranktest[adj.p.value < 0.001]$feats
        if (length(spat_genes) < 20) spat_genes <- head(ranktest$feats, 50)
        cor_dt <- step(section, paste(tag, "detectSpatialCorFeats"), b,
            detectSpatialCorFeats(g, spat_unit = su, method = "network",
                spatial_network_name = "kNN_network", subset_feats = spat_genes))
        if (!is.null(cor_dt)) {
            k <- min(16L, max(2L, length(spat_genes) %/% 5L))
            cor_dt <- step(section, paste(tag, "clusterSpatialCorFeats"), b,
                clusterSpatialCorFeats(cor_dt, name = "spat_netw_clus", k = k))
        }
        if (!is.null(cor_dt)) {
            show("heatmSpatialCorFeats", function() heatmSpatialCorFeats(g,
                spatCorObject = cor_dt, use_clus_name = "spat_netw_clus",
                show_plot = FALSE, return_plot = TRUE, save_plot = FALSE))
            cg <- showSpatialCorFeats(cor_dt, use_clus_name = "spat_netw_clus",
                                      show_top_feats = 1)
            clus <- stats::setNames(cg$clus, cg$feat_ID)
            run("createMetafeats", function(g) createMetafeats(g, spat_unit = su,
                expression_values = "normalized", feat_clusters = clus,
                name = "cluster_metagene"))
            show("spatCellPlot(cluster_metagene)", function() spatCellPlot(g,
                spat_unit = su, spat_enr_names = "cluster_metagene",
                cell_annotation_values = as.character(head(sort(unique(clus)), 4)),
                gradient_style = "sequential", point_size = 1, cow_n_col = 2,
                return_plot = TRUE, show_plot = FALSE, save_plot = FALSE))
        }
    }
    g
}

# established Giotto functions beyond the tutorial's pipeline
extra_steps <- function(g, su, b, section) {
    tag <- sprintf("[%s]", su)
    mk <- step(section, paste(tag, "findMarkers_one_vs_all (scran)"), b,
        findMarkers_one_vs_all(g, spat_unit = su, cluster_column = "leiden_clus",
            method = "scran", expression_values = "normalized", min_feats = 5),
        function(x) sprintf("%d rows", nrow(x)))
    if (!is.null(mk)) {
        top <- mk[, head(feats, 10), by = cluster]
        sm <- makeSignMatrixPAGE(sign_names = paste0("c", unique(top$cluster)),
                                 sign_list = split(top$V1, top$cluster))
        g2 <- step(section, paste(tag, "runPAGEEnrich"), b, runPAGEEnrich(g,
            spat_unit = su, sign_matrix = sm, expression_values = "normalized",
            verbose = FALSE))
        if (!is.null(g2)) g <- g2
    }
    if (requireNamespace("spatialDWLS", quietly = TRUE)) {
        sm_dwls <- step(section, paste(tag, "makeSignMatrixDWLS"), b,
            makeSignMatrixDWLS(g, spat_unit = su, sign_gene = unique(top$V1),
                cell_type_vector = pDataDT(g, spat_unit = su)$leiden_clus))
        if (!is.null(sm_dwls)) step(section, paste(tag, "runDWLSDeconv"), b,
            runDWLSDeconv(g, spat_unit = su, sign_matrix = sm_dwls,
                cluster_column = "leiden_clus"))
    } else {
        skip_step(section, paste(tag, "runDWLSDeconv"), b, "spatialDWLS not installed")
    }
    step(section, paste(tag, "cellProximityEnrichment"), b, cellProximityEnrichment(g,
        spat_unit = su, cluster_column = "leiden_clus",
        spatial_network_name = "kNN_network", number_of_simulations = 50))
    step(section, paste(tag, "spatialAutoCorGlobal"), b, spatialAutoCorGlobal(g,
        spat_unit = su, feats = head(feats(g, su), 20),
        spatial_network_to_use = "kNN_network"))
    if (has_python()) {
        step(section, paste(tag, "doLouvainCluster"), b, doLouvainCluster(g,
            spat_unit = su, resolution = 0.5))
    } else {
        skip_step(section, paste(tag, "doLouvainCluster"), b,
                  "needs python + community (not set up here)")
    }
    pl <- step(section, paste(tag, "violinPlot"), b, violinPlot(g, spat_unit = su,
        feats = head(feats(g, su), 3), cluster_column = "leiden_clus",
        return_plot = TRUE, show_plot = FALSE, save_plot = FALSE))
    if (!is.null(pl)) try(print(pl), silent = TRUE)
    sl <- getSpatialLocations(g, spat_unit = su, output = "data.table")
    step(section, paste(tag, "subsetGiotto(half the cells)"), b, subsetGiotto(g,
        spat_unit = su, cell_ids = head(spatIDs(g, su), nbins(g, su) %/% 2)),
        function(x) sprintf("%d", nbins(x, su)))
    step(section, paste(tag, "subsetGiottoLocs(left half)"), b, subsetGiottoLocs(g,
        spat_unit = su, x_max = stats::median(sl$sdimx)),
        function(x) sprintf("%d", nbins(x, su)))
    g
}

save_load <- function(g, su, b, section, proj) {
    tag <- sprintf("[%s]", su)
    if (b == "memory") {
        f <- file.path(WORK, paste0("save_", su))
        step(section, paste(tag, "saveGiotto + loadGiotto"), b, {
            saveGiotto(g, dir = f, overwrite = TRUE, verbose = FALSE)
            loadGiotto(file.path(f, "saveGiottoDir"), verbose = FALSE)
        }, function(x) sprintf("%d", nbins(x, su)))
    } else {
        step(section, paste(tag, "saveGiotto(name) + loadGiotto"), b, {
            saveGiotto(g, name = "workflow")
            loadGiotto(proj, name = "workflow")
        }, function(x) sprintf("%d", nbins(x, su)))
    }
}

compare_backends <- function(gm, gd, su, section) {
    tag <- sprintf("[%s]", su)
    agree(section, paste(tag, "same cells and features after filter"),
          setequal(spatIDs(gd, su), spatIDs(gm, su)) && setequal(feats(gd, su), feats(gm, su)),
          sprintf("%d x %d", length(feats(gm, su)), nbins(gm, su)))
    agree(section, paste(tag, "same normalized values"),
          same_mat(mat(gd, su, "normalized"), mat(gm, su, "normalized")))
    hv <- function(g) { f <- fDataDT(g, spat_unit = su); sort(f$feat_ID[f$hvf == "yes"]) }
    agree(section, paste(tag, "same HVFs"), identical(hv(gd), hv(gm)), sprintf("%d", length(hv(gm))))
    pc <- function(g) { d <- getDimReduction(g, spat_unit = su, output = "matrix"); d[sort(rownames(d)), 1:3] }
    cc <- abs(diag(stats::cor(pc(gd), pc(gm))))
    agree(section, paste(tag, "PC1-3 match up to sign (|r| > 0.99)"), all(cc > 0.99),
          paste(round(cc, 4), collapse = " "))
    k <- function(g) length(unique(pDataDT(g, spat_unit = su)$leiden_clus))
    agree(section, paste(tag, "Leiden: >= 3 clusters on both"), k(gm) >= 3 && k(gd) >= 3,
          sprintf("memory %d, disk %d", k(gm), k(gd)))
}

build <- list(
    bin008 = function(b) createGiottoVisiumHDObjectBin(BO, bin = 8, tissue_only = TRUE,
        backend = be(b, "E8"), verbose = FALSE),
    cell = function(b) createGiottoVisiumHDObjectCell(SO, backend = be(b, "Ec"), verbose = FALSE)
)
E <- list()
for (su in names(build)) for (b in BACKENDS) {
    g <- step("E", sprintf("[%s] build", su), b, build[[su]](b))
    if (is.null(g)) next
    g <- tutorial_pipeline(g, su, b, "E")
    g <- extra_steps(g, su, b, "E")
    save_load(g, su, b, "E", PROJ(if (su == "cell") "Ec" else "E8"))
    E[[su]][[b]] <- g
}
for (su in names(E)) {
    if (!is.null(E[[su]]$memory) && !is.null(E[[su]]$disk)) {
        compare_backends(E[[su]]$memory, E[[su]]$disk, su, "E")
    }
}

# ================================================================ F ==========
# The tutorial's own aggregation: 2 um bin points into polygons with
# calculateOverlap() / overlapToMatrix() on the object, then the same pipeline
# on the new unit. Hexagons are tessellated over the points' extent, as in the
# tutorial.
cat("\n== F. aggregation (tutorial path) ==\n")
F <- list(); EXT <- list()
for (b in BACKENDS) {
    g <- step("F", "cells + 2 um points", b, createGiottoVisiumHDObjectCell(SO,
        load_transcripts = TRUE, backend = be(b, "F"), verbose = FALSE))
    if (is.null(g)) next
    gp <- getFeatureInfo(g, return_giottoPoints = TRUE)
    EXT[[b]] <- as.vector(terra::ext(getPolygonInfo(g, "cell", return_giottoPolygon = TRUE)[]))
    g2 <- step("F", "calculateOverlap(cell, rna)", b,
               calculateOverlap(g, spatial_info = "cell", feat_info = "rna"))
    if (!is.null(g2)) g <- g2
    g2 <- step("F", "overlapToMatrix(cell, rna, name = from_points)", b,
               overlapToMatrix(g, poly_info = "cell", feat_info = "rna", name = "from_points"))
    if (!is.null(g2)) {
        g <- g2
        F[[b]]$cell <- mat(g, "cell", "from_points")
    }
    hex <- step("F", "tessellate(hexagon 50 px over the points)", b, tessellate(
        extent = terra::ext(gp), shape = "hexagon", shape_size = 50, name = "hex50"),
        function(x) sprintf("%d hexagons", length(spatIDs(x))))
    if (is.null(hex)) next
    g2 <- step("F", "add hexagons + addSpatialCentroidLocations", b, {
        gh <- setGiotto(g, hex, verbose = FALSE)
        addSpatialCentroidLocations(gh, poly_info = "hex50")
    })
    if (is.null(g2)) next
    g <- g2
    g2 <- step("F", "calculateOverlap(hex50, rna)", b,
               calculateOverlap(g, spatial_info = "hex50", feat_info = "rna"))
    if (!is.null(g2)) g <- g2
    g2 <- step("F", "overlapToMatrix(hex50, rna, name = raw)", b,
               overlapToMatrix(g, poly_info = "hex50", feat_info = "rna", name = "raw"))
    if (is.null(g2)) next
    g <- g2
    F[[b]]$hex <- mat(g, "hex50")
    F[[b]]$g <- tutorial_pipeline(g, "hex50", b, "F")
}
dev.off()

shipped <- local({
    h <- hdf5r::H5File$new(file.path(SO, "raw_feature_cell_matrix.h5"), "r"); on.exit(h$close_all())
    Matrix::sparseMatrix(i = h[["matrix/indices"]][] + 1, p = as.numeric(h[["matrix/indptr"]][]),
        x = as.numeric(h[["matrix/data"]][]), dims = h[["matrix/shape"]][],
        dimnames = list(h[["matrix/features/name"]][], h[["matrix/barcodes"]][]))
})
for (b in names(F)) if (!is.null(F[[b]]$cell)) {
    agree("F", sprintf("[%s] points in cells: total == shipped cell matrix total", b),
          sum(F[[b]]$cell) == sum(shipped),
          sprintf("%d vs %d", sum(F[[b]]$cell), sum(shipped)))
}
if (!is.null(F$memory$cell) && !is.null(F$disk$cell)) {
    agree("F", "cells from points: disk == memory", same_mat(F$disk$cell, F$memory$cell))
}
if (!is.null(F$memory$hex) && !is.null(F$disk$hex)) {
    agree("F", "hexagon matrix: disk == memory", same_mat(F$disk$hex, F$memory$hex),
          sprintf("%d hexagons", ncol(F$memory$hex)))
}
if (!is.null(F$memory$g) && !is.null(F$disk$g)) {
    compare_backends(F$memory$g, F$disk$g, "hex50", "F")
}
if (length(EXT) == 2L) {
    # known, deferred: ext() of a disk polygon store spans the centroids
    agree("F", "ext(cell polygons): disk == memory [known issue]",
          isTRUE(all.equal(EXT$disk, EXT$memory)),
          sprintf("memory %s | disk %s", paste(round(EXT$memory), collapse = " "),
                  paste(round(EXT$disk), collapse = " ")))
}

# ============================================================== report ======
rep <- rbindlist(REPORT)
fwrite(rep, paste0(OUT, "_workflow_report.csv"))
cat("\n== summary ==\n")
print(rep[, .N, by = .(section, status)][order(section, status)])
bad <- rep[status %in% c("ERROR", "DIFFER")]
if (nrow(bad)) {
    cat("\nnot ok:\n")
    print(bad[, .(section, backend, step, status, detail)], nrows = 200)
}
cat(sprintf("\nreport: %s_workflow_report.csv\nplots:  %s_workflow_plots.pdf\n", OUT, OUT))
