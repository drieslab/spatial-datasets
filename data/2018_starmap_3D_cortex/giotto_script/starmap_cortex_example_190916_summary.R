library(Giotto)

STARMAP_results_folder = '/path/to/results/folder/'


## 1. PREPARE DATA ####
STARMAP_data_folder = '/path/to/data/folder/'
expr = read.table(paste0(STARMAP_data_folder, '/', 'STARmap_3D_data_expression.txt'))
cell_loc = read.table(paste0(STARMAP_data_folder, '/', 'STARmap_3D_data_cell_locations.txt'))

## 2. Create & Process Giotto ####
gobject_folder = paste0(STARMAP_results_folder,'/','2_Gobject/')
if(!file.exists(gobject_folder)) dir.create(gobject_folder, recursive = T)

STAR_test <- createGiottoObject(raw_exprs = expr, spatial_locs = cell_loc)
filterDistributions(STAR_test, detection = 'genes')
filterDistributions(STAR_test, detection = 'cells')
filterCombinations(STAR_test, expression_thresholds = c(1, 1,2), gene_det_in_min_cells = c(20000, 20000, 30000), min_det_genes_per_cell = c(10, 20, 25))
STAR_test <- filterGiotto(gobject = STAR_test,
                          gene_det_in_min_cells = 20000,
                          min_det_genes_per_cell = 20)

## normalize & adjust
STAR_test <- normalizeGiotto(gobject = STAR_test, scalefactor = 10000, verbose = T)
STAR_test <- addStatistics(gobject = STAR_test)
STAR_test <- adjustGiottoMatrix(gobject = STAR_test, expression_values = c('normalized'),
                                batch_columns = NULL, covariate_columns = c('nr_genes', 'total_expr'),
                                return_gobject = TRUE,
                                update_slot = c('custom'))

real_3D = visPlot(gobject = STAR_test,
                  sdimx = "sdimx", sdimy = "sdimy", sdimz = "sdimz",
                  point_size = 1, axis_scale = "real", z_ticks = 2)
htmlwidgets::saveWidget(plotly::as_widget(real_3D), file = paste0(gobject_folder,'/', 'starmap_real_3D.html'))

cube_3D = visPlot(gobject = STAR_test,
                  sdimx = "sdimx", sdimy = "sdimy", sdimz = "sdimz",
                  point_size = 1, axis_scale = "cube", z_ticks = 2)
htmlwidgets::saveWidget(plotly::as_widget(cube_3D), file = paste0(gobject_folder,'/', 'starmap_cube_3D.html'))


## 3. Dimension reduction ####
STAR_test <- calculateHVG(gobject = STAR_test, method = 'cov_groups', zscore_threshold = 0.5, nr_expression_groups = 3)
STAR_test <- runPCA(gobject = STAR_test, genes_to_use = NULL, scale_unit = F)
signPCA(STAR_test)
STAR_test <- runUMAP(STAR_test, dimensions_to_use = 1:8, n_components = 3, n_threads = 4)

## 4. Cluster ####
cluster_folder = paste0(STARMAP_results_folder,'/','4_Cluster/')
if(!file.exists(cluster_folder)) dir.create(cluster_folder, recursive = T)

## sNN network (default)
STAR_test <- createNearestNetwork(gobject = STAR_test, dimensions_to_use = 1:8, k = 15)
## Leiden clustering
STAR_test <- doLeidenCluster(gobject = STAR_test, resolution = 0.2, n_iterations = 1000,
                             name = 'leiden_0.2',
                             python_path = "/Users/rubendries/Bin/anaconda3/envs/py36/bin/python")
STAR_UMAP <- plotUMAP(gobject = STAR_test, cell_color = 'leiden_0.2', 
                      point_size = 1.5,
                      plot_method = "plotly",
                      show_NN_network = T, 
                      edge_alpha = 0.05,
                      dim1_to_use = 1,
                      dim2_to_use = 2,
                      dim3_to_use = 3)
htmlwidgets::saveWidget(plotly::as_widget(STAR_UMAP), file = paste0(cluster_folder,'/', 'cluster_UMAP.html'))



### 5. co-visualize ####
covis_folder = paste0(STARMAP_results_folder,'/','5_Covisuals/')
if(!file.exists(covis_folder)) dir.create(covis_folder, recursive = T)

coPlot = visSpatDimPlot(gobject = STAR_test,
                        cell_color = 'leiden_0.2',
                        dim3_to_use = 3,
                        sdimz = "sdimz",
                        axis_scale = "real",
                        z_ticks = 2,
                        dim_point_size = 1,
                        spatial_point_size = 1,
                        show_NN_network = F)
htmlwidgets::saveWidget(plotly::as_widget(coPlot), file = paste0(covis_folder,'/', 'coPlot.html'))

showClusterHeatmap(gobject = STAR_test, cluster_column = 'leiden_0.2')





## 6. differential expression ####
DEG_folder = paste0(STARMAP_results_folder,'/','6_DEG/')
if(!file.exists(DEG_folder)) dir.create(DEG_folder, recursive = T)

markers = findMarkers_one_vs_all(gobject = STAR_test,
                                 method = 'gini',
                                 expression_values = 'normalized',
                                 cluster_column = 'leiden_0.2',
                                 min_genes = 5, rank_score = 2)
markers[, head(.SD, 4), by = 'cluster']

# violinplot
png(file = paste0(DEG_folder,'/','DEG_violinplot.png'), height = 920, width = 480)
violinPlot(STAR_test, genes = unique(markers$genes), cluster_column = 'leiden_0.2')
dev.off()

# genes heatmap
png(file = paste0(DEG_folder,'/','DEG_heatmap_cells.png'), height = 720, width = 720)
plotHeatmap(STAR_test, genes = STAR_test@gene_ID, cluster_column = 'leiden_0.2',
            legend_nrows = 2, expression_values = 'scaled',
            cluster_order = 'correlation', gene_order = 'correlation',
            show_plot = F)
dev.off()

# cluster heatmap
png(file = paste0(DEG_folder,'/','DEG_heatmap_clusters.png'), height = 720, width = 720)
plotMetaDataHeatmap(STAR_test, expression_values = 'scaled',
                    metadata_cols = c('leiden_0.2'))
dev.off()





## 7. cell type annotation ####
# --------------------------- #
annotation_folder = paste0(STARMAP_results_folder,'/','7_annotation/')
if(!file.exists(annotation_folder)) dir.create(annotation_folder, recursive = T)

## general cell types
clusters_cell_types_cortex = c('excit','excit','excit', 'inh', 'excit',
                               'other', 'other', 'other', 'inh', 'inh')
names(clusters_cell_types_cortex) = c(1:10)
STAR_test = annotateGiotto(gobject = STAR_test, annotation_vector = clusters_cell_types_cortex,
                           cluster_column = 'leiden_0.2', name = 'general_cell_types')

png(file = paste0(annotation_folder,'/', 'umap_general_cell_type.png'))
plotUMAP(STAR_test, plot_method = 'ggplot', cell_color = 'general_cell_types', point_size = 1.5)
dev.off()

png(file = paste0(annotation_folder,'/', 'cluster_heatmap_general_cell_type.png'))
plotMetaDataHeatmap(STAR_test, expression_values = 'scaled',
                    metadata_cols = c('general_cell_types'))
dev.off()



## detailed cell types
clusters_cell_types_cortex = c('L5','L4','L2/3', 'PV', 'L6',
                               'Astro', 'Olig1', 'Olig2', 'Calretinin', 'SST')
names(clusters_cell_types_cortex) = c(1:10)
STAR_test = annotateGiotto(gobject = STAR_test, annotation_vector = clusters_cell_types_cortex,
                           cluster_column = 'leiden_0.2', name = 'cell_types')


png(file = paste0(annotation_folder,'/', 'umap_cell_type.png'))
plotUMAP(STAR_test, plot_method = 'ggplot', cell_color = 'cell_types', point_size = 1.5)
dev.off()

png(file = paste0(annotation_folder,'/', 'cluster_heatmap_cell_type.png'))
plotMetaDataHeatmap(STAR_test, expression_values = 'scaled',
                    metadata_cols = c('cell_types'))
dev.off()


# create consistent color code
mynames = unique(pDataDT(STAR_test)$cell_types)
mycolorcode = Giotto:::getDistinctColors(n = 10)
names(mycolorcode) = mynames

realCellTypes =  visPlot(STAR_test, cell_color = 'cell_types', axis_scale = 'real',
                         sdimx = 'sdimx', sdimy = 'sdimy', sdimz = 'sdimz',
                         show_grid = F, cell_color_code = mycolorcode)
htmlwidgets::saveWidget(plotly::as_widget(realCellTypes), file = paste0(annotation_folder,'/', 'realCellTypes.html'))

## subsets
excit = visPlot(STAR_test, cell_color = 'cell_types', plot_method = 'plotly',
                sdimx = 'sdimx', sdimy = 'sdimy', sdimz = 'sdimz', axis_scale = 'real',
                select_cell_groups = c('L6','L5','L4','L2/3'),
                show_grid = F, cell_color_code = mycolorcode)
htmlwidgets::saveWidget(plotly::as_widget(excit), file = paste0(annotation_folder,'/', 'realCellTypes_excit.html'))

inhib = visPlot(STAR_test, cell_color = 'cell_types', plot_method = 'plotly',
                sdimx = 'sdimx', sdimy = 'sdimy', sdimz = 'sdimz', axis_scale = 'real',
                select_cell_groups = c('PV','Calretinin', 'SST'),
                show_grid = F, cell_color_code = mycolorcode)
htmlwidgets::saveWidget(plotly::as_widget(inhib), file = paste0(annotation_folder,'/', 'realCellTypes_inh.html'))

other = visPlot(STAR_test, cell_color = 'cell_types', plot_method = 'plotly',
                sdimx = 'sdimx', sdimy = 'sdimy', sdimz = 'sdimz', axis_scale = 'real',
                select_cell_groups = c('Astro', 'Olig1', 'Olig2'),
                show_grid = F, cell_color_code = mycolorcode)
htmlwidgets::saveWidget(plotly::as_widget(other), file = paste0(annotation_folder,'/', 'realCellTypes_other.html'))


## 8. spatial grid ####
# ------------------- #
grid_folder = paste0(STARMAP_results_folder,'/','8_grid/')
if(!file.exists(grid_folder)) dir.create(grid_folder, recursive = T)

STAR_test <- createSpatialGrid(gobject = STAR_test,
                               sdimx_stepsize = 100,
                               sdimy_stepsize = 100,
                               sdimz_stepsize = 20,
                               minimum_padding = 0)

mycolorcode = c('red', 'blue')
names(mycolorcode) = c("L2/3", "L6")

png(paste0(grid_folder,'/','gridplot.png'))
visPlot(STAR_test, cell_color = 'cell_types', sdimx = 'sdimx', sdimy = 'sdimy',
        show_grid = T, grid_color = 'green', spatial_grid_name = 'spatial_grid',
        point_size = 1.5, plot_method = 'ggplot',
        select_cell_groups = c("L2/3", "L6"), other_cells_alpha = 0.5, cell_color_code = mycolorcode)
dev.off()

#### spatial patterns ####
pattern_VC = detectSpatialPatterns(gobject = STAR_test, 
                                   expression_values = 'normalized',
                                   spatial_grid_name = 'spatial_grid',
                                   min_cells_per_grid = 5, 
                                   scale_unit = T, 
                                   PC_zscore = 1, 
                                   show_plot = T)

dim3_pattern = showPattern(pattern_VC,  plot_dim = 3, point_size = 4)
htmlwidgets::saveWidget(plotly::as_widget(dim3_pattern), file = paste0(grid_folder,'/', 'dim3_pattern.html'))

png(paste0(grid_folder,'/','patterngenes.png'))
showPatternGenes(pattern_VC, dimension = 1)
dev.off()

## 9. spatial network ####
## ---------------------##
spatnet_folder = paste0(STARMAP_results_folder,'/','9_spatial_network/')
if(!file.exists(spatnet_folder)) dir.create(spatnet_folder, recursive = T)

STAR_test <- createSpatialNetwork(gobject = STAR_test, k = 3)

networkplot = visPlot(gobject = STAR_test, show_network = T,
                      sdimx = "sdimx",sdimy = "sdimy",sdimz = "sdimz",
                      network_color = 'blue', spatial_network_name = 'spatial_network',axis_scale = "real",z_ticks = 2,
                      point_size = 4, cell_color = 'cell_types')
htmlwidgets::saveWidget(plotly::as_widget(networkplot), file = paste0(spatnet_folder,'/', 'networkplot.html'))



## 10. spatial genes ####
## ---------------------##
spatgenes_folder = paste0(STARMAP_results_folder,'/','10_spatial_genes/')
if(!file.exists(spatgenes_folder)) dir.create(spatgenes_folder, recursive = T)

kmtest = binGetSpatialGenes(STAR_test, bin_method = 'kmeans',
                            do_fisher_test = F, community_expectation = 5,
                            spatial_network_name = 'spatial_network', verbose = T)

ranktest = binGetSpatialGenes(STAR_test, bin_method = 'rank',
                              do_fisher_test = F, community_expectation = 5,
                              spatial_network_name = 'spatial_network', verbose = T)

png(file = paste0(spatgenes_folder,'/','spatial_genes.png'), width = 820, height = 360)
visSpatDimGenePlot(STAR_test, plot_method = 'ggplot',
                   genes = c('Cux2', 'Pcp4', 'Gja1', 'Mbp'), plot_alignment = 'vertical', cow_n_col = 4,
                   genes_high_color = 'red', genes_mid_color = 'white', genes_low_color = 'darkblue', midpoint = 8)
dev.off()


## 11. spatial HMRF domains ####
## ---------------------------##

# not available


## 12. cell-cell preferential proximity ####
## ---------------------------------------##
cellproxim_folder = paste0(STARMAP_results_folder,'/','12_cell_proxim/')
if(!file.exists(cellproxim_folder)) dir.create(cellproxim_folder, recursive = T)

## calculate frequently seen proximities
cell_proximities = cellProximityEnrichment(gobject = STAR_test,
                                           cluster_column = 'cell_types',
                                           spatial_network_name = 'spatial_network',
                                           number_of_simulations = 400)

png(file = paste0(cellproxim_folder,'/','cell_proximity_barplot.png'))
cellProximityBarplot(CPscore = cell_proximities, min_orig_ints = 5, min_sim_ints = 5)
dev.off()

png(file = paste0(cellproxim_folder,'/','cell_proximity_heatmap.png'))
cellProximityHeatmap(CPscore = cell_proximities, order_cell_types = T, scale = T)
dev.off()

STAR_astro_pv <- cellProximityVisPlot(gobject = STAR_test, interaction_name = "Astro-PV", spatial_network_name = 'spatial_network',
                                      axis_scale = 'real', 
                                      cluster_column = 'cell_types',
                                      sdimx = "sdimx",sdimy = "sdimy",sdimz = "sdimz",
                                      show_other_cells = F,
                                      cell_color = 'cell_types', 
                                      show_network = T,
                                      network_color = 'blue', 
                                      point_size_select = 4)
htmlwidgets::saveWidget(plotly::as_widget(STAR_astro_pv), file = paste0(cellproxim_folder,'/', 'STAR_astro_pv.html'))






