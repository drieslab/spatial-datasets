
## osmFISH SS cortex dataset example
## works with Giotto version v0.1.2
library(Giotto)

## folder to save results to ##
osmFISH_results_folder = '/path/to/results/folder/'


## 1. PREPARE DATA ####
data_dir = '/Volumes/Ruben_Seagate/Dropbox/Projects/GC_lab/Ruben_Dries/190225_spatial_package/Data/osmFISH_data/'
## ss cortex expression DATA ##
osm_exprs = read.table(file = paste0(data_dir,'/','osmFISH_prep_expression.txt'))
## prepare cell locations
osm_locs = read.table(file = paste0(data_dir,'/','osmFISH_prep_cell_coordinates.txt'))
osm_locs = osm_locs[rownames(osm_locs) %in% colnames(osm_exprs),]



## 2. Create & Process Giotto ####
gobject_folder = paste0(osmFISH_results_folder,'/','2_Gobject/')
if(!file.exists(gobject_folder)) dir.create(gobject_folder, recursive = T)

osm_test <- createGiottoObject(raw_exprs = osm_exprs, spatial_locs = osm_locs)
metadata = fread(file = paste0(data_dir,'/','osmFISH_prep_cell_metadata.txt'))
osm_test = addCellMetadata(osm_test, new_metadata = metadata,
                           by_column = T, column_cell_ID = 'CellID')
osm_test <- filterGiotto(gobject = osm_test,
                         expression_threshold = 1,
                         gene_det_in_min_cells = 10,
                         min_det_genes_per_cell = 10,
                         expression_values = c('raw'),
                         verbose = T)

## normalize
# 1. standard z-score way
osm_test <- normalizeGiotto(gobject = osm_test)

# 2. osmFISH way
raw_expr_matrix = osm_test@raw_exprs
norm_genes = (raw_expr_matrix/rowSums(raw_expr_matrix)) * nrow(raw_expr_matrix)
norm_genes_cells = t((t(norm_genes)/colSums(norm_genes)) * ncol(raw_expr_matrix))
osm_test@custom_expr = norm_genes_cells

## add gene & cell statistics
osm_test <- addStatistics(gobject = osm_test)

## visualize original annotations ##
png(file = paste0(gobject_folder,'/','original_clusters.png'), height = 920, width = 920)
visPlot(gobject = osm_test, sdimx = 'sdimx', sdimy = 'sdimy', cell_color = 'ClusterName')
dev.off()

png(file = paste0(gobject_folder,'/','original_regions.png'), height = 920, width = 920)
visPlot(gobject = osm_test, sdimx = 'sdimx', sdimy = 'sdimy', cell_color = 'Region')
dev.off()



## 3. Dimension reduction ####
dimred_folder = paste0(osmFISH_results_folder,'/','3_DimRed/')
if(!file.exists(dimred_folder)) dir.create(dimred_folder, recursive = T)

## highly variable genes (HVG)
# only 33 genes so use all genes

## run PCA on expression values (default)
osm_test <- runPCA(gobject = osm_test, expression_values = 'custom', scale_unit = F)
signPCA(gobject = osm_test, expression_values = 'custom', scale_unit = F)

png(file = paste0(dimred_folder,'/','PCA_screeplot.png'), height = 480, width = 480)
signPCA(gobject = osm_test, expression_values = 'custom', scale_unit = F)
dev.off()

png(file = paste0(dimred_folder,'/','PCA_reduction.png'), height = 920, width = 920)
plotPCA(osm_test)
dev.off()

## run UMAP and tSNE on PCA space (default)
osm_test <- runUMAP(osm_test, dimensions_to_use = 1:31, expression_values = 'custom')
png(file = paste0(dimred_folder,'/','UMAP_reduction.png'), height = 920, width = 920)
plotUMAP(gobject = osm_test)
dev.off()

osm_test <- runtSNE(osm_test, dimensions_to_use = 1:31, perplexity = 70, check_duplicates = F)
png(file = paste0(dimred_folder,'/','tSNE_reduction.png'), height = 920, width = 920)
plotTSNE(gobject = osm_test)
dev.off()




## 4. Cluster ####
# -------------- #
cluster_folder = paste0(osmFISH_results_folder,'/','4_Cluster/')
if(!file.exists(cluster_folder)) dir.create(cluster_folder, recursive = T)

## hierarchical clustering
osm_test = doHclust(gobject = osm_test, expression_values = 'custom', k = 34)
png(file = paste0(cluster_folder,'/','UMAP_hclust.png'), height = 920, width = 920)
plotUMAP(gobject = osm_test, cell_color = 'hclust', point_size = 2.5,
         show_NN_network = F, edge_alpha = 0.05, plot_method = 'ggplot')
dev.off()

## kmeans clustering
osm_test = doKmeans(gobject = osm_test, expression_values = 'custom', centers = 32, nstart = 2000)
png(file = paste0(cluster_folder,'/','UMAP_kmeans.png'), height = 920, width = 920)
plotUMAP(gobject = osm_test, cell_color = 'kmeans',
         point_size = 2.5, show_NN_network = F, edge_alpha = 0.05, plot_method = 'ggplot')
dev.off()

## Leiden clustering
# sNN network (default)
osm_test <- createNearestNetwork(gobject = osm_test, dimensions_to_use = 1:31, k = 15)
osm_test <- doLeidenCluster(gobject = osm_test, resolution = 0.05, n_iterations = 1000,
                            python_path = "/Users/rubendries/Bin/anaconda3/envs/py36/bin/python")
png(file = paste0(cluster_folder,'/','UMAP_leiden.png'), height = 920, width = 920)
plotUMAP(gobject = osm_test, cell_color = 'leiden_clus', point_size = 2.5,
         show_NN_network = F, edge_alpha = 0.05, plot_method = 'ggplot')
dev.off()

# merge small groups based on similarity
leiden_similarities = getClusterSimilarity(osm_test,
                                           expression_values = 'custom',
                                           cluster_column = 'leiden_clus')
osm_test = mergeClusters(osm_test, expression_values = 'custom',
                         cluster_column = 'leiden_clus',
                         new_cluster_name = 'leiden_clus_m',
                         max_group_size = 30, force_min_group_size = 20,
                         return_gobject = T)
png(file = paste0(cluster_folder,'/','UMAP_leiden_merged.png'), height = 920, width = 920)
plotUMAP(gobject = osm_test, cell_color = 'leiden_clus_m', point_size = 2.5,
         show_NN_network = F, edge_alpha = 0.05, plot_method = 'ggplot')
dev.off()

## show cluster relationships
png(file = paste0(cluster_folder,'/','leiden_merged_heatmap.png'), height = 920, width = 920)
showClusterHeatmap(gobject = osm_test, expression_values = 'custom', cluster_column = 'leiden_clus_m')
dev.off()

png(file = paste0(cluster_folder,'/','leiden_merged_dendrogram.png'), height = 920, width = 920)
showClusterDendrogram(gobject = osm_test, expression_values = 'custom', cluster_column = 'leiden_clus_m')
dev.off()


## 5. co-visualize ####
# ------------------- #
# expression and spatial
covis_folder = paste0(osmFISH_results_folder,'/','5_Covisuals/')
if(!file.exists(covis_folder)) dir.create(covis_folder, recursive = T)

png(file = paste0(covis_folder,'/','covis_leiden_merged.png'), height = 920, width = 920)
visSpatDimPlot(gobject = osm_test, cell_color = 'leiden_clus_m', sdimx = 'sdimx', sdimy = 'sdimy',
               dim_point_size = 2, spatial_point_size = 2)
dev.off()

png(file = paste0(covis_folder,'/','covis_leiden_merged_selected.png'), height = 920, width = 920)
visSpatDimPlot(gobject = osm_test, cell_color = 'leiden_clus_m', sdimx = 'sdimx', sdimy = 'sdimy',
               dim_point_size = 2, spatial_point_size = 2, select_cell_groups = 'm_8')
dev.off()





## 6. differential expression ####
# ------------------------------ #
DEG_folder = paste0(osmFISH_results_folder,'/','6_DEG/')
if(!file.exists(DEG_folder)) dir.create(DEG_folder, recursive = T)

## split dendrogram nodes ##
dendsplits = getDendrogramSplits(gobject = osm_test,
                                 expression_values = 'custom',
                                 cluster_column = 'leiden_clus_m')
split_3_markers = findGiniMarkers(gobject = osm_test, expression_values = 'custom', cluster_column = 'leiden_clus_m',
                                  group_1 = unlist(dendsplits[3]$tree_1), group_2 = unlist(dendsplits[3]$tree_2))

## Individual populations ##
markers = findMarkers_one_vs_all(gobject = osm_test,
                                 method = 'scran',
                                 expression_values = 'custom',
                                 cluster_column = 'leiden_clus_m',
                                 min_genes = 2, rank_score = 2)

## violinplot
topgenes = markers[, head(.SD, 1), by = 'cluster_ID']$gene_ID
png(file = paste0(DEG_folder,'/','violinplot_leiden_merged.png'), height = 920, width = 480)
violinPlot(osm_test, genes = unique(topgenes), cluster_column = 'leiden_clus_m', expression_values = 'custom', strip_text = 5)
dev.off()

## cluster heatmap
ranked_genes = c('Bmp4', 'Itpr2', 'Tmem2', 'Ctps', 'Plp1',
                 'Sox10','Foxj1', 'Aldoc', 'Gfap', 'Acta2',
                 'Mrc1', 'Vtn', 'Crhbp', 'Slc32a1', 'Gad2',
                 'Syt6', 'Serpinf1', 'Cpne5', 'Lamp5', 'Hexb',
                 'Kcnip2', 'Tbr1', 'Ttr', 'Apln', 'Anln',
                 'Crh', 'Vip', 'Cnr1', 'Pthlh', 'Rorb',
                 'Flt1', 'Mfge8', 'Pdgfra')

png(file = paste0(DEG_folder,'/','cluster_heatmap_leiden_merged.png'), height = 480, width = 480)
plotMetaDataHeatmap(osm_test, expression_values = 'custom',
                    metadata_cols = c('leiden_clus_m'), custom_gene_order = ranked_genes)
dev.off()


## 7. cell type annotation ####
# --------------------------- #
annotation_folder = paste0(osmFISH_results_folder,'/','7_annotation/')
if(!file.exists(annotation_folder)) dir.create(annotation_folder, recursive = T)

## create vector with names
clusters_SS_cortex = c('OOP', 'OL1', 'OL2', 'OL3', 'OL4',
                       'Ependymal', 'unknown', 'Astro_Gfap', 'vSMC', 'Pericytes',
                       'IN1', 'IN2', 'Pyr1', 'Astro', 'IN3',
                       'IN4', 'Pyr2', 'Miglia1', 'IN5', 'Pyr3',
                       'Choroid', 'Vend1', 'OL5', 'IN6', 'IN7',
                       'IN8', 'IN9', 'Pyr4', 'Pyr5', 'Pyr6',
                       'Vend2', 'Astro_Mfge8', 'OPC')
names(clusters_SS_cortex) = c('m_1', '18', 'm_2', 'm_5', 'm_8',
                              'm_10', 'm_21', '9', 'm_17', 'm_19',
                              'm_11', 'm_14', 'm_6', '30', 'm_3',
                              'm_16', 'm_7', 'm_12', '11', '13',
                              'm_15', 'm_18', '27', 'm_20', '20',
                              '17', '31', '33', '22', 'm_4',
                              'm_13', '8', 'm_9')
osm_test = annotateGiotto(gobject = osm_test, annotation_vector = clusters_SS_cortex,
                          cluster_column = 'leiden_clus_m', name = 'leiden_clus_m_types')
png(file = paste0(annotation_folder,'/','annotation_leiden_merged_first.png'), height = 920, width = 920)
visSpatDimPlot(gobject = osm_test, cell_color = 'leiden_clus_m_types', sdimx = 'sdimx', sdimy = 'sdimy',
               dim_point_size = 2, spatial_point_size = 2)
dev.off()



## compare clusters with osmFISH paper
clusters_det_SS_cortex = c('Olig_COP', 'Olig_NF', 'Olig_MF', 'Olig_mat', 'Olig_mat',
                           'Ependymal', 'unknown', 'Astro_Gfap', 'vSMC', 'Pericytes',
                           'Inh_Crhbp', 'Inh_IC', 'Pyr_L6', 'Periv_Macro', 'Pyr_Cpne5',
                           'unknown', 'Pyr_L2/3', 'Microglia', 'Hippocampus', 'Pyr_L5',
                           'Choroid', 'vEnd', 'unknown', 'Inh_Anln', 'Inh_Crh',
                           'Inh_Vip', 'Inh_Pthlh', 'Pyr_Apln', 'Pyr_Kcnip2', 'Pyr_L4',
                           'vEnd', 'Astro_Mfge8', 'Olig_precursor')
names(clusters_det_SS_cortex) = c('m_1', '18', 'm_2', 'm_5', 'm_8',
                                  'm_10', 'm_21', '9', 'm_17', 'm_19',
                                  'm_11', 'm_14', 'm_6', '30', 'm_3',
                                  'm_16', 'm_7', 'm_12', '11', '13',
                                  'm_15', 'm_18', '27', 'm_20', '20',
                                  '17', '31', '33', '22', 'm_4',
                                  'm_13', '8', 'm_9')
osm_test = annotateGiotto(gobject = osm_test, annotation_vector = clusters_det_SS_cortex,
                          cluster_column = 'leiden_clus_m', name = 'det_cell_types')
png(file = paste0(annotation_folder,'/','annotation_leiden_merged_detailed.png'), height = 920, width = 920)
visSpatDimPlot(gobject = osm_test, cell_color = 'det_cell_types', sdimx = 'sdimx', sdimy = 'sdimy',
               dim_point_size = 2, spatial_point_size = 2)
dev.off()

## coarse cell types
clusters_coarse_SS_cortex = c('Olig', 'Olig', 'Olig', 'Olig', 'Olig',
                              'Ependymal', 'unknown', 'Astro', 'vSMC', 'Pericytes',
                              'Inh', 'Inh', 'Pyr', 'Periv_Macro', 'Pyr',
                              'unknown', 'Pyr', 'Microglia', 'Hippocampus', 'Pyr',
                              'Choroid', 'vEnd', 'unknown', 'Inh', 'Inh',
                              'Inh', 'Inh', 'Pyr', 'Pyr', 'Pyr',
                              'vEnd', 'Astro', 'Olig')
names(clusters_coarse_SS_cortex) = c('Olig_COP', 'Olig_NF', 'Olig_MF', 'Olig_mat', 'Olig_mat',
                                     'Ependymal', 'unknown', 'Astro_Gfap', 'vSMC', 'Pericytes',
                                     'Inh_Crhbp', 'Inh_IC', 'Pyr_L6', 'Periv_Macro', 'Pyr_Cpne5',
                                     'unknown', 'Pyr_L2/3', 'Microglia', 'Hippocampus', 'Pyr_L5',
                                     'Choroid', 'vEnd', 'unknown', 'Inh_Anln', 'Inh_Crh',
                                     'Inh_Vip', 'Inh_Pthlh', 'Pyr_Apln', 'Pyr_Kcnip2', 'Pyr_L4',
                                     'vEnd', 'Astro_Mfge8', 'Olig_precursor')
osm_test = annotateGiotto(gobject = osm_test, annotation_vector = clusters_coarse_SS_cortex,
                          cluster_column = 'det_cell_types', name = 'coarse_cell_types')
png(file = paste0(annotation_folder,'/','annotation_leiden_merged_coarse.png'), height = 920, width = 920)
visSpatDimPlot(gobject = osm_test, cell_color = 'coarse_cell_types', sdimx = 'sdimx', sdimy = 'sdimy',
               dim_point_size = 2, spatial_point_size = 2)
dev.off()


## 8. spatial grid ####
# ------------------- #
grid_folder = paste0(osmFISH_results_folder,'/','8_grid/')
if(!file.exists(grid_folder)) dir.create(grid_folder, recursive = T)

osm_test <- createSpatialGrid(gobject = osm_test,
                              sdimx_stepsize = 2000,
                              sdimy_stepsize = 2000,
                              minimum_padding = 0)
png(file = paste0(grid_folder,'/','grid_det_cell_types.png'), height = 920, width = 920)
visPlot(osm_test, cell_color = 'det_cell_types', sdimx = 'sdimx', sdimy = 'sdimy',
        show_grid = T, grid_color = 'lightblue', spatial_grid_name = 'spatial_grid',
        point_size = 1.5, plot_method = 'ggplot')
dev.off()



#### spatial patterns ####
pattern_osm = detectSpatialPatterns(gobject = osm_test, 
                                    expression_values = 'custom',
                                    spatial_grid_name = 'spatial_grid',
                                    min_cells_per_grid = 5, 
                                    scale_unit = T, 
                                    PC_zscore = 1, 
                                    show_plot = T)
png(file = paste0(grid_folder,'/','pattern1_pca.png'), height = 480, width = 480)
showPattern(pattern_osm, dimension = 1,  plot_dim = 2, point_size = 4)
dev.off()
png(file = paste0(grid_folder,'/','pattern1_pca_genes.png'), height = 320, width = 200)
showPatternGenes(pattern_osm, dimension = 1)
dev.off()

png(file = paste0(grid_folder,'/','pattern3_pca.png'), height = 480, width = 480)
showPattern(pattern_osm, dimension = 3,  plot_dim = 2, point_size = 4)
dev.off()
png(file = paste0(grid_folder,'/','pattern3_pca_genes.png'), height = 320, width = 200)
showPatternGenes(pattern_osm, dimension = 3)
dev.off()



## 9. spatial network ####
## ---------------------##
spatnet_folder = paste0(osmFISH_results_folder,'/','9_spatial_network/')
if(!file.exists(spatnet_folder)) dir.create(spatnet_folder, recursive = T)

osm_test <- createSpatialNetwork(gobject = osm_test, k = 5)
png(file = paste0(spatnet_folder,'/','spatial_network_k5.png'), height = 920, width = 920)
visPlot(gobject = osm_test, show_network = T,
        sdimx = "sdimx",sdimy = "sdimy",
        network_color = 'blue', spatial_network_name = 'spatial_network',
        point_size = 1, cell_color = 'det_cell_types')
dev.off()


## 10. spatial genes ####
## ---------------------##
spatgenes_folder = paste0(osmFISH_results_folder,'/','10_spatial_genes/')
if(!file.exists(spatgenes_folder)) dir.create(spatgenes_folder, recursive = T)

kmtest = binGetSpatialGenes(osm_test, bin_method = 'kmeans',
                            do_fisher_test = T, community_expectation = 5,
                            spatial_network_name = 'spatial_network', verbose = T)
ranktest = binGetSpatialGenes(osm_test, bin_method = 'rank',
                              do_fisher_test = T, community_expectation = 5,
                              spatial_network_name = 'spatial_network', verbose = T)

visSpatDimGenePlot(osm_test, plot_method = 'ggplot', expression_values = 'normalized',
                   genes = c('Rorb', 'Syt6', 'Gfap', 'Kcnip2'),
                   plot_alignment = 'vertical', cow_n_col = 4,
                   genes_high_color = 'red', genes_mid_color = 'white', genes_low_color = 'darkblue', midpoint = 4)

png(file = paste0(spatgenes_folder,'/','spatial_network_k5_genes.png'), height = 400, width = 920)
visSpatDimGenePlot(osm_test, plot_method = 'ggplot', expression_values = 'scaled',
                   genes = c('Rorb', 'Syt6', 'Gfap', 'Kcnip2'),
                   plot_alignment = 'vertical', cow_n_col = 4,
                   genes_high_color = 'red', genes_mid_color = 'white', genes_low_color = 'darkblue', midpoint = 0)
dev.off()



## 11. spatial HMRF domains ####
## ---------------------------##
hmrf_folder = paste0(osmFISH_results_folder,'/','11_HMRF/')
if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = T)

spatial_genes = unique(c(as.character(kmtest[1:15]$genes),
                         as.character(ranktest[1:15]$genes)))

# do HMRF with different betas
HMRF_spatial_genes = doHMRF(gobject = osm_test, expression_values = 'scaled',
                            spatial_genes = spatial_genes,
                            k = 9,
                            betas = c(5,2,3), 
                            output_folder = paste0(hmrf_folder, '/', 'Spatial_genes/SG_top15_k9_scaled'),
                            python_path = "/Users/rubendries/Bin/anaconda3/envs/py36/bin/pythonw")

## view results of HMRF
viewHMRFresults(gobject = osm_test,
                HMRFoutput = HMRF_spatial_genes,
                k = 9, betas_to_view = seq(28, 34, by = 2),
                point_size = 2)

## add HMRF of interest to giotto object
osm_test = addHMRF(gobject = osm_test,
                   HMRFoutput = HMRF_spatial_genes,
                   k = 9, betas_to_add = 30,
                   hmrf_name = 'HMRF')

## visualize
visPlot(gobject = VC_test, cell_color = 'HMRF_k9_b.30', point_size = 2)



## 12. cell-cell preferential proximity ####
## ---------------------------------------##
cellproxim_folder = paste0(osmFISH_results_folder,'/','12_cell_proxim/')
if(!file.exists(cellproxim_folder)) dir.create(cellproxim_folder, recursive = T)

## calculate frequently seen proximities
cell_proximities = cellProximityEnrichment(gobject = osm_test,
                                           cluster_column = 'det_cell_types',
                                           spatial_network_name = 'spatial_network',
                                           number_of_simulations = 400)

## barplot
png(file = paste0(cellproxim_folder,'/','barplot_cell_cell_enrichment.png'), height = 920, width = 920)
cellProximityBarplot(CPscore = cell_proximities, min_orig_ints = 25, min_sim_ints = 25)
dev.off()

## heatmap
png(file = paste0(cellproxim_folder,'/','heatmap_cell_cell_enrichment.png'), height = 920, width = 920)
cellProximityHeatmap(CPscore = cell_proximities, order_cell_types = T, scale = T,
                     color_breaks = c(-1.5, 0, 1.5), color_names = c('blue', 'white', 'red'))
dev.off()

## network
png(file = paste0(cellproxim_folder,'/','network_cell_cell_enrichment.png'), height = 920, width = 920)
cellProximityNetwork(CPscore = cell_proximities)
dev.off()

cellProximityNetwork(CPscore = cell_proximities, remove_self_edges = T)

## visualization of selected cell-cell interaction
spec_interaction = "Astro_Gfap--Olig_mat"
png(file = paste0(cellproxim_folder,'/','cell_cell_enrichment_selected.png'), height = 920, width = 920)
cellProximityVisPlot(gobject = osm_test,
                     interaction_name = spec_interaction,
                     cluster_column = 'det_cell_types',
                     cell_color = 'det_cell_types', coord_fix_ratio = 0.5,
                     point_size_select = 4, point_size_other = 2)
dev.off()



