
library(data.table)

## 2016 ##

# spatial transcriptomics
ST_OB_data_1 = list(
  dataset = 'ST_OB1',
  spatial_locs = "https://raw.githubusercontent.com/RubD/spatial-datasets/master/data/2016_ST_olfactory_bulb/cell_locations/Rep11_MOB_0_location.txt",
  expr_matrix = "https://raw.githubusercontent.com/RubD/spatial-datasets/master/data/2016_ST_olfactory_bulb/count_matrix/Rep11_MOB_0_expr.txt",
  metadata = c(NA)
)

ST_OB_data_2 = list(
  dataset = 'ST_OB2',
  spatial_locs = "https://raw.githubusercontent.com/RubD/spatial-datasets/master/data/2016_ST_olfactory_bulb/cell_locations/Rep12_MOB_0_location.txt",
  expr_matrix = "https://raw.githubusercontent.com/RubD/spatial-datasets/master/data/2016_ST_olfactory_bulb/count_matrix/Rep12_MOB_0_expr.txt",
  metadata = c(NA)
)



## 2018 ##

# codex
codex_spleen_data = list(
  dataset = 'codex_spleen',
  spatial_locs = "https://raw.githubusercontent.com/RubD/spatial-datasets/master/data/2018_codex_spleen/cell_locations/codex_BALBc_3_coord.txt",
  expr_matrix = "https://raw.githubusercontent.com/RubD/spatial-datasets/master/data/2018_codex_spleen/count_matrix/codex_BALBc_3_expression.txt.gz",
  metadata = c("https://raw.githubusercontent.com/RubD/spatial-datasets/master/data/2018_codex_spleen/cell_locations/codex_BALBc_3_annotation.txt",
               "https://raw.githubusercontent.com/RubD/spatial-datasets/master/data/2018_codex_spleen/cell_locations/cell_type_annotation.csv")
)


# cyCIF
cyCIF_PDAC_data = list(
  dataset = 'cycif_PDAC',
  spatial_locs = "https://raw.githubusercontent.com/RubD/spatial-datasets/master/data/2018_CyCIF_PDAC/cell_locations/cyCIF_PDAC_coord.txt",
  expr_matrix = "https://raw.githubusercontent.com/RubD/spatial-datasets/master/data/2018_CyCIF_PDAC/count_matrix/cyCIF_PDAC_expression.txt.gz",
  metadata = c("https://raw.githubusercontent.com/RubD/spatial-datasets/master/data/2018_CyCIF_PDAC/cell_locations/cyCIF_PDAC_annot.txt")
)


# merfish
merfish_preoptic_data = list(
  dataset = 'merfish_preoptic',
  spatial_locs = "https://raw.githubusercontent.com/RubD/spatial-datasets/master/data/2018_merFISH_science_hypo_preoptic/cell_locations/merFISH_3D_data_cell_locations.txt",
  expr_matrix = "https://raw.githubusercontent.com/RubD/spatial-datasets/master/data/2018_merFISH_science_hypo_preoptic/count_matrix/merFISH_3D_data_expression.txt.gz",
  metadata = c("https://raw.githubusercontent.com/RubD/spatial-datasets/master/data/2018_merFISH_science_hypo_preoptic/cell_locations/merFISH_3D_metadata.txt")
)


# osmfish
osmfish_SS_data = list(
  dataset = 'osmfish_SS_cortex',
  spatial_locs = "https://raw.githubusercontent.com/RubD/spatial-datasets/master/data/2018_osmFISH_SScortex/cell_locations/osmFISH_prep_cell_coordinates.txt",
  expr_matrix = "https://raw.githubusercontent.com/RubD/spatial-datasets/master/data/2018_osmFISH_SScortex/count_matrix/osmFISH_prep_expression.txt",
  metadata = c(NA)
)


# starmap
starmap_cortex_data = list(
  dataset = 'starmap_3D_cortex',
  spatial_locs = "https://raw.githubusercontent.com/RubD/spatial-datasets/master/data/2018_starmap_3D_cortex/cell_locations/STARmap_3D_data_cell_locations.txt",
  expr_matrix = "https://raw.githubusercontent.com/RubD/spatial-datasets/master/data/2018_starmap_3D_cortex/count_matrix/STARmap_3D_data_expression.txt",
  metadata = c(NA)
)




## 2019

# seqfish ss cortex
seqfish_SS_data = list(
  dataset = 'seqfish_SS_cortex',
  spatial_locs = "https://raw.githubusercontent.com/RubD/spatial-datasets/master/data/2019_seqfish_plus_SScortex/cell_locations/cortex_svz_centroids_coord.txt",
  expr_matrix = "https://raw.githubusercontent.com/RubD/spatial-datasets/master/data/2019_seqfish_plus_SScortex/count_matrix/cortex_svz_expression.txt",
  metadata = c("https://raw.githubusercontent.com/RubD/spatial-datasets/master/data/2019_seqfish_plus_SScortex/cell_locations/cortex_svz_centroids_annot.txt")
)

# seqfish OB
seqfish_OB_data = list(
  dataset = 'seqfish_OB',
  spatial_locs = "https://raw.githubusercontent.com/RubD/spatial-datasets/master/data/2019_seqfish_plus_olfactory_bulb/cell_locations/OB_centroids_coord.txt",
  expr_matrix = "https://raw.githubusercontent.com/RubD/spatial-datasets/master/data/2019_seqfish_plus_olfactory_bulb/count_matrix/OB_expression.txt",
  metadata = c("https://raw.githubusercontent.com/RubD/spatial-datasets/master/data/2019_seqfish_plus_olfactory_bulb/cell_locations/OB_centroids_annot.txt")
)






# create datasets.txt
# need to be placed in /extdata directory of Giotto package
datasets = as.data.table(rbind(ST_OB_data_1,
                               ST_OB_data_2, 
                               codex_spleen_data,
                               cyCIF_PDAC_data,
                               starmap_cortex_data,
                               osmfish_SS_data,
                               merfish_preoptic_data, 
                               seqfish_SS_data,
                               seqfish_OB_data), row.names = F)
fwrite(datasets, './datasets.txt', sep = '\t')




