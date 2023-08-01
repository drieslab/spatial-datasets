
library(data.table)

## 2016 ####


# spatial transcriptomics ####
ST_OB_data_1 = list(
  dataset = 'ST_OB1',
  spatial_locs = "https://raw.githubusercontent.com/drieslab/spatial-datasets/master/data/2016_ST_olfactory_bulb/cell_locations/Rep11_MOB_0_location.txt",
  expr_matrix = "https://raw.githubusercontent.com/drieslab/spatial-datasets/master/data/2016_ST_olfactory_bulb/count_matrix/Rep11_MOB_0_expr.txt",
  metadata = c(NA),
  segmentations = c(NA)
)

ST_OB_data_2 = list(
  dataset = 'ST_OB2',
  spatial_locs = "https://raw.githubusercontent.com/drieslab/spatial-datasets/master/data/2016_ST_olfactory_bulb/cell_locations/Rep12_MOB_0_location.txt",
  expr_matrix = "https://raw.githubusercontent.com/drieslab/spatial-datasets/master/data/2016_ST_olfactory_bulb/count_matrix/Rep12_MOB_0_expr.txt",
  metadata = c(NA),
  segmentations = c(NA)
)



## 2018 ####

# codex ####
codex_spleen_data = list(
  dataset = 'codex_spleen',
  spatial_locs = "https://raw.githubusercontent.com/drieslab/spatial-datasets/master/data/2018_codex_spleen/cell_locations/codex_BALBc_3_coord.txt",
  expr_matrix = "https://raw.githubusercontent.com/drieslab/spatial-datasets/master/data/2018_codex_spleen/count_matrix/codex_BALBc_3_expression.txt.gz",
  metadata = c("https://raw.githubusercontent.com/drieslab/spatial-datasets/master/data/2018_codex_spleen/cell_locations/codex_BALBc_3_annotation.txt",
               "https://raw.githubusercontent.com/drieslab/spatial-datasets/master/data/2018_codex_spleen/cell_locations/cell_type_annotation.csv"),
  segmentations = c(NA)
)


# cyCIF ####
cyCIF_PDAC_data = list(
  dataset = 'cycif_PDAC',
  spatial_locs = "https://raw.githubusercontent.com/drieslab/spatial-datasets/master/data/2018_CyCIF_PDAC/cell_locations/cyCIF_PDAC_coord.txt",
  expr_matrix = "https://raw.githubusercontent.com/drieslab/spatial-datasets/master/data/2018_CyCIF_PDAC/count_matrix/cyCIF_PDAC_expression.txt.gz",
  metadata = c("https://raw.githubusercontent.com/drieslab/spatial-datasets/master/data/2018_CyCIF_PDAC/cell_locations/cyCIF_PDAC_annot.txt",
               "https://raw.githubusercontent.com/drieslab/spatial-datasets/master/data/2018_CyCIF_PDAC/raw_data/canvas.png"),
  segmentations = c(NA)
)


# merfish ####
merfish_preoptic_data = list(
  dataset = 'merfish_preoptic',
  spatial_locs = "https://raw.githubusercontent.com/drieslab/spatial-datasets/master/data/2018_merFISH_science_hypo_preoptic/cell_locations/merFISH_3D_data_cell_locations.txt",
  expr_matrix = "https://raw.githubusercontent.com/drieslab/spatial-datasets/master/data/2018_merFISH_science_hypo_preoptic/count_matrix/merFISH_3D_data_expression.txt.gz",
  metadata = c("https://raw.githubusercontent.com/drieslab/spatial-datasets/master/data/2018_merFISH_science_hypo_preoptic/cell_locations/merFISH_3D_metadata.txt"),
  segmentations = c(NA)
)


# osmfish ####
osmfish_SS_data = list(
  dataset = 'osmfish_SS_cortex',
  spatial_locs = "https://raw.githubusercontent.com/drieslab/spatial-datasets/master/data/2018_osmFISH_SScortex/cell_locations/osmFISH_prep_cell_coordinates.txt",
  expr_matrix = "https://raw.githubusercontent.com/drieslab/spatial-datasets/master/data/2018_osmFISH_SScortex/count_matrix/osmFISH_prep_expression.txt",
  metadata = c("https://raw.githubusercontent.com/drieslab/spatial-datasets/master/data/2018_osmFISH_SScortex/raw_data/osmFISH_prep_cell_metadata.txt"),
  segmentations = c(NA)
)


# starmap ####
starmap_cortex_data = list(
  dataset = 'starmap_3D_cortex',
  spatial_locs = "https://raw.githubusercontent.com/drieslab/spatial-datasets/master/data/2018_starmap_3D_cortex/cell_locations/STARmap_3D_data_cell_locations.txt",
  expr_matrix = "https://raw.githubusercontent.com/drieslab/spatial-datasets/master/data/2018_starmap_3D_cortex/count_matrix/STARmap_3D_data_expression.txt",
  metadata = c(NA),
  segmentations = c(NA)
)




## 2019 ####

# seqfish ss cortex ####
seqfish_SS_data = list(
  dataset = 'seqfish_SS_cortex',
  spatial_locs = "https://raw.githubusercontent.com/drieslab/spatial-datasets/master/data/2019_seqfish_plus_SScortex/cell_locations/cortex_svz_centroids_coord.txt",
  expr_matrix = "https://raw.githubusercontent.com/drieslab/spatial-datasets/master/data/2019_seqfish_plus_SScortex/count_matrix/cortex_svz_expression.txt",
  metadata = c("https://raw.githubusercontent.com/drieslab/spatial-datasets/master/data/2019_seqfish_plus_SScortex/cell_locations/cortex_svz_centroids_annot.txt"),
  segmentations = c(NA)
)

# seqfish OB ####
seqfish_OB_data = list(
  dataset = 'seqfish_OB',
  spatial_locs = "https://raw.githubusercontent.com/drieslab/spatial-datasets/master/data/2019_seqfish_plus_olfactory_bulb/cell_locations/OB_centroids_coord.txt",
  expr_matrix = "https://raw.githubusercontent.com/drieslab/spatial-datasets/master/data/2019_seqfish_plus_olfactory_bulb/count_matrix/OB_expression.txt",
  metadata = c("https://raw.githubusercontent.com/drieslab/spatial-datasets/master/data/2019_seqfish_plus_olfactory_bulb/cell_locations/OB_centroids_annot.txt"),
  segmentations = c(NA)
)


# slideseq cerebellum ####
slideseq_cerebellum_data = list(
  dataset = 'slideseq_cerebellum',
  spatial_locs = "https://zenodo.org/record/4034228/files/BeadLocationsForR.csv",
  expr_matrix = "https://zenodo.org/record/4034228/files/MappedDGEForR.csv",
  metadata = c("https://zenodo.org/record/4034228/files/l1.cerebellum.cellID.txt",
               "https://zenodo.org/record/4034228/files/l1.cerebellum.cells.nonblood.txt",
               "https://zenodo.org/record/4034228/files/l1.cerebellum.cells.txt",
               "https://zenodo.org/record/4034228/files/l1.cerebellum.class.txt",
               "https://zenodo.org/record/4034228/files/l1.cerebellum.genes.txt",
               "https://zenodo.org/record/4034228/files/l1.cerebellum.txt"),
  segmentations = c(NA)
)



## 2020 ####

# ST SCC ####
ST_SCC_data = list(
  dataset = 'ST_SCC',
  spatial_locs = c("https://github.com/drieslab/spatial-datasets/raw/master/data/2020_ST_SCC/cell_locations/P2_1_spatial_locs.csv",
                   "https://github.com/drieslab/spatial-datasets/raw/master/data/2020_ST_SCC/cell_locations/P2_2_spatial_locs.csv",
                   "https://github.com/drieslab/spatial-datasets/raw/master/data/2020_ST_SCC/cell_locations/P2_3_spatial_locs.csv"),
  expr_matrix = c("https://github.com/drieslab/spatial-datasets/raw/master/data/2020_ST_SCC/count_matrix/P2_1_expression.csv",
                  "https://github.com/drieslab/spatial-datasets/raw/master/data/2020_ST_SCC/count_matrix/P2_2_expression.csv",
                  "https://github.com/drieslab/spatial-datasets/raw/master/data/2020_ST_SCC/count_matrix/P2_3_expression.csv"),
  metadata = c("https://github.com/drieslab/spatial-datasets/raw/master/data/2020_ST_SCC/scRNAseq/p2_metadata.rds",
               "https://github.com/drieslab/spatial-datasets/raw/master/data/2020_ST_SCC/scRNAseq/scrna_expr.rds",
               "https://github.com/drieslab/spatial-datasets/raw/master/data/2020_ST_SCC/scRNAseq/normalized_sc_matrix.RDS",
               "https://github.com/drieslab/spatial-datasets/raw/master/data/2020_ST_SCC/scRNAseq/cell_type_vector.RDS",
               "https://github.com/drieslab/spatial-datasets/raw/master/data/2020_ST_SCC/scRNAseq/sign_list.RDS",
               "https://github.com/drieslab/spatial-datasets/raw/master/data/2020_ST_SCC/lowres_images/P2_1_0.0625.jpg",
               "https://github.com/drieslab/spatial-datasets/raw/master/data/2020_ST_SCC/lowres_images/P2_2_0.0625.jpg",
               "https://github.com/drieslab/spatial-datasets/raw/master/data/2020_ST_SCC/lowres_images/P2_3_0.0625.jpg",
               "https://github.com/drieslab/spatial-datasets/raw/master/data/2020_ST_SCC/imgReg.zip",
               "https://github.com/drieslab/spatial-datasets/raw/master/data/2020_ST_SCC/PairsLigRec.txt"),
  segmentations = c(NA)
)


## 2022 ####

# scRNA-seq prostate data
scRNAseq_prostate_data = list(
  dataset = 'scRNA_prostate',
  spatial_locs = c(NA),
  expr_matrix = c("https://github.com/drieslab/spatial-datasets/raw/master/data/2022_scRNAseq_human_prostate/count_matrix/prostate_sc_expression_matrix.csv.gz"),
  metadata = c("https://github.com/drieslab/spatial-datasets/raw/master/data/2022_scRNAseq_human_prostate/cell_metadata/prostate_sc_metadata.csv"),
  segmentations = c(NA)
)


# scRNA-seq mouse visium brain data
scRNAseq_mouse_brain_data = list(
  dataset = 'scRNA_mouse_brain',
  spatial_locs = c(NA),
  expr_matrix = c("https://github.com/drieslab/spatial-datasets/raw/master/data/2022_scRNAseq_mouse_brain/count_matrix/brain_sc_expression_matrix.txt.gz"),
  metadata = c("https://github.com/drieslab/spatial-datasets/raw/master/data/2022_scRNAseq_mouse_brain/cell_metadata/brain_sc_metadata.csv"),
  segmentations = c(NA)
)

# resolve biosciences molecular cartography human lung 873_C1
mol_cart_lung_873_C1 = list(
  dataset = 'mol_cart_lung_873_C1',
  spatial_locs = c(NA),
  expr_matrix = c(NA),
  metadata = c(NA),
  segmentations = c("https://github.com/drieslab/spatial-datasets/raw/master/data/2022_mol_cart_human_lung/segmentations/Resolve_hLung_stardist.geojson")
)

## 2023 ####

# spatial genomics mini kidney data
sg_mini_kidney_data = list(
  dataset = 'sg_mini_kidney',
  spatial_locs = c("https://github.com/drieslab/spatial-datasets/data/2023_spatial_genomics_mouse_kidney/Raw.zip"),
  expr_matrix = c(NA),
  metadata = c(NA),
  segmentations = c(NA)
)


# create datasets.txt
# need to be placed in /extdata directory of Giotto package
datasets = data.table::as.data.table(rbind(ST_OB_data_1,
                                           ST_OB_data_2, 
                                           codex_spleen_data,
                                           cyCIF_PDAC_data,
                                           starmap_cortex_data,
                                           osmfish_SS_data,
                                           merfish_preoptic_data, 
                                           seqfish_SS_data,
                                           seqfish_OB_data,
                                           slideseq_cerebellum_data,
                                           ST_SCC_data,
                                           scRNAseq_prostate_data,
                                           scRNAseq_mouse_brain_data,
                                           mol_cart_lung_873_C1,
                                           sg_mini_kidney_data
                                           ),
                                     row.names = FALSE)
data.table::fwrite(datasets, './datasets.txt', sep = '\t')


