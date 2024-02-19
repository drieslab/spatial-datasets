# spatial datasets
Repository of spatial datasets analyzed by [Giotto](https://rubd.github.io/Giotto/). Each dataset has its own folder with subfolders for:  
- count matrix  
- locations  
- giotto R-script with secondary analyis (when available)   
- raw data  (optional)  

These datasets can be downloaded directly from within the Giotto package using the function **getSpatialDataset** as illustrated in the examples. 

Visium datasets are not available, but the original directory can be downloaded from the 10X Visium website and directly used as input for the wrapper
**createGiottoVisiumObject** to create a Giotto object.  




## Available datasets

### 2023

#### Spatial Genomics Mouse Kidney Mini Dataset
- [**website**]([https://resolvebiosciences.com/open-dataset/?dataset=human-lung-dataset](https://spatialgenomics.com/data/#downloads))

### 2022

#### Molecular Cartography Human Lung
- [**website**](https://resolvebiosciences.com/open-dataset/?dataset=human-lung-dataset)
- [directory](./data/2022_mol_cart_human_lung/)

#### single-cell Human Prostate
- [**paper**]()
- [directory](./data/2022_scRNAseq_human_prostate)


#### single-cell Mouse Brain
- [**paper**]() 
- [directory](./data/2022_scRNAseq_mouse_brain)


-------------------------------------------------------------------------


### 2020

#### Spatial Transcriptomics squamous cell carcinoma
- [**paper**](https://www.sciencedirect.com/science/article/pii/S0092867420306723)  
- [directory](./data/2020_ST_SCC)

-------------------------------------------------------------------------

### 2019

#### Visium kidney
- [**website**](https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Mouse_Kidney)  

#### Visium brain
- [**website**](https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Adult_Mouse_Brain)  

#### merFISH 
- [**paper**](https://www.pnas.org/content/116/39/19490)  
- [directory](./data/2019_merFISH_PNAS)

#### seqFISH+ somatosensory cortex
- [**paper**](https://www.nature.com/articles/s41586-019-1049-y)  
- [directory](./data/2019_seqfish_plus_SScortex/)

#### seqFISH+ olfactory bulb
- [**paper**](https://www.nature.com/articles/s41586-019-1049-y)  
- [directory](./data/2019_seqfish_plus_olfactory_bulb/)

#### slide-seq cerebellum
- [**paper**](https://science.sciencemag.org/content/363/6434/1463)  
- [directory](./data/2019_slideseq_cerebellum/)

-------------------------------------------------------------------------


### 2018

#### STARmap 3D cortex
- [**paper**](https://science.sciencemag.org/content/361/6400/eaat5691)  
- [directory](./data/2018_starmap_3D_cortex/)

#### osmFISH somatosensory cortex
- [**paper**](https://www.nature.com/articles/s41592-018-0175-z)  
- [directory](./data/2018_osmFISH_SScortex/)

#### merFISH 3D hypothalamic preoptic region
- [**paper**](https://science.sciencemag.org/content/362/6416/eaau5324)
- [directory](./data/2018_merFISH_science_hypo_preoptic/)

#### CODEX spleen
- [**paper**](http://www.sciencedirect.com/science/article/pii/S0092867418309048)
- [directory](./data/2018_codex_spleen/)

#### CyCIF PDAC
- [**paper**](https://doi.org/10.7554/eLife.31657)
- [directory](./data/2018_CyCIF_PDAC/)

#### MIBI TNBC
- [**paper**](https://www.cell.com/fulltext/S0092-8674(18)31100-0)
- [directory](./data/2019_slideseq_cerebellum/)
- not yet available

--------------------------------------------------------------------------

### 2016

#### Spatial Transcriptomics olfactory bulb
- [**paper**](https://science.sciencemag.org/content/353/6294/78)
- [directory](./data/2016_ST_olfactory_bulb/)

## Soon available datasets  

1. slide-seq
2. MIBI
3. ...
