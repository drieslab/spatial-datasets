# atera_mini_tiny

10x Atera v1, human, FFPE (tissue not recorded): a tiny complete run made by 10x Genomics.

Released as [`atera_mini_tiny-v1`](https://github.com/drieslab/spatial-datasets/releases/tag/atera_mini_tiny-v1).

**Source:** 10x Genomics' tiny Atera test run (`run_01`, `slide_01`, `region_01`), Atera v1 chemistry, instrument software 3.2.1.2, analysis software `atera-1.0.0.0`, output format 1.2. The panel is a 53-target R&D panel (`RD`). The run is not otherwise public.

**Crop:** none. This is the whole run as 10x produced it.

**Counts:** 13 cells, 53 targets, 1,011 decoded transcripts (Q20).

**Differences from the vendor layout:** none. All output files are unchanged: zarr archives, `experiment.spatial`, `morphology_2d` / `morphology_3d` OME-TIFFs, metrics and the HTML summaries.

**Provenance:** made by 10x Genomics, so there is no crop script. Read it with `GiottoDisk::importAteraDisk()` or `Giotto::createGiottoAteraObject(backend = )`, and download it with `GiottoData::getSpatialDataset("atera_mini_tiny")`.

**Size:** 3,643,979 bytes. sha256 `acb3a6557d76a25eba42ce4d0056efb841659ab9e899857ff844aaa802d06e28`
