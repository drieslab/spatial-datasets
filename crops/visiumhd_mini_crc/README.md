# visiumhd_mini_crc

A 400 µm window of 10x's Visium HD human colorectal cancer FFPE run (Space
Ranger 4.0.1), across the tissue edge (55 % of 2 µm bins in tissue): 2, 8 and
16 µm binned outputs (h5 and matrix-market), 806 segmented cells with their
cell and nucleus geojson, `barcode_mappings.parquet` and the cropped hires and
lowres images. The window sits on the 16 µm grid and cells are kept only when
all their 2 µm bins are inside, so the 8 / 16 µm and cell matrices are exact
sums of the kept 2 µm bins. Fullres pixel coordinates are translated so the
cropped images start at (0, 0); barcodes and array indices are unchanged.

Source: `binned_outputs`, `segmented_outputs`, `spatial` and
`barcode_mappings.parquet` from
`https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_<file>`,
extracted into one `outs/`.

```sh
VHD_SOURCE_DIR=/path/to/outs VHD_OUT_DIR=/path/to/visiumhd_mini_crc \
  Rscript make_visiumhd_mini_crc.R
VHD_SOURCE_DIR=/path/to/outs VHD_OUT_DIR=/path/to/visiumhd_mini_crc \
  Rscript test_visiumhd_mini_crc.R
VHD_MINI_DIR=/path/to/visiumhd_mini_crc Rscript workflow_visiumhd_mini_crc.R
```

- `test_visiumhd_mini_crc.R` accepts the crop: internal consistency, fidelity
  to the source, both readers, image alignment, size.
- `workflow_visiumhd_mini_crc.R` runs every reader option and the website
  tutorial's pipeline in memory and on disk, and writes a report (one line per
  step, plus memory-vs-disk comparisons).

Published as
[`visiumhd_mini_crc-v1`](https://github.com/drieslab/spatial-datasets/releases/tag/visiumhd_mini_crc-v1).
`make_visiumhd_mini_crc.R` and `visiumhd_crop_geojson.py` are the exact files
that cut it; their sha256 is in its `crop_manifest.json`. The manifest they
write predates the schema, so the published one was converted to schema v1
without rerunning the crop. A v2 rerun updates the script to write schema v1
first.
