# xenium_mini_breast_ffpe

A 1200 µm square window (6,764 cells, 280 genes in a 541-row feature axis,
3,197,806 transcripts) of 10x's
[Xenium FFPE Human Breast Biomarkers](https://www.10xgenomics.com/datasets/xenium-ffpe-human-breast-biomarkers)
run, section `S1_Bot`, panel `Xenium Human Breast Cancer HK Gene Panel (280g
Custom)`, Xenium Onboard Analysis 4.0.1.1.

The source is a **zarr export** (`.zarr.zip`), not a legacy `_outs` bundle, so
the crop is rebuilt rather than sliced: every table is written in 10x's legacy
parquet/MEX layout and the fixture loads under stock
`createGiottoXeniumObject(xenium_dir = ...)` with no GiottoDisk and no zarr
reader. It therefore does not exercise the zarr path.

Coordinates are translated to the origin, so the fixture spans `[0, 1200]` µm
in both axes and lines up with the cropped DAPI image, which every reader
anchors at (0, 0). `crop_manifest.json` records the shift under
`crop.translation`; add it back to compare against the source run.

Two things a reader should know before relying on it:

- **Per-transcript `cell_id` is recovered, not copied.** `transcripts.zarr`
  carries no cell assignment, so it is read back out of the vendor's own
  segmentation label rasters (`cells.zarr` `masks/1` and `masks/0`). Measured
  against the vendor matrix: ratio 0.99623, 72.4 % of cells exact, median
  per-cell ratio 1.00000, r 0.998620. The lookup is 2D where the vendor's is
  z-aware, so a small residual remains.
- **The nucleus layer does not join to `cell_id`.** Some cells own more than
  one nucleus ring, and Giotto then names the polygons `nucleus_1`,
  `nucleus_2`, so `subsetGiotto()` drops the layer with a warning. This is
  upstream behaviour, not a property of this crop: the published
  `xenium_mini_lung_ffpe` fixture does the same.

`crop_manifest.json` (schema v1) records the source, window, translation,
per-file actions, counts and every deviation from the vendor layout.

## Rebuilding

```sh
XENIUM_ZARR_DIR=/path/to/source_outs \
XENIUM_MINI_OUTDIR=/path/to/out/xenium_mini_breast_ffpe \
  Rscript make_xenium_mini.R 1200

./package_mini.sh /path/to/out /path/to/xenium_mini_breast_ffpe.zip

XENIUM_MINI_DIR=/path/to/out/xenium_mini_breast_ffpe \
  Rscript test_xenium_mini_breast_ffpe.R
```

`make_xenium_mini.R` refuses to overwrite an existing fixture; set
`XENIUM_MINI_OVERWRITE=1` to replace one. Set `SOURCE_DATE_EPOCH` and
`MINI_ZIP_MTIME` to the same instant for a reproducible archive — two builds
in different directories then produce the same sha256.

R needs arrow, GiottoDisk, data.table, Matrix and jsonlite. The three python
helpers need tifffile, imagecodecs and pyarrow; point `XENIUM_MINI_PYTHON` at
an interpreter that has them. Python does the image and mask work because the
source morphology image is JPEG-2000 inside BigTIFF, which GDAL cannot open.

- `recover_assignment.py` recovers `cell_id` / `overlaps_nucleus` from the masks.
- `extract_nucleus.py` rebuilds `nucleus_boundaries.parquet` from
  `polygon_sets/0`, which the zarr→parquet conversion omits.
- `crop_morphology.py` crops one morphology channel and re-encodes it to
  deflate so terra can read it.
- `patch_experiment.py` applies the `experiment.xenium` edits; it refuses to
  add keys the vendor file does not already have.

`test_xenium_mini_breast_ffpe.R` accepts the crop: the importer under default
arguments, a downstream pipeline, every count in `crop_manifest.json` against
the files, image registration, size, and internal consistency. Its pixel-exact
comparison against the source image is opt-in
(`XENIUM_TEST_PIXEL_EXACT=1`) because it reads ~4.4 GB.
