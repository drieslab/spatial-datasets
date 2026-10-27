# Cropped vendor-format minis

Small spatial crops of public vendor runs, published as assets on the
[`drieslab/spatial-datasets` releases](https://github.com/drieslab/spatial-datasets/releases)
and downloaded with `GiottoData::getSpatialDataset("<name>")`. Each
folder holds the script that cut the crop from the full run and the check that
accepts it, where the script is available.

Each crop gets its own `<dataset>-v<N>` release and ships a `crop_manifest.json`
valid against GiottoData's
[crop manifest schema](https://github.com/giotto-suite/GiottoData/blob/dev/.claude/skills/add-spatial-dataset/crop_manifest.schema.json),
recording the source run, the window and the script's sha256. The lung,
Stereo-seq and CosMx crops predate this: they share the `crops-2026-09`
release, and their manifests, where present, follow an earlier layout.

| dataset | technology | made by | checked by |
|---|---|---|---|
| `xenium_mini_lung_ffpe` | 10x Xenium | `make_xenium_mini_lung_cancer.R` | `test_xenium_mini_lung_cancer.R` |
| `stereoseq_mini_mouse_eyeball` | Stereo-seq | `make_stereoseq_mini_mouse_eyeball.R` + `stereoseq_crop.py` | `test_stereoseq_mini_mouse_eyeball.R` + `audit_fixture.py` |
| `cosmx_mini_colon_wtx` | Bruker CosMx | see its README | see its README |
| `visiumhd_mini_crc` | 10x Visium HD | `make_visiumhd_mini_crc.R` + `visiumhd_crop_geojson.py` | `test_visiumhd_mini_crc.R`, `workflow_visiumhd_mini_crc.R` |

To rerun a crop, download the full run (see each README), point the script at
it and run it with `Rscript`. Rerunning gives the same content; HDF5 files
differ byte-wise because they embed write times.
