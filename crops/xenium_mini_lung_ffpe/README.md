# xenium_mini_lung_ffpe

A 425 µm window (1531 cells, 541 features) of 10x's
`Xenium_V1_humanLung_Cancer_FFPE_outs` (Xenium 2.0.0.6). It was made under the
name `xenium_mini_lung_cancer`, which is still its `fixture_name`.

The scripts are the exact files the published crop was made with: their sha256
matches the one recorded in its `crop_manifest.json`. They therefore keep
the paths they ran with; edit `SOURCE_DIR` / `OUT_DIR` in
`make_xenium_mini_lung_cancer.R` and `OUT_DIR` / `PNG_OUT` in
`test_xenium_mini_lung_cancer.R` to rerun elsewhere.
