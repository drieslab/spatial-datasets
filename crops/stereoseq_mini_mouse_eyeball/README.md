# stereoseq_mini_mouse_eyeball

A 1.2 mm window of Stereo-seq run C04687E314 (SAW 8.2.0, mouse eyeball), on a
600-unit lattice so every bin level keeps its own grid.

The scripts are the exact files the published crop was made with: their sha256
matches the one recorded in its `crop_manifest.json`. They therefore keep
the paths they ran with; edit `SOURCE_DIR` / `BASE_DIR` in
`make_stereoseq_mini_mouse_eyeball.R` and `BASE` in
`test_stereoseq_mini_mouse_eyeball.R` to rerun elsewhere. Python dependencies
are provisioned through reticulate.
