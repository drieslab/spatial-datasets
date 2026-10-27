#!/usr/bin/env python3
"""Keep the features of a Space Ranger segmentation GeoJSON whose integer
`cell_id` is in a list, and translate every vertex by (-dx, -dy) fullres px.

usage: visiumhd_crop_geojson.py SRC DST IDS_TXT DX DY
IDS_TXT holds one integer cell_id per line. Standard library only.
"""
import json
import sys


def shift(coords, dx, dy):
    if coords and isinstance(coords[0], (int, float)):
        return [round(coords[0] - dx, 6), round(coords[1] - dy, 6)] + coords[2:]
    return [shift(c, dx, dy) for c in coords]


def main(src, dst, ids_txt, dx, dy):
    dx, dy = float(dx), float(dy)
    with open(ids_txt) as fh:
        keep = {int(line) for line in fh if line.strip()}
    with open(src) as fh:
        g = json.load(fh)
    n_src = len(g["features"])
    feats = [f for f in g["features"] if f["properties"]["cell_id"] in keep]
    for f in feats:
        f["geometry"]["coordinates"] = shift(f["geometry"]["coordinates"], dx, dy)
    g["features"] = feats
    with open(dst, "w") as fh:
        json.dump(g, fh, separators=(", ", ": "))
    print(json.dumps({"features_source": n_src, "features_kept": len(feats)}))


if __name__ == "__main__":
    main(*sys.argv[1:6])
