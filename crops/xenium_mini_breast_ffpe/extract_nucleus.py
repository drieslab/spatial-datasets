#!/usr/bin/env python3
"""Extract nucleus boundaries for the cropped cells from the vendor cells.zarr.

The zarr->parquet conversion the builder consumes emits cell boundaries only,
so nucleus_boundaries.parquet has to come straight from cells.zarr
polygon_sets/0 ("DAPI-based nuclei segmentation"). The published
xenium_mini_lung_ffpe fixture ships this file, which is why the Xenium mini
does too.

Also reports the nucleus segmentation-method split for the kept cells, so
experiment.xenium's segmented_nucleus_* fields can be recomputed for the crop
instead of carrying the section's values.

NOTE: a cell may own more than one nucleus ring -- 6,232 of 254,480 do on this
section. The rings are written one row per vertex, keyed by cell_id, exactly
as the vendor's legacy nucleus_boundaries.parquet is.
"""
import argparse, json, zipfile
import numpy as np, imagecodecs
import pyarrow as pa, pyarrow.parquet as pq

ALPH = "abcdefghijklmnop"


def full(z, name):
    m = json.loads(z.read(name + "/.zarray"))
    ch, sh = m["chunks"], m["shape"]
    out = np.zeros(sh, dtype=m["dtype"])
    nr = -(-sh[0] // ch[0])
    nc = -(-sh[1] // ch[1]) if len(sh) > 1 else 1
    for i in range(nr):
        for j in range(nc):
            key = "%s/%d.%d" % (name, i, j) if len(sh) > 1 else "%s/%d" % (name, i)
            blk = np.frombuffer(imagecodecs.blosc_decode(z.read(key)), dtype=m["dtype"]).reshape(ch)
            r0 = i * ch[0]
            if len(sh) > 1:
                c0 = j * ch[1]
                sub = blk[:min(ch[0], sh[0] - r0), :min(ch[1], sh[1] - c0)]
                out[r0:r0 + sub.shape[0], c0:c0 + sub.shape[1]] = sub
            else:
                sub = blk[:min(ch[0], sh[0] - r0)]
                out[r0:r0 + sub.shape[0]] = sub
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells-zarr", required=True)
    ap.add_argument("--keep-cells", required=True)
    ap.add_argument("--xmin", type=float, required=True)
    ap.add_argument("--ymin", type=float, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--stats-out", required=True)
    a = ap.parse_args()

    z = zipfile.ZipFile(a.cells_zarr)
    cid = full(z, "cell_id")
    nib = np.stack([(cid[:, 0] >> (28 - 4 * i)) & 0xF for i in range(8)], axis=1)
    names = np.array(["".join(r) for r in np.array(list(ALPH))[nib]])
    names = np.char.add(names, np.char.add("-", cid[:, 1].astype(str)))

    idx = full(z, "polygon_sets/0/cell_index")
    nv = full(z, "polygon_sets/0/num_vertices").ravel()
    meth = full(z, "polygon_sets/0/method").ravel()
    verts = full(z, "polygon_sets/0/vertices")
    idx = idx.ravel()

    keep = pq.read_table(a.keep_cells, columns=["cell_id"]).column("cell_id").to_pylist()
    keep_set = set(keep)
    ring_name = names[idx]
    sel = np.fromiter((n in keep_set for n in ring_name), dtype=bool, count=len(ring_name))
    print("nucleus rings in fixture: %d over %d cells (%d cells have >1 ring)"
          % (sel.sum(), len(set(ring_name[sel])), sel.sum() - len(set(ring_name[sel]))))

    desc = json.loads(z.read(".zattrs"))["segmentation_methods"]
    mm = {}
    for m_ in np.unique(meth[sel]):
        mm[desc[int(m_)]] = int((meth[sel] == int(m_)).sum())
    print("nucleus methods:", json.dumps(mm))

    rows_id, rows_x, rows_y = [], [], []
    for r in np.nonzero(sel)[0]:
        k = int(nv[r])
        v = verts[r, :2 * k]
        rows_id.extend([ring_name[r]] * k)
        rows_x.append(v[0::2]); rows_y.append(v[1::2])
    vx = np.concatenate(rows_x).astype(np.float64) - a.xmin
    vy = np.concatenate(rows_y).astype(np.float64) - a.ymin

    tbl = pa.table({"cell_id": pa.array(rows_id, type=pa.string()),
                    "vertex_x": pa.array(vx), "vertex_y": pa.array(vy)})
    # same ordering rule as every other file the builder writes, so the
    # fixture's sha256 does not depend on iteration order
    tbl = tbl.sort_by([("cell_id", "ascending"),
                       ("vertex_y", "ascending"), ("vertex_x", "ascending")])
    pq.write_table(tbl, a.out, compression="zstd")
    print("wrote %s: %d vertex rows" % (a.out, tbl.num_rows))

    json.dump({"rings": int(sel.sum()),
               "cells_with_nucleus": len(set(ring_name[sel])),
               "methods": mm}, open(a.stats_out, "w"), indent=1)


if __name__ == "__main__":
    main()
