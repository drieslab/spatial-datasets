#!/usr/bin/env python3
"""Recover per-transcript cell assignment for the Xenium mini from the vendor
segmentation masks.

The transcripts zarr carries no per-transcript cell_id -- its tile arrays are
codeword_identity, gene_identity, gene_offset, id, location, quality_score,
status, uuid, valid.  The assignment lives in cells.zarr as masks/1, 10x's own
cell-segmentation label raster, whose value at a pixel is the 1-based row index
into the cells array.  masks/0 is the matching nucleus raster.

Validated before use: the 32-bit cell_id encoding reproduces all 254,480
staged cell_id strings exactly, and the mask label at a cell centroid equals
its row index + 1 for 400/400 sampled cells.
"""
import argparse, json, zipfile
import numpy as np, imagecodecs
import pyarrow as pa, pyarrow.parquet as pq

ALPH = "abcdefghijklmnop"
PX_UM = 0.2125


def _zarr_chunk(z, name, meta, i, j):
    return np.frombuffer(imagecodecs.blosc_decode(z.read("%s/%d.%d" % (name, i, j))),
                         dtype=meta["dtype"]).reshape(meta["chunks"])


def read_window(z, name, y0, y1, x0, x1):
    """Rows [y0,y1) and cols [x0,x1) of a 2-D blosc zarr array."""
    meta = json.loads(z.read(name + "/.zarray"))
    ch, sh = meta["chunks"], meta["shape"]
    out = np.zeros((y1 - y0, x1 - x0), dtype=meta["dtype"])
    for ci in range(y0 // ch[0], (y1 - 1) // ch[0] + 1):
        for cj in range(x0 // ch[1], (x1 - 1) // ch[1] + 1):
            blk = _zarr_chunk(z, name, meta, ci, cj)
            by0, bx0 = ci * ch[0], cj * ch[1]
            ys, ye = max(y0, by0), min(y1, by0 + ch[0], sh[0])
            xs, xe = max(x0, bx0), min(x1, bx0 + ch[1], sh[1])
            if ys >= ye or xs >= xe:
                continue
            out[ys - y0:ye - y0, xs - x0:xe - x0] = blk[ys - by0:ye - by0, xs - bx0:xe - bx0]
    return out


def read_cell_ids(z):
    meta = json.loads(z.read("cell_id/.zarray"))
    ch, sh = meta["chunks"], meta["shape"]
    out = np.zeros(sh, dtype=meta["dtype"])
    for i in range(-(-sh[0] // ch[0])):
        for j in range(-(-sh[1] // ch[1])):
            blk = _zarr_chunk(z, "cell_id", meta, i, j)
            r0, c0 = i * ch[0], j * ch[1]
            sub = blk[:min(ch[0], sh[0] - r0), :min(ch[1], sh[1] - c0)]
            out[r0:r0 + sub.shape[0], c0:c0 + sub.shape[1]] = sub
    return out


def encode_ids(pref, suf):
    """10x's cell_id string: 8 nibbles of the uint32 over 'a'..'p', then -suffix."""
    nib = np.stack([(pref >> (28 - 4 * i)) & 0xF for i in range(8)], axis=1)
    chars = np.array(list(ALPH))[nib]
    return np.char.add(np.array(["".join(r) for r in chars]),
                       np.char.add("-", suf.astype(str)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells-zarr", required=True)
    ap.add_argument("--transcripts", required=True, help="the fixture's transcripts.parquet")
    ap.add_argument("--keep-cells", required=True, help="the fixture's cells.parquet")
    ap.add_argument("--xmin", type=float, required=True, help="crop origin in SOURCE um")
    ap.add_argument("--ymin", type=float, required=True)
    ap.add_argument("--window", type=float, required=True)
    ap.add_argument("--out", help="write a patched transcripts.parquet here")
    ap.add_argument("--stats-out", help="write the counts below as json here")
    a = ap.parse_args()

    z = zipfile.ZipFile(a.cells_zarr)
    s = 1.0 / PX_UM
    x0, y0 = int(np.floor(a.xmin * s)), int(np.floor(a.ymin * s))
    x1, y1 = int(np.ceil((a.xmin + a.window) * s)), int(np.ceil((a.ymin + a.window) * s))
    print("crop pixels: x [%d,%d) y [%d,%d)  -> %d x %d" % (x0, x1, y0, y1, x1 - x0, y1 - y0))

    cell_mask = read_window(z, "masks/1", y0, y1, x0, x1)
    nuc_mask = read_window(z, "masks/0", y0, y1, x0, x1)
    print("masks read: cell labels %d distinct, nucleus labels %d distinct"
          % (len(np.unique(cell_mask)), len(np.unique(nuc_mask))))

    cid = read_cell_ids(z)
    all_ids = encode_ids(cid[:, 0], cid[:, 1])

    tx = pq.read_table(a.transcripts)
    # the fixture is translated to the origin; put it back in the source frame
    sx = tx.column("x_location").to_numpy() + a.xmin
    sy = tx.column("y_location").to_numpy() + a.ymin
    px = np.clip((sx * s).astype(np.int64) - x0, 0, x1 - x0 - 1)
    py = np.clip((sy * s).astype(np.int64) - y0, 0, y1 - y0 - 1)

    lab = cell_mask[py, px]                      # 0 = background
    nlab = nuc_mask[py, px]

    # Referential integrity: a transcript may fall in a cell whose centroid is
    # outside the window, and that cell is not in cells.parquet. Pointing at it
    # would be a dangling reference, so it stays UNASSIGNED and is counted.
    keep = set(pq.read_table(a.keep_cells, columns=["cell_id"]).column("cell_id").to_pylist())
    name = np.where(lab > 0, all_ids[np.clip(lab, 1, len(all_ids)) - 1], "UNASSIGNED")
    in_fixture = np.fromiter((n in keep for n in name), dtype=bool, count=len(name))
    dangling = int((lab > 0).sum() - in_fixture.sum())
    final = np.where(in_fixture, name, "UNASSIGNED")

    n = len(final)
    assigned = int(in_fixture.sum())
    print()
    print("transcripts            : %d" % n)
    print("in a cell per the mask : %d (%.4f)" % (int((lab > 0).sum()), (lab > 0).mean()))
    print("  of those, cell outside the window (left UNASSIGNED): %d" % dangling)
    print("ASSIGNED in fixture    : %d (%.4f)" % (assigned, assigned / n))
    print("distinct cells hit     : %d" % len(set(final[in_fixture])))
    print("overlaps_nucleus       : %d (%.4f)" % (int((nlab > 0).sum()), (nlab > 0).mean()))

    if a.stats_out:
        # the builder writes these into crop_manifest.json, so they are
        # measured here rather than restated by hand from this log
        with open(a.stats_out, "w") as fh:
            json.dump({"transcripts": n,
                       "in_cell_per_mask": int((lab > 0).sum()),
                       "cell_outside_window": dangling,
                       "assigned": assigned,
                       "distinct_cells": len(set(final[in_fixture])),
                       "overlaps_nucleus": int((nlab > 0).sum())}, fh)

    if a.out:
        cols = {c: tx.column(c) for c in tx.column_names}
        cols["cell_id"] = pa.array(final, type=pa.string())
        cols["overlaps_nucleus"] = pa.array((nlab > 0).astype(np.int32), type=pa.int32())
        pq.write_table(pa.table(cols).select(tx.column_names), a.out, compression="zstd")
        print("\nwrote %s" % a.out)


if __name__ == "__main__":
    main()
