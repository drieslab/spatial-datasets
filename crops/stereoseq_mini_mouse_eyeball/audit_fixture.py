"""Structural audit of the cropped Stereo-seq fixture (called from the test script)."""
import os
import h5py
import numpy as np

BORDER_PAD = 32767


def audit(out_dir, source_dir, sn, x0, y0, size):
    res = []
    def ck(name, cond, detail=""):
        res.append({"name": name, "pass": bool(cond), "detail": str(detail)})

    fe = os.path.join(out_dir, "feature_expression")

    # ---------------- bin GEF ----------------
    with h5py.File(os.path.join(fe, "%s.tissue.gef" % sn), "r") as f:
        mids = {}
        for b in sorted(f["geneExp"], key=lambda s: int(s[3:])):
            bs = int(b[3:])
            e = f["geneExp"][b]["expression"][:]
            g = f["geneExp"][b]["gene"][:]
            cnt = g["count"].astype(np.int64)
            ck("bin GEF %s: offset == exclusive cumsum" % b,
               np.array_equal(g["offset"].astype(np.int64),
                              np.concatenate(([0], np.cumsum(cnt)[:-1]))))
            ck("bin GEF %s: sum(count) == nrow(expression)" % b,
               cnt.sum() == len(e), "%d vs %d" % (cnt.sum(), len(e)))
            ck("bin GEF %s: coords inside [0,%d)" % (b, size),
               e["x"].min() >= 0 and e["y"].min() >= 0 and
               e["x"].max() < size and e["y"].max() < size,
               "x[%d,%d] y[%d,%d]" % (e["x"].min(), e["x"].max(), e["y"].min(), e["y"].max()))
            ck("bin GEF %s: bin origins are multiples of %d" % (b, bs),
               bool((e["x"] % bs == 0).all() and (e["y"] % bs == 0).all()))
            a = dict(f["geneExp"][b]["expression"].attrs)
            ck("bin GEF %s: extent attrs match data" % b,
               int(a["minX"][0]) == 0 and int(a["minY"][0]) == 0 and
               int(a["maxX"][0]) == int(e["x"].max()) and
               int(a["maxY"][0]) == int(e["y"].max()) and
               int(a["maxExp"][0]) == int(e["count"].max()))
            # wholeExp grid must reproduce from the expression rows
            w = f["wholeExp"][b]; W = w[:]; wa = dict(w.attrs)
            nx = ny = size // bs
            flat = (e["x"] // bs).astype(np.int64) * ny + (e["y"] // bs).astype(np.int64)
            mid = np.bincount(flat, weights=e["count"].astype(np.int64),
                              minlength=nx * ny).reshape(nx, ny)
            gct = np.bincount(flat, minlength=nx * ny).reshape(nx, ny)
            ck("bin GEF %s: wholeExp grid reproduces" % b,
               W.shape == (nx, ny) and
               np.array_equal(W["MIDcount"], mid.astype(W.dtype["MIDcount"])) and
               np.array_equal(W["genecount"], gct.astype(W.dtype["genecount"])))
            ck("bin GEF %s: wholeExp attrs match" % b,
               int(wa["lenX"][0]) == size and int(wa["lenY"][0]) == size and
               int(wa["number"][0]) == int((gct > 0).sum()) and
               int(wa["maxMID"][0]) == int(mid.max()) and
               int(wa["maxGene"][0]) == int(gct.max()))
            mids[b] = int(e["count"].astype(np.int64).sum())
        ck("bin GEF: total MID identical across all 8 bin levels",
           len(set(mids.values())) == 1, sorted(set(mids.values())))
        ck("bin GEF: stat/ and contour/ dropped",
           "stat" not in f and "contour" not in f)

    # ---------------- cell-bin GEFs ----------------
    for fn in ["%s.cellbin.gef" % sn, "%s.adjusted.cellbin.gef" % sn]:
        tag = fn.replace("%s." % sn, "")
        with h5py.File(os.path.join(fe, fn), "r") as f:
            cb = f["cellBin"]
            cell = cb["cell"][:]; gene = cb["gene"][:]
            ce = cb["cellExp"][:]; ge = cb["geneExp"][:]
            n = len(cell)
            ck("%s: cell ids are 0..n-1" % tag,
               np.array_equal(cell["id"], np.arange(n, dtype=cell.dtype["id"])))
            ck("%s: cell offset == cumsum(geneCount)" % tag,
               np.array_equal(cell["offset"].astype(np.int64),
                              np.concatenate(([0], np.cumsum(cell["geneCount"].astype(np.int64))[:-1]))))
            ck("%s: sum(geneCount) == nrow(cellExp)" % tag,
               cell["geneCount"].astype(np.int64).sum() == len(ce))
            cc = gene["cellCount"].astype(np.int64)
            ck("%s: sum(cellCount) == nrow(geneExp)" % tag,
               cc.sum() == len(ge), "%d vs %d" % (cc.sum(), len(ge)))
            nz = cc > 0
            off = gene["offset"].astype(np.int64)
            ck("%s: non-zero gene offsets contiguous" % tag,
               np.array_equal(off[nz], np.concatenate(([0], np.cumsum(cc[nz])[:-1]))))
            ck("%s: zero-count genes have offset 0" % tag, bool((off[~nz] == 0).all()))
            # the two indexes must encode the same triplets
            ci = np.repeat(np.arange(n, dtype=np.int64), cell["geneCount"].astype(np.int64))
            gi = np.repeat(np.arange(len(gene), dtype=np.int64), cc)
            A = np.stack([ci, ce["geneID"].astype(np.int64), ce["count"].astype(np.int64)], 1)
            B = np.stack([ge["cellID"].astype(np.int64), gi, ge["count"].astype(np.int64)], 1)
            A = A[np.lexsort((A[:, 1], A[:, 0]))]
            B = B[np.lexsort((B[:, 1], B[:, 0]))]
            ck("%s: cellExp and geneExp encode identical triplets" % tag,
               np.array_equal(A, B))
            ck("%s: every geneExp cellID in range" % tag,
               bool((ge["cellID"] < n).all()))
            ck("%s: every cellExp geneID in range" % tag,
               bool((ce["geneID"] < len(gene)).all()))
            ck("%s: centroids inside [0,%d)" % (tag, size),
               cell["x"].min() >= 0 and cell["y"].min() >= 0 and
               cell["x"].max() < size and cell["y"].max() < size,
               "x[%d,%d] y[%d,%d]" % (cell["x"].min(), cell["x"].max(),
                                      cell["y"].min(), cell["y"].max()))
            # attrs
            a = dict(cb["cell"].attrs)
            ar = cell["area"].astype(np.float64)
            ck("%s: cell attrs recomputed (bbox + area stats)" % tag,
               int(a["minX"][0]) == int(cell["x"].min()) and
               int(a["maxX"][0]) == int(cell["x"].max()) and
               int(a["minY"][0]) == int(cell["y"].min()) and
               int(a["maxY"][0]) == int(cell["y"].max()) and
               abs(float(a["averageArea"][0]) - ar.mean()) < 1e-3 and
               int(a["maxArea"][0]) == int(ar.max()))
            r = dict(f.attrs)
            ck("%s: root maxX/maxY/offset updated" % tag,
               int(r["maxX"][0]) == size - 1 and int(r["maxY"][0]) == size - 1 and
               int(r["offsetX"][0]) == 0 and int(r["offsetY"][0]) == 0)
            # blockIndex
            bs_ = cb["blockSize"][:]; bi = cb["blockIndex"][:]
            nbx, nby = int(bs_[2]), int(bs_[3])
            blk = (cell["y"] // int(bs_[1])).astype(np.int64) * nbx + \
                  (cell["x"] // int(bs_[0])).astype(np.int64)
            ck("%s: cells sorted by block, blockIndex is the CSR prefix" % tag,
               bool(np.all(np.diff(blk) >= 0)) and len(bi) == nbx * nby + 1 and
               int(bi[-1]) == n and
               np.array_equal(bi.astype(np.int64),
                              np.concatenate(([0], np.cumsum(np.bincount(blk, minlength=nbx * nby))))))
            # borders
            bd = cb["cellBorder"][:]
            ck("%s: cellBorder rows == cells, sentinel preserved" % tag,
               bd.shape[0] == n and bool((bd == BORDER_PAD).any()))
            v = bd[..., 0] != BORDER_PAD
            ax = bd[..., 0].astype(np.int64) + cell["x"][:, None].astype(np.int64)
            ay = bd[..., 1].astype(np.int64) + cell["y"][:, None].astype(np.int64)
            ba = dict(cb["cellBorder"].attrs)
            ck("%s: cellBorder bbox attrs match absolute vertices" % tag,
               int(ba["minX"][0]) == int(ax[v].min()) and int(ba["maxX"][0]) == int(ax[v].max()) and
               int(ba["minY"][0]) == int(ay[v].min()) and int(ba["maxY"][0]) == int(ay[v].max()))
            ck("%s: codedCellBlock dropped" % tag, "codedCellBlock" not in f)
            # exon roll-ups
            ck("%s: cellExon == sum(cellExpExon) per cell" % tag,
               np.array_equal(cb["cellExon"][:].astype(np.int64),
                              np.bincount(ci, weights=cb["cellExpExon"][:].astype(np.int64),
                                          minlength=n).astype(np.int64)))
            ck("%s: geneExon == sum(geneExpExon) per gene" % tag,
               np.array_equal(cb["geneExon"][:].astype(np.int64),
                              np.bincount(gi, weights=cb["geneExpExon"][:].astype(np.int64),
                                          minlength=len(gene)).astype(np.int64)))

    # ---------------- image registration vs the untouched source ----------------
    import tifffile, zarr
    out = {}
    for fn in ["%s_HE_regist.tif" % sn, "%s_HE_mask.tif" % sn]:
        crop = tifffile.imread(os.path.join(out_dir, "image", fn))
        with tifffile.TiffFile(os.path.join(source_dir, "image", fn)) as tf:
            st = tf.pages[0].aszarr(); z = zarr.open(st, mode="r")
            tile = np.asarray(z[y0:y0 + size, x0:x0 + size]); st.close()
        ck("image %s: pixel-exact vs source[%d:%d, %d:%d]" % (fn, y0, y0 + size, x0, x0 + size),
           np.array_equal(crop, tile), "%s" % (crop.shape,))
        out[fn] = crop
    return res, out


def registration(out_dir, source_dir, sn, x0, y0, size):
    """Prove the coordinate translation and the image crop used the same window.

    Exact set equality against the source beats any hit-rate heuristic: the
    _HE_mask is a cell/nucleus segmentation covering ~1.4% of the chip, so its
    enrichment collapses inside a mask-dense window even when alignment is perfect.
    """
    import tifffile
    res = []
    def ck(name, cond, detail=""):
        res.append({"name": name, "pass": bool(cond), "detail": str(detail)})

    def occupied(path, b="bin100"):
        with h5py.File(path, "r") as f:
            e = f["geneExp"][b]["expression"][:]
        return np.unique(np.stack([e["x"], e["y"]], 1), axis=0)

    src_gef = os.path.join(source_dir, "feature_expression", "%s.tissue.gef" % sn)
    dst_gef = os.path.join(out_dir, "feature_expression", "%s.tissue.gef" % sn)
    S = occupied(src_gef)
    S = S[(S[:, 0] >= x0) & (S[:, 0] < x0 + size) &
          (S[:, 1] >= y0) & (S[:, 1] < y0 + size)] - np.array([x0, y0])
    F = occupied(dst_gef)
    ck("registration: bin100 origins == source origins - (x0,y0)",
       np.array_equal(np.unique(S, axis=0), np.unique(F, axis=0)),
       "%d fixture vs %d source-in-window" % (len(F), len(S)))

    with h5py.File(os.path.join(source_dir, "feature_expression",
                                "%s.adjusted.cellbin.gef" % sn), "r") as f:
        c = f["cellBin/cell"][:]
    m = ((c["x"] >= x0) & (c["x"] < x0 + size) & (c["y"] >= y0) & (c["y"] < y0 + size))
    SC = np.stack([c["x"][m] - x0, c["y"][m] - y0], 1)
    with h5py.File(os.path.join(out_dir, "feature_expression",
                                "%s.adjusted.cellbin.gef" % sn), "r") as f:
        fc = f["cellBin/cell"][:]
    FC = np.stack([fc["x"], fc["y"]], 1)
    ck("registration: cell centroids == source centroids - (x0,y0)",
       np.array_equal(np.unique(SC, axis=0), np.unique(FC, axis=0)),
       "%d fixture vs %d source-in-window" % (len(FC), len(SC)))

    # the image must carry the same window: mask value under every bin origin
    fm = tifffile.imread(os.path.join(out_dir, "image", "%s_HE_mask.tif" % sn))
    sm = tifffile.imread(os.path.join(source_dir, "image", "%s_HE_mask.tif" % sn))
    fy, fx = F[:, 1], F[:, 0]
    ck("registration: mask under every bin origin matches the source",
       np.array_equal(fm[fy, fx], sm[fy + y0, fx + x0]), "%d origins" % len(F))

    # context: chip-wide, unflipped indexing is 15x enriched over random, which is
    # what fixes the row-order convention (mask[y, x], no vertical flip)
    A = occupied(src_gef)
    ay = np.clip(A[:, 1] + 50, 0, sm.shape[0] - 1)
    ax = np.clip(A[:, 0] + 50, 0, sm.shape[1] - 1)
    rng = np.random.default_rng(0)
    r = (sm[rng.integers(0, sm.shape[0], len(ax)), rng.integers(0, sm.shape[1], len(ax))] > 0).mean()
    up = (sm[ay, ax] > 0).mean()
    fl = (sm[sm.shape[0] - 1 - ay, ax] > 0).mean()
    ck("registration: source convention is mask[y, x] (no flip)",
       up > fl and up > 5 * r,
       "unflipped %.1f%%, flipped %.1f%%, random %.1f%% (%.1fx enrichment)"
       % (100 * up, 100 * fl, 100 * r, up / r))
    return res
