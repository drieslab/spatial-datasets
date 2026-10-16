"""
stereoseq_crop.py -- spatially crop a SAW 8.2.0 Stereo-seq ``outs/`` bundle.

Driven by make_stereoseq_mini_mouse_eyeball.R, which provisions the python
dependencies (h5py, numpy, tifffile, imagecodecs, zarr, anndata) through
reticulate/uv and then calls ``main()``.

Everything is defined by ONE rectangular window in chip units (1 unit = 500 nm):

    keep bins whose ORIGIN is inside the window
    keep cells whose CENTROID is inside the window

and every kept coordinate is translated so the fixture's origin is (0, 0).
The window origin and size are multiples of 600 = lcm(1,5,10,20,50,100,150,200)
so that, after translation, every bin level's origins remain multiples of its
own bin size.

GEF index structures are rebuilt, not copied: ``geneExp/<bin>/gene`` offsets,
``cellBin/cell`` offsets, both cell-major and gene-major cell-bin indexes,
``wholeExp`` grids, ``blockIndex``, and every dataset attribute that records an
extent, a min/max/median/average or a count.
"""

import json
import os
import shutil

import h5py
import numpy as np

GZIP = 6
BORDER_PAD = 32767          # cellBorder padding sentinel
BLOCK = 256                 # cellBin block size in chip units


# --------------------------------------------------------------------------- helpers
def _ds(grp, name, data, gzip_level=GZIP):
    """Create a dataset, compressing whenever HDF5 allows chunking."""
    if getattr(data, "ndim", 0) >= 1 and data.size > 0:
        return grp.create_dataset(name, data=data, compression="gzip",
                                  compression_opts=gzip_level)
    return grp.create_dataset(name, data=data)


def _set_attrs(obj, src_attrs, updates=None):
    """Copy attributes preserving each one's original dtype, applying updates."""
    updates = updates or {}
    for k, v in src_attrs.items():
        if k in updates:
            nv = updates[k]
            if isinstance(v, np.ndarray):
                v = np.array(np.atleast_1d(nv), dtype=v.dtype)
            elif isinstance(v, (np.generic,)):
                v = np.array(nv, dtype=v.dtype)[()]
            else:
                v = type(v)(nv)
        obj.attrs[k] = v


def _copy_verbatim(src, dst, name):
    _ds(dst, name, src[name][:])


def _cumsum0(counts):
    """Exclusive prefix sum as int64."""
    c = np.asarray(counts, dtype=np.int64)
    return np.concatenate(([0], np.cumsum(c)[:-1]))


# --------------------------------------------------------------------------- bin GEF
def crop_bin_gef(src_path, dst_path, x0, y0, size):
    """Crop a square-bin GEF (``*.tissue.gef``) to [x0,x0+size) x [y0,y0+size)."""
    report = {"bins": {}}
    with h5py.File(src_path, "r") as f, h5py.File(dst_path, "w") as o:
        src_root = dict(f.attrs)
        og = o.create_group("geneExp")
        ow = o.create_group("wholeExp")

        occ_bin1_total = occ_bin1_kept = None

        for b in sorted(f["geneExp"], key=lambda s: int(s[3:])):
            bs = int(b[3:])
            assert x0 % bs == 0 and y0 % bs == 0 and size % bs == 0, \
                f"window not aligned to {b}"
            sb = f["geneExp"][b]
            expr = sb["expression"][:]
            gene = sb["gene"][:]

            cnt = gene["count"].astype(np.int64)
            assert np.array_equal(gene["offset"].astype(np.int64), _cumsum0(cnt)), \
                f"{b}: gene table is not gene-major contiguous"
            assert cnt.sum() == len(expr)

            m = ((expr["x"] >= x0) & (expr["x"] < x0 + size) &
                 (expr["y"] >= y0) & (expr["y"] < y0 + size))
            gidx = np.repeat(np.arange(len(gene), dtype=np.int64), cnt)

            new_expr = expr[m].copy()
            new_expr["x"] -= x0
            new_expr["y"] -= y0

            new_cnt = np.bincount(gidx[m], minlength=len(gene))
            new_gene = gene.copy()
            new_gene["count"] = new_cnt.astype(gene.dtype["count"])
            new_gene["offset"] = _cumsum0(new_cnt).astype(gene.dtype["offset"])

            g = og.create_group(b)
            d = _ds(g, "expression", new_expr)
            _set_attrs(d, dict(sb["expression"].attrs), {
                "minX": 0, "minY": 0,
                "maxX": int(new_expr["x"].max()) if len(new_expr) else 0,
                "maxY": int(new_expr["y"].max()) if len(new_expr) else 0,
                "maxExp": int(new_expr["count"].max()) if len(new_expr) else 0,
            })
            dg = _ds(g, "gene", new_gene)
            _set_attrs(dg, dict(sb["gene"].attrs))
            if "exon" in sb:
                _ds(g, "exon", sb["exon"][:][m])

            # wholeExp: dense (lenX/bs, lenY/bs) grid of (MIDcount, genecount).
            # One expression row == one (gene, bin) pair, so row counts per bin
            # are exactly the number of distinct genes in that bin.
            sw = f["wholeExp"][b]
            nx = ny = size // bs
            flat = (new_expr["x"] // bs).astype(np.int64) * ny + \
                   (new_expr["y"] // bs).astype(np.int64)
            mid = np.bincount(flat, weights=new_expr["count"].astype(np.int64),
                              minlength=nx * ny).reshape(nx, ny)
            gct = np.bincount(flat, minlength=nx * ny).reshape(nx, ny)
            grid = np.zeros((nx, ny), dtype=sw.dtype)
            grid["MIDcount"] = mid.astype(sw.dtype["MIDcount"])
            grid["genecount"] = gct.astype(sw.dtype["genecount"])
            dw = _ds(ow, b, grid)
            _set_attrs(dw, dict(sw.attrs), {
                "lenX": size, "lenY": size, "minX": 0, "minY": 0,
                "maxGene": int(gct.max()), "maxMID": int(mid.max()),
                "number": int((gct > 0).sum()),
            })

            if bs == 1:
                occ_bin1_total = int(dict(sw.attrs)["number"][0])
                occ_bin1_kept = int((gct > 0).sum())

            report["bins"][b] = {
                "rows_source": int(len(expr)), "rows_kept": int(len(new_expr)),
                "genes": int(len(gene)),
                "genes_nonzero": int((new_cnt > 0).sum()),
                "occupied_bins": int((gct > 0).sum()),
                "mid_total": int(mid.sum()),
            }

        # gef_area is an area measure whose exact SAW definition is not
        # recoverable from the file; areas are additive, so scale it by the
        # fraction of occupied bin1 positions retained.
        frac = occ_bin1_kept / occ_bin1_total
        _set_attrs(o, src_root, {"gef_area": float(src_root["gef_area"][0]) * frac})
        report["gef_area_scale"] = frac
        report["dropped"] = ["stat/", "contour/"]
    return report


# ----------------------------------------------------------------------- cellbin GEF
def crop_cellbin_gef(src_path, dst_path, x0, y0, size):
    """Crop a cell-bin GEF (``*.cellbin.gef`` / ``*.adjusted.cellbin.gef``)."""
    with h5py.File(src_path, "r") as f, h5py.File(dst_path, "w") as o:
        cb = f["cellBin"]
        cell = cb["cell"][:]
        gene = cb["gene"][:]
        cellExp = cb["cellExp"][:]
        geneExp = cb["geneExp"][:]
        border = cb["cellBorder"][:]
        cellExpExon = cb["cellExpExon"][:]
        geneExpExon = cb["geneExpExon"][:]

        assert np.array_equal(cell["id"], np.arange(len(cell), dtype=cell.dtype["id"])), \
            "cell$id is expected to be the row index"

        keep = ((cell["x"] >= x0) & (cell["x"] < x0 + size) &
                (cell["y"] >= y0) & (cell["y"] < y0 + size))
        idx = np.flatnonzero(keep)

        # Cells are stored sorted by 256-unit block (row-major). The block grid
        # changes under translation, so re-sort into the new block order.
        nbx = -(-size // BLOCK)
        nby = -(-size // BLOCK)
        bx = ((cell["x"][idx] - x0) // BLOCK).astype(np.int64)
        by = ((cell["y"][idx] - y0) // BLOCK).astype(np.int64)
        blk = by * nbx + bx
        idx = idx[np.argsort(blk, kind="stable")]
        blk = np.sort(blk, kind="stable")
        n = len(idx)

        # ---- cell-major side -------------------------------------------------
        gcount = cell["geneCount"][idx].astype(np.int64)
        coff = cell["offset"][idx].astype(np.int64)
        take = (np.concatenate([np.arange(a, a + b) for a, b in zip(coff, gcount)])
                if n else np.zeros(0, dtype=np.int64))
        new_cellExp = cellExp[take]
        new_cellExpExon = cellExpExon[take]

        new_cell = cell[idx].copy()
        new_cell["id"] = np.arange(n, dtype=cell.dtype["id"])
        new_cell["offset"] = _cumsum0(gcount).astype(cell.dtype["offset"])
        new_cell["x"] -= x0
        new_cell["y"] -= y0

        # ---- gene-major side, rebuilt from the same triplets -----------------
        cell_of_row = np.repeat(np.arange(n, dtype=np.int64), gcount)
        gene_of_row = new_cellExp["geneID"].astype(np.int64)
        cnt_of_row = new_cellExp["count"]
        order = np.lexsort((cell_of_row, gene_of_row))   # gene-major, cellID ascending

        new_geneExp = np.empty(len(order), dtype=geneExp.dtype)
        new_geneExp["cellID"] = cell_of_row[order].astype(geneExp.dtype["cellID"])
        new_geneExp["count"] = cnt_of_row[order]
        new_geneExpExon = new_cellExpExon[order]

        cellCount = np.bincount(gene_of_row, minlength=len(gene))
        expCount = np.bincount(gene_of_row, weights=cnt_of_row.astype(np.int64),
                               minlength=len(gene))
        maxMID = np.zeros(len(gene), dtype=np.int64)
        if len(gene_of_row):
            np.maximum.at(maxMID, gene_of_row, cnt_of_row.astype(np.int64))

        new_gene = gene.copy()
        new_gene["cellCount"] = cellCount.astype(gene.dtype["cellCount"])
        new_gene["expCount"] = expCount.astype(gene.dtype["expCount"])
        new_gene["maxMIDcount"] = maxMID.astype(gene.dtype["maxMIDcount"])
        # zero-count genes carry offset 0 in the source; non-zero genes are contiguous
        off = np.zeros(len(gene), dtype=np.int64)
        nz = cellCount > 0
        off[nz] = _cumsum0(cellCount[nz])
        new_gene["offset"] = off.astype(gene.dtype["offset"])

        # ---- exon roll-ups ---------------------------------------------------
        # cells are kept whole, so cellExon is unchanged for kept cells; geneExon
        # must be recomputed (genes lose cells). NB the source geneExon overflows
        # a uint16 accumulator for 5 genes; the recomputed values are correct.
        new_cellExon = cb["cellExon"][:][idx]
        new_geneExon = np.bincount(gene_of_row,
                                   weights=new_cellExpExon.astype(np.int64),
                                   minlength=len(gene)).astype(cb["geneExon"].dtype)

        # ---- borders (cell-relative, sentinel-padded: never translated) ------
        new_border = border[idx]
        valid = new_border[..., 0] != BORDER_PAD
        if valid.any():
            ax = new_border[..., 0].astype(np.int64) + new_cell["x"][:, None].astype(np.int64)
            ay = new_border[..., 1].astype(np.int64) + new_cell["y"][:, None].astype(np.int64)
            bminx, bmaxx = int(ax[valid].min()), int(ax[valid].max())
            bminy, bmaxy = int(ay[valid].min()), int(ay[valid].max())
        else:
            bminx = bmaxx = bminy = bmaxy = 0

        # ---- write -----------------------------------------------------------
        _set_attrs(o, dict(f.attrs), {
            "maxX": size - 1, "maxY": size - 1, "offsetX": 0, "offsetY": 0,
        })
        g = o.create_group("cellBin")
        _set_attrs(g, dict(cb.attrs))

        dc = _ds(g, "cell", new_cell)
        area = new_cell["area"].astype(np.float64)
        dnb = new_cell["dnbCount"].astype(np.float64)
        exp = new_cell["expCount"].astype(np.float64)
        gcf = new_cell["geneCount"].astype(np.float64)
        _set_attrs(dc, dict(cb["cell"].attrs), {
            "minX": int(new_cell["x"].min()), "maxX": int(new_cell["x"].max()),
            "minY": int(new_cell["y"].min()), "maxY": int(new_cell["y"].max()),
            "minArea": int(area.min()), "maxArea": int(area.max()),
            "medianArea": float(np.median(area)), "averageArea": float(area.mean()),
            "minDnbCount": int(dnb.min()), "maxDnbCount": int(dnb.max()),
            "medianDnbCount": float(np.median(dnb)), "averageDnbCount": float(dnb.mean()),
            "minExpCount": int(exp.min()), "maxExpCount": int(exp.max()),
            "medianExpCount": float(np.median(exp)), "averageExpCount": float(exp.mean()),
            "minGeneCount": int(gcf.min()), "maxGeneCount": int(gcf.max()),
            "medianGeneCount": float(np.median(gcf)), "averageGeneCount": float(gcf.mean()),
        })

        dg = _ds(g, "gene", new_gene)
        _set_attrs(dg, dict(cb["gene"].attrs), {
            "minCellCount": int(cellCount[nz].min()) if nz.any() else 0,
            "maxCellCount": int(cellCount.max()),
            "minExpCount": int(expCount[nz].min()) if nz.any() else 0,
            "maxExpCount": int(expCount.max()),
        })

        dge = _ds(g, "geneExp", new_geneExp)
        _set_attrs(dge, dict(cb["geneExp"].attrs),
                   {"maxCount": int(new_geneExp["count"].max()) if len(new_geneExp) else 0})
        dce = _ds(g, "cellExp", new_cellExp)
        _set_attrs(dce, dict(cb["cellExp"].attrs),
                   {"maxCount": int(new_cellExp["count"].max()) if len(new_cellExp) else 0})

        db = _ds(g, "cellBorder", new_border)
        # ax/ay were built from the ALREADY-translated centroids, so the bbox is
        # already in fixture coordinates -- do not subtract the origin again.
        _set_attrs(db, dict(cb["cellBorder"].attrs), {
            "minX": bminx, "maxX": bmaxx, "minY": bminy, "maxY": bmaxy,
        })

        _ds(g, "cellExon", new_cellExon)
        _ds(g, "cellExpExon", new_cellExpExon)
        _ds(g, "geneExon", new_geneExon)
        _ds(g, "geneExpExon", new_geneExpExon)

        blockIndex = np.concatenate(
            ([0], np.cumsum(np.bincount(blk, minlength=nbx * nby)))
        ).astype(cb["blockIndex"].dtype)
        _ds(g, "blockSize", np.array([BLOCK, BLOCK, nbx, nby],
                                     dtype=cb["blockSize"].dtype))
        _ds(g, "blockIndex", blockIndex)
        if "cellTypeList" in cb:
            _copy_verbatim(cb, g, "cellTypeList")

    return {
        "cells_source": int(len(cell)), "cells_kept": int(n),
        "genes": int(len(gene)), "genes_nonzero": int(nz.sum()),
        "nnz_source": int(len(cellExp)), "nnz_kept": int(len(new_cellExp)),
        "blocks": [int(nbx), int(nby)],
        "dropped": ["codedCellBlock/"],
    }


# ------------------------------------------------------------------------------ images
def crop_images(src_dir, dst_dir, x0, y0, size, keep=("_HE_regist.tif", "_HE_mask.tif")):
    import tifffile
    import zarr

    os.makedirs(dst_dir, exist_ok=True)
    out = {}
    for fn in sorted(os.listdir(src_dir)):
        if not any(fn.endswith(k) for k in keep):
            continue
        src = os.path.join(src_dir, fn)
        with tifffile.TiffFile(src) as tf:
            page = tf.pages[0]
            store = page.aszarr()
            z = zarr.open(store, mode="r")
            sub = np.asarray(z[y0:y0 + size, x0:x0 + size])
            store.close()
            photometric = "rgb" if sub.ndim == 3 else "minisblack"
        dst = os.path.join(dst_dir, fn)
        tifffile.imwrite(dst, sub, photometric=photometric, compression="deflate")
        out[fn] = {"shape": list(sub.shape), "dtype": str(sub.dtype),
                   "bytes": os.path.getsize(dst)}
    return out


# ------------------------------------------------------------------------------- h5ad
def crop_h5ad(src_path, dst_path, x0, y0, size):
    import anndata as ad

    a = ad.read_h5ad(src_path)
    if "x" in a.obs and "y" in a.obs:
        X = a.obs["x"].to_numpy().astype(np.int64)
        Y = a.obs["y"].to_numpy().astype(np.int64)
    else:
        sp = np.asarray(a.obsm["spatial"])
        X, Y = sp[:, 0].astype(np.int64), sp[:, 1].astype(np.int64)
    m = (X >= x0) & (X < x0 + size) & (Y >= y0) & (Y < y0 + size)
    b = a[m].copy()

    if "x" in b.obs and "y" in b.obs:
        b.obs["x"] = b.obs["x"].to_numpy().astype(np.int64) - x0
        b.obs["y"] = b.obs["y"].to_numpy().astype(np.int64) - y0
    if "spatial" in b.obsm:
        sp = np.asarray(b.obsm["spatial"]).copy()
        sp[:, 0] -= x0
        sp[:, 1] -= y0
        b.obsm["spatial"] = sp

    # per-observation fields survive subsetting; run-level aggregates do not
    counts = b.raw.X if b.raw is not None else b.X
    n_cells = np.asarray((counts > 0).sum(axis=0)).ravel()
    n_counts = np.asarray(counts.sum(axis=0)).ravel()
    for frame in (b.var, (b.raw.var if b.raw is not None else None)):
        if frame is None:
            continue
        if "n_cells" in frame:
            frame["n_cells"] = n_cells.astype(frame["n_cells"].dtype, copy=False)
        if "n_counts" in frame:
            frame["n_counts"] = n_counts.astype(frame["n_counts"].dtype, copy=False)

    # run-level results that do not survive a crop, plus the opaque coded cell
    # block (the same structure dropped from the cell-bin GEF -- anndata cannot
    # even write it back, and it describes the full-section block grid)
    dropped = [k for k in ("rank_genes_groups", "gene_exp_leiden",
                           "pca_variance_ratio", "coded_cell_block")
               if k in b.uns]
    for k in dropped:
        del b.uns[k]

    # defensively drop anything else anndata cannot round-trip (raw bytes)
    def _prune(d, prefix=""):
        for k in list(d.keys()):
            v = d[k]
            if isinstance(v, bytes):
                del d[k]
                dropped.append(prefix + k)
            elif hasattr(v, "keys"):
                _prune(v, prefix + k + "/")
    _prune(b.uns)

    b.write_h5ad(dst_path, compression="gzip")
    return {"obs_source": int(a.n_obs), "obs_kept": int(b.n_obs),
            "vars": int(b.n_vars), "uns_dropped": dropped,
            "bytes": os.path.getsize(dst_path)}


# ------------------------------------------------------------------------------- main
def main(source_dir, out_dir, x0, y0, size,
         bin_gef, cellbin_gefs, h5ads):
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    fe_out = os.path.join(out_dir, "feature_expression")
    an_out = os.path.join(out_dir, "analysis")
    os.makedirs(fe_out)
    os.makedirs(an_out)

    rep = {"bin_gef": {}, "cellbin_gef": {}, "images": {}, "h5ad": {}}

    src = os.path.join(source_dir, "feature_expression", bin_gef)
    rep["bin_gef"][bin_gef] = crop_bin_gef(src, os.path.join(fe_out, bin_gef), x0, y0, size)

    for fn in cellbin_gefs:
        src = os.path.join(source_dir, "feature_expression", fn)
        rep["cellbin_gef"][fn] = crop_cellbin_gef(src, os.path.join(fe_out, fn), x0, y0, size)

    rep["images"] = crop_images(os.path.join(source_dir, "image"),
                                os.path.join(out_dir, "image"), x0, y0, size)

    for fn in h5ads:
        rep["h5ad"][fn] = crop_h5ad(os.path.join(source_dir, "analysis", fn),
                                    os.path.join(an_out, fn), x0, y0, size)
    return rep
