#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Crop one Xenium morphology_focus channel to the mini window and re-encode it.
#
#   python3 crop_morphology.py --src <ome.tif> --out <ome.tif> \
#       --xmin 2400 --xmax 3600 --ymin 3600 --ymax 4800 --pixel-size 0.2125
#
# Python, not R: Xenium ships morphology_focus as JPEG-2000 inside BigTIFF
# (compression 34712), which GDAL 3.6.4 cannot open, so terra cannot read the
# source. The output is DEFLATE, which GDAL reads natively.
#
# ACCEPTED TRADE-OFF: the shipped image is NOT byte-identical to the vendor's
# and does not exercise the JPEG-2000 read path. The source series carries one
# resolution level, so no pyramid is written.
# ---------------------------------------------------------------------------
import argparse, json, os, sys, uuid
import numpy as np
import tifffile

ap = argparse.ArgumentParser()
ap.add_argument("--src", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--xmin", type=float, required=True)
ap.add_argument("--xmax", type=float, required=True)
ap.add_argument("--ymin", type=float, required=True)
ap.add_argument("--ymax", type=float, required=True)
ap.add_argument("--pixel-size", type=float, required=True)
ap.add_argument("--geom-out", help="write the crop geometry as json here")
a = ap.parse_args()

px = a.pixel_size
x0, x1 = int(a.xmin / px), int(a.xmax / px)
y0, y1 = int(a.ymin / px), int(a.ymax / px)

with tifffile.TiffFile(a.src) as tf:
    page = tf.pages[0]
    print("[img] source %d x %d %s, compression %s"
          % (page.shape[1], page.shape[0], page.dtype, page.compression),
          flush=True)
    if y1 > page.shape[0] or x1 > page.shape[1]:
        sys.exit("window falls outside the image: need %d x %d, have %d x %d"
                 % (x1, y1, page.shape[1], page.shape[0]))
    sub = page.asarray()[y0:y1, x0:x1].copy()

print("[img] crop %d x %d %s, px [%d:%d, %d:%d]"
      % (sub.shape[1], sub.shape[0], sub.dtype, x0, x1, y0, y1), flush=True)

# tifffile stamps a fresh uuid1 into the OME-XML on every write, so identical
# pixels hashed differently. Derive it from the name to keep rebuilds stable.
ome_uuid = "urn:uuid:" + str(
    uuid.uuid5(uuid.NAMESPACE_URL, "xenium_mini/" + os.path.basename(a.out))
)

os.makedirs(os.path.dirname(a.out), exist_ok=True)
tifffile.imwrite(
    a.out, sub,
    photometric="minisblack",
    tile=(256, 256),
    compression="deflate",
    metadata={"axes": "YX", "UUID": ome_uuid,
              "PhysicalSizeX": px, "PhysicalSizeXUnit": "um",
              "PhysicalSizeY": px, "PhysicalSizeYUnit": "um"},
)
if a.geom_out:
    # the manifest must report the window that was cut, not one recomputed
    # with different rounding: floor/ceil on the same um bounds gives 5648 px
    # where the truncation above gives 5647
    with open(a.geom_out, "w") as fh:
        json.dump({"x_min": x0, "x_max": x1, "y_min": y0, "y_max": y1,
                   "width": int(sub.shape[1]), "height": int(sub.shape[0])}, fh)

print("[img] wrote %s  %s bytes" % (a.out, format(os.path.getsize(a.out), ",")),
      flush=True)
