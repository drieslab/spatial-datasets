#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Apply a minimal patch to a vendor experiment.xenium.
#
#   python3 patch_experiment.py --src <in> --out <out> --patch <patch.json>
#
# <patch.json> is {"set": {...}, "delete": [...]}. Keys in "set" must already
# exist in the source: this edits the vendor schema, it never extends it.
# Crop provenance belongs in crop_manifest.json, not here.
#
# Python, not R: jsonlite's digits = NA is 15 significant digits, so a
# read/write round trip silently altered nine vendor floats nothing had asked
# to change. Python's json emits the shortest text that round-trips, which is
# what the vendor wrote, so unnamed fields come out byte-identical.
# ---------------------------------------------------------------------------
import argparse, json, sys

ap = argparse.ArgumentParser()
ap.add_argument("--src", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--patch", required=True)
a = ap.parse_args()

with open(a.src) as fh:
    xe = json.load(fh)
with open(a.patch) as fh:
    patch = json.load(fh)

dele = patch.get("delete", [])
if not isinstance(dele, list):
    # a bare string would be iterated character by character, deleting nothing
    sys.exit("patch 'delete' must be a list, got %s" % type(dele).__name__)

unknown = [k for k in patch.get("set", {}) if k not in xe]
if unknown:
    sys.exit("refusing to add non-vendor keys to experiment.xenium: %s"
             % ", ".join(sorted(unknown)))

for k, v in patch.get("set", {}).items():
    xe[k] = v
for k in dele:
    xe.pop(k, None)

with open(a.out, "w") as fh:
    json.dump(xe, fh, indent=2)
    fh.write("\n")

print("[xenium] patched %d field(s), dropped %d"
      % (len(patch.get("set", {})), len(dele)), flush=True)
