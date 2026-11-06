#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Zip a built mini fixture reproducibly and print its sha256.
#
#   ./package_mini.sh <build_dir> <output.zip>
#
# <build_dir> is the directory that CONTAINS xenium_mini_breast_ffpe/, i.e.
# whatever XENIUM_MINI_OUTDIR's parent was.
#
# Not plain `zip -r`: a zip records each entry's mtime and stores entries in
# the order given, so the same files zip to different bytes. TZ is fixed (zip
# writes local time with no zone), every mtime is pinned, and the entry list
# is sorted. Build determinism is separate and lives in make_xenium_mini.R;
# set SOURCE_DATE_EPOCH when building to compare two builds.
# ---------------------------------------------------------------------------
set -euo pipefail

SRC=${1:?usage: package_mini.sh <build_dir> <output.zip>}
OUT=${2:?usage: package_mini.sh <build_dir> <output.zip>}
NAME=xenium_mini_breast_ffpe
STAMP=${MINI_ZIP_MTIME:-2026-10-27 00:00:00}

[ -d "$SRC/$NAME" ] || { echo "no $NAME/ under $SRC" >&2; exit 1; }
OUT=$(readlink -f "$OUT")

export TZ=UTC
stage=$(mktemp -d)
trap 'rm -rf "$stage"' EXIT

cp -r "$SRC/$NAME" "$stage/"
find "$stage" -exec touch -d "$STAMP" {} +
rm -f "$OUT"
( cd "$stage" && find "$NAME" -print | LC_ALL=C sort | zip -X -q "$OUT" -@ )

printf '%s\n' "$OUT"
printf 'size   : %s bytes\n' "$(stat -c%s "$OUT")"
printf 'sha256 : %s\n' "$(sha256sum "$OUT" | cut -d' ' -f1)"
