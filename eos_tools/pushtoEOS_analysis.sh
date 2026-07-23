#!/usr/bin/env bash
# pushtoEOS_analysis.sh
# Stage the eft_bkg_fullsyst_v9 results and push to EOS www as RECO3D_samebins.
# Run after cmsenv.

set -euo pipefail

SCRIPT="$(dirname "$0")/pushtoEOS_parallel.py"
BASE=/grid_mnt/data__data.polcms/cms/adufour/spritz_fabian/configs/eft_bkg_fullsyst_v9
EOS_DEST=RECO3D_samebins
WORKERS=${1:-8}

STAGE="$(dirname "$0")/stage_${EOS_DEST}"

echo "Staging into ${STAGE} ..."
mkdir -p "$STAGE"

cp "$BASE"/eft_summary_multivar* "$STAGE"/

for obs in mll rapll_abs costhetastar triple_diff; do
  for sub in scans nuis_combine nuis_lin_quad; do
    mkdir -p "$STAGE/$obs/$sub"
    cp -r "$BASE/datacards/inc_mm/$obs/$sub/." "$STAGE/$obs/$sub/"
  done
done

echo "Pushing to EOS as ${EOS_DEST} (workers=${WORKERS}) ..."
python3 "$SCRIPT" -j "$WORKERS" "$STAGE" "$EOS_DEST"

echo "Cleaning up staging dir ..."
rm -rf "$STAGE"
