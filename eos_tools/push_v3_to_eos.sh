#!/bin/bash
# Push propcorr v3 (PDF soft cut) results to EOS www as RECO3D_propcorr_v3:
#   top level : multivar summaries
#   <obs>/    : eft_summary_two_panel, scan/ (1D scans), dist/ (distributions)
#   (2D scans are pushed by run_2d_perpair.sh into 2D_mll / 2D_triple_diff / 2D_mll_vs_triple_diff)
# Needs a valid proxy/ticket for EOS; sources cmsenv itself.
# Usage: bash push_v3_to_eos.sh [workers]
set -e

C=/grid_mnt/data__data.polcms/cms/adufour/spritz_fabian/configs
V3=$C/propcorr_v3_pdfcut
EOS_DEST=RECO3D_propcorr_v3
WORKERS=${1:-8}
PUSH="$(cd "$(dirname "$0")" && pwd)/pushtoEOS_parallel.py"
STAGE=/grid_mnt/data__data.polcms/cms/adufour/stage_${EOS_DEST}

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /grid_mnt/data__data.polcms/cms/adufour/CMSSW_spritz/CMSSW_14_1_0_pre4/src && eval "$(scramv1 runtime -sh)" && cd - >/dev/null

[ -e "$STAGE" ] && { echo "staging dir $STAGE exists, remove it first"; exit 1; }
mkdir -p "$STAGE"

cp "$V3"/summary_plots/*.png "$V3"/summary_plots/*.pdf "$STAGE"/
for obs in mll rapll_abs costhetastar triple_diff; do
  d=$V3/datacards/inc_mm/$obs
  mkdir -p "$STAGE/$obs/scan" "$STAGE/$obs/dist"
  cp "$d"/eft_summary_two_panel.* "$STAGE/$obs/" 2>/dev/null || true
  cp "$d"/scan/scan_*.png "$d"/scan/scan_*.pdf "$STAGE/$obs/scan/"
  cp -r "$V3/plots/inc_mm/$obs/." "$STAGE/$obs/dist/"
done

du -sh "$STAGE"
python3 "$PUSH" -j "$WORKERS" --overwrite "$STAGE" "$EOS_DEST"

rm -rf "$STAGE"
