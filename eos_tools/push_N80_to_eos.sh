#!/bin/bash
# Push the same-binning study (N80) to EOS www as RECO3D_same_binning_N80 (notes/binning_study.md):
#   top level : eft_summary_multivar_{N80,v3} (1D rescans, ranges sized to the intervals), widths_N80_vs_v3.txt
#   compare/  : eft_summary_compare_<obs> (v3 vs N80, rescans)
#   <obs>/scan/ : N80 1D scans (rescan, stat+syst with stat-only overlaid)
#   lhe/      : LHE-level binning checks (lhe_binning_check.py --set N80, lhe_binning_saturation.py)
# Needs a valid proxy/ticket for EOS; sources cmsenv itself.
# Usage: bash push_N80_to_eos.sh [workers]
set -e

C=/grid_mnt/data__data.polcms/cms/adufour/spritz_fabian/configs
N80=$C/same_binning_N80
CMP=$N80/compare_ref_rescan
LHE=/grid_mnt/data__data.polcms/cms/adufour/binning_study
EOS_DEST=RECO3D_same_binning_N80
WORKERS=${1:-8}
PUSH="$(cd "$(dirname "$0")" && pwd)/pushtoEOS_parallel.py"
STAGE=/grid_mnt/data__data.polcms/cms/adufour/stage_${EOS_DEST}

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /grid_mnt/data__data.polcms/cms/adufour/CMSSW_spritz/CMSSW_14_1_0_pre4/src && eval "$(scramv1 runtime -sh)" && cd - >/dev/null

[ -e "$STAGE" ] && { echo "staging dir $STAGE exists, remove it first"; exit 1; }
mkdir -p "$STAGE/compare" "$STAGE/lhe"

cp "$CMP"/eft_summary_multivar_*.png "$CMP"/eft_summary_multivar_*.pdf "$STAGE"/
cp "$CMP"/widths_N80_vs_v3.md "$STAGE"/widths_N80_vs_v3.txt
cp "$CMP"/eft_summary_compare_*.png "$CMP"/eft_summary_compare_*.pdf "$STAGE/compare/"
for obs in mll rapll_abs costhetastar triple_diff; do
  mkdir -p "$STAGE/$obs/scan"
  cp "$N80"/datacards/inc_mm/$obs/rescan/scan_*.png "$N80"/datacards/inc_mm/$obs/rescan/scan_*.pdf "$STAGE/$obs/scan/"
  cp "$N80"/datacards/inc_mm/$obs/rescan/eft_summary_two_panel.* "$STAGE/$obs/"
done
cp "$LHE"/lhe_N80/*.png "$LHE"/lhe_saturation/saturation.png "$STAGE/lhe/"
cp "$LHE"/lhe_N80/sigma_table.md "$STAGE/lhe/sigma_table_N80.txt"
cp "$LHE"/lhe_saturation/saturation.md "$STAGE/lhe/saturation.txt"

du -sh "$STAGE"
python3 "$PUSH" -j "$WORKERS" "$STAGE" "$EOS_DEST"

rm -rf "$STAGE"   # staging copies only
