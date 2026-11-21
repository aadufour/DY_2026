#!/bin/bash
# 1D rescans with per-operator ranges sized to the 95% interval (adapt_scan_ranges.py), in <vardir>/rescan/,
# leaving the original scans untouched. Two passes: ranges from the original scans, then refined from pass 1.
# Then summaries (--logscale), A vs B comparison, multi-variable summaries and the width table, all from rescan/.
#
# Usage: rescan_1d.sh <config_dir_A> <label_A> <config_dir_B> <label_B> <outdir for the comparison>
# e.g.   rescan_1d.sh .../propcorr_v3_pdfcut v3 .../same_binning_N80 N80 .../same_binning_N80/compare_ref_rescan
A=$1; LA=$2; B=$3; LB=$4; OUT=$5
CT=/grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/combine_tools
VARS="mll rapll_abs costhetastar triple_diff"
log() { echo "[$(date '+%m-%d %H:%M')] $*"; }
source ~/.bashrc >/dev/null 2>&1
dy_combine_morphing >/dev/null 2>&1
for D in $A $B; do
  for v in $VARS; do
    V=$D/datacards/inc_mm/$v
    for pass in 1 2; do
      src=$V; [ $pass -eq 2 ] && src=$V/rescan
      python3 $CT/adapt_scan_ranges.py $V --src $src > $V/rescan_ranges_pass$pass.log 2>&1 || { log "ranges failed $V"; exit 1; }
      cd $V/rescan
      runScans.py 1 scan --points 100 --metadata metadata.json > scan_pass$pass.log 2>&1
      runScans.py 1 scan --stat --points 100 --metadata metadata_stat.json > scan_stat_pass$pass.log 2>&1
    done
    log "$(basename $D) $v: $(ls $V/rescan/higgsCombine.*.individual.MultiDimFit.mH125.root | grep -vc _stat) + $(ls $V/rescan/higgsCombine.*_stat.individual.MultiDimFit.mH125.root | wc -l) stat scans"
  done
done
mkdir -p $OUT
for v in $VARS; do
  for D in $A $B; do (cd $D/datacards/inc_mm/$v/rescan && python3 $CT/makeSummary.py --indir . --logscale > summary.log 2>&1); done
  python3 $CT/makeSummaryCompare.py --indir1 $A/datacards/inc_mm/$v/rescan --indir2 $B/datacards/inc_mm/$v/rescan --label1 "$LA" --label2 "$LB" \
      --logscale --verbose -o $OUT/eft_summary_compare_$v > $OUT/compare_$v.log 2>&1
done
for D in $A $B; do
  L=$LA; [ $D = $B ] && L=$LB
  (cd $OUT && python3 $CT/makeSummaryMultiVar.py --mll $D/datacards/inc_mm/mll/rescan --rapll $D/datacards/inc_mm/rapll_abs/rescan \
      --costhetastar $D/datacards/inc_mm/costhetastar/rescan --triple_diff $D/datacards/inc_mm/triple_diff/rescan --logscale > multivar_$L.log 2>&1 \
   && mv eft_summary_multivar.png eft_summary_multivar_$L.png && mv eft_summary_multivar.pdf eft_summary_multivar_$L.pdf)
done
(cd $OUT && python3 $CT/compare_1d_widths.py $A $B --labels $LA $LB --subdir rescan -o widths_${LB}_vs_${LA}.md > widths.log 2>&1)
log "ALL DONE"
