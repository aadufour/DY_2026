source ~/.bashrc
D=/grid_mnt/data__data.polcms/cms/adufour/spritz_fabian/configs/propcorr_v3_pdfcut
V1=/grid_mnt/data__data.polcms/cms/adufour/spritz_fabian/configs/propcorr_v1_minnlo/datacards
CT=/grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/combine_tools
VARS="mll rapll_abs costhetastar triple_diff"
( dy_combine_morphing
  for v in $VARS; do cd $D/datacards/inc_mm/$v; echo "== $v $(date)"
    createCombineJson.py --datacard datacard.txt --binname w1_ --output jsonComb.json > combjson.log 2>&1; echo "json: $(cat jsonComb.json | head -c 80)"
    python3 $CT/createWS.py 1 > createWS.log 2>&1; echo "ws rc=$? $(ls model_*.root | wc -l) models"
    runScans.py 1 initial > initial.log 2>&1; echo "initial rc=$?"
    runScans.py 1 scan --points 100 > scan.log 2>&1; echo "scan rc=$? $(ls higgsCombine.*.individual.MultiDimFit.mH125.root 2>/dev/null | wc -l) scans"
    runScans.py 1 initial --stat > initial_stat.log 2>&1
    runScans.py 1 scan --stat --points 100 > scan_stat.log 2>&1; echo "stat rc=$?"
    runPlots_compare.py 1 --label "Stat + Syst" --compare-stat > plots.log 2>&1; echo "plots rc=$?"
  done )
( dy_analysis
  for v in $VARS; do cd $D/datacards/inc_mm/$v
    python3 $CT/makeSummary.py --indir . > summary.log 2>&1; echo "summary $v rc=$?"
    python3 $CT/makeSummaryCompare.py --indir1 $V1/inc_mm/$v --indir2 . --label1 "MiNNLO v1 (no PDF cut)" --label2 "MiNNLO v3 (PDF cut)" --verbose -o $D/compare_v1minnlo/eft_summary_compare_$v > $D/compare_v1minnlo/compare_$v.log 2>&1; echo "cmp $v rc=$?"
  done
  cd $D/summary_plots; python3 $CT/makeSummaryMultiVar.py --mll ../datacards/inc_mm/mll --rapll ../datacards/inc_mm/rapll_abs --costhetastar ../datacards/inc_mm/costhetastar --triple_diff ../datacards/inc_mm/triple_diff > multivar.log 2>&1; echo "multivar rc=$?" )
echo ALL1DDONE
