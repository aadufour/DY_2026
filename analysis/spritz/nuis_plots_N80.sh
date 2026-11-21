#!/bin/bash
# N80: nuis_combine (drawNuisances.py) + nuis_eft (drawNuisancesEFT.py) per variable, push to RECO3D_same_binning_N80/<var>/
# (as propcorr_v3_pdfcut/nuis_plots.sh; the 1D scans are already on EOS from the rescans, eos_tools/push_N80_to_eos.sh)
C=/grid_mnt/data__data.polcms/cms/adufour/spritz_fabian/configs/same_binning_N80
W=/grid_mnt/data__data.polcms/cms/adufour; CT=$W/DY_2026/analysis/combine_tools
export X509_USER_PROXY=$W/proxy.pem
source ~/.bashrc >/dev/null 2>&1; dy_combine_morphing >/dev/null 2>&1
for v in mll rapll_abs costhetastar triple_diff; do
  cd $C/datacards/inc_mm/$v || continue
  [ -d nuis_combine ] || nice python3 $CT/drawNuisances.py --datacard datacard.txt --output nuis_combine --ncores 8 > nuis_combine.log 2>&1
  [ -d nuis_eft ] || nice python3 $CT/drawNuisancesEFT.py --shapes shapes.root --outdir nuis_eft --ncores 8 > nuis_eft.log 2>&1
  echo "[$(date +%H:%M)] $v: nuis_combine $(ls nuis_combine | wc -l), nuis_eft $(ls nuis_eft | wc -l)"
  for d in nuis_combine nuis_eft; do
    python3 $W/DY_2026/eos_tools/pushtoEOS_parallel.py -j 8 $PWD/$d RECO3D_same_binning_N80/$v/$d 2>&1 | grep -E "Done|ERROR: some|Web URL"
  done
done
echo "ALL DONE $(date)"
