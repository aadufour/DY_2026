#!/bin/bash
# Profiled 1D scans: for each operator of OPS, scan it with all other OPS floating (--floatOtherPOIs=1),
# all nuisances profiled, Asimov. Run inside a card dir that has model_<OPS joined by _>.root
# (createWS.py 3 --ops OPS) and the metadata file (ranges of the scanned operators; profiled ones get +-300).
# Usage: run_profiled_scans.sh OPS [POINTS=200] [METADATA=metadata_profiled.json]
# Output: higgsCombine.<OPS_>__<op>_profiled.individual.MultiDimFit.mH125.root (+ --savePOIs: other POIs saved)
OPS=$1; POINTS=${2:-200}; META=${3:-metadata_profiled.json}
CT=/grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/combine_tools
source ~/.bashrc >/dev/null 2>&1; dy_combine_morphing >/dev/null 2>&1
for op in ${OPS//,/ }; do
  t0=$(date +%s)
  nice python3 $CT/runScans.py 3 initial --doOnly $OPS --signalPOIs $op --metadata $META > initial_$op.log 2>&1
  nice python3 $CT/runScans.py 3 scan --doOnly $OPS --signalPOIs $op --metadata $META --points $POINTS --savePOIs > scan_$op.log 2>&1
  echo "[$(date +%H:%M)] $op rc=$? $(( $(date +%s) - t0 )) s"
done
echo "ALL DONE $(date)"
