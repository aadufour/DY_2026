#!/bin/bash
# Unattended spritz chain on LLR for a config dir whose condor jobs are already submitted (run in tmux):
#   1. wait for the condor cluster; retry once the jobs that did not run or logged "Error for chunk"
#   2. merge (spritz merge, 4 cpus)          -> condor/results_merged_new.pkl
#   3. MiNNLO postproc (pp_minnlo.py, 4 cpus) -> histos.root
#   4. cards (spritz-cards-eft)              -> datacards/inc_mm/<var>/
#   5. 1D fits per variable (as propcorr_v3_pdfcut run_1d.sh: initial + 100-point scans, stat-only, plots)
#   6. summaries (--logscale) + comparison with a reference production
#
# Usage: run_chain_llr.sh <config_dir> <cluster> <reference_config_dir> <label> <reference_label>
# e.g.   run_chain_llr.sh .../configs/same_binning_N80 823600 .../configs/propcorr_v3_pdfcut "N80" "v3"
D=$1; CLUSTER=$2; REF=$3; LAB=$4; REFLAB=$5
RD=/grid_mnt/data__data.polcms/cms/adufour
SIF=$RD/spritz-env.sif
SP=$RD/spritz_fabian
CT=$RD/DY_2026/analysis/combine_tools
SPR=$RD/DY_2026/analysis/spritz
VARS="mll rapll_abs costhetastar triple_diff"
APPT="apptainer exec -B /cvmfs -B /grid_mnt -B /etc/grid-security/certificates -B $SP/data/Full2018v9/samples/samples.json:/opt/spritz/data/Full2018v9/samples/samples.json --env SPRITZ_PATH=$SP --env PYTHONPATH=$SP/src:$SP $SIF"
log() { echo "[$(date '+%m-%d %H:%M')] $*"; }
wait_cluster() { while condor_q "$1" -af ProcId 2>/dev/null | grep -q .; do sleep 300; done; }
cd "$D" || exit 1

# ---- 1. condor
if [ ! -f condor/results_merged_new.pkl ]; then
  log "waiting for cluster $CLUSTER"; wait_cluster "$CLUSTER"
  bad() { cd condor; { find job_* -name chunks_job.pkl -size -100k -printf "%h\n"; grep -l "Error for chunk" job_*/err.txt 2>/dev/null | xargs -r -n1 dirname; } | sort -u; cd ..; }
  BAD=$(bad | paste -sd,)
  if [ -n "$BAD" ]; then
    log "retrying $(echo $BAD | tr , '\n' | wc -l) jobs"
    (cd condor && sed "s/^queue 1 Folder in .*/queue 1 Folder in $BAD/" submit.jdl > submit_retry.jdl && \
       c=$(condor_submit submit_retry.jdl | awk '/submitted to cluster/{gsub("\\.","",$NF); print $NF}') && echo $c > retry_cluster.txt)
    wait_cluster "$(cat condor/retry_cluster.txt)"
    BAD=$(bad | paste -sd,)
    [ -n "$BAD" ] && { log "STILL BAD after retry: $BAD"; exit 1; }
  fi
  log "condor done, all 3000 jobs complete"
  # ---- 2. merge
  bash merge.sh "$D" 4; tail -1 merge.log | grep -q "MERGE DONE" || { log "merge failed"; exit 1; }
  log "merged: $(ls -la condor/results_merged_new.pkl | awk '{print $5}') bytes"
fi

# ---- 3. postproc
if [ ! -f histos.root ]; then
  $APPT python3 pp_minnlo.py > postproc.log 2>&1; rc=$?
  [ $rc -eq 0 ] && [ -f histos.root ] || { log "postproc failed rc=$rc"; exit 1; }
  log "postproc done: histos.root $(ls -la histos.root | awk '{print $5}') bytes"
fi

# ---- 4. cards
if [ ! -f datacards/inc_mm/triple_diff/datacard.txt ]; then
  $APPT bash $SPR/spritz-cards-eft > cards.log 2>&1 || { log "cards failed"; exit 1; }
  log "cards: $(ls -d datacards/inc_mm/*/ | wc -l) variables"
fi

# ---- 5. 1D fits (as propcorr_v3_pdfcut/run_1d.sh)
source ~/.bashrc >/dev/null 2>&1
( dy_combine_morphing >/dev/null 2>&1
  for v in $VARS; do cd $D/datacards/inc_mm/$v
    [ -f higgsCombine.cbl_stat.individual.MultiDimFit.mH125.root ] && continue
    [ -f metadata.json ] || cp $REF/datacards/inc_mm/$v/metadata.json .   # same scan ranges as the reference
    createCombineJson.py --datacard datacard.txt --binname w1_ --output jsonComb.json > combjson.log 2>&1
    python3 $CT/createWS.py 1 > createWS.log 2>&1
    runScans.py 1 initial > initial.log 2>&1
    runScans.py 1 scan --points 100 > scan.log 2>&1
    runScans.py 1 initial --stat > initial_stat.log 2>&1
    runScans.py 1 scan --stat --points 100 > scan_stat.log 2>&1
    runPlots_compare.py 1 --label "Stat + Syst" --compare-stat > plots.log 2>&1
    log "1D $v: $(ls higgsCombine.*.individual.MultiDimFit.mH125.root 2>/dev/null | grep -vc _stat) scans, $(ls higgsCombine.*_stat.individual.MultiDimFit.mH125.root 2>/dev/null | wc -l) stat-only"
  done )
# ---- 6. summaries (analysis env, as propcorr_v3_pdfcut/run_1d.sh)
( dy_analysis >/dev/null 2>&1
  mkdir -p $D/summary_plots $D/compare_ref
  for v in $VARS; do cd $D/datacards/inc_mm/$v
    python3 $CT/makeSummary.py --indir . --logscale > summary.log 2>&1
    python3 $CT/makeSummaryCompare.py --indir1 $REF/datacards/inc_mm/$v --indir2 . --label1 "$REFLAB" --label2 "$LAB" \
        --verbose -o $D/compare_ref/eft_summary_compare_$v > $D/compare_ref/compare_$v.log 2>&1
  done
  cd $D/summary_plots
  python3 $CT/makeSummaryMultiVar.py --mll ../datacards/inc_mm/mll --rapll ../datacards/inc_mm/rapll_abs \
      --costhetastar ../datacards/inc_mm/costhetastar --triple_diff ../datacards/inc_mm/triple_diff --logscale > multivar.log 2>&1
)
log "ALL DONE"
