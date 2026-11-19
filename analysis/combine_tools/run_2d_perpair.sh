#!/bin/bash
# Full 2D scan chain with per-pair boxes for one variable, unattended (run in tmux on LLR).
#
#   0. initial fits (text2workspace + initial) for pairs that lack them, copied from SRC's condor_initial
#   1. coarse pass: 50x50 grid, box = 5x the 1D 95% interval per operator (from V1D 1D scans)
#      + up to 2 rounds re-running "open"/missing pairs with the doubled box
#   2. fine pass: 100x100 grid in each pair's box (95% extent + 35%), 25 pieces per short job (2 jobs per pair)
#   3. one more fine pass for pairs still open/tight/missing
#   4. final/ (symlinks to the best scan per pair), plots (plot_2d_scans.py), xrdcp to EOS
#
# Usage: run_2d_perpair.sh <VARDIR> <SRC_VARDIR_with_condor_initial> <V1D_dir_with_1D_scans> <keep_pairs_csv> <eos_subdir>
# e.g.   run_2d_perpair.sh .../propcorr_v3_2d/datacards/inc_mm/triple_diff .../propcorr_v3_2d/datacards/inc_mm/mll \
#            .../propcorr_v3_pdfcut/datacards/inc_mm/triple_diff cHQ1_cHQ3,cQl1_cQl3 RECO3D_propcorr_v3/2D_triple_diff
VARDIR=$1; SRC=$2; V1D=$3; KEEP=$4; EOSDIR=$5
CT=$(cd "$(dirname "$0")" && pwd)
FULL=$VARDIR/full
log() { echo "[$(date '+%m-%d %H:%M')] $*"; }
wait_cluster() { while condor_q "$1" -af ProcId 2>/dev/null | grep -q .; do sleep 120; done; }
submit() { condor_submit "$1" | awk '/submitted to cluster/{gsub("\\.","",$NF); print $NF}'; }

source ~/.bashrc >/dev/null 2>&1; dy_combine_morphing >/dev/null 2>&1
mkdir -p "$FULL" || exit 1
cd "$VARDIR" || exit 1
ALL=$(python3 -c "import json,itertools;o=list(json.load(open('metadata.json'))['operators']);print(','.join('_'.join(c) for c in itertools.combinations(o,2)))")
log "pairs: $(echo $ALL | tr , '\n' | wc -l)"

# ---------------- 0. initial fits ----------------
if [ ! -d "$FULL/initial" ]; then
  mkdir -p "$FULL/initial/scripts" "$FULL/initial/logs"
  n=0; i=0; cur=""
  for s in "$SRC"/condor_initial/scripts/*.sh; do
    p=$(grep -o "model_[A-Za-z0-9]*_[A-Za-z0-9]*\.root" "$s" | head -1 | sed 's/model_//;s/\.root//')
    [ -f "$VARDIR/higgsCombine.initialFit_$p.MultiDimFit.mH125.root" ] && continue
    if [ $((n % 15)) -eq 0 ]; then
      cur=$FULL/initial/scripts/b_$i.sh; i=$((i+1))
      grep -v -E "^(combine|text2workspace|set -e)" "$s" | sed "s|cd $SRC\$|cd $VARDIR|" > "$cur"; chmod +x "$cur"
      echo "$cur" >> "$FULL/initial/joblist.txt"
    fi
    grep -E "^(text2workspace|combine)" "$s" >> "$cur"; n=$((n+1))
  done
  if [ $n -gt 0 ]; then
    sed -e "s|$SRC/condor_scan_bundled|$FULL/initial|g" -e "s/request_memory = .*/request_memory = 4G/" \
        "$SRC/condor_scan_bundled/scan.sub" > "$FULL/initial/scan.sub"
    c=$(submit "$FULL/initial/scan.sub"); log "initial: $n pairs in $i jobs, cluster $c"; wait_cluster $c
  fi
fi
log "initial fits present: $(ls $VARDIR/higgsCombine.initialFit_*.root | wc -l)"

# ---------------- scan helper ----------------
# run_scan <dir> <metadata> <pairs_csv> <points> <split> <bundle>
run_scan() {
  local D=$1 M=$2 P=$3 PTS=$4 SPL=$5 BUN=$6
  mkdir -p "$D"; cd "$D" || return 1
  for p in $(echo $P | tr , ' '); do ln -sfn "$VARDIR/higgsCombine.initialFit_$p.MultiDimFit.mH125.root" .; done
  cp "$M" metadata_pairs.json
  python3 "$CT/runScans.py" 2 scan --points=$PTS --doSplitPoints=$SPL --metadata metadata_pairs.json --pairs "$P" \
      --condor --condor-dir condor_scan > runscans.log 2>&1
  python3 "$CT/bundle_scan_jobs.py" $BUN condor_scan condor_scan_bundled
  local c; c=$(submit condor_scan_bundled/scan.sub); log "$(basename $D): $(echo $P | tr , '\n' | wc -l) pairs, cluster $c"
  wait_cluster $c
  for p in $(echo $P | tr , ' '); do
    ls higgsCombine.$p.individual.POINTS.*.root >/dev/null 2>&1 && \
      hadd -f higgsCombine.$p.individual.MultiDimFit.mH125.root higgsCombine.$p.individual.POINTS.*.MultiDimFit.mH125.root >/dev/null 2>&1
  done
  log "$(basename $D): merged $(ls higgsCombine.*.individual.MultiDimFit.mH125.root 2>/dev/null | wc -l)"
  cd "$VARDIR"
}
# boxes <metadata_in> <scan_dir> <pairs_csv> <out> <listfile> : updates boxes, writes pairs needing another pass
boxes() {
  python3 "$CT/make_pair_boxes.py" --metadata "$1" --scan-dir "$2" --pairs "$3" --keep "$KEEP" --pad 0.35 \
      -o "$4" --list-out "$5" > "$5.log" 2>&1
  tail -1 "$5.log"
}

# ---------------- 1. coarse ----------------
if [ ! -f "$FULL/meta_coarse.json" ]; then
  python3 - "$V1D" "$FULL/meta_coarse.json" <<'EOF'
import json, sys, ROOT
v1d, out = sys.argv[1:]
meta = json.load(open("metadata.json"))
for op in meta["operators"]:
    d = ROOT.RDataFrame("limit", f"{v1d}/higgsCombine.{op}.individual.MultiDimFit.mH125.root").AsNumpy([f"k_{op}", "deltaNLL"])
    x, z = d[f"k_{op}"], 2 * (d["deltaNLL"] - d["deltaNLL"].min())
    r = 5 * max(abs(x[z < 3.84].min()), abs(x[z < 3.84].max()))
    meta["operators"][op] = [-float(f"{r:.3g}"), float(f"{r:.3g}")]
    print(op, meta["operators"][op])
json.dump(meta, open(out, "w"), indent=4)
EOF
fi
META=$FULL/meta_coarse.json
if [ ! -f "$FULL/meta_pass1.json" ]; then
  run_scan "$FULL/pass1" "$META" "$ALL" 2500 5 5
  boxes "$META" "$FULL/pass1" "$ALL" "$FULL/meta_pass1.json" "$FULL/redo_pass1.txt"
fi
META=$FULL/meta_pass1.json
for r in 1 2; do
  [ -f "$FULL/meta_pass1_r$r.json" ] && { META=$FULL/meta_pass1_r$r.json; continue; }
  prev=$([ $r = 1 ] && echo redo_pass1 || echo redo_pass1_r$((r-1)))
  OPEN=$(grep -E " open |MISSING" "$FULL/$prev.txt.log" | awk '{print $1}' | paste -sd,)
  [ -z "$OPEN" ] && break
  run_scan "$FULL/pass1_r$r" "$META" "$OPEN" 2500 5 5
  boxes "$META" "$FULL/pass1_r$r" "$OPEN" "$FULL/meta_pass1_r$r.json" "$FULL/redo_pass1_r$r.txt"
  META=$FULL/meta_pass1_r$r.json
done
log "coarse done, boxes from $META"

# ---------------- 2. fine ----------------
if [ ! -f "$FULL/meta_pass2.json" ]; then
  run_scan "$FULL/pass2" "$META" "$ALL" 10000 50 25
  boxes "$META" "$FULL/pass2" "$ALL" "$FULL/meta_pass2.json" "$FULL/redo_pass2.txt"
fi
REDO=$(grep -v -x -F -f <(echo $KEEP | tr , '\n') "$FULL/redo_pass2.txt" | paste -sd,)
if [ -n "$REDO" ] && [ ! -f "$FULL/meta_pass3.json" ]; then
  run_scan "$FULL/pass3" "$FULL/meta_pass2.json" "$REDO" 10000 50 25
  boxes "$FULL/meta_pass2.json" "$FULL/pass3" "$REDO" "$FULL/meta_pass3.json" "$FULL/redo_pass3.txt"
fi

# ---------------- 4. final, plots, EOS ----------------
mkdir -p "$FULL/final"; cd "$FULL/final"
for p in $(echo $ALL | tr , ' '); do
  f=higgsCombine.$p.individual.MultiDimFit.mH125.root
  if [ -f "$FULL/pass3/$f" ]; then ln -sfn "$FULL/pass3/$f" .; elif [ -f "$FULL/pass2/$f" ]; then ln -sfn "$FULL/pass2/$f" .; fi
done
log "final: $(ls | wc -l) pairs"
# plots with runPlots.py (mkEFTScan, no interpolation; plot_2d_scans.py upsampled with a cubic spline and drew
# fake islands around isolated outliers, 2026-11-19)
mkdir -p "$FULL/plots"; cd "$FULL/plots"
for f in "$FULL"/final/higgsCombine.*.root; do ln -sfn "$f" .; done
cp "$VARDIR/metadata.json" .
python3 "$CT/runPlots.py" 2 --cores 8 > "$FULL/plots.log" 2>&1
rm -f higgsCombine.*.root metadata.json
log "plots: $(ls $FULL/plots/*.png | wc -l)"
python3 "$CT/../../eos_tools/pushtoEOS_parallel.py" -j 8 --overwrite "$FULL/plots" "$EOSDIR" > "$FULL/eos.log" 2>&1 || echo "EOS FAIL, see $FULL/eos.log"
log "ALL DONE"
