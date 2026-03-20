#!/usr/bin/env bash
# make_all_histograms.sh
#
# Runs make_combine_histograms.py for every operator, writing per-operator
# ROOT files, then merges them into a single histograms.root.
#
# Usage:
#   bash make_all_histograms.sh [--C 1.0] [--lumi 59740] [--unrolled]
#
# All extra flags are forwarded to make_combine_histograms.py.

set -euo pipefail

SCRIPT="/Users/albertodufour/code/DY2026/datacard_test/make_combine_histograms.py"
OUTDIR="/Users/albertodufour/code/DY2026/restriction_table"
FINAL_ROOT="${OUTDIR}/histograms.root"

OPERATORS=(
    cHDD cHWB cHj1
    cHj3 cHu  cHd
    cHl1 cHl3 cHe
    cll1 clj1 clj3
    ceu  ced  cje
    clu  cld
)

EXTRA_ARGS="$@"

TMPDIR_HIST="${OUTDIR}/tmp_hists"
mkdir -p "${TMPDIR_HIST}"

echo "=============================================="
echo " Running make_combine_histograms.py for each operator"
echo " Extra args: ${EXTRA_ARGS:-<none>}"
echo "=============================================="

for OP in "${OPERATORS[@]}"; do
    OUTFILE="${TMPDIR_HIST}/histograms_${OP}.root"
    echo ""
    echo "--- Operator: ${OP} ---"
    python3 "${SCRIPT}" --op "${OP}" --output "${OUTFILE}" --no-plot ${EXTRA_ARGS}
done

echo ""
echo "=============================================="
echo " Merging into ${FINAL_ROOT}"
echo "=============================================="

python3 - <<'PYEOF'
import os
import uproot
import numpy as np

OUTDIR   = "/Users/albertodufour/code/DY2026/datacard_test"
TMPDIR   = os.path.join(OUTDIR, "tmp_hists")
FINAL    = os.path.join(OUTDIR, "histograms.root")
CHANNEL  = "triple_DY"

OPERATORS = [
    'cHDD', 'cHWB', 'cHj1',
    'cHj3', 'cHu',  'cHd',
    'cHl1', 'cHl3', 'cHe',
    'cll1', 'clj1', 'clj3',
    'ceu',  'ced',  'cje',
    'clu',  'cld',
]

with uproot.recreate(FINAL) as out:
    # Write sm and data_obs once, from the first operator's file
    first_file = os.path.join(TMPDIR, f"histograms_{OPERATORS[0]}.root")
    with uproot.open(first_file) as f:
        for name in ("sm", "data_obs"):
            key = f"{CHANNEL}/{name}"
            out[key] = f[key]
            print(f"  wrote {key}  (from {OPERATORS[0]})")

    # Write quad_<op> and sm_lin_quad_<op> from each operator's file
    for op in OPERATORS:
        src = os.path.join(TMPDIR, f"histograms_{op}.root")
        with uproot.open(src) as f:
            for kind in (f"quad_{op}", f"sm_lin_quad_{op}"):
                key = f"{CHANNEL}/{kind}"
                out[key] = f[key]
                print(f"  wrote {key}")

print(f"\nDone: {FINAL}")

# Readback summary
print("\nReadback check:")
with uproot.open(FINAL) as f:
    for name in f[CHANNEL].keys():
        h = f[f"{CHANNEL}/{name}"]
        print(f"  {name:<30}  nbins={h.member('fXaxis').member('fNbins')}  integral={h.values().sum():.4e}")
PYEOF

echo ""
echo "All done. Output: ${FINAL_ROOT}"
