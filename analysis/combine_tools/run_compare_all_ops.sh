#!/bin/bash
# run_compare_all_ops.sh
# Runs compare_ops.py for every SMEFT operator, linear and quadratic
# components only (no mixed/cross terms). Plots land in ./compare_ops_plots
# by default (or wherever --outdir points, forwarded via $1 if given).

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="${1:-./compare_ops_plots}"

OPERATORS=(
    cHDD cHWB cbWRe cbBRe cHj1 cHQ1
    cHj3 cHQ3 cHu  cHd   cHbq cHl1
    cHl3 cHe  cll1 clj1  clj3 cQl1
    cQl3 ceu  ced  cbe   cje  cQe
    clu  cld  cbl
)

for op in "${OPERATORS[@]}"; do
    for component in lin quad; do
        echo "=== $op ($component) ==="
        python3 "$SCRIPT_DIR/compare_ops.py" --op1 "$op" --component "$component" --outdir "$OUTDIR"
    done
done

echo "Done. Plots in $OUTDIR"
