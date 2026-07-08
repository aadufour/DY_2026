#!/bin/bash
# run_compare_all_ops.sh
# Runs compare_ops.py for every SMEFT operator, linear and quadratic
# components only (no mixed/cross terms).
#
# Usage:
#   ./run_compare_all_ops.sh [--outdir DIR] [--no-normalize]

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="./compare_ops_plots"
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --outdir)
            OUTDIR="$2"
            shift 2
            ;;
        --no-normalize)
            EXTRA_ARGS+=("--no-normalize")
            shift
            ;;
        *)
            echo "Unknown argument: $1"
            exit 1
            ;;
    esac
done

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
        python3 "$SCRIPT_DIR/compare_ops.py" --op1 "$op" --component "$component" --outdir "$OUTDIR" "${EXTRA_ARGS[@]}"
    done
done

echo "Done. Plots in $OUTDIR"
