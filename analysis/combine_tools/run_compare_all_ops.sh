#!/bin/bash
# run_compare_all_ops.sh
# Runs compare_ops.py for every SMEFT operator, linear and quadratic
# components only (no mixed/cross terms).
#
# Defaults to the high-stats caches (propcorr's 50_120 top-up + baseline's
# matching 1M-event top-up), not the original lower-stats ones -- override
# with --prop-cache/--base-cache if you want the originals.
#
# Usage:
#   ./run_compare_all_ops.sh [--outdir DIR] [--no-normalize] [--prop-cache PATH] [--base-cache PATH]

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="./compare_ops_plots"
PROP_CACHE="/grid_mnt/data__data.polcms/cms/adufour/LHE/propcorr/CACHE/lhe_cache_propcorr_parallel.pkl"
BASE_CACHE="/grid_mnt/data__data.polcms/cms/adufour/LHE/baseline/CACHE/lhe_cache_baseline_highstats.pkl"
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
        --prop-cache)
            PROP_CACHE="$2"
            shift 2
            ;;
        --base-cache)
            BASE_CACHE="$2"
            shift 2
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
        python3 "$SCRIPT_DIR/compare_ops.py" --op1 "$op" --component "$component" --outdir "$OUTDIR" \
            --prop-cache "$PROP_CACHE" --base-cache "$BASE_CACHE" "${EXTRA_ARGS[@]}"
    done
done

echo "Done. Plots in $OUTDIR"
