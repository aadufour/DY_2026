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
#   ./run_compare_all_ops.sh [--outdir DIR] [--no-normalize] [--prop-cache PATH] [--base-cache PATH] [--lims VAL]
#
# Operators that actually enter the W/Z/H propagator-width corrections
# (SMEFTsim_practical_guide.pdf, Appendix A, eqs. A.10/A.12/A.19 -- cHQ1/cHQ3
# only via the Z width, since W can't decay to t+b) get their plots routed
# into $OUTDIR/propcorr_relevant/ instead of $OUTDIR/ directly.

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
        --lims)
            EXTRA_ARGS+=("--lims" "$2")
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

# enters delta-Gamma_W and/or delta-Gamma_Z per Appendix A
PROPCORR_OPS=(cHDD cHWB cHl1 cHl3 cll1 cHj1 cHj3 cHQ1 cHQ3)

for op in "${OPERATORS[@]}"; do
    op_outdir="$OUTDIR"
    for pc_op in "${PROPCORR_OPS[@]}"; do
        if [[ "$op" == "$pc_op" ]]; then
            op_outdir="$OUTDIR/propcorr_relevant"
            break
        fi
    done
    for component in lin quad; do
        echo "=== $op ($component) -> $op_outdir ==="
        python3 "$SCRIPT_DIR/compare_ops.py" --op1 "$op" --component "$component" --outdir "$op_outdir" \
            --prop-cache "$PROP_CACHE" --base-cache "$BASE_CACHE" "${EXTRA_ARGS[@]}"
    done
done

echo "Done. Plots in $OUTDIR (propcorr-relevant operators in $OUTDIR/propcorr_relevant)"
