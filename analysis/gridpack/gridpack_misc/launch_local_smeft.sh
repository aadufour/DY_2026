#!/bin/bash
# Generate 10k events locally for each DYSMEFTMll bin
# Run this script on llruicms01 from the mg5amcnlo directory:
#   cd /grid_mnt/data__data.polcms/cms/adufour/MG5/mg5amcnlo
#   bash /path/to/launch_local_smeft.sh

set -e

MG5_DIR="/grid_mnt/data__data.polcms/cms/adufour/MG5/mg5amcnlo"
REWEIGHT_CARD="${MG5_DIR}/reweight_card.dat"

# Check reweight card is present
if [[ ! -f "$REWEIGHT_CARD" ]]; then
    echo "ERROR: reweight card not found at $REWEIGHT_CARD"
    echo "Copy it first:  scp <local>:/Users/albertodufour/code/DY2026/analysis/gridpack/gridpack_misc/cards/reweight_card.dat ${REWEIGHT_CARD}"
    exit 1
fi

# Bin name -> "mmll_min mmll_max"
declare -A BINS
BINS["DYSMEFTMll50_120"]="50.0 120.0"
BINS["DYSMEFTMll120_200"]="120.0 200.0"
BINS["DYSMEFTMll200_400"]="200.0 400.0"
BINS["DYSMEFTMll400_600"]="400.0 600.0"
BINS["DYSMEFTMll600_800"]="600.0 800.0"
BINS["DYSMEFTMll800_1000"]="800.0 1000.0"
BINS["DYSMEFTMll1000_3000"]="1000.0 3000.0"

for FOLDER in "${!BINS[@]}"; do
    PROC_DIR="${MG5_DIR}/${FOLDER}"

    if [[ ! -d "$PROC_DIR" ]]; then
        echo "WARNING: $PROC_DIR not found, skipping."
        continue
    fi

    echo "=============================="
    echo "Processing: $FOLDER"
    echo "=============================="

    read -r MLL_MIN MLL_MAX <<< "${BINS[$FOLDER]}"
    CARDS_DIR="${PROC_DIR}/Cards"

    # --- run_card.dat ---
    # Set nevents, mmll, mmllmax
    RUN_CARD="${CARDS_DIR}/run_card.dat"
    if [[ ! -f "$RUN_CARD" ]]; then
        echo "ERROR: $RUN_CARD not found, skipping $FOLDER."
        continue
    fi

    sed -i \
        -e "s|^\s*[0-9]*\s*=\s*nevents.*|  10000 = nevents ! Number of unweighted events requested|" \
        -e "s|^\s*[0-9.-]*\s*=\s*mmll\s.*|  ${MLL_MIN}  = mmll    ! min invariant mass of l+l- (same flavour) lepton pair|" \
        -e "s|^\s*[0-9.-]*\s*=\s*mmllmax\s.*|  ${MLL_MAX}  = mmllmax ! max invariant mass of l+l- (same flavour) lepton pair|" \
        "$RUN_CARD"

    echo "  Run card updated: nevents=10000, mmll=${MLL_MIN}, mmllmax=${MLL_MAX}"

    # --- reweight_card.dat ---
    cp "$REWEIGHT_CARD" "${CARDS_DIR}/reweight_card.dat"
    echo "  Reweight card installed."

    # --- launch generate_events ---
    echo "  Launching generate_events for $FOLDER ..."
    cd "$PROC_DIR"
    ./bin/generate_events -f 2>&1 | tee "${PROC_DIR}/generate_events.log"
    cd "$MG5_DIR"

    echo "  Done: $FOLDER"
done

echo ""
echo "All bins completed."
