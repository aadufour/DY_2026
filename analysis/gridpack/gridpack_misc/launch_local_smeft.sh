#!/bin/bash
# Generate 10k SMEFT events locally for each DYSMEFTMll bin.
#
# Step 1: run mg5_aMC with the proc card to generate+output DYSMEFTMll50_120
# Step 2: copy that compiled process to all other bins (MEs are identical)
# Step 3: for each bin, patch the run_card and call generate_events
#
# Usage on llruicms01:
#   bash /grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/gridpack/gridpack_misc/launch_local_smeft.sh

set -e

MG5="/grid_mnt/data__data.polcms/cms/adufour/MG5/mg5amcnlo"
REPO="/grid_mnt/data__data.polcms/cms/adufour/DY_2026"
CARDS="${REPO}/analysis/gridpack/gridpack_misc/cards"

REFERENCE_BIN="DYSMEFTMll50_120"
REFERENCE_DIR="${MG5}/${REFERENCE_BIN}"

BINS=(
    "DYSMEFTMll50_120     50.0   120.0"
    "DYSMEFTMll120_200   120.0   200.0"
    "DYSMEFTMll200_400   200.0   400.0"
    "DYSMEFTMll400_600   400.0   600.0"
    "DYSMEFTMll600_800   600.0   800.0"
    "DYSMEFTMll800_1000  800.0  1000.0"
    "DYSMEFTMll1000_3000 1000.0 3000.0"
)

# -------------------------------------------------------
# Step 1: generate+output the reference bin if not ready
# -------------------------------------------------------
if [[ ! -f "${REFERENCE_DIR}/bin/generate_events" ]]; then
    echo "Initializing process in ${REFERENCE_BIN} via mg5_aMC ..."

    # Write a temporary proc card that outputs to the reference bin
    TMP_PROC=$(mktemp /tmp/smeft_proc_XXXX.mg5)
    cat > "$TMP_PROC" <<EOF
import model SMEFTsim_topU3l_MwScheme_UFO-all_massless
define p = g u c d s u~ c~ d~ s~
define j = g u c d s u~ c~ d~ s~
define l+ = e+ mu+
define l- = e- mu-
define vl = ve vm vt
define vl~ = ve~ vm~ vt~
define p = p b b~
generate p p > mu+ mu- QCD=0 SMHLOOP=0 NP<=1
output ${REFERENCE_DIR}
EOF

    cd "$MG5"
    ./bin/mg5_aMC "$TMP_PROC"
    rm "$TMP_PROC"
else
    echo "Reference process directory already initialized, skipping generate step."
fi

# -------------------------------------------------------
# Step 2: copy reference to all other bins
# -------------------------------------------------------
for BIN in "${BINS[@]}"; do
    read -r TAG MLL_MIN MLL_MAX <<< "$BIN"
    TARGET="${MG5}/${TAG}"

    if [[ "$TAG" == "$REFERENCE_BIN" ]]; then
        continue
    fi

    if [[ ! -f "${TARGET}/bin/generate_events" ]]; then
        echo "Copying compiled process to ${TAG} ..."
        cp -r "$REFERENCE_DIR" "$TARGET"
    else
        echo "${TAG} already initialized, skipping copy."
    fi
done

# -------------------------------------------------------
# Step 3: for each bin, patch cards and generate events
# -------------------------------------------------------
for BIN in "${BINS[@]}"; do
    read -r TAG MLL_MIN MLL_MAX <<< "$BIN"
    PROC_DIR="${MG5}/${TAG}"
    CARDS_DIR="${PROC_DIR}/Cards"

    echo "=============================="
    echo "Bin: $TAG  [${MLL_MIN}, ${MLL_MAX}]"
    echo "=============================="

    # Install cards
    cp "${CARDS}/run_card.dat"      "${CARDS_DIR}/run_card.dat"
    cp "${CARDS}/reweight_card.dat" "${CARDS_DIR}/reweight_card.dat"

    # Patch run_card for this bin
    sed -i \
        -e "s|^\s*[0-9]*\s*=\s*nevents.*|  10000 = nevents ! Number of unweighted events requested|" \
        -e "s|^\s*[0-9.-]*\s*=\s*mmll\b.*|  ${MLL_MIN}  = mmll    ! min invariant mass of l+l- (same flavour) lepton pair|" \
        -e "s|^\s*[0-9.-]*\s*=\s*mmllmax\b.*|  ${MLL_MAX}  = mmllmax ! max invariant mass of l+l- (same flavour) lepton pair|" \
        -e "s|^\s*\S*\s*=\s*run_tag.*|  ${TAG} = run_tag ! name of the run|" \
        "${CARDS_DIR}/run_card.dat"

    echo "  run_card: nevents=10000, mmll=${MLL_MIN}, mmllmax=${MLL_MAX}"

    # Generate events
    cd "$PROC_DIR"
    ./bin/generate_events -f 2>&1 | tee "${PROC_DIR}/generate_events_${TAG}.log"
    cd - > /dev/null

    echo "  Done: output in ${PROC_DIR}/Events/"
done

echo ""
echo "All bins completed."
