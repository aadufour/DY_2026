#!/bin/bash
# Set up DY_pdf_LO_HI directories by copying a base MG5 process folder
# and patching the mmll/mmllmax cuts in the run card.
#
# Usage:
#   bash setup_DY_bins.sh
#
# Assumes this is run on llrcms (Linux — uses GNU sed, no '' after -i).

BIN_EDGES=(50 120 200 400 600 800 1000 3000)

BASE_SRC="/home/llr/cms/adufour/MG5/mg5amcnlo/DY_pdf"
BASE_DST="/home/llr/cms/adufour/MG5/mg5amcnlo"

for (( i=0; i<${#BIN_EDGES[@]}-1; i++ )); do
    MIN=${BIN_EDGES[$i]}
    MAX=${BIN_EDGES[$((i+1))]}

    DST="${BASE_DST}/DY_pdf_${MIN}_${MAX}"

    echo "Setting up ${DST} ..."
    cp -r "$BASE_SRC" "$DST"

    RUN_CARD="${DST}/Cards/run_card.dat"

    # GNU sed (Linux): no '' after -i
    sed -i "s|^[[:space:]]*[0-9.-]*[[:space:]]*= mmll[[:space:]]*!.*| ${MIN}.0   = mmll    ! min invariant mass of l+l- (same flavour) lepton pair|" "$RUN_CARD"
    sed -i "s|^[[:space:]]*[0-9.-]*[[:space:]]*= mmllmax[[:space:]]*!.*| ${MAX}.0  = mmllmax ! max invariant mass of l+l- (same flavour) lepton pair|" "$RUN_CARD"

    echo "  mmll lines:"
    grep "mmll" "$RUN_CARD"
done

echo "Done."
