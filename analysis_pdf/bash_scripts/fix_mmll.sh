#!/bin/bash
# Fix mmll/mmllmax cuts in already-existing DY_pdf_LO_HI run cards.
# Run this on llrcms before re-generating events.

BIN_EDGES=(50 120 200 400 600 800 1000 3000)
BASE="/home/llr/cms/adufour/MG5/mg5amcnlo"

for (( i=0; i<${#BIN_EDGES[@]}-1; i++ )); do
    MIN=${BIN_EDGES[$i]}
    MAX=${BIN_EDGES[$((i+1))]}
    CARD="${BASE}/DY_pdf_${MIN}_${MAX}/Cards/run_card.dat"

    sed -i "s|^[[:space:]]*[0-9.-]*[[:space:]]*= mmll[[:space:]]*!.*| ${MIN}.0   = mmll    ! min invariant mass of l+l- (same flavour) lepton pair|" "$CARD"
    sed -i "s|^[[:space:]]*[0-9.-]*[[:space:]]*= mmllmax[[:space:]]*!.*| ${MAX}.0  = mmllmax ! max invariant mass of l+l- (same flavour) lepton pair|" "$CARD"

    echo "DY_pdf_${MIN}_${MAX}:"
    grep "mmll" "$CARD"
done
