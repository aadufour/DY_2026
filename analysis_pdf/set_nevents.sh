#!/bin/bash
# Usage: bash set_nevents.sh 50000

NEVENTS=$1
BASE="/home/llr/cms/adufour/MG5/mg5amcnlo"

if [ -z "$NEVENTS" ]; then
    echo "Usage: $0 <nevents>"
    exit 1
fi

for RUN_CARD in ${BASE}/DY_pdf_*/Cards/run_card.dat; do
    sed -i "s|^[[:space:]]*[0-9]*[[:space:]]*= nevents.*|  ${NEVENTS} = nevents ! Number of unweighted events requested|" "$RUN_CARD"
    echo "$(dirname $(dirname $RUN_CARD) | xargs basename): $(grep nevents $RUN_CARD)"
done
