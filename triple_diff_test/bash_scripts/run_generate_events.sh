#!/bin/bash

BIN_EDGES=(50 120 200 400 600 800 1000 3000)
BASE="/Users/albertodufour/MG5_2_9_18/mg5amcnlo"

for (( i=0; i<${#BIN_EDGES[@]}-1; i++ )); do
    MIN=${BIN_EDGES[$i]}
    MAX=${BIN_EDGES[$((i+1))]}

    DIR="${BASE}/DY_all_${MIN}_${MAX}"

    echo "============================================"
    echo "Running generate_events in ${DIR} ..."
    echo "============================================"

    cd "$DIR" && ./bin/generate_events

    echo "Done with DY_all_${MIN}_${MAX}."
done

echo "All bins completed."
