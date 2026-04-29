#!/bin/sh
set -e

SEED=$(( $1 + 0 ))
EVENTS="${2#*=}"
GRIDPACK=$3

echo "PRINTING PWD run_gen_only"
pwd
echo "SEED: ${SEED}, EVENTS: ${EVENTS}, GRIDPACK: ${GRIDPACK}"

python3 copy_gridpack.py -i ${GRIDPACK} --use-xrootd
FILENAME=$(basename "$GRIDPACK")

echo "----->", ${FILENAME}

cmssw-el7 -- bash chain_step_0_test.sh ${SEED} ${EVENTS} ${FILENAME}
