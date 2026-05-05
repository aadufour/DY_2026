#!/bin/sh
set -e

SEED=$(( $1 + 0 ))
EVENTS="${2#*=}"
GRIDPACK=$3

echo "PRINTING PWD run_chain_test"
pwd
echo "SEED: ${SEED}, EVENTS: ${EVENTS}, GRIDPACK: ${GRIDPACK}"

python3 copy_gridpack.py -i ${GRIDPACK} --use-xrootd
FILENAME=$(basename "$GRIDPACK")

echo "----->", ${FILENAME}

cmssw-el7 -- bash run_chain_inner.sh ${SEED} ${EVENTS} ${FILENAME}
