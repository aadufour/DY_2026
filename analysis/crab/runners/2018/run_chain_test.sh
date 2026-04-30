#!/bin/sh
set -e
#set -x

SEED=$(( $1 + 0 ))
EVENTS="${2#*=}"
GRIDPACK=$3

echo "PRINTING PWD run_test"
pwd

python3 copy_gridpack.py -i ${GRIDPACK} --use-xrootd
FILENAME=$(basename "$GRIDPACK")

echo "----->", ${FILENAME}

cmssw-el7 -- bash chain_step_0_test.sh ${SEED} ${EVENTS} ${FILENAME}
bash chain_step_1_test.sh ${SEED} ${EVENTS} 
bash chain_step_2_test.sh ${SEED} ${EVENTS}
bash chain_step_3_test.sh ${SEED} ${EVENTS}
bash chain_step_4_test.sh ${SEED} ${EVENTS}
bash chain_step_5_test.sh ${SEED} ${EVENTS}

