#!/bin/sh
set -e

SEED=$1
EVENTS=$2
FILENAME=$3

bash chain_step_0_test.sh ${SEED} ${EVENTS} ${FILENAME}
bash chain_step_1_test.sh ${SEED} ${EVENTS}
bash chain_step_2_test.sh ${SEED} ${EVENTS}
bash chain_step_3_test.sh ${SEED} ${EVENTS}
bash chain_step_4_test.sh ${SEED} ${EVENTS}
bash chain_step_5_test.sh ${SEED} ${EVENTS}
