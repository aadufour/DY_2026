#!/bin/sh
set -e
#set -x

SEED=$(( $1 + 0 ))

echo "PRINTING PWD run_test"
pwd



bash chain_step_0_test.sh ${SEED}
bash chain_step_1_test.sh ${SEED}
bash chain_step_2_test.sh ${SEED}
bash chain_step_3_test.sh ${SEED}
bash chain_step_4_test.sh ${SEED}
bash chain_step_5_test.sh ${SEED}

