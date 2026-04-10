#!/bin/bash
# run_gridpack.sh
#
# Condor job script: generate events for one mll bin from a CMS gridpack.
#
# Arguments (passed by condor via Arguments = ...):
#   $1  mll bin string, e.g. "50_120"
#   $2  number of events, e.g. 10000
#   $3  random seed, e.g. 1
#   $4  number of cores, e.g. 1
#
# Gridpack location  : /home/adufour/gridpack_tests/DYSMEFTMll${BIN}/
# Output location    : /grid_mnt/data__data.polcms/cms/adufour/MG5/mg5amcnlo/DYSMEFTMll${BIN}/Events/run_01/

set -e

BIN=$1
NEVENTS=$2
SEED=$3
NCORES=$4

GRIDPACK_DIR="/home/adufour/gridpack_tests/DYSMEFTMll${BIN}"
OUTPUT_DIR="/grid_mnt/data__data.polcms/cms/adufour/MG5/mg5amcnlo/DYSMEFTMll${BIN}/Events/run_01"

echo "============================="
echo "Bin      : ${BIN}"
echo "N events : ${NEVENTS}"
echo "Seed     : ${SEED}"
echo "Cores    : ${NCORES}"
echo "Gridpack : ${GRIDPACK_DIR}"
echo "Output   : ${OUTPUT_DIR}"
echo "============================="

# Sanity checks
if [ ! -d "${GRIDPACK_DIR}" ]; then
    echo "ERROR: gridpack directory not found: ${GRIDPACK_DIR}"
    exit 1
fi
if [ ! -f "${GRIDPACK_DIR}/runcmsgrid_LO.sh" ]; then
    echo "ERROR: runcmsgrid_LO.sh not found in ${GRIDPACK_DIR}"
    exit 1
fi

# Clean and recreate output directory
rm -rf "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"

# Run the gridpack from the output directory so all output lands there
cd "${OUTPUT_DIR}"
bash "${GRIDPACK_DIR}/runcmsgrid_LO.sh" "${NEVENTS}" "${SEED}" "${NCORES}"

# Rename to match the naming convention expected by build_cache.py
if [ -f "cmsgrid_final.lhe.gz" ]; then
    mv cmsgrid_final.lhe.gz unweighted_events.lhe.gz
    echo "Renamed cmsgrid_final.lhe.gz -> unweighted_events.lhe.gz"
fi

echo ""
echo "Done. Output in ${OUTPUT_DIR}"
ls -lh "${OUTPUT_DIR}"
