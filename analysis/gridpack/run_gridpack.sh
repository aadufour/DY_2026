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

set -e

export X509_USER_PROXY=/grid_mnt/data__data.polcms/cms/adufour/.t3/proxy.cert


BIN=$1
NEVENTS=$2
SEED=$3
NCORES=$4

GRID_BASE="/grid_mnt/data__data.polcms/cms/adufour"
GRIDPACK_DIR="${GRID_BASE}/gridpacks/DYSMEFTMll${BIN}"
OUTPUT_DIR="${GRID_BASE}/LHE/DYSMEFTMll${BIN}"

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
if [ ! -f "${GRIDPACK_DIR}/runcmsgrid.sh" ]; then
    echo "ERROR: runcmsgrid.sh not found in ${GRIDPACK_DIR}"
    exit 1
fi

# Remove stale lock file if present
rm -f "${GRIDPACK_DIR}/process/madevent/RunWeb"

# Clean and recreate output directory
rm -rf "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"

# Run from inside the gridpack directory (runcmsgrid.sh uses relative paths)
cd "${GRIDPACK_DIR}"
bash runcmsgrid.sh "${NEVENTS}" "${SEED}" "${NCORES}"

# Move output to destination
mv "${GRIDPACK_DIR}/cmsgrid_final.lhe" "${OUTPUT_DIR}/unweighted_events.lhe"

echo "Done. Output in ${OUTPUT_DIR}"
ls -lh "${OUTPUT_DIR}"
