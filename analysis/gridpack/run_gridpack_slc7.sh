#!/bin/bash
# run_gridpack_slc7.sh
# Production event generation from slc7 gridpacks with 325300/325500 PDF systematics.

set -e

export X509_USER_PROXY=/grid_mnt/data__data.polcms/cms/adufour/.t3/proxy.cert

BIN=$1
NEVENTS=$2
SEED=$3
NCORES=$4

GRID_BASE="/grid_mnt/data__data.polcms/cms/adufour"
GRIDPACK_DIR="${GRID_BASE}/gridpacks/SYST_slc7/DYSMEFTMll${BIN}"
OUTPUT_DIR="${GRID_BASE}/LHE/SYST_slc7/DYSMEFTMll${BIN}"

echo "============================="
echo "Bin      : ${BIN}"
echo "N events : ${NEVENTS}"
echo "Seed     : ${SEED}"
echo "Cores    : ${NCORES}"
echo "Gridpack : ${GRIDPACK_DIR}"
echo "Output   : ${OUTPUT_DIR}"
echo "============================="

if [ ! -d "${GRIDPACK_DIR}" ]; then
    echo "ERROR: gridpack directory not found: ${GRIDPACK_DIR}"
    exit 1
fi
if [ ! -f "${GRIDPACK_DIR}/runcmsgrid.sh" ]; then
    echo "ERROR: runcmsgrid.sh not found in ${GRIDPACK_DIR}"
    exit 1
fi

rm -f "${GRIDPACK_DIR}/process/madevent/RunWeb"
rm -rf "${GRIDPACK_DIR}/process/madevent/Events/GridRun_*"
rm -rf "${GRIDPACK_DIR}/process/madevent/Events/GridRun_PostProc_*"

rm -rf "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"

cd "${GRIDPACK_DIR}"
bash runcmsgrid.sh "${NEVENTS}" "${SEED}" "${NCORES}"

mv "${GRIDPACK_DIR}/cmsgrid_final.lhe" "${OUTPUT_DIR}/unweighted_events.lhe"

echo "Done. Output in ${OUTPUT_DIR}"
ls -lh "${OUTPUT_DIR}"
