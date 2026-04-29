#!/bin/bash
# run_gridpack_slc7_val.sh
#
# Like run_gridpack_syst.sh but for slc7 gridpacks — runs runcmsgrid.sh
# inside the cmssw-cc7 Singularity container (required for CMSSW_10_6_19).
#
# Arguments (passed by condor):
#   $1  mll bin string, e.g. "50_120"
#   $2  number of events, e.g. 5000
#   $3  random seed, e.g. 1
#   $4  number of cores, e.g. 8

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

rm -f  "${GRIDPACK_DIR}/process/madevent/RunWeb"
rm -rf "${GRIDPACK_DIR}/process/madevent/Events/GridRun_*"
rm -rf "${GRIDPACK_DIR}/process/madevent/Events/GridRun_PostProc_*"

rm -rf "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"

/cvmfs/cms.cern.ch/common/cmssw-cc7 -B /grid_mnt -B /home -- \
    "source /cvmfs/cms.cern.ch/cmsset_default.sh && cd ${GRIDPACK_DIR} && bash runcmsgrid.sh ${NEVENTS} ${SEED} ${NCORES}"

mv "${GRIDPACK_DIR}/cmsgrid_final.lhe" "${OUTPUT_DIR}/unweighted_events.lhe"

echo "Done. Output in ${OUTPUT_DIR}"
ls -lh "${OUTPUT_DIR}"
