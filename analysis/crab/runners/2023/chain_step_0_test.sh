#!/bin/sh
set -e
#set -x

SEED=$1
EVENTS=$2
GRIDPACK=$3

RUN_DIR=${PWD}



#== CMSSW: gridpack=zee_dim6_mll100-200
#== CMSSW: 1
#== CMSSW: events=1000

echo "ARGUMENTS IN ORDER"
echo ${SEED}
echo ${EVENTS}
echo ${GRIDPACK}
echo "------------"

# mv SMP-RunIISummer20UL18wmLHEGEN-00061*.py SMP-RunIISummer20UL18wmLHEGEN-00061_1_cfg.py 

ls 

echo ">> Setting RUN_DIR to ${RUN_DIR}"

CMSSW_RELEASE=CMSSW_13_0_14
SCRAM_ARCH=el8_amd64_gcc11

if [ "${CMSSW_RELEASE}" != "local" ]; then
    if [ -d ${CMSSW_RELEASE} ]; then
      echo ">> Cleaning up existing ${CMSSW_RELEASE} directory"
      rm -r ${CMSSW_RELEASE}
    fi
    echo ">> Setting up release area for ${CMSSW_RELEASE} and arch ${SCRAM_ARCH}"
    if [ ! -d ${CMSSW_RELEASE} ]; then
      scram project CMSSW ${CMSSW_RELEASE}
    fi

    cd ${CMSSW_RELEASE}/src
    eval `scramv1 runtime -sh`
    cd -

fi



python3 ${RUN_DIR}/modifyCfg.py ${RUN_DIR}/GEN-Run3Summer23wmLHEGS-00327_1_cfg.py ${RUN_DIR}/step_0_cfg.py --randomSeeds=${SEED} --events=${EVENTS} --gridpack="${RUN_DIR}/${GRIDPACK}" --strategy=0 

echo "PRINTING PWD chain, where FrameworkJobReport.xml will be"
pwd

cmsRun -e -j FrameworkJobReport.xml ${RUN_DIR}/step_0_cfg.py jobNum=$1 ${EVENTS}
