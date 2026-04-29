#!/bin/sh
set -e
#set -x

SEED=$1

RUN_DIR=${PWD}
echo ">> Setting RUN_DIR to ${RUN_DIR}"

CMSSW_RELEASE=CMSSW_10_6_19_patch3
SCRAM_ARCH=slc7_amd64_gcc700

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


xrdcp -f root://eosuser.cern.ch//eos/user/g/gboldrin/Zee_dim6_LHE/mll_binned/gridpacks_v2_2025_02_07/zee_dim6_mll50-100_slc7_amd64_gcc700_CMSSW_10_6_19_tarball.tar.xz .


python ${RUN_DIR}/modifyCfg.py ${RUN_DIR}/SMP-RunIISummer20UL17wmLHEGEN-00065_1_cfg.py ${RUN_DIR}/step_0_cfg.py --events=20 --randomSeeds=${SEED}

echo "PRINTING PWD chain, where FrameworkJobReport.xml will be"
pwd

cmsRun -e -j FrameworkJobReport.xml ${RUN_DIR}/step_0_cfg.py
