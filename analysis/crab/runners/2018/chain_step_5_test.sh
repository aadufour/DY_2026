#!/bin/sh
set -e
#set -x

SEED=$1
EVENTS=$2

RUN_DIR=${PWD}
echo ">> Setting RUN_DIR to ${RUN_DIR}"

CMSSW_RELEASE=CMSSW_10_6_26
SCRAM_ARCH=slc7_amd64_gcc700

# Check if a tarred version exists
if [ -f ${CMSSW_RELEASE}.tar.gz ]; then

  if [ -d ${CMSSW_RELEASE} ]; then
      echo ">> Cleaning up existing ${CMSSW_RELEASE} directory"
      rm -r ${CMSSW_RELEASE}
    fi

    echo ">> Extracting ${CMSSW_RELEASE} from tarball"
    tar -xzf ${CMSSW_RELEASE}.tar.gz
    cd ${CMSSW_RELEASE}/src 
    eval `scramv1 runtime -sh`
    scramv1 b ProjectRename
    scramv1 b -j 8
    if [ -d PhysicsTools ]; then
        cd PhysicsTools/NanoAOD; 
        scramv1 b -j 8
        cd -
	eval `scramv1 runtime -sh`
    fi
    # redundancy does not harm 
    eval `scramv1 runtime -sh`
    cd ../..

elif [ "${CMSSW_RELEASE}" != "local" ]; then
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
    # patching GenWeightTable to retrieve correctly the reweighting weights 
    git-cms-addpkg PhysicsTools/NanoAOD

    awk 'NR == 578 {
    print "      std::regex weightgroupRwgt(\"<weightgroup\\\\s+(?:name)=\\\"(.*)\\\"\\\\s+(?:weight_name_strategy)=\\\"(.*)\\\"\\\\s*>\");"
    next
    }
    { print }
    ' PhysicsTools/NanoAOD/plugins/GenWeightsTableProducer.cc > tmp && mv tmp PhysicsTools/NanoAOD/plugins/GenWeightsTableProducer.cc
    
    eval `scramv1 runtime -sh`
    scramv1 b -j 8
    eval `scramv1 runtime -sh`
    cd -

fi

echo ${PWD} 

python ${RUN_DIR}/modifyCfg.py ${RUN_DIR}/SMP-RunIISummer20UL18NanoAODv9-00051_1_cfg.py ${RUN_DIR}/step_5_cfg.py --randomSeeds=${SEED} --events=${EVENTS}

cmsRun -e -j FrameworkJobReport.xml ${RUN_DIR}/step_5_cfg.py
