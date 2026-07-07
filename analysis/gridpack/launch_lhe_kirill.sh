#! /bin/bash

export X509_USER_PROXY=/home/llr/cms/biriukov/.globus/user_proxy.pem

source /cvmfs/cms.cern.ch/cmsset_default.sh

tar xavf WZto3LNu-1Jets-4FS_amcatnloFXFX-pythia8_el8_amd64_gcc10_CMSSW_12_4_8_tarball.tar.xz 

./runcmsgrid.sh 10000 $SEEDID 8

mkdir -p $SAVEFILE/lhe_output
mv cmsgrid_final.lhe $SAVEFILE/lhe_output/output_$SEEDID.lhe