#! /bin/bash

export X509_USER_PROXY=/home/llr/cms/adufour/.t3/proxy.cert

source /cvmfs/cms.cern.ch/cmsset_default.sh

tar xavf $TARBALL

./runcmsgrid.sh 10000 $SEEDID 8

mkdir -p $SAVEFILE/lhe_output
mv cmsgrid_final.lhe $SAVEFILE/lhe_output/output_$SEEDID.lhe
