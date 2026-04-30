#!/bin/bash

mkdir -p chain_production
mkdir -p chain_production/runners

cp -r runners/2018 chain_production/runners
cp -r 2018_LO chain_production
cp do_nothing_cfg.py modifyCfg.py copy_gridpack.py chain_production
cp CMSSW_10_6_26.tar.gz chain_production

sed '6i mll_bin = "mll_50_120" ' crab_sub_2018_chain.py > chain_production/crab_sub_2018_chain_mll_50_120.py
sed '6i mll_bin = "mll_120_200" ' crab_sub_2018_chain.py > chain_production/crab_sub_2018_chain_mll_120_200.py
sed '6i mll_bin = "mll_200_400" ' crab_sub_2018_chain.py > chain_production/crab_sub_2018_chain_mll_200_400.py
sed '6i mll_bin = "mll_400_600" ' crab_sub_2018_chain.py > chain_production/crab_sub_2018_chain_mll_400_600.py
sed '6i mll_bin = "mll_600_800" ' crab_sub_2018_chain.py > chain_production/crab_sub_2018_chain_mll_600_800.py
sed '6i mll_bin = "mll_800_1000" ' crab_sub_2018_chain.py > chain_production/crab_sub_2018_chain_mll_800_1000.py
sed '6i mll_bin = "mll_1000_3000" ' crab_sub_2018_chain.py > chain_production/crab_sub_2018_chain_mll_1000_3000.py

cd chain_production

if [ $# -ne 0 ]; then
    if [ "$1" -eq 1 ]; then
        crab submit crab_sub_2018_chain_mll_50_120.py
        crab submit crab_sub_2018_chain_mll_120_200.py
        crab submit crab_sub_2018_chain_mll_200_400.py
        crab submit crab_sub_2018_chain_mll_400_600.py
        crab submit crab_sub_2018_chain_mll_600_800.py
        crab submit crab_sub_2018_chain_mll_800_1000.py
        crab submit crab_sub_2018_chain_mll_1000_3000.py
    else
        echo "Pass 1 if you want to crab submit"
    fi
fi
