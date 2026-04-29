#!/bin/bash

mkdir -p LO_generation
mkdir -p LO_generation/runners

cp -r runners/2018 LO_generation/runners
cp -r 2018_LO LO_generation
cp do_nothing_cfg.py modifyCfg.py copy_gridpack.py LO_generation

sed '6i mll_bin = "mll_50_120" ' crab_sub_2018_LO.py > LO_generation/crab_sub_2018_LO_mll_50_120.py
sed '6i mll_bin = "mll_120_200" ' crab_sub_2018_LO.py > LO_generation/crab_sub_2018_LO_mll_120_200.py
sed '6i mll_bin = "mll_200_400" ' crab_sub_2018_LO.py > LO_generation/crab_sub_2018_LO_mll_200_400.py
sed '6i mll_bin = "mll_400_600" ' crab_sub_2018_LO.py > LO_generation/crab_sub_2018_LO_mll_400_600.py
sed '6i mll_bin = "mll_600_800" ' crab_sub_2018_LO.py > LO_generation/crab_sub_2018_LO_mll_600_800.py
sed '6i mll_bin = "mll_800_1000" ' crab_sub_2018_LO.py > LO_generation/crab_sub_2018_LO_mll_800_1000.py
sed '6i mll_bin = "mll_1000_3000" ' crab_sub_2018_LO.py > LO_generation/crab_sub_2018_LO_mll_1000_3000.py

cd LO_generation

if [ $# -ne 0 ]; then
   if [ "$1" -eq 1 ]; then
     crab submit crab_sub_2018_LO_mll_50_120.py
     crab submit crab_sub_2018_LO_mll_120_200.py
     crab submit crab_sub_2018_LO_mll_200_400.py
     crab submit crab_sub_2018_LO_mll_400_600.py
     crab submit crab_sub_2018_LO_mll_600_800.py
     crab submit crab_sub_2018_LO_mll_800_1000.py
     crab submit crab_sub_2018_LO_mll_1000_3000.py
   else
     echo "Pass 1 if you want to crab submit"
   fi
fi
