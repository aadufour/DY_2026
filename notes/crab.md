# SLC7 gridpack rebuilding (23/04/2026)
=============================================================
REBUILDING GRIDPACKS FOR SLC7 (required for 2018 UL GEN chain)
=============================================================

REASON
------
El9 gridpacks (el9_amd64_gcc11_CMSSW_13_2_9) are incompatible with
the 2018 UL GEN cfg which runs in CMSSW_10_6_19 (SLC7).
Must rebuild as slc7_amd64_gcc700_CMSSW_10_6_19.

KEY CHANGES FROM EL9 WORKFLOW
------------------------------
1. genproductions branch: must use mg265UL (not master)
   cd /grid_mnt/data__data.polcms/cms/adufour/genproductions
   git checkout mg265UL

2. gridpack_job.sh: run inside cmssw-cc7 container with /grid_mnt bind:
   /cvmfs/cms.cern.ch/common/cmssw-cc7 -B /grid_mnt -B /home -- \
   "source /cvmfs/cms.cern.ch/cmsset_default.sh && cd /grid_mnt/.../MadGraph5_aMCatNLO/ && bash gridpack_generation.sh $1 cards/DY_SMEFT_Gridpacks/SYST/$1 local ALL"

3. gridpacks.submit: keep WNTag=el9, keep +SingularityCmd=""
   (container handled inside gridpack_job.sh, not by condor)

PROC CARD FIXES (MG5 2.9.18 → 2.6.5 compatibility)
-----------------------------------------------------
Remove these lines from all *_proc_card.dat (not in MG5 2.6.5):
   set max_t_for_channel
   set zerowidth_tchannel
   set low_mem_multicore_nlo_generation
   set loop_color_flows

RUN CARD FIX
------------
Remove from all *_run_card.dat:
   sde_strategy   (causes Fortran compile error: no IMPLICIT type)

OUTPUT
------
Tarballs: DYSMEFTMll{bin}_slc7_amd64_gcc700_CMSSW_10_6_19_tarball.tar.xz
Location: /grid_mnt/data__data.polcms/cms/adufour/genproductions/bin/MadGraph5_aMCatNLO/
Size: ~48M each (vs 35M el9) — 700 process files, reweight_card.dat present
=============================================================

crab_sub_2018.py

file sh e py di config

/grid_mnt/data__data.polcms/cms/adufour/3DYGeneration/2018_LO


generate p p > l+ l-
add process p p > l+ l- j 

genera entrambi: eventi con 0j e 1j

pythia probabilmente aggiunge un jet (ISR/FSR) aggiunge radiazione qcd, adronizzazione,... tutto quello che non è perturbativa (MG5 è perturbativo).
Mathing and merging pythia, sottrarre contributi che sono doppi



#non tocco
SMP-RunIISummer20UL18RECO-00035_1_cfg.py
SMP-RunIISummer20UL18SIM-00035_1_cfg.py

#da giocare con SMP-RunIISummer20UL18wmLHEGEN-00061_1_cfg.py


Occhio a CMSSW release!
Mi adatto alla versione hardcoded in 2018!!

girare in slc7 su qualxsiasi macchina cmssw-cc7

event oper job nei cfg files per crab


# GEN-only CRAB submission (28/04/2026)
=============================================================
FIRST CRAB STEP: LHE+GEN only (step 0 of full chain)
=============================================================

GOAL
----
Validate CRAB submission before attempting full GEN->NanoAOD chain.
500 events/job x 1000 jobs = 500k events per mll bin, all 7 bins.

OUTPUT
------
SMP-RunIISummer20UL18wmLHEGEN-00061.root — GEN-level ROOT file with
showered particles + all EFT and PDF weights embedded.
Stored at: /store/user/aldufour/3DY_SMEFTsim_LO/ on T2_FR_GRIF_LLR
(physically: eos.grif.fr:11000/eos/grif/cms/llr/store/user/aldufour/)

JOB FLOW
--------
1. run_gen_only.sh       — CRAB entry point. Downloads gridpack from EOS
                           via xrootd, then calls chain_step_0_test.sh.
2. copy_gridpack.py      — xrdcp root://eosuser.cern.ch//eos/user/a/aldufour/gridpacks/...
3. chain_step_0_test.sh  — sets up CMSSW_10_6_19_patch3 (slc7), calls
                           modifyCfg.py to inject seed/events/gridpack into
                           the cfg, then runs cmsRun.
4. cmsRun                — MadGraph5 (gridpack) -> Pythia8 -> Photos QED FSR

KEY FILES (analysis/crab/)
--------------------------
crab_sub_2018_LO.py          — CRAB config (scriptExe, inputFiles, storage, etc.)
submit_2018_LO.sh            — generates 7 per-bin configs and submits
runners/2018/run_gen_only.sh — new: step 0 only (stripped from run_chain_test.sh)

CHANGES VS GIACOMO'S FULL CHAIN
---------------------------------
- scriptExe: run_chain_test.sh -> run_gen_only.sh
- outputFiles: NanoAOD root -> SMP-RunIISummer20UL18wmLHEGEN-00061.root
- inputFiles: stripped to step 0 only (no CMSSW_10_6_26.tar.gz, no steps 1-5)
- maxMemoryMB: 5000 -> 2500 (dry run measured only 591 MB actual usage)
- storageSite: T2_CH_CERN -> T2_FR_GRIF_LLR (no write access to CMS EOS yet)

LXPLUS SETUP (needed every new session)
----------------------------------------
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd ~/CMSSW_13_3_3/src && cmsenv && cd -   # any el9 CMSSW works, just for crab
source /cvmfs/cms.cern.ch/common/crab-setup.sh
voms-proxy-init --voms cms --valid 192:00

CHECK STATUS
------------
cd ~/crab_gen/LO_generation
for bin in 50_120 120_200 200_400 400_600 600_800 800_1000 1000_3000; do
    echo "=== mll_${bin} ==="
    crab status -d DYSMEFTMll-nanoaod18_SMEFTsim_mll_${bin}_newPythia_plus_Photos/crab_DYSMEFTMll-nanoaod18_SMEFTsim_mll_${bin} 2>&1 | grep -E "Status|Warning|jobs"
done

NOTES
-----
- T2_CH_CERN refused: /store/user/aldufour not provisioned there yet
- CMSSW_10_6_19_patch3 (slc7) not available on lxplus el9 — irrelevant,
  jobs set it up themselves on the grid worker nodes
- modifyCfg.py must be in inputFiles — patches cfg at runtime
- Gridpacks on EOS: /eos/user/a/aldufour/gridpacks/DYSMEFTMll<bin>_slc7_...tar.xz
