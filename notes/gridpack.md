# GRIDPACK (8/04/2026)
modo compatto per girare mg5 con tutte le restriction card, reweight,...
più veloce, riproducibile, necesario per reco level
outputs una tarball che posso dare in pasto all'analisi completa cmssw (pythia e compagnia fino a nanoaod)


ssh llrcms
git clone https://github.com/cms-sw/genproductions.git
-> cartella genproductions


# condor (9/04/2026)
preso da kirill e reinterpretato
cd /grid_mnt/data__data.polcms/cms/adufour/genproductions/bin/MadGraph5_aMCatNLO/
gridpack_job.sh
gridpacks.submit



=============================================================
CONDOR GRIDPACK SUBMISSION AT LLR T3 - QUICK REFERENCE
=============================================================

PREREQUISITES
-------------
Create a long-lived VOMS proxy:
    voms-proxy-init --voms cms --valid 168:00 --out ~/.t3/proxy.cert

Verify:
    voms-proxy-info --file ~/.t3/proxy.cert


FILES
-----
All files live in:
    /grid_mnt/data__data.polcms/cms/adufour/genproductions/bin/MadGraph5_aMCatNLO/


1) gridpack_job.sh  (the script that runs on the worker node)
-------------------------------------------------------------
#!/bin/bash
export X509_USER_PROXY=/home/llr/cms/adufour/.t3/proxy.cert
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /grid_mnt/data__data.polcms/cms/adufour/genproductions/bin/MadGraph5_aMCatNLO/
./gridpack_generation.sh $1 cards/DY_SMEFT_Gridpacks/$1 local ALL


2) gridpacks.submit  (the condor configuration)
------------------------------------------------
Universe = vanilla
Executable = gridpack_job.sh
Arguments = $(bin)
input = /dev/null
output = logs/$(bin).out
error = logs/$(bin).err
log = logs/$(bin).log
getenv = true

request_memory = 20G
request_cpus = 8

T3Queue = long
WNTag = el9
+SingularityCmd = ""
include : /opt/exp_soft/cms/t3/t3queue |
requirements = regexp("llrgrwnvm[0-9]+.in2p3.fr", Machine) == FALSE

max_retries = 0

queue bin from (
    DYSMEFTMll50_120
    DYSMEFTMll120_200
    DYSMEFTMll200_400
    DYSMEFTMll400_600
    DYSMEFTMll600_800
    DYSMEFTMll800_1000
    DYSMEFTMll1000_3000
)


BEFORE SUBMITTING
-----------------
- Make sure logs/ directory exists:
    mkdir -p logs

- Make sure each proc card has the correct output name:
    for bin in DYSMEFTMll50_120 ...; do
        grep "^output" cards/DY_SMEFT_Gridpacks/$bin/${bin}_proc_card.dat
    done
  Output must be exactly: output <binname>

- Make sure no leftover directories from previous runs:
    for bin in DYSMEFTMll50_120 ...; do
        rm -rf $bin/ ${bin}.log
    done


CONDOR COMMANDS
---------------
Submit jobs:
    condor_submit -name llrt3condor gridpacks.submit

Monitor jobs:
    /opt/exp_soft/cms/t3/t3stat | grep adufour

Monitor continuously:
    watch -n 30 '/opt/exp_soft/cms/t3/t3stat | grep adufour'
    (exit with Ctrl+C)

Cancel all your jobs:
    condor_rm -name llrt3condor adufour

Check job queue:
    condor_q -name llrt3condor


JOB STATUS CODES
----------------
Q = Queued (waiting for a node)
R = Running
C = Completed
H = Held (something went wrong)


MONITORING OUTPUT
-----------------
Follow a specific job live:
    tail -f logs/DYSMEFTMll50_120.out
    tail -f logs/DYSMEFTMll50_120.err

Check tarballs when done:
    ls -lh *tarball*.tar.xz


OUTPUT
------
Tarballs appear in:
    /grid_mnt/data__data.polcms/cms/adufour/genproductions/bin/MadGraph5_aMCatNLO/
Named as:
    <binname>_el9_amd64_gcc11_CMSSW_13_2_9_tarball.tar.xz


NOTES
-----
- The proxy is valid for 7 days. Renew before submitting long jobs.
- Jobs run independently of your terminal. Safe to close SSH.
- Each job uses 8 cores and 20GB RAM on the T3 worker node.
- T3 long queue allows up to 10 days runtime.
- gridpack_generation.sh downloads a fresh MG5 v2.9.18 each time.
- The SMEFTsim model is injected from the local tarball at:
    cards/SMEFTsim_topU3l_MwScheme_UFO.tar.gz
  via the custom patch in gridpack_generation.sh.
=============================================================