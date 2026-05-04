import CRABClient
from WMCore.Configuration import Configuration
from multiprocessing import Process
config = Configuration()

gp_path = "/eos/user/a/aldufour/gridpacks/DYSMEFTMll{}_slc7_amd64_gcc700_CMSSW_10_6_19_tarball.tar.xz".format(mll_bin.split("mll_")[1])

events_per_job = 500
PROD = 'DYSMEFTMll-nanoaod18_SMEFTsim_' + mll_bin

config.section_('General')
config.General.workArea = PROD + '_chain'
config.General.requestName = PROD

config.section_('JobType')
config.JobType.scriptExe = 'runners/2018/run_chain_test.sh'
config.JobType.psetName = 'do_nothing_cfg.py'
config.JobType.pluginName = 'PrivateMC'
config.JobType.outputFiles = ['SMP-RunIISummer20UL18NanoAODv9-00051.root']
config.JobType.inputFiles = [
    'CMSSW_10_6_26.tar.gz',  # Patched NanoAOD: fixes GenWeightsTableProducer for EFT reweighting weights
    'copy_gridpack.py',
    'modifyCfg.py',
    'runners/2018/run_chain_test.sh',
    'runners/2018/chain_step_0_test.sh',
    'runners/2018/chain_step_1_test.sh',
    'runners/2018/chain_step_2_test.sh',
    'runners/2018/chain_step_3_test.sh',
    'runners/2018/chain_step_4_test.sh',
    'runners/2018/chain_step_5_test.sh',
    '2018_LO/SMP-RunIISummer20UL18wmLHEGEN-00061_1_cfg.py',
    '2018_LO/SMP-RunIISummer20UL18SIM-00035_1_cfg.py',
    '2018_LO/SMP-RunIISummer20UL18DIGIPremix-00035_1_cfg.py',
    '2018_LO/SMP-RunIISummer20UL18HLT-00035_1_cfg.py',
    '2018_LO/SMP-RunIISummer20UL18RECO-00035_1_cfg.py',
    '2018_LO/SMP-RunIISummer20UL18MiniAODv2-00051_1_cfg.py',
    '2018_LO/SMP-RunIISummer20UL18NanoAODv9-00051_1_cfg.py',
]
config.JobType.disableAutomaticOutputCollection = False
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxMemoryMB = 2500
config.JobType.numCores = 1

config.section_('Data')
config.Data.unitsPerJob = events_per_job
NJOBS = 1000
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.splitting = 'EventBased'
config.Data.publication = False
config.Data.outputPrimaryDataset = PROD
config.Data.outputDatasetTag = PROD
config.Data.outLFNDirBase = '/store/user/aldufour/3DY_SMEFTsim_LO/'

config.section_('User')

config.section_('Site')
config.Site.storageSite = 'T2_FR_GRIF_LLR'

config.JobType.scriptArgs = ['nEvents=' + str(config.Data.unitsPerJob)]
config.JobType.scriptArgs.append('inputGridpack=' + gp_path)
print('Submitting jobs with script args --> ' + ' '.join(config.JobType.scriptArgs))
print('Submitting jobs with unitsPerJob --> ' + str(config.Data.unitsPerJob) + ' totalUnits --> ' + str(config.Data.totalUnits), ' primary dataset --> ', str(config.Data.outputPrimaryDataset))
