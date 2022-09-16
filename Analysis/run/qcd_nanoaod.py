from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferLogs = False
config.General.transferOutputs = True
config.General.workArea = 'crab_projects'
config.General.requestName = 'QCD_nanoaod'
config.section_('JobType')
config.JobType.numCores = 1
config.JobType.sendExternalFolder = True
config.JobType.pyCfgParams = ['isMC=True', 'GlobalTagMC=123X_mcRun3_2021_realistic_v15', 'trigProcess=HLTX']
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.JobType.psetName = '../test/ScoutingNanoAOD_cfg.py'
config.JobType.maxMemoryMB = 2000
config.section_('Data')
config.Data.outputDatasetTag = 'ScoutingNanoAOD'
config.Data.outputPrimaryDataset = 'QCD_Pt15to7000_TuneCP5_14TeV-pythia8'
config.Data.userInputFiles = open('files_qcd.txt').readlines()
config.Data.publication = False
config.Data.unitsPerJob = 12
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.allowNonValidInputDataset = True
#config.Data.outLFNDirBase = '/eos/user/a/adlintul/Run3/samples/MC'
config.section_('Site')
config.Site.storageSite = 'T3_CH_CERNBOX'
config.section_('User')
config.section_('Debug')
