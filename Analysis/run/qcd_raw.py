from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferLogs = False
config.General.transferOutputs = True
config.General.workArea = 'crab_projects'
config.General.requestName = 'QCD_Pt15to7000_TuneCP5_14TeV-pythia8'
config.section_('JobType')
config.JobType.numCores = 1
config.JobType.sendExternalFolder = True
#config.JobType.pyCfgParams = ['outputFile=nano.root', 'inputDataset=/QCD_Pt300to7000_TuneCP5_14TeV-pythia8']
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.JobType.psetName = './scoutingPF.py'
config.JobType.maxMemoryMB = 2000
config.section_('Data')
config.Data.outputDatasetTag = 'QCD_Pt15to7000_TuneCP5_14TeV-pythia8'
config.Data.inputDataset = "/QCD_Pt15to7000_TuneCP5_14TeV-pythia8/Run3Summer21DR-FlatPU0to80FEVT_castor_120X_mcRun3_2021_realistic_v6-v1/GEN-SIM-DIGI-RAW"
config.Data.publication = False
config.Data.unitsPerJob = 5
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.allowNonValidInputDataset = True
#config.Data.outLFNDirBase = '/eos/user/a/adlintul/Run3/samples/MC'
config.section_('Site')
config.Site.storageSite = 'T3_CH_CERNBOX'
config.section_('User')
config.section_('Debug')
