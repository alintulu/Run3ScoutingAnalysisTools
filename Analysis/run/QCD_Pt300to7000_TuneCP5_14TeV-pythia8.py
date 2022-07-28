from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferLogs = False
config.General.transferOutputs = True
config.General.workArea = 'crab_projects'
config.General.requestName = 'QCD_Pt300to7000_TuneCP5_14TeV-pythia8'
config.section_('JobType')
config.JobType.numCores = 1
config.JobType.sendExternalFolder = True
config.JobType.pyCfgParams = ['outputFile=nano.root', 'inputDataset=/QCD_Pt300to7000_TuneCP5_14TeV-pythia8']
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.JobType.psetName = '../test/ScoutingNanoAOD_cfg.py'
config.JobType.maxMemoryMB = 2000
config.section_('Data')
config.Data.outputPrimaryDataset = 'QCD_Pt300to7000_TuneCP5_14TeV-pythia8'
config.Data.outputDatasetTag = 'QCD_Pt300to7000_TuneCP5_14TeV-pythia8'
config.Data.userInputFiles = open('files_QCD.txt').readlines()
config.Data.publication = False
config.Data.unitsPerJob = 5
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.allowNonValidInputDataset = True
config.Data.outLFNDirBase = '/store/group/ml/Tagging4ScoutingHackathon/Adelina/DeepNtuples/12_3_0/ScoutingAK8-v00/Ntuples'
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
config.section_('User')
config.section_('Debug')
