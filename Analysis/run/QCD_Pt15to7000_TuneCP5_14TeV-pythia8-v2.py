from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferLogs = False
config.General.transferOutputs = True
config.General.workArea = 'crab_projects_ak8scouting_QCD_Pt15to7000'
config.General.requestName = 'QCD_Pt15to7000_TuneCP5_14TeV-pythia8-v2'
config.section_('JobType')
config.JobType.numCores = 1
config.JobType.sendExternalFolder = True
config.JobType.pyCfgParams = ['inputDataset=/QCD_Pt15to7000_TuneCP5_14TeV-pythia8/Run3Winter21DRMiniAOD-FlatPU0to80_Scouting_Patatrack_castor_112X_mcRun3_2021_realistic_v16-v2/MINIAODSIM']
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.JobType.psetName = '../test/reHLT_2in1_crab.py'
config.JobType.maxMemoryMB = 2000
config.section_('Data')
config.Data.inputDataset = '/QCD_Pt15to7000_TuneCP5_14TeV-pythia8/Run3Winter21DRMiniAOD-FlatPU0to80_Scouting_Patatrack_castor_112X_mcRun3_2021_realistic_v16-v2/MINIAODSIM'
config.Data.outputDatasetTag = 'DeepNtuplesAK8-v00'
config.Data.publication = False
config.Data.unitsPerJob = 5
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.allowNonValidInputDataset = True
config.Data.outLFNDirBase = '/store/group/ml/Tagging4ScoutingHackathon/Adelina/DeepNtuples/ScoutingAK8-v00'
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
config.section_('User')
config.section_('Debug')
