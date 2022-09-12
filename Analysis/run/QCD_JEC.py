from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferLogs = False
config.General.transferOutputs = True
config.General.workArea = 'crab_projects'
config.General.requestName = 'QCD_JEC'
config.section_('JobType')
config.JobType.numCores = 1
config.JobType.sendExternalFolder = True
#config.JobType.pyCfgParams = ['inputDataset=/TTbar_TuneCP5_14TeV-pythia8/mkomm-ML_210512-d1606bbcd4ad268d18a3191de15a9732/USER']
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.JobType.psetName = '../test/ScoutingNanoAOD_cfg.py'
config.JobType.maxMemoryMB = 2000
config.section_('Data')
config.Data.inputDataset = '/QCD_Pt15to7000_TuneCP5_14TeV-pythia8/Run3Winter21DRMiniAOD-FlatPU0to80FEVT_castor_112X_mcRun3_2021_realistic_v16-v2/GEN-SIM-DIGI-RAW'
config.Data.outputDatasetTag = 'QCD_JEC'
config.Data.publication = False
config.Data.unitsPerJob = 50
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.allowNonValidInputDataset = True
config.Data.outLFNDirBase = '/store/group/ml/Tagging4ScoutingHackathon/Adelina'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
config.section_('User')
config.section_('Debug')
