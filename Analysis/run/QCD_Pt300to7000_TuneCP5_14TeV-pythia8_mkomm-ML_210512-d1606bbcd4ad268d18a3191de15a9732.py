from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferLogs = False
config.General.transferOutputs = True
config.General.workArea = 'crab_projects_ak8scouting_QCD_Pt300to7000'
config.General.requestName = 'QCD_Pt300to7000_TuneCP5_14TeV-pythia8'
config.section_('JobType')
config.JobType.numCores = 1
config.JobType.sendExternalFolder = True
config.JobType.pyCfgParams = ['inputDataset=/QCD_Pt300to7000_TuneCP5_14TeV-pythia8/mkomm-ML_210512-d1606bbcd4ad268d18a3191de15a9732/USER']
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.JobType.psetName = '../test/reHLT_2in1_crab.py'
config.JobType.maxMemoryMB = 2000
config.section_('Data')
config.Data.inputDataset = '/QCD_Pt300to7000_TuneCP5_14TeV-pythia8/mkomm-ML_210512-d1606bbcd4ad268d18a3191de15a9732/USER'
config.Data.outputDatasetTag = 'DeepNtuplesAK8-v00'
config.Data.publication = False
config.Data.unitsPerJob = 5
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.allowNonValidInputDataset = True
config.Data.outLFNDirBase = '/store/group/ml/Tagging4ScoutingHackathon/Adelina/DeepNtuples/ScoutingAK8-v00'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
config.section_('User')
config.section_('Debug')
