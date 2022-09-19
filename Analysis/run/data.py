from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'ScoutingPFRun3'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../test/ScoutingNanoAOD_cfg.py'
config.JobType.pyCfgParams = ['outputFile=scouting.root', 'isMC=False', 'GlobalTagData=123X_dataRun3_HLT_v14']
config.Data.inputDataset = '/ScoutingPFRun3/Run2022B-v1/RAW'
config.Data.inputDBS = 'global'
config.Data.splitting = 'Automatic'
config.Data.LumiMask = 'Cert_Collisions2022_355100_357101_Golden.json'
config.Data.publication = False
config.Data.outputDatasetTag = '2022B'
#config.Data.ignoreLocality = True
config.Data.outLFNDirBase = '/store/group/ml/Tagging4ScoutingHackathon/Adelina/Run3/samples'

#config.Site.whitelist = ['T2_US*','T2_CH*']
config.Site.storageSite = 'T2_CH_CERN'
