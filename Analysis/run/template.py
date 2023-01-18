from CRABClient.UserUtilities import config
config = config()

config.General.transferLogs = False
config.General.transferOutputs = True
config.General.workArea = 'crab_projects'
config.General.requestName = 'NAME'

config.JobType.pluginName = 'Analysis'
config.JobType.numCores = 1
config.JobType.sendExternalFolder = True
config.JobType.pyCfgParams = ['outputFile=ntuple.root', 'inputDataset=DATASET']
config.JobType.allowUndistributedCMSSW = True
config.JobType.psetName = '../test/ScoutingNanoAOD_cfg.py'
#config.JobType.maxMemoryMB = 2000

config.Data.inputDataset = 'DATASET'
config.Data.outputDatasetTag = 'CAMPAIGN'
config.Data.publication = False
config.Data.inputDBS = 'global'
config.Data.splitting = 'Automatic'
config.Data.outLFNDirBase = '/store/group/ml/Tagging4ScoutingHackathon/Adelina/hbb/ak8'

config.Site.storageSite = 'T2_CH_CERN'
