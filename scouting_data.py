# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: scouting --filein root://cmseos.fnal.gov//store/user/dmahon/condor/RPVSingleStopRun3MC/MINIAOD/MINIAOD-2000_100-1.root --eventcontent NANOAODGEN --datatier NANOAOD --fileout file:nanogen.root --conditions auto:phase1_2022_realistic --step NANOGEN --geometry DB:Extended --era Run3 --no_exec --mc --customise Run3ScoutingAnalysisTools/Analysis/Run3Scouting_cff.customise -n 1
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process('NANOGEN',Run3)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('PhysicsTools.NanoAOD.nanogen_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
     #'root://cmseos.fnal.gov//store/user/dmahon/condor/RPVSingleStopRun3MC/MINIAOD/MINIAOD-2000_100-1.root'
     '/store/data/Run2022F/ScoutingPFRun3/RAW/v1/000/361/303/00001/4f57c471-20d2-44db-b4f1-20e89b7c9fda.root'
     #'/store/mc/Run3Summer22DRPremix/MinBias_TuneCP5_14TeV-pythia8/AODSIM/Pilot_124X_mcRun3_2022_realistic_v11-v2/50000/b255101d-13ab-4068-a46c-a5ea9e1a8a6c.root'
     #'file:miniaodscouting.root'
    ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    makeTriggerResults = cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('scouting nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.NANOAODGENoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAOD'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:nanogen.root'),
    outputCommands = process.NANOAODGENEventContent.outputCommands
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2022_realistic', '')

# L1 RAW to Digi to PATTrigger
process.load("PhysicsTools.NanoAOD.triggerObjects_cff")
process.load("L1Trigger.Configuration.L1TRawToDigi_cff")
process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi")
process.load("PhysicsTools.PatAlgos.slimming.selectedPatTrigger_cfi")
process.load("PhysicsTools.PatAlgos.slimming.slimmedPatTrigger_cfi")
process.gtStage2Digis.InputLabel = cms.InputTag("hltFEDSelectorL1")
process.nanoAOD_seq = cms.Sequence(process.L1TRawToDigi+process.patTrigger+process.selectedPatTrigger+process.slimmedPatTrigger)
process.trigger_task = cms.Task(process.unpackedPatTrigger,process.triggerObjectTable,process.l1bits)
# Path and EndPath definitions
process.nanoAOD_step = cms.Path(process.nanoAOD_seq,process.trigger_task)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODGENoutput_step = cms.EndPath(process.NANOAODGENoutput)

# Schedule definition
process.schedule = cms.Schedule(process.nanoAOD_step,process.endjob_step,process.NANOAODGENoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.NanoAOD.nanogen_cff
from PhysicsTools.NanoAOD.nanogen_cff import customizeNanoGENFromMini 

#call to customisation function customizeNanoGENFromMini imported from PhysicsTools.NanoAOD.nanogen_cff
process = customizeNanoGENFromMini(process)

# Automatic addition of the customisation function from Run3ScoutingAnalysisTools.Analysis.Run3Scouting_cff
from Run3ScoutingAnalysisTools.Analysis.Run3Scouting_cff import customise 

#call to customisation function customise imported from Run3ScoutingAnalysisTools.Analysis.Run3Scouting_cff
process = customise(process)

# End of customisation functions

#process.schedule.remove(process.nanoAOD_step)
process.NANOAODGENoutput.outputCommands.extend(['keep edmTriggerResults_*_*_*'])
process.NANOAODGENoutput.isPFScouting = cms.untracked.bool(True)

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
