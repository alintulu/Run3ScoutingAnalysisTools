# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: reHLT --python_filename scouting_reHLT.py --eventcontent RAWMINIAODSIM --customise HLTrigger/Configuration/customizeHLTforPatatrack.customizeHLTforPatatrackTriplets --filein file:/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/Scouting/Run3/ML_210512/QCD_Pt15to7000_TuneCP5_14TeV-pythia8/ML_210512/210514_170522/0000/scouting_1.root --fileout file:scouting_1.root --conditions auto:phase1_2021_realistic --step HLT:GRun --geometry DB:Extended --era Run3 --no_exec --mc -n 10
import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_cff import Run3

# Set parameters externally
from FWCore.ParameterSet.VarParsing import VarParsing
params = VarParsing('analysis')

params.maxEvents = -1
params.outputFile = "nano.root"

params.register('inputDataset',
    '',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Input dataset"
)

params.register('outputFileMini',
    'mini.root',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output path for MiniAOD"
)

params.parseArguments()

process = cms.Process('reHLT',Run3)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('HLTrigger.Configuration.HLT_GRun_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(params.maxEvents),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(params.inputFiles),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
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
    annotation = cms.untracked.string('reHLT nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RAWMINIAODSIMoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    ),
    dropMetaData = cms.untracked.string('ALL'),
    eventAutoFlushCompressedSize = cms.untracked.int32(20971520),
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string(params.outputFileMini),
    outputCommands = cms.untracked.vstring( 'drop *',
      'keep *_hltFEDSelectorL1_*_*',
      'keep *_hltScoutingEgammaPacker_*_*',
      'keep *_hltScoutingMuonPacker_*_*',
      'keep *_hltScoutingPFPacker_*_*',
      'keep *_hltScoutingPrimaryVertexPacker_*_*',
      'keep *_hltScoutingTrackPacker_*_*',
      'keep edmTriggerResults_*_*_*',
      'keep *_prunedGenParticles_*_*',
      'keep *_gtStage2Digis_*_*',
      'keep *_TriggerResults_*_*',
      'keep *_packedGenParticles_*_*',
      'keep *_packedPFCandidates_*_*',
      'keep *_slimmedJets_*_*',
      'keep *_slimmedGenJets_*_*',
      'keep *_slimmedAddPileupInfo_*_*',
      'keep *_slimmedSecondaryVertices_*_*',
      'keep *_offlineSlimmedPrimaryVertices_*_*',
      'keep *_fixedGridRhoFastjetAll_*_*',
    ),
    overrideBranchesSplitLevel = cms.untracked.VPSet(
        cms.untracked.PSet(
            branch = cms.untracked.string('patPackedCandidates_packedPFCandidates__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('recoGenParticles_prunedGenParticles__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('patTriggerObjectStandAlones_slimmedPatTrigger__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('patPackedGenParticles_packedGenParticles__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('patJets_slimmedJets__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('recoVertexs_offlineSlimmedPrimaryVertices__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('recoVertexs_offlineSlimmedPrimaryVerticesWithBS__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('recoCaloClusters_reducedEgamma_reducedESClusters_*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEBRecHits_*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEERecHits_*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('recoGenJets_slimmedGenJets__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('patJets_slimmedJetsPuppi__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedESRecHits_*'),
            splitLevel = cms.untracked.int32(99)
        )
    ),
    overrideInputFileSplitLevels = cms.untracked.bool(True),
    splitLevel = cms.untracked.int32(0)
)

L1Info = ['L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7', 'L1_DoubleMu_12_5','L1_DoubleMu_15_7','L1_TripleMu_5_3_3','L1_TripleMu_5_5_3','L1_QuadMu0','L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4','L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18','L1_DoubleMu4_SQ_OS_dR_Max1p2','L1_SingleMu22','L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4','L1_DoubleMu4p5_SQ_OS_dR_Max1p2','L1_DoubleMu4p5_SQ_OS','L1_DoubleMu0er1p5_SQ_dR_Max1p4','L1_DoubleMu0er2p0_SQ_dR_Max1p4','L1_DoubleMu0_SQ']

# Define the services needed for the treemaker
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(params.outputFile)
)

process.mmtree = cms.EDAnalyzer('ScoutingNanoAOD',

        doL1 = cms.bool(False),
        triggerConfiguration = cms.PSet(
                hltResults            = cms.InputTag('TriggerResults','','HLT'),
                l1tResults            = cms.InputTag(''),
                daqPartitions         = cms.uint32(1),
                l1tIgnoreMaskAndPrescale = cms.bool(False),
                throw                 = cms.bool(False)
        ),
        ReadPrescalesFromFile = cms.bool( False ),
        AlgInputTag       = cms.InputTag("gtStage2Digis"),
        l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis"),
        l1tExtBlkInputTag = cms.InputTag("gtStage2Digis"),
        l1Seeds           = cms.vstring(L1Info),
        pfcandsParticleNet = cms.InputTag("hltScoutingPFPacker"),
        genpart          = cms.InputTag("prunedGenParticles"),
        isQCD            = cms.bool( '/QCD_' in params.inputDataset ),
        slimjet          = cms.InputTag("slimmedJets")
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2021_realistic', '')

# Path and EndPath definitions
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWMINIAODSIMoutput_step = cms.EndPath(process.RAWMINIAODSIMoutput)
process.mmtree_step = cms.EndPath(process.mmtree)

# Schedule definition
# process.schedule imported from cff in HLTrigger.Configuration
process.schedule.extend([process.endjob_step,process.RAWMINIAODSIMoutput_step,process.mmtree_step])
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforPatatrack
from HLTrigger.Configuration.customizeHLTforPatatrack import customizeHLTforPatatrackTriplets 

#call to customisation function customizeHLTforPatatrackTriplets imported from HLTrigger.Configuration.customizeHLTforPatatrack
process = customizeHLTforPatatrackTriplets(process)

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC 

#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforMC(process)

# End of customisation functions


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
