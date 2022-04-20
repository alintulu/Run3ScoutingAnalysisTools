import FWCore.ParameterSet.Config as cms

# Set parameters externally
from FWCore.ParameterSet.VarParsing import VarParsing
params = VarParsing('analysis')

params.register(
    'GlobalTagMC',
    '122X_mcRun3_2021_realistic_v5',
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'Process name for the HLT paths'
)

params.maxEvents = -1
params.outputFile = "nano.root"

params.register('inputDataset',
    '',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Input dataset"
)

# Define the process
process = cms.Process("LL")

# Parse command line arguments
params.parseArguments()

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 5

# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary      = cms.untracked.bool(True),
    #SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# How many events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# Input EDM files
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(params.inputFiles)
)

# Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

##--- l1 stage2 digis ---
process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
process.gtStage2Digis.InputLabel = cms.InputTag( "hltFEDSelectorL1" )
process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')

# Load the global tag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = params.GlobalTagMC

# Define the services needed for the treemaker
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(params.outputFile)
)

#from DarkPhotonAnalysis.DimuonAnalysis2018.TriggerPaths_cfi import getL1Conf
L1Info = ['L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7', 'L1_DoubleMu_12_5','L1_DoubleMu_15_7','L1_TripleMu_5_3_3','L1_TripleMu_5_5_3','L1_QuadMu0','L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4','L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18','L1_DoubleMu4_SQ_OS_dR_Max1p2','L1_SingleMu22','L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4','L1_DoubleMu4p5_SQ_OS_dR_Max1p2','L1_DoubleMu4p5_SQ_OS','L1_DoubleMu0er1p5_SQ_dR_Max1p4','L1_DoubleMu0er2p0_SQ_dR_Max1p4','L1_DoubleMu0_SQ']

process.genParticlesMerged = cms.EDProducer("MergedGenParticleProducer",
    inputPruned = cms.InputTag("prunedGenParticles"),
    inputPacked = cms.InputTag("packedGenParticles")
)
# Make tree
process.mmtree = cms.EDAnalyzer('ScoutingNanoAOD',

    	triggerresults   = cms.InputTag("TriggerResults", "", "HLT"),
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
	    #vertices         = cms.InputTag("hltScoutingMuonPacker","displacedVtx"),
        muons            = cms.InputTag("hltScoutingMuonPacker"),
        electrons        = cms.InputTag("hltScoutingEgammaPacker"),
        photons          = cms.InputTag("hltScoutingEgammaPacker"),
        pfcands          = cms.InputTag("hltScoutingPFPacker"),
        pfjets           = cms.InputTag("hltScoutingPFPacker"),
        tracks           = cms.InputTag("hltScoutingTrackPacker"),
        pfcandsParticleNet = cms.InputTag("hltScoutingPFPacker"),
        genpart          = cms.InputTag("prunedGenParticles"),
        # NEW
        slimjet          = cms.InputTag("slimmedJets"),
        isQCD            = cms.bool( '/QCD_' in params.inputDataset )
        # genParticles = cms.InputTag("genParticlesMerged")
    	#pileupinfo       = cms.InputTag("addPileupInfo"),
    	#geneventinfo     = cms.InputTag("generator"),

)
process.p = cms.Path(                  process.mmtree)
