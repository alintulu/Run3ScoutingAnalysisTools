import FWCore.ParameterSet.Config as cms

# Set parameters externally
from FWCore.ParameterSet.VarParsing import VarParsing
params = VarParsing('analysis')

# Define the process
process = cms.Process("LL")

# Parse command line arguments
params.parseArguments()

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5

# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary      = cms.untracked.bool(True),
    #SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# How many events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(params.maxEvents))

# Input EDM files
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(params.inputFiles)
)

# Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

# Define the services needed for the treemaker
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(params.outputFile)
)

# Make tree
process.mmtree = cms.EDAnalyzer('ScoutingNanoAOD',
        pfcandsParticleNet = cms.InputTag("hltScoutingPFPacker"),
        genpart          = cms.InputTag("prunedGenParticles")
)
process.p = cms.Path(process.mmtree)
