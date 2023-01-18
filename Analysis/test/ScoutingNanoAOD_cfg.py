import FWCore.ParameterSet.Config as cms

# Set parameters externally
from FWCore.ParameterSet.VarParsing import VarParsing
params = VarParsing('analysis')

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

# Create softdrop groomed GEN jets
from RecoJets.JetProducers.ak8GenJets_cfi import ak8GenJets
process.ak8GenJetsWithNu = ak8GenJets.clone(
    src='genParticles',
    rParam=cms.double(0.8),
    jetPtMin=100.0
)

process.ak8GenJetsWithNuSoftDrop = process.ak8GenJetsWithNu.clone(
    useSoftDrop=cms.bool(True),
    zcut=cms.double(0.1),
    beta=cms.double(0.0),
    R0=cms.double(0.8),
    useExplicitGhosts=cms.bool(True)
)

process.genJetSeq = cms.Sequence(
    process.ak8GenJetsWithNu+
    process.ak8GenJetsWithNuSoftDrop
)

# Get flavour association
process.load('PhysicsTools.NanoAOD.jetMC_cff')
process.patJetPartonsNano.particles = cms.InputTag("genParticles")
process.genJetFlavourAssociation.jets = cms.InputTag("ak8GenJetsWithNuSoftDrop")
process.flavourSeq = cms.Sequence(
   process.patJetPartonsNano+
   process.genJetFlavourAssociation
)

# Make tree
process.mmtree = cms.EDAnalyzer('ScoutingNanoAOD',
        pfcandsParticleNet = cms.InputTag("hltScoutingPFPacker"),
        genpart          = cms.InputTag("genParticles"),
        isQCD            = cms.bool( '/QCD_' in params.inputDataset ),
        ak8genjet        = cms.InputTag('ak8GenJetsWithNuSoftDrop'),
)
process.p = cms.Path(process.genJetSeq*process.flavourSeq*process.mmtree)
