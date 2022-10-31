import FWCore.ParameterSet.Config as cms
from RecoJets.Configuration.RecoPFJets_cff import ak4PFJetsCHS

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

#process.TFileService = cms.Service("TFileService",
#    fileName = cms.string(params.outputFile)
#)
#
#process.hltFixedGridRhoFastjetAll = cms.EDProducer("FixedGridRhoProducerFastjet",
#    gridSpacing = cms.double(0.55),
#    maxRapidity = cms.double(5.0),
#    pfCandidatesTag = cms.InputTag("hltParticleFlow")
#)
#
#process.mmtree = cms.EDAnalyzer('ScoutingNanoAOD',
#        scoutpart = cms.InputTag("hltScoutingPFPacker"),
#        ak8genjet = cms.InputTag('ak8GenJets'),
#        ak8hltjet = cms.InputTag("hltAK8PFJets"),
#        ak8hltcorrjet = cms.InputTag("hltAK8PFJetsCorrected"),
#        scoutrho = cms.InputTag("hltScoutingPFPacker","rho"),
#        hltrho = cms.InputTag("hltFixedGridRhoFastjetAll"),
#)
#
#process.recoJets = cms.EDProducer(
#  "Run3ScoutingToRecoJetProducer",
#  scoutingjet=cms.InputTag("hltScoutingPFPacker"),
#  scoutingparticle=cms.InputTag("hltScoutingPFPacker"),
#)
#
#process.pfcands = cms.EDProducer(
#  "Run3ScoutingToPFCandidateProducer",
#  scoutingparticle=cms.InputTag("hltScoutingPFPacker"),
#)

process.packs = cms.EDProducer(
  "Run3ScoutingToPackedCandidateProducer",
  scoutingparticle=cms.InputTag("hltScoutingPFPacker"),
  scoutingvertex=cms.InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx")
)

process.chs = cms.EDFilter(
  'CandPtrSelector',
  src = cms.InputTag('packs'),
  cut = cms.string('fromPV')
)

process.ak4 = ak4PFJetsCHS.clone( 
   src = cms.InputTag('chs'),
   doAreaFastjet = True,
   rParam = 0.4,
   jetAlgorithm = 'AntiKt'
)

process.p = cms.Path(process.packs*process.chs*process.ak4)

process.outputFile = cms.OutputModule('PoolOutputModule',
  fileName = cms.untracked.string('output.root'),
 )

process.endPath = cms.EndPath(process.outputFile)
