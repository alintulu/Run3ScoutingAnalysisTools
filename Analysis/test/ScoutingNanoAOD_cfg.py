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
	fileNames = cms.untracked.vstring(params.inputFiles),
        inputCommands=cms.untracked.vstring(['drop *', 'keep *_slimmedGenJets_*_*'])
)

# Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(params.outputFile)
)

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
#
#process.vertexs = cms.EDProducer(
#  "Run3ScoutingToVertexReco",
#  scoutingvertex=cms.InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx")
#)
#
#process.packs = cms.EDProducer(
#  "Run3ScoutingToPackedCandidateProducer",
#  scoutingparticle=cms.InputTag("hltScoutingPFPacker"),
#  scoutingvertex=cms.InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx")
#)
#
#process.pv = cms.EDFilter(
#  'CandPtrSelector',
#  src = cms.InputTag('packs'),
#  cut = cms.string('fromPV')
#)
#
#process.chs = ak4PFJetsCHS.clone( 
#   src = cms.InputTag('pv'),
#   doAreaFastjet = True,
#   rParam = 0.4,
#   jetAlgorithm = 'AntiKt'
#)
#
#process.p = cms.Path(process.packs*process.pv*process.chs)

process.load('PhysicsTools.NanoAOD.nanogen_cff')
process.load("PhysicsTools.NanoAOD.NanoAODEDMEventContent_cff")
from PhysicsTools.NanoAOD.taus_cff import *
from PhysicsTools.NanoAOD.jetMC_cff import *
from PhysicsTools.NanoAOD.globals_cff import genTable,genFilterTable
from PhysicsTools.NanoAOD.met_cff import metMCTable
from PhysicsTools.NanoAOD.genparticles_cff import *
from PhysicsTools.NanoAOD.particlelevel_cff import *
from PhysicsTools.NanoAOD.genWeightsTable_cfi import *
from PhysicsTools.NanoAOD.genVertex_cff import *
from PhysicsTools.NanoAOD.common_cff import Var,CandVars

nanoMetadata = cms.EDProducer("UniqueStringProducer",
    strings = cms.PSet(
        tag = cms.string("untagged"),
    )
)

nanogenSequence = cms.Sequence(
    #nanoMetadata+
    #cms.Sequence(particleLevelTask)+
    genJetTable+
    patJetPartonsNano+
    genJetFlavourAssociation+
    genJetFlavourTable+
    genJetAK8Table+
    genJetAK8FlavourAssociation+
    genJetAK8FlavourTable
    #cms.Sequence(genTauTask)+
    #genTable+
    #genFilterTable+
    #cms.Sequence(genParticleTablesTask)+
    #cms.Sequence(genVertexTablesTask)+
    #tautagger+
    #rivetProducerHTXS+
    #cms.Sequence(particleLevelTablesTask)+
    #metMCTable+
    #genWeightsTable
)

cms.Path(process.nanogenSequence)

process.out = cms.OutputModule("NanoAODOutputModule",
    fileName=cms.untracked.string("nano.root"),
    outputCommands = process.NanoAODEDMEventContent.outputCommands,
)

process.end = cms.EndPath(process.out)

#process.outputFile = cms.OutputModule('PoolOutputModule',
#  fileName = cms.untracked.string('output.root'),
# )
#
#process.endPath = cms.EndPath(process.outputFile)
