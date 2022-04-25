# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: Configuration/WMassNanoGen/python/ZJ_MiNNLO_mWPilotTest_CP5_Photos_cff.py --fileout file:ZJ_MiNNLO_mWPilotTest_CP5_Photos.root --mc --eventcontent NANOAODSIM --datatier NANOAOD --conditions auto:mc --step LHE,GEN,NANOGEN --python_filename configs/ZJ_MiNNLO_mWPilotTest_CP5_Photos_cfg.py --customise_commands process.RandomNumberGeneratorService.externalLHEProducer.initialSeed=999 -n 10 --no_exec
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.jetMC_cff import *

from FWCore.ParameterSet.VarParsing import VarParsing
params = VarParsing('analysis')

params.maxEvents = -1
params.outputFile = "nano.root"

process = cms.Process('NANOGEN')

# Parse command line arguments
params.parseArguments()

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic50ns13TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('PhysicsTools.NanoAOD.nanogen_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(params.maxEvents),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(params.outputFile)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    params.inputFiles 
    )
)

process.NANOAODSIMoutput = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('output.root')
)

process.nanogenSequence = cms.Sequence(
   genJetTable+
   patJetPartonsNano+
   genJetFlavourAssociation
)

process.mmtree = cms.EDAnalyzer('ScoutingNanoAOD',
  pfcandsParticleNet = cms.InputTag("hltScoutingPFPacker"),
  jetFlavourInfos    = cms.InputTag("genJetFlavourAssociation")
)

process.content = cms.EDAnalyzer("EventContentAnalyzer")

# Path and EndPath definitions
process.nanoAOD_step = cms.Path(process.nanogenSequence*process.content*process.mmtree)
#process.NANOAODSIMoutput_step = cms.EndPath(process.NANOAODSIMoutput)

#process.schedule = cms.Schedule(process.nanoAOD_step,process.NANOAODSIMoutput_step)
