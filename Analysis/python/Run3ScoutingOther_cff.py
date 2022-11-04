import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

def addScouting(process):

   process.run3ScoutingTable = cms.EDProducer("Run3ScoutingTableProducer",
       vertex = cms.InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx"),
       photon = cms.InputTag("hltScoutingEgammaPacker"),
       muon = cms.InputTag("hltScoutingMuonPacker"),
       electron = cms.InputTag("hltScoutingEgammaPacker"),
       track = cms.InputTag("hltScoutingTrackPacker"),
       metpt = cms.InputTag("hltScoutingPFPacker","pfMetPt"),
       metphi = cms.InputTag("hltScoutingPFPacker","pfMetPhi"),
       rho = cms.InputTag("hltScoutingPFPacker", "rho"),
   )
