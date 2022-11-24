import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

def addParticles(process):

   process.particleTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
       src = cms.InputTag("pfcands"),
       name = cms.string("ScoutingParticle"),
       cut = cms.string(""),
       doc = cms.string("ScoutingParticle"),
       singleton = cms.bool(False),
       extension = cms.bool(False), # this is the main table
       externalVariables = cms.PSet(
          trkNormchi2 = ExtVar(cms.InputTag("pfcands", "normchi2"), float, doc="normchi2 of best track", precision=6),
          trkDz = ExtVar(cms.InputTag("pfcands", "dz"), float, doc="dz of best track", precision=6),
          trkDxy = ExtVar(cms.InputTag("pfcands", "dxy"), float, doc="dxy of best track", precision=6),
          trkDzsig = ExtVar(cms.InputTag("pfcands", "dzsig"), float, doc="dzsig of best track", precision=6),
          trkDxysig = ExtVar(cms.InputTag("pfcands", "dxysig"), float, doc="dxysig of best track", precision=6),
          trkLostInnerHits = ExtVar(cms.InputTag("pfcands", "lostInnerHits"), int, doc="lostInnerHits of best track"),
          trkQuality = ExtVar(cms.InputTag("pfcands", "quality"), int, doc="quality of best track"),
          trkPt = ExtVar(cms.InputTag("pfcands", "trkPt"), float, doc="pt of best track", precision=6),
          trkEta = ExtVar(cms.InputTag("pfcands", "trkEta"), float, doc="eta of best track", precision=6),
          trkPhi = ExtVar(cms.InputTag("pfcands", "trkPhi"), float, doc="phi of best track", precision=6),
       ),
       variables = cms.PSet(
          CandVars,
       ),
   )

   process.particleTask = cms.Task(process.particleTable)
   process.schedule.associate(process.particleTask)
