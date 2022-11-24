import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

def scoutingToReco(process):

   process.pfcands = cms.EDProducer(
     "Run3ScoutingToPFCandidateProducer2",
     scoutingparticle=cms.InputTag("hltScoutingPFPacker"),
   )

   process.pfcandsCHS = cms.EDProducer(
     "Run3ScoutingToPFCandidateProducer2",
     scoutingparticle=cms.InputTag("hltScoutingPFPacker"),
     CHS = cms.bool(True)
   )

   process.pfcandsSK = cms.EDProducer(
     "Run3ScoutingToPFCandidateProducer2",
     scoutingparticle=cms.InputTag("hltScoutingPFPacker"),
     softKiller = cms.bool(True)
   )

   #process.ak4Jets = cms.EDProducer(
   #  "Run3ScoutingToRecoJetProducer",
   #  scoutingjet=cms.InputTag("hltScoutingPFPacker"),
   #  recopfcand=cms.InputTag("pfcands"),
   #)

   process.recoTask = cms.Task(
      process.pfcands,
      process.pfcandsCHS,
      process.pfcandsSK,
      #process.ak4Jets,
   )
   process.schedule.associate(process.recoTask)
