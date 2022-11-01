import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

def scoutingToReco(process):

   process.pfcands = cms.EDProducer(
     "Run3ScoutingToPFCandidateProducer",
     scoutingparticle=cms.InputTag("hltScoutingPFPacker"),
   )

   process.recojets = cms.EDProducer(
     "Run3ScoutingToRecoJetProducer",
     scoutingjet=cms.InputTag("hltScoutingPFPacker"),
     scoutingparticle=cms.InputTag("pfcands"),
   )

def addParticles(process):

   process.particleTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
       src = cms.InputTag("pfcands"),
       name = cms.string("Run3ScoutingParticle"),
       cut = cms.string(""),
       doc = cms.string("Run3ScoutingParticle"),
       singleton = cms.bool(False),
       extension = cms.bool(False), # this is the main table
       variables = cms.PSet(
          CandVars,
       ),
   )

def addJets(process):

   process.jetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
       src = cms.InputTag("recojets"),
       name = cms.string("Run3ScoutingJet"),
       cut = cms.string(""),
       doc = cms.string("Run3ScoutingJet"),
       singleton = cms.bool(False),
       extension = cms.bool(False), # this is the main table
       variables = cms.PSet(
          P3Vars,
          chHEF = Var("chargedHadronEnergy()/(chargedHadronEnergy()+neutralHadronEnergy()+photonEnergy()+electronEnergy()+muonEnergy())", float, doc="charged Hadron Energy Fraction", precision= 6),
          neHEF = Var("neutralHadronEnergy()/(chargedHadronEnergy()+neutralHadronEnergy()+photonEnergy()+electronEnergy()+muonEnergy())", float, doc="neutral Hadron Energy Fraction", precision= 6),
          chEmEF = Var("(electronEnergy()+muonEnergy())/(chargedHadronEnergy()+neutralHadronEnergy()+photonEnergy()+electronEnergy()+muonEnergy())", float, doc="charged Electromagnetic Energy Fraction", precision= 6),
          neEmEF = Var("(photonEnergy())/(chargedHadronEnergy()+neutralHadronEnergy()+photonEnergy()+electronEnergy()+muonEnergy())", float, doc="neutral Electromagnetic Energy Fraction", precision= 6),
          muEmEF = Var("(muonEnergy())/(chargedHadronEnergy()+neutralHadronEnergy()+photonEnergy()+electronEnergy()+muonEnergy())", float, doc="muon Energy Fraction", precision= 6),
          nCh = Var("chargedHadronMultiplicity()", int, doc="number of charged hadrons in the jet"),
          nNh = Var("neutralHadronMultiplicity()", int, doc="number of neutral hadrons in the jet"),
          nMuons = Var("muonMultiplicity()", int, doc="number of muons in the jet"),
          nElectrons = Var("electronMultiplicity()", int, doc="number of electrons in the jet"),
          nPhotons = Var("photonMultiplicity()", int, doc="number of photons in the jet"),
          nConstituents = Var("numberOfDaughters()", "uint8", doc="Number of particles in the jet")
       ),
   )
