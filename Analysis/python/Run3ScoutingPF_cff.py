import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

def scoutingToReco(process):

   process.pfcands = cms.EDProducer(
     "Run3ScoutingToPFCandidateProducer",
     scoutingparticle=cms.InputTag("hltScoutingPFPacker"),
   )

   process.ak4Jets = cms.EDProducer(
     "Run3ScoutingToRecoJetProducer",
     scoutingjet=cms.InputTag("hltScoutingPFPacker"),
     scoutingparticle=cms.InputTag("pfcands"),
   )

def addParticles(process):

   process.particleTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
       src = cms.InputTag("pfcands"),
       name = cms.string("ScoutingParticle"),
       cut = cms.string(""),
       doc = cms.string("ScoutingParticle"),
       singleton = cms.bool(False),
       extension = cms.bool(False), # this is the main table
       variables = cms.PSet(
          CandVars,
       ),
   )

def addAK4Jets(process):

   process.ak4MatchGen = cms.EDProducer("MatchJetToGenJetProducer",
       jets = cms.InputTag("ak4Jets"),
       genjets = cms.InputTag("slimmedGenJets"),
       nameTable = cms.string("ScoutingJet"),
   ) 

   process.ak4JetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
       src = cms.InputTag("ak4Jets"),
       name = cms.string("ScoutingJet"),
       cut = cms.string(""),
       doc = cms.string("ScoutingJet"),
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

def addAK8Jets(process):

   from RecoJets.JetProducers.ak8PFJets_cfi import ak8PFJets

   process.ak8Jets = ak8PFJets.clone(
      src = ("pfcands"),
   )

   process.ak8MatchGen = cms.EDProducer("MatchJetToGenJetProducer",
       jets = cms.InputTag("ak8Jets"),
       genjets = cms.InputTag("slimmedAK8GenJets"),
       nameTable = cms.string("ScoutingFatJet"),
   ) 
 
   process.ak8JetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
       src = cms.InputTag("ak8Jets"),
       name = cms.string("ScoutingFatJet"),
       cut = cms.string(""),
       doc = cms.string("ScoutingFatJet"),
       singleton = cms.bool(False),
       extension = cms.bool(False), # this is the main table
       variables = cms.PSet(
          P3Vars,
       ),
   )

