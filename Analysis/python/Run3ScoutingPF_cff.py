import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

def scoutingToReco(process):

   process.pfcands = cms.EDProducer(
     "Run3ScoutingToPFCandidateProducer",
     #"Run3ScoutingToPackedCandidateProducer",
     scoutingparticle=cms.InputTag("hltScoutingPFPacker"),
     #scoutingvertex=cms.InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx"),
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
       externalVariables = cms.PSet(
          normchi2 = ExtVar(cms.InputTag("pfcands", "normchi2"), float, doc="normchi2 of best track", precision=6),
          dz = ExtVar(cms.InputTag("pfcands", "dz"), float, doc="dz of best track", precision=6),
          dxy = ExtVar(cms.InputTag("pfcands", "dxy"), float, doc="dxy of best track", precision=6),
          dzsig = ExtVar(cms.InputTag("pfcands", "dzsig"), float, doc="dzsig of best track", precision=6),
          dxysig = ExtVar(cms.InputTag("pfcands", "dxysig"), float, doc="dxysig of best track", precision=6),
          lostInnerHits = ExtVar(cms.InputTag("pfcands", "lostInnerHits"), int, doc="lostInnerHits of best track"),
          quality = ExtVar(cms.InputTag("pfcands", "quality"), int, doc="quality of best track"),
          trkPt = ExtVar(cms.InputTag("pfcands", "trkPt"), float, doc="pt of best track", precision=6),
          trkEta = ExtVar(cms.InputTag("pfcands", "trkEta"), float, doc="eta of best track", precision=6),
          trkPhi = ExtVar(cms.InputTag("pfcands", "trkPhi"), float, doc="phi of best track", precision=6),
       ),
       variables = cms.PSet(
          CandVars,
       ),
   )

def addAK4Jets(process):

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

   process.ak4MatchGenTable = cms.EDProducer("MatchJetToGenJetTableProducer",
       jets = cms.InputTag("ak4Jets"),
       genjets = cms.InputTag("slimmedGenJets"),
       nameTable = cms.string("ScoutingJet"),
   )

def addAK8Jets(process):

   from RecoJets.JetProducers.ak8PFJets_cfi import ak8PFJets

   process.ak8Jets = ak8PFJets.clone(
      src = ("pfcands"),
      useSoftDrop = cms.bool(True),
      zcut = cms.double(0.1),
      beta = cms.double(0.0),
      R0   = cms.double(0.8),
      useExplicitGhosts = cms.bool(True),
      writeCompound = cms.bool(True),
      jetCollInstanceName=cms.string("SubJets"),
      jetPtMin = 170.0
   )
 
   process.ak8JetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
       src = cms.InputTag("ak8Jets", "SubJets"),
       name = cms.string("ScoutingFatJet"),
       cut = cms.string(""),
       doc = cms.string("ScoutingFatJet"),
       singleton = cms.bool(False),
       extension = cms.bool(False), # this is the main table
       variables = cms.PSet(
          P4Vars,
       ),
   )

   process.AK8ParticleNetJetTagInfos = cms.EDProducer("DeepBoostedJetTagInfoScoutingProducer",
       jet_radius = cms.double( 0.8 ),
       min_jet_pt = cms.double( 170.0 ),
       max_jet_eta = cms.double( 2.5 ),
       min_pt_for_track_properties = cms.double( 0.95 ),
       min_pt_for_pfcandidates = cms.double( 0.1 ),
       use_puppiP4 = cms.bool( False ),
       include_neutrals = cms.bool( True ),
       sort_by_sip2dsig = cms.bool( False ),
       min_puppi_wgt = cms.double( -1.0 ),
       flip_ip_sign = cms.bool( False ),
       sip3dSigMax = cms.double( -1.0 ),
       use_hlt_features = cms.bool( False ),
       pf_candidates = cms.InputTag( "pfcands" ),
       jets = cms.InputTag( "ak8Jets" ),
       puppi_value_map = cms.InputTag( "" ),
       normchi2_value_map = cms.InputTag("pfcands", "normchi2"),
       dz_value_map = cms.InputTag("pfcands", "dz"),
       dxy_value_map = cms.InputTag("pfcands", "dxy"),
       dzsig_value_map = cms.InputTag("pfcands", "dzsig"),
       dxysig_value_map = cms.InputTag("pfcands", "dxysig"),
       lostInnerHits_value_map = cms.InputTag("pfcands", "lostInnerHits"),
       quality_value_map = cms.InputTag("pfcands", "quality"),
       trkPt_value_map = cms.InputTag("pfcands", "trkPt"),
       trkEta_value_map = cms.InputTag("pfcands", "trkEta"),
       trkPhi_value_map = cms.InputTag("pfcands", "trkPhi"),
   )

   process.AK8ParticleNetONNXJetTags = cms.EDProducer( "BoostedJetONNXJetTagsProducer",
    src = cms.InputTag( "AK8ParticleNetJetTagInfos" ),
    preprocess_json = cms.string( "RecoBTag/Combined/data/HLT/ParticleNetAK4/V00/preprocess-with-tauh.json" ),
    preprocessParams = cms.PSet(  ),
    model_path = cms.FileInPath( "RecoBTag/Combined/data/HLT/ParticleNetAK4/V00/particle-net-with-tauh.onnx" ),
    flav_names = cms.vstring( 'probtauh',
      'probb',
      'probc',
      'probuds',
      'probg' ),
    debugMode = cms.untracked.bool( False )
)

   process.ak8MatchGenTable = cms.EDProducer("MatchJetToGenJetTableProducer",
       jets = cms.InputTag("ak8Jets"),
       genjets = cms.InputTag("slimmedAK8GenJets"),
       nameTable = cms.string("ScoutingFatJet"),
   ) 
