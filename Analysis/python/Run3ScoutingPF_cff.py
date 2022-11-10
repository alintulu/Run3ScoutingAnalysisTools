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
     recopfcand=cms.InputTag("pfcands"),
   )

   process.recoTask = cms.Task(
      process.pfcands,
      process.ak4Jets,
   )
   process.schedule.associate(process.recoTask)


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

def addAK4Jets(process):

   process.ak4JetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
       src = cms.InputTag("ak4Jets"),
       name = cms.string("ScoutingJet"),
       cut = cms.string(""),
       doc = cms.string("ScoutingJet"),
       singleton = cms.bool(False),
       extension = cms.bool(False), # this is the main table
       externalVariables = cms.PSet(
          mass = ExtVar(cms.InputTag("ak4Jets", "mass"), float, doc="mass", precision=6),
          area = ExtVar(cms.InputTag("ak4Jets", "jetArea"), float, doc="area", precision=6),
       ),
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

   process.ak4JetTask = cms.Task(process.ak4JetTable)
   process.schedule.associate(process.ak4JetTask)

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
      jetCollInstanceName=cms.string("SoftDrop"),
      jetPtMin = 170.0
   )
 
   process.ak8JetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
       src = cms.InputTag("ak8Jets", "SoftDrop"),
       name = cms.string("ScoutingFatJet"),
       cut = cms.string(""),
       doc = cms.string("ScoutingFatJet"),
       singleton = cms.bool(False),
       extension = cms.bool(False), # this is the main table
       variables = cms.PSet(
          P4Vars,
       ),
   )

   process.pfParticleNetMassRegressionJetTagInfos = cms.EDProducer("DeepBoostedJetTagInfoScoutingProducer",
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

   from RecoBTag.ONNXRuntime.boostedJetONNXJetTagsProducer_cfi import boostedJetONNXJetTagsProducer

   process.pfParticleNetJetTags = boostedJetONNXJetTagsProducer.clone(
       src = cms.InputTag("pfParticleNetMassRegressionJetTagInfos"),
       preprocess_json = cms.string("Run3ScoutingAnalysisTools/Models/preprocess_doublebtag.json"),
       model_path = cms.FileInPath("Run3ScoutingAnalysisTools/Models/doublebtag.onnx"),
       flav_names = ["probHbb", "probHcc","probHqq", "probQCDall"],
       debugMode = cms.untracked.bool(True),
   )

   process.pfParticleNetDiscriminatorsJetTags = cms.EDProducer("BTagProbabilityToDiscriminator",
       discriminators = cms.VPSet(
           cms.PSet(name = cms.string("HbbVsQCD"),
             numerator = cms.VInputTag('pfParticleNetMassRegressionJetTagInfos:probHbb'),
             denominator = cms.VInputTag('pfParticleNetMassRegressionJetTagInfos:probHbb','pfParticleNetMassRegressionJetTagInfos:probQCDall')
           ),
       )
   )

   #setattr(process.ak8JetTable.variables, "ParticleNetHbbvsQCD", Var("bDiscriminator('pfParticleNetDiscriminatorsJetTags:HbbVsQCD')", float, doc="ParticleNet tagger H(->bb) vs QCD discriminator", precision=10),)
   #setattr(process.ak8JetTable.variables, "ParticleNetHbbvsQCD", Var(cms.InputTag("pfParticleNetDiscriminatorsJetTags:HbbVsQCD"), float, doc="ParticleNet tagger H(->bb) vs QCD discriminator", precision=10),)

   process.pfParticleNetMassRegressionJetTags = boostedJetONNXJetTagsProducer.clone(
       src = cms.InputTag("pfParticleNetMassRegressionJetTagInfos"),
       preprocess_json = cms.string("Run3ScoutingAnalysisTools/Models/preprocess_massreg.json"),
       model_path = cms.FileInPath("Run3ScoutingAnalysisTools/Models/massreg.onnx"),
       flav_names = cms.vstring( ["mass"]),
       debugMode = cms.untracked.bool(True),
   )

   #setattr(process.ak8JetTable.variables, "ParticleNetMass",  Var("bDiscriminator('pfParticleNetMassRegressionJetTags:mass')", float, doc="ParticleNet mass regression",precision=10))

   process.ak8MatchGenTable = cms.EDProducer("MatchJetToGenJetTableProducer",
       jets = cms.InputTag("ak8Jets"),
       genjets = cms.InputTag("slimmedAK8GenJets"),
       nameTable = cms.string("ScoutingFatJet"),
   )
 
   process.ak8JetTask = cms.Task(
       process.ak8Jets,
       process.ak8JetTable,
       #process.pfParticleNetMassRegressionJetTagInfos,
       #process.pfParticleNetJetTags,
       #process.pfParticleNetDiscriminatorsJetTags,
       #process.pfParticleNetMassRegressionJetTags
   )
   process.schedule.associate(process.ak8JetTask)
