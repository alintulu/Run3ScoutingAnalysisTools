import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

def addAK4Jets(process, isMC):
 
   from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
   process.ak4Jets = ak4PFJets.clone(
      src = ("pfcands"),
   )

   process.ak4ParticleNetJetTagInfos = cms.EDProducer("DeepBoostedJetTagInfoProducer",
       jet_radius = cms.double( 0.4 ),
       min_jet_pt = cms.double( 5.0 ),
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
       jets = cms.InputTag( "ak4Jets" ),
       puppi_value_map = cms.InputTag( "" ),
       use_scouting_features = cms.bool( True ),
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
   process.ak4ParticleNetJetTags = cms.EDProducer("BoostedJetONNXValueMapProducer",
       jets = cms.InputTag("ak4Jets"),
       src = cms.InputTag("ak4ParticleNetJetTagInfos"),
       preprocess_json = cms.string("Run3ScoutingAnalysisTools/Models/preprocess_flavourtag.json"),
       model_path = cms.FileInPath("Run3ScoutingAnalysisTools/Models/flavourtag.onnx"),
       flav_names = cms.vstring(["probb", "probbb","probc", "probcc", "probuds", "probg", "probundef"]),
       debugMode = cms.untracked.bool(False),
   )

   process.ak4JetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
       src = cms.InputTag("ak4Jets"),
       name = cms.string("ScoutingJet"),
       cut = cms.string(""),
       doc = cms.string("ScoutingJet"),
       singleton = cms.bool(False),
       extension = cms.bool(False), # this is the main table
       externalVariables = cms.PSet(
          particleNet_prob_b = ExtVar(cms.InputTag('ak4ParticleNetJetTags:probb'), float, doc="ParticleNet prob b", precision=10),
          particleNet_prob_bb = ExtVar(cms.InputTag('ak4ParticleNetJetTags:probbb'), float, doc="ParticleNet prob b", precision=10),
          particleNet_prob_c = ExtVar(cms.InputTag('ak4ParticleNetJetTags:probc'), float, doc="ParticleNet prob c", precision=10),
          particleNet_prob_cc = ExtVar(cms.InputTag('ak4ParticleNetJetTags:probcc'), float, doc="ParticleNet prob cc", precision=10),
          particlenet_prob_uds = ExtVar(cms.InputTag('ak4ParticleNetJetTags:probuds'), float, doc="particlenet prob uds", precision=10),
          particleNet_prob_g = ExtVar(cms.InputTag('ak4ParticleNetJetTags:probg'), float, doc="ParticleNet prob g", precision=10),
          particleNet_prob_undef = ExtVar(cms.InputTag('ak4ParticleNetJetTags:probundef'), float, doc="ParticleNet prob undef", precision=10),
       ),
       variables = cms.PSet(
          P4Vars,
          area = Var("jetArea()", float, doc="jet catchment area, for JECs",precision=10),
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

   process.ak4JetTask = cms.Task(
      process.ak4Jets,
      process.ak4ParticleNetJetTagInfos,
      process.ak4ParticleNetJetTags,
      process.ak4JetTable,
   )

   process.schedule.associate(process.ak4JetTask)

   if (isMC):
      process.ak4MatchGen = cms.EDProducer("RecoJetToGenJetDeltaRValueMapProducer2",
          src = cms.InputTag("ak4Jets"),
          matched = cms.InputTag("slimmedGenJets"),
          distMax = cms.double(0.4),
          value = cms.string("index"),
      )
      process.ak4MatchGenTask = cms.Task(process.ak4MatchGen)
      externalVariables = getattr(process.ak4JetTable, 'externalVariables', cms.PSet())
      externalVariables.genJetIdx = ExtVar(cms.InputTag("ak4MatchGen"), int, doc="gen jet idx")
      process.ak4JetTable.externalVariables = externalVariables
      process.schedule.associate(process.ak4MatchGenTask)

def addAK4CHSJets(process, isMC):
 
   from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
   process.ak4JetsCHS = ak4PFJets.clone(
      src = ("pfcandsCHS"),
   )

   process.ak4CHSParticleNetJetTagInfos = cms.EDProducer("DeepBoostedJetTagInfoProducer",
       jet_radius = cms.double( 0.4 ),
       min_jet_pt = cms.double( 5.0 ),
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
       pf_candidates = cms.InputTag( "pfcandsCHS" ),
       jets = cms.InputTag( "ak4JetsCHS" ),
       puppi_value_map = cms.InputTag( "" ),
       use_scouting_features = cms.bool( True ),
       normchi2_value_map = cms.InputTag("pfcandsCHS", "normchi2"),
       dz_value_map = cms.InputTag("pfcandsCHS", "dz"),
       dxy_value_map = cms.InputTag("pfcandsCHS", "dxy"),
       dzsig_value_map = cms.InputTag("pfcandsCHS", "dzsig"),
       dxysig_value_map = cms.InputTag("pfcandsCHS", "dxysig"),
       lostInnerHits_value_map = cms.InputTag("pfcandsCHS", "lostInnerHits"),
       quality_value_map = cms.InputTag("pfcandsCHS", "quality"),
       trkPt_value_map = cms.InputTag("pfcandsCHS", "trkPt"),
       trkEta_value_map = cms.InputTag("pfcandsCHS", "trkEta"),
       trkPhi_value_map = cms.InputTag("pfcandsCHS", "trkPhi"),
   )

   from RecoBTag.ONNXRuntime.boostedJetONNXJetTagsProducer_cfi import boostedJetONNXJetTagsProducer
   process.ak4CHSParticleNetJetTags = cms.EDProducer("BoostedJetONNXValueMapProducer",
       jets = cms.InputTag("ak4JetsCHS"),
       src = cms.InputTag("ak4CHSParticleNetJetTagInfos"),
       preprocess_json = cms.string("Run3ScoutingAnalysisTools/Models/preprocess_flavourtag.json"),
       model_path = cms.FileInPath("Run3ScoutingAnalysisTools/Models/flavourtag.onnx"),
       flav_names = cms.vstring(["probb", "probbb","probc", "probcc", "probuds", "probg", "probundef"]),
       debugMode = cms.untracked.bool(False),
   )

   process.ak4CHSJetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
       src = cms.InputTag("ak4JetsCHS"),
       name = cms.string("ScoutingJetCHS"),
       cut = cms.string(""),
       doc = cms.string("ScoutingJetCHS"),
       singleton = cms.bool(False),
       extension = cms.bool(False), # this is the main table
       externalVariables = cms.PSet(
          particleNet_prob_b = ExtVar(cms.InputTag('ak4CHSParticleNetJetTags:probb'), float, doc="ParticleNet prob b", precision=10),
          particleNet_prob_bb = ExtVar(cms.InputTag('ak4CHSParticleNetJetTags:probbb'), float, doc="ParticleNet prob b", precision=10),
          particleNet_prob_c = ExtVar(cms.InputTag('ak4CHSParticleNetJetTags:probc'), float, doc="ParticleNet prob c", precision=10),
          particleNet_prob_cc = ExtVar(cms.InputTag('ak4CHSParticleNetJetTags:probcc'), float, doc="ParticleNet prob cc", precision=10),
          particlenet_prob_uds = ExtVar(cms.InputTag('ak4CHSParticleNetJetTags:probuds'), float, doc="particlenet prob uds", precision=10),
          particleNet_prob_g = ExtVar(cms.InputTag('ak4CHSParticleNetJetTags:probg'), float, doc="ParticleNet prob g", precision=10),
          particleNet_prob_undef = ExtVar(cms.InputTag('ak4CHSParticleNetJetTags:probundef'), float, doc="ParticleNet prob undef", precision=10),
       ),
       variables = cms.PSet(
          P4Vars,
          area = Var("jetArea()", float, doc="jet catchment area, for JECs",precision=10),
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

   process.ak4CHSJetTask = cms.Task(
      process.pfcandsCHS,
      process.ak4JetsCHS,
      process.ak4CHSParticleNetJetTagInfos,
      process.ak4CHSParticleNetJetTags,
      process.ak4CHSJetTable,
   )

   process.schedule.associate(process.ak4CHSJetTask)

   if (isMC):
      process.ak4CHSMatchGen = cms.EDProducer("RecoJetToGenJetDeltaRValueMapProducer2",
          src = cms.InputTag("ak4JetsCHS"),
          matched = cms.InputTag("slimmedGenJets"),
          distMax = cms.double(0.4),
          value = cms.string("index"),
      )
      process.ak4CHSMatchGenTask = cms.Task(process.ak4CHSMatchGen)
      externalVariables = getattr(process.ak4CHSJetTable, 'externalVariables', cms.PSet())
      externalVariables.genJetIdx = ExtVar(cms.InputTag("ak4CHSMatchGen"), int, doc="gen jet idx")
      process.ak4CHSJetTable.externalVariables = externalVariables
      process.schedule.associate(process.ak4CHSMatchGenTask)

def addAK4SoftKillerJets(process, isMC):
 
   from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
   process.ak4JetsSK = ak4PFJets.clone(
      src = ("pfcandsSK"),
      useExplicitGhosts = cms.bool(True)
   )

   process.ak4SKParticleNetJetTagInfos = cms.EDProducer("DeepBoostedJetTagInfoProducer",
       jet_radius = cms.double( 0.4 ),
       min_jet_pt = cms.double( 5.0 ),
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
       pf_candidates = cms.InputTag( "pfcandsSK" ),
       jets = cms.InputTag( "ak4JetsSK" ),
       puppi_value_map = cms.InputTag( "" ),
       use_scouting_features = cms.bool( True ),
       normchi2_value_map = cms.InputTag("pfcandsSK", "normchi2"),
       dz_value_map = cms.InputTag("pfcandsSK", "dz"),
       dxy_value_map = cms.InputTag("pfcandsSK", "dxy"),
       dzsig_value_map = cms.InputTag("pfcandsSK", "dzsig"),
       dxysig_value_map = cms.InputTag("pfcandsSK", "dxysig"),
       lostInnerHits_value_map = cms.InputTag("pfcandsSK", "lostInnerHits"),
       quality_value_map = cms.InputTag("pfcandsSK", "quality"),
       trkPt_value_map = cms.InputTag("pfcandsSK", "trkPt"),
       trkEta_value_map = cms.InputTag("pfcandsSK", "trkEta"),
       trkPhi_value_map = cms.InputTag("pfcandsSK", "trkPhi"),
   )

   from RecoBTag.ONNXRuntime.boostedJetONNXJetTagsProducer_cfi import boostedJetONNXJetTagsProducer
   process.ak4SKParticleNetJetTags = cms.EDProducer("BoostedJetONNXValueMapProducer",
       jets = cms.InputTag("ak4JetsSK"),
       src = cms.InputTag("ak4SKParticleNetJetTagInfos"),
       preprocess_json = cms.string("Run3ScoutingAnalysisTools/Models/preprocess_flavourtag.json"),
       model_path = cms.FileInPath("Run3ScoutingAnalysisTools/Models/flavourtag.onnx"),
       flav_names = cms.vstring(["probb", "probbb","probc", "probcc", "probuds", "probg", "probundef"]),
       debugMode = cms.untracked.bool(False),
   )

   process.ak4SKJetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
       src = cms.InputTag("ak4JetsSK"),
       name = cms.string("ScoutingJetSK"),
       cut = cms.string(""),
       doc = cms.string("ScoutingJetSK"),
       singleton = cms.bool(False),
       extension = cms.bool(False), # this is the main table
       externalVariables = cms.PSet(
          particleNet_prob_b = ExtVar(cms.InputTag('ak4SKParticleNetJetTags:probb'), float, doc="ParticleNet prob b", precision=10),
          particleNet_prob_bb = ExtVar(cms.InputTag('ak4SKParticleNetJetTags:probbb'), float, doc="ParticleNet prob b", precision=10),
          particleNet_prob_c = ExtVar(cms.InputTag('ak4SKParticleNetJetTags:probc'), float, doc="ParticleNet prob c", precision=10),
          particleNet_prob_cc = ExtVar(cms.InputTag('ak4SKParticleNetJetTags:probcc'), float, doc="ParticleNet prob cc", precision=10),
          particlenet_prob_uds = ExtVar(cms.InputTag('ak4SKParticleNetJetTags:probuds'), float, doc="particlenet prob uds", precision=10),
          particleNet_prob_g = ExtVar(cms.InputTag('ak4SKParticleNetJetTags:probg'), float, doc="ParticleNet prob g", precision=10),
          particleNet_prob_undef = ExtVar(cms.InputTag('ak4SKParticleNetJetTags:probundef'), float, doc="ParticleNet prob undef", precision=10),
       ),
       variables = cms.PSet(
          P4Vars,
          area = Var("jetArea()", float, doc="jet catchment area, for JECs",precision=10),
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

   process.ak4SKJetTask = cms.Task(
      process.pfcandsSK,
      process.ak4JetsSK,
      process.ak4SKParticleNetJetTagInfos,
      process.ak4SKParticleNetJetTags,
      process.ak4SKJetTable,
   )

   process.schedule.associate(process.ak4SKJetTask)

   if (isMC):
      process.ak4SKMatchGen = cms.EDProducer("RecoJetToGenJetDeltaRValueMapProducer2",
          src = cms.InputTag("ak4JetsSK"),
          matched = cms.InputTag("slimmedGenJets"),
          distMax = cms.double(0.4),
          value = cms.string("index"),
      )
      process.ak4SKMatchGenTask = cms.Task(process.ak4SKMatchGen)
      externalVariables = getattr(process.ak4SKJetTable, 'externalVariables', cms.PSet())
      externalVariables.genJetIdx = ExtVar(cms.InputTag("ak4SKMatchGen"), int, doc="gen jet idx")
      process.ak4SKJetTable.externalVariables = externalVariables
      process.schedule.associate(process.ak4SKMatchGenTask)
