import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

def addAK8Jets(process, isMC):

   from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
   process.ak8Jets = ak4PFJets.clone(
      src = ("pfcands"),
      rParam   = 0.8,
      jetPtMin = 50.0,
   )

   process.ak8JetsSoftDrop = ak4PFJets.clone(
      src = ("pfcands"),
      rParam   = 0.8,
      jetPtMin = 50.0,
      useSoftDrop = cms.bool(True),
      zcut = cms.double(0.1),
      beta = cms.double(0.0),
      R0   = cms.double(0.8),
      useExplicitGhosts = cms.bool(True),
      writeCompound = cms.bool(True),
      jetCollInstanceName=cms.string("SubJets"),
   )

   process.ak8JetsSoftDropMass = cms.EDProducer("RecoJetDeltaRValueMapProducer",
      src = cms.InputTag("ak8Jets"),
      matched = cms.InputTag("ak8JetsSoftDrop"),                                         
      distMax = cms.double(0.8),
      value = cms.string('mass')  
   )

   from RecoJets.JetProducers.ECF_cff import ecfNbeta1
   process.ecfNbeta1 = ecfNbeta1.clone(src = cms.InputTag("ak8Jets"), srcWeights="")

   from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness
   process.Njettiness = Njettiness.clone(src = cms.InputTag("ak8Jets"), srcWeights="")

   process.ak8ParticleNetJetTagInfos = cms.EDProducer("DeepBoostedJetTagInfoProducer",
       jet_radius = cms.double( 0.8 ),
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
       jets = cms.InputTag( "ak8Jets" ),
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
   process.ak8ParticleNetJetTags = cms.EDProducer("BoostedJetONNXValueMapProducer",
       jets = cms.InputTag("ak8Jets"),
       src = cms.InputTag("ak8ParticleNetJetTagInfos"),
       preprocess_json = cms.string("Run3ScoutingAnalysisTools/Models/preprocess_doublebtag.json"),
       model_path = cms.FileInPath("Run3ScoutingAnalysisTools/Models/doublebtag.onnx"),
       flav_names = cms.vstring(["probHbb", "probHcc","probHqq", "probQCDall"]),
       debugMode = cms.untracked.bool(False),
   )

   process.ak8ParticleNetMassRegressionJetTags = cms.EDProducer("BoostedJetONNXValueMapProducer",
       jets = cms.InputTag("ak8Jets"), 
       src = cms.InputTag("ak8ParticleNetJetTagInfos"),
       preprocess_json = cms.string("Run3ScoutingAnalysisTools/Models/preprocess_massreg.json"),
       model_path = cms.FileInPath("Run3ScoutingAnalysisTools/Models/massreg.onnx"),
       flav_names = cms.vstring(["mass"]),
       debugMode = cms.untracked.bool(False),
   )

   process.ak8JetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
       src = cms.InputTag("ak8Jets"),
       name = cms.string("ScoutingFatJet"),
       cut = cms.string(""),
       doc = cms.string("ScoutingFatJet"),
       singleton = cms.bool(False),
       extension = cms.bool(False), # this is the main table
       externalVariables = cms.PSet(
          msoftdrop = ExtVar(cms.InputTag('ak8JetsSoftDropMass'), float, doc="Softdrop mass", precision=10),
          n2b1 = ExtVar(cms.InputTag('ecfNbeta1:ecfN2'), float, doc="N2 with beta=1", precision=10),
          n3b1 = ExtVar(cms.InputTag('ecfNbeta1:ecfN3'), float, doc="N3 with beta=1", precision=10),
          tau1 = ExtVar(cms.InputTag('Njettiness:tau1'), float, doc="Nsubjettiness (1 axis)", precision=10),
          tau2 = ExtVar(cms.InputTag('Njettiness:tau2'), float, doc="Nsubjettiness (2 axis)", precision=10),
          tau3 = ExtVar(cms.InputTag('Njettiness:tau3'), float, doc="Nsubjettiness (3 axis)", precision=10),
          tau4 = ExtVar(cms.InputTag('Njettiness:tau4'), float, doc="Nsubjettiness (4 axis)", precision=10),
          particleNet_mass = ExtVar(cms.InputTag('ak8ParticleNetMassRegressionJetTags:mass'), float, doc="ParticleNet regress mass", precision=10),
          particleNet_prob_Hbb = ExtVar(cms.InputTag('ak8ParticleNetJetTags:probHbb'), float, doc="ParticleNet prob Hbb", precision=10),
          particleNet_prob_Hcc = ExtVar(cms.InputTag('ak8ParticleNetJetTags:probHcc'), float, doc="ParticleNet prob Hcc", precision=10),
          particleNet_prob_Hqq = ExtVar(cms.InputTag('ak8ParticleNetJetTags:probHqq'), float, doc="ParticleNet prob Hqq", precision=10),
          particleNet_prob_QCD = ExtVar(cms.InputTag('ak8ParticleNetJetTags:probQCDall'), float, doc="ParticleNet probbQCD", precision=10),
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
 
   process.ak8JetTask = cms.Task(
       process.ak8Jets,
       process.ak8JetsSoftDrop,
       process.ak8JetsSoftDropMass,
       process.ecfNbeta1,
       process.Njettiness,
       process.ak8ParticleNetJetTagInfos,
       process.ak8ParticleNetJetTags,
       process.ak8ParticleNetMassRegressionJetTags,
       process.ak8JetTable,
   )

   process.schedule.associate(process.ak8JetTask)

   if (isMC):
      process.ak8MatchGen = cms.EDProducer("RecoJetToGenJetDeltaRValueMapProducer2",
          src = cms.InputTag("ak8Jets"),
          matched = cms.InputTag("slimmedGenJetsAK8"),
          distMax = cms.double(0.8),
          value = cms.string("index"),
      )
      process.ak8MatchGenTask = cms.Task(process.ak8MatchGen)
      externalVariables = getattr(process.ak8JetTable, 'externalVariables', cms.PSet())
      externalVariables.genJetAK8Idx = ExtVar(cms.InputTag("ak8MatchGen"), int, doc="gen jet idx")
      process.ak8JetTable.externalVariables = externalVariables
      process.schedule.associate(process.ak8MatchGenTask)

def addAK8CHSJets(process, isMC):

   from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
   process.ak8CHSJets = ak4PFJets.clone(
      src = ("pfcandsCHS"),
      rParam   = 0.8,
      jetPtMin = 50.0,
   )

   process.ak8CHSJetsSoftDrop = ak4PFJets.clone(
      src = ("pfcandsCHS"),
      rParam   = 0.8,
      jetPtMin = 50.0,
      useSoftDrop = cms.bool(True),
      zcut = cms.double(0.1),
      beta = cms.double(0.0),
      R0   = cms.double(0.8),
      useExplicitGhosts = cms.bool(True),
      writeCompound = cms.bool(True),
      jetCollInstanceName=cms.string("SubJets"),
   )

   process.ak8CHSJetsSoftDropMass = cms.EDProducer("RecoJetDeltaRValueMapProducer",
      src = cms.InputTag("ak8CHSJets"),
      matched = cms.InputTag("ak8CHSJetsSoftDrop"),                                         
      distMax = cms.double(0.8),
      value = cms.string('mass')  
   )

   from RecoJets.JetProducers.ECF_cff import ecfNbeta1
   process.ecfNbeta1CHS = ecfNbeta1.clone(src = cms.InputTag("ak8CHSJets"), srcWeights="")

   from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness
   process.NjettinessCHS = Njettiness.clone(src = cms.InputTag("ak8CHSJets"), srcWeights="")

   process.ak8CHSParticleNetJetTagInfos = cms.EDProducer("DeepBoostedJetTagInfoProducer",
       jet_radius = cms.double( 0.8 ),
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
       jets = cms.InputTag( "ak8CHSJets" ),
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
   process.ak8CHSParticleNetJetTags = cms.EDProducer("BoostedJetONNXValueMapProducer",
       jets = cms.InputTag("ak8CHSJets"),
       src = cms.InputTag("ak8CHSParticleNetJetTagInfos"),
       preprocess_json = cms.string("Run3ScoutingAnalysisTools/Models/preprocess_doublebtag.json"),
       model_path = cms.FileInPath("Run3ScoutingAnalysisTools/Models/doublebtag.onnx"),
       flav_names = cms.vstring(["probHbb", "probHcc","probHqq", "probQCDall"]),
       debugMode = cms.untracked.bool(False),
   )

   process.ak8CHSParticleNetMassRegressionJetTags = cms.EDProducer("BoostedJetONNXValueMapProducer",
       jets = cms.InputTag("ak8CHSJets"), 
       src = cms.InputTag("ak8CHSParticleNetJetTagInfos"),
       preprocess_json = cms.string("Run3ScoutingAnalysisTools/Models/preprocess_massreg.json"),
       model_path = cms.FileInPath("Run3ScoutingAnalysisTools/Models/massreg.onnx"),
       flav_names = cms.vstring(["mass"]),
       debugMode = cms.untracked.bool(False),
   )

   process.ak8CHSJetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
       src = cms.InputTag("ak8CHSJets"),
       name = cms.string("ScoutingFatJetCHS"),
       cut = cms.string(""),
       doc = cms.string("ScoutingFatJetCHS"),
       singleton = cms.bool(False),
       extension = cms.bool(False), # this is the main table
       externalVariables = cms.PSet(
          msoftdrop = ExtVar(cms.InputTag('ak8CHSJetsSoftDropMass'), float, doc="Softdrop mass", precision=10),
          n2b1 = ExtVar(cms.InputTag('ecfNbeta1CHS:ecfN2'), float, doc="N2 with beta=1", precision=10),
          n3b1 = ExtVar(cms.InputTag('ecfNbeta1CHS:ecfN3'), float, doc="N3 with beta=1", precision=10),
          tau1 = ExtVar(cms.InputTag('NjettinessCHS:tau1'), float, doc="Nsubjettiness (1 axis)", precision=10),
          tau2 = ExtVar(cms.InputTag('NjettinessCHS:tau2'), float, doc="Nsubjettiness (2 axis)", precision=10),
          tau3 = ExtVar(cms.InputTag('NjettinessCHS:tau3'), float, doc="Nsubjettiness (3 axis)", precision=10),
          tau4 = ExtVar(cms.InputTag('NjettinessCHS:tau4'), float, doc="Nsubjettiness (4 axis)", precision=10),
          particleNet_mass = ExtVar(cms.InputTag('ak8CHSParticleNetMassRegressionJetTags:mass'), float, doc="ParticleNet regress mass", precision=10),
          particleNet_prob_Hbb = ExtVar(cms.InputTag('ak8CHSParticleNetJetTags:probHbb'), float, doc="ParticleNet prob Hbb", precision=10),
          particleNet_prob_Hcc = ExtVar(cms.InputTag('ak8CHSParticleNetJetTags:probHcc'), float, doc="ParticleNet prob Hcc", precision=10),
          particleNet_prob_Hqq = ExtVar(cms.InputTag('ak8CHSParticleNetJetTags:probHqq'), float, doc="ParticleNet prob Hqq", precision=10),
          particleNet_prob_QCD = ExtVar(cms.InputTag('ak8CHSParticleNetJetTags:probQCDall'), float, doc="ParticleNet probbQCD", precision=10),
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
 
   process.ak8CHSJetTask = cms.Task(
       process.ak8CHSJets,
       process.ak8CHSJetsSoftDrop,
       process.ak8CHSJetsSoftDropMass,
       process.ecfNbeta1CHS,
       process.NjettinessCHS,
       process.ak8CHSParticleNetJetTagInfos,
       process.ak8CHSParticleNetJetTags,
       process.ak8CHSParticleNetMassRegressionJetTags,
       process.ak8CHSJetTable,
   )

   process.schedule.associate(process.ak8CHSJetTask)

   if (isMC):
      process.ak8CHSMatchGen = cms.EDProducer("RecoJetToGenJetDeltaRValueMapProducer2",
          src = cms.InputTag("ak8CHSJets"),
          matched = cms.InputTag("slimmedGenJetsAK8"),
          distMax = cms.double(0.8),
          value = cms.string("index"),
      )
      process.ak8CHSMatchGenTask = cms.Task(process.ak8CHSMatchGen)
      externalVariables = getattr(process.ak8CHSJetTable, 'externalVariables', cms.PSet())
      externalVariables.genJetAK8Idx = ExtVar(cms.InputTag("ak8CHSMatchGen"), int, doc="gen jet idx")
      process.ak8CHSJetTable.externalVariables = externalVariables
      process.schedule.associate(process.ak8CHSMatchGenTask)

def addAK8SoftKillerJets(process, isMC):

   from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
   process.ak8SKJets = ak4PFJets.clone(
      src = ("pfcandsSK"),
      rParam   = 0.8,
      jetPtMin = 50.0,
      useExplicitGhosts = cms.bool(True)
   )

   process.ak8SKJetsSoftDrop = ak4PFJets.clone(
      src = ("pfcandsSK"),
      rParam   = 0.8,
      jetPtMin = 50.0,
      useSoftDrop = cms.bool(True),
      zcut = cms.double(0.1),
      beta = cms.double(0.0),
      R0   = cms.double(0.8),
      useExplicitGhosts = cms.bool(True),
      writeCompound = cms.bool(True),
      jetCollInstanceName=cms.string("SubJets"),
   )

   process.ak8SKJetsSoftDropMass = cms.EDProducer("RecoJetDeltaRValueMapProducer",
      src = cms.InputTag("ak8SKJets"),
      matched = cms.InputTag("ak8SKJetsSoftDrop"),                                         
      distMax = cms.double(0.8),
      value = cms.string('mass')  
   )

   from RecoJets.JetProducers.ECF_cff import ecfNbeta1
   process.ecfNbeta1SK = ecfNbeta1.clone(src = cms.InputTag("ak8SKJets"), srcWeights="")

   from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness
   process.NjettinessSK = Njettiness.clone(src = cms.InputTag("ak8SKJets"), srcWeights="")

   process.ak8SKParticleNetJetTagInfos = cms.EDProducer("DeepBoostedJetTagInfoProducer",
       jet_radius = cms.double( 0.8 ),
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
       jets = cms.InputTag( "ak8SKJets" ),
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
   process.ak8SKParticleNetJetTags = cms.EDProducer("BoostedJetONNXValueMapProducer",
       jets = cms.InputTag("ak8SKJets"),
       src = cms.InputTag("ak8SKParticleNetJetTagInfos"),
       preprocess_json = cms.string("Run3ScoutingAnalysisTools/Models/preprocess_doublebtag.json"),
       model_path = cms.FileInPath("Run3ScoutingAnalysisTools/Models/doublebtag.onnx"),
       flav_names = cms.vstring(["probHbb", "probHcc","probHqq", "probQCDall"]),
       debugMode = cms.untracked.bool(False),
   )

   process.ak8SKParticleNetMassRegressionJetTags = cms.EDProducer("BoostedJetONNXValueMapProducer",
       jets = cms.InputTag("ak8SKJets"), 
       src = cms.InputTag("ak8SKParticleNetJetTagInfos"),
       preprocess_json = cms.string("Run3ScoutingAnalysisTools/Models/preprocess_massreg.json"),
       model_path = cms.FileInPath("Run3ScoutingAnalysisTools/Models/massreg.onnx"),
       flav_names = cms.vstring(["mass"]),
       debugMode = cms.untracked.bool(False),
   )

   process.ak8SKJetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
       src = cms.InputTag("ak8SKJets"),
       name = cms.string("ScoutingFatJetSK"),
       cut = cms.string(""),
       doc = cms.string("ScoutingFatJetSK"),
       singleton = cms.bool(False),
       extension = cms.bool(False), # this is the main table
       externalVariables = cms.PSet(
          msoftdrop = ExtVar(cms.InputTag('ak8SKJetsSoftDropMass'), float, doc="Softdrop mass", precision=10),
          n2b1 = ExtVar(cms.InputTag('ecfNbeta1SK:ecfN2'), float, doc="N2 with beta=1", precision=10),
          n3b1 = ExtVar(cms.InputTag('ecfNbeta1SK:ecfN3'), float, doc="N3 with beta=1", precision=10),
          tau1 = ExtVar(cms.InputTag('NjettinessSK:tau1'), float, doc="Nsubjettiness (1 axis)", precision=10),
          tau2 = ExtVar(cms.InputTag('NjettinessSK:tau2'), float, doc="Nsubjettiness (2 axis)", precision=10),
          tau3 = ExtVar(cms.InputTag('NjettinessSK:tau3'), float, doc="Nsubjettiness (3 axis)", precision=10),
          tau4 = ExtVar(cms.InputTag('NjettinessSK:tau4'), float, doc="Nsubjettiness (4 axis)", precision=10),
          particleNet_mass = ExtVar(cms.InputTag('ak8SKParticleNetMassRegressionJetTags:mass'), float, doc="ParticleNet regress mass", precision=10),
          particleNet_prob_Hbb = ExtVar(cms.InputTag('ak8SKParticleNetJetTags:probHbb'), float, doc="ParticleNet prob Hbb", precision=10),
          particleNet_prob_Hcc = ExtVar(cms.InputTag('ak8SKParticleNetJetTags:probHcc'), float, doc="ParticleNet prob Hcc", precision=10),
          particleNet_prob_Hqq = ExtVar(cms.InputTag('ak8SKParticleNetJetTags:probHqq'), float, doc="ParticleNet prob Hqq", precision=10),
          particleNet_prob_QCD = ExtVar(cms.InputTag('ak8SKParticleNetJetTags:probQCDall'), float, doc="ParticleNet probbQCD", precision=10),
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
 
   process.ak8SKJetTask = cms.Task(
       process.ak8SKJets,
       process.ak8SKJetsSoftDrop,
       process.ak8SKJetsSoftDropMass,
       process.ecfNbeta1SK,
       process.NjettinessSK,
       process.ak8SKParticleNetJetTagInfos,
       process.ak8SKParticleNetJetTags,
       process.ak8SKParticleNetMassRegressionJetTags,
       process.ak8SKJetTable,
   )

   process.schedule.associate(process.ak8SKJetTask)

   if (isMC):
      process.ak8SKMatchGen = cms.EDProducer("RecoJetToGenJetDeltaRValueMapProducer2",
          src = cms.InputTag("ak8SKJets"),
          matched = cms.InputTag("slimmedGenJetsAK8"),
          distMax = cms.double(0.8),
          value = cms.string("index"),
      )
      process.ak8SKMatchGenTask = cms.Task(process.ak8SKMatchGen)
      externalVariables = getattr(process.ak8SKJetTable, 'externalVariables', cms.PSet())
      externalVariables.genJetAK8Idx = ExtVar(cms.InputTag("ak8SKMatchGen"), int, doc="gen jet idx")
      process.ak8SKJetTable.externalVariables = externalVariables
      process.schedule.associate(process.ak8SKMatchGenTask)

