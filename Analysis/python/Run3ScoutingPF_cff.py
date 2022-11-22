import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

def scoutingToReco(process):

   process.pfcands = cms.EDProducer(
     "Run3ScoutingToPFCandidateProducer",
     scoutingparticle=cms.InputTag("hltScoutingPFPacker"),
   )

   #process.ak4Jets = cms.EDProducer(
   #  "Run3ScoutingToRecoJetProducer",
   #  scoutingjet=cms.InputTag("hltScoutingPFPacker"),
   #  recopfcand=cms.InputTag("pfcands"),
   #)

   process.recoTask = cms.Task(
      process.pfcands,
      #process.ak4Jets,
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

