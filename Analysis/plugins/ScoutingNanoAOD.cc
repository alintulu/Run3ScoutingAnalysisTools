// Standard C++ includes
#include <memory>
#include <vector>
#include <iostream>

// ROOT includes
#include <TTree.h>
#include <TLorentzVector.h>
#include <TPRegexp.h>

// CMSSW framework includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// CMSSW data formats
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"

// Other relevant CMSSW includes
#include "CommonTools/UtilAlgos/interface/TFileService.h" 
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"


#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/Scouting/interface/Run3ScoutingElectron.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPhoton.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/Scouting/interface/Run3ScoutingTrack.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticleParticleNet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"

#include <DataFormats/TrackReco/interface/TrackBase.h>

#include "DataFormats/Math/interface/libminifloat.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DataFormats/TrackReco/interface/fillCovariance.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/RefToBase.h"

// Root include files
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

// User include files

#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/XConePlugin.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/contrib/RecursiveSoftDrop.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"
#include "fastjet/JadePlugin.hh"
#include "fastjet/contrib/SoftKiller.hh"

#include "DataFormats/BTauReco/interface/TaggingVariable.h"

#include "Run3ScoutingAnalysisTools/Analysis/interface/FatJetMatching.h"

using namespace std;
using namespace deepntuples;

FatJetMatching ak8_match;

class ScoutingNanoAOD : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns, edm::one::WatchLuminosityBlocks> {
public:
  explicit ScoutingNanoAOD(const edm::ParameterSet&);
  ~ScoutingNanoAOD();
		
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	
	
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void clearVars();
  const edm::InputTag triggerResultsTag;
  const edm::EDGetTokenT<std::vector<Run3ScoutingParticleParticleNet> >  	pfcandsParticleNetToken;
  const edm::EDGetTokenT<reco::GenParticleCollection>      genpartsToken;

  std::vector<std::string> triggerPathsVector;
  std::map<std::string, int> triggerPathsMap;

  bool doL1;       
  triggerExpression::Data triggerCache_;
      
  // Generator-level information
  // Flags for the different types of triggers used in the analysis
  // For now we are interested in events passing either the single or double lepton triggers
  unsigned char                trig;
       
  edm::InputTag                algInputTag_;       
  edm::InputTag                extInputTag_;       
  edm::EDGetToken              algToken_;
  //l1t::L1TGlobalUtil          *l1GtUtils_;
  std::unique_ptr<l1t::L1TGlobalUtil> l1GtUtils_;
  std::vector<std::string>     l1Seeds_;
  std::vector<bool>            l1Result_;
        
  // TTree carrying the event weight information
  TTree* tree;

  //PFCand ParticleNet
  vector<Float16_t> pfcand_pt_log_nopuppi;
  vector<Float16_t> pfcand_e_log_nopuppi;
  vector<Float16_t> pfcand_etarel;
  vector<Float16_t> pfcand_phirel;
  vector<Float16_t> pfcand_abseta;
  vector<Float16_t> pfcand_charge;
  vector<Float16_t> pfcand_isEl;
  vector<Float16_t> pfcand_isMu;
  vector<Float16_t> pfcand_isGamma;
  vector<Float16_t> pfcand_isChargedHad;
  vector<Float16_t> pfcand_isNeutralHad;
  vector<Float16_t> pfcand_lostInnerHits;
  vector<Float16_t> pfcand_normchi2;
  vector<Float16_t> pfcand_quality;
  vector<Float16_t> pfcand_dz;
  vector<Float16_t> pfcand_dzsig;
  vector<Float16_t> pfcand_dxy;
  vector<Float16_t> pfcand_dxysig;
  vector<Float16_t> pfcand_btagEtaRel;
  vector<Float16_t> pfcand_btagPtRatio;
  vector<Float16_t> pfcand_btagPParRatio;

  //Jet kinematics
  float fj_pt;
  float fj_eta;
  float fj_phi;
  float fj_mass;

  //ParticleNet Jet label
  int label_Top_bcq;
  int label_Top_bqq;
  int label_Top_bc;
  int label_Top_bq;
  int label_W_cq;
  int label_W_qq;
  int label_Z_bb;
  int label_Z_cc;
  int label_Z_qq;
  int label_H_bb;
  int label_H_cc;
  int label_H_qqqq;
  int label_H_tautau;
  int label_H_qq;
  int label_QCD_all;

  bool isQCD;
};

ScoutingNanoAOD::ScoutingNanoAOD(const edm::ParameterSet& iConfig):
  pfcandsParticleNetToken  (consumes<std::vector<Run3ScoutingParticleParticleNet> > (iConfig.getParameter<edm::InputTag>("pfcandsParticleNet"))),
  genpartsToken            (consumes<reco::GenParticleCollection> (iConfig.getParameter<edm::InputTag>("genpart"))),
  isQCD                    (iConfig.existsAs<bool>("isQCD") ? iConfig.getParameter<bool>("isQCD") : false)
{
  usesResource("TFileService");
  if (doL1) {
   algInputTag_ = iConfig.getParameter<edm::InputTag>("AlgInputTag");
   extInputTag_ = iConfig.getParameter<edm::InputTag>("l1tExtBlkInputTag");
   algToken_ = consumes<BXVector<GlobalAlgBlk>>(algInputTag_);
   l1Seeds_ = iConfig.getParameter<std::vector<std::string> >("l1Seeds");
    /* l1GtUtils_ = new l1t::L1TGlobalUtil(iConfig,consumesCollector());*/	
   l1GtUtils_ = std::make_unique<l1t::L1TGlobalUtil>(
    iConfig, consumesCollector(), *this, algInputTag_, extInputTag_, l1t::UseEventSetupIn::Event);
  }
  else {
    l1Seeds_ = std::vector<std::string>();
    l1GtUtils_ = 0;
  }

 // Access the TFileService
  edm::Service<TFileService> fs;

  // Create the TTree
  tree = fs->make<TTree>("tree", "tree");

  tree->Branch("pfcand_pt_log_nopuppi", &pfcand_pt_log_nopuppi);
  tree->Branch("pfcand_e_log_nopuppi", &pfcand_e_log_nopuppi);
  tree->Branch("pfcand_etarel", &pfcand_etarel);
  tree->Branch("pfcand_phirel", &pfcand_phirel);
  tree->Branch("pfcand_abseta", &pfcand_abseta);
  tree->Branch("pfcand_charge", &pfcand_charge);
  tree->Branch("pfcand_isEl", &pfcand_isEl);
  tree->Branch("pfcand_isMu", &pfcand_isMu);
  tree->Branch("pfcand_isGamma", &pfcand_isGamma);
  tree->Branch("pfcand_isChargedHad", &pfcand_isChargedHad);
  tree->Branch("pfcand_isNeutralHad", &pfcand_isNeutralHad);
  tree->Branch("pfcand_lostInnerHits", &pfcand_lostInnerHits);
  tree->Branch("pfcand_normchi2", &pfcand_normchi2);
  tree->Branch("pfcand_quality", &pfcand_quality);
  tree->Branch("pfcand_dz", &pfcand_dz);
  tree->Branch("pfcand_dzsig", &pfcand_dzsig);
  tree->Branch("pfcand_dxy", &pfcand_dxy);
  tree->Branch("pfcand_dxysig", &pfcand_dxysig);
  tree->Branch("pfcand_btagEtaRel", &pfcand_btagEtaRel);
  tree->Branch("pfcand_btagPtRatio", &pfcand_btagPtRatio);
  tree->Branch("pfcand_btagPParRatio", &pfcand_btagPParRatio);

  tree->Branch("fj_pt", &fj_pt);
  tree->Branch("fj_eta", &fj_eta);
  tree->Branch("fj_phi", &fj_phi);
  tree->Branch("fj_mass", &fj_mass);
  tree->Branch("fj_msd", &fj_msd);
  tree->Branch("fj_n2b1", &fj_n2b1);

  tree->Branch("label_Top_bcq", &label_Top_bcq);
  tree->Branch("label_Top_bqq", &label_Top_bqq);
  tree->Branch("label_Top_bc", &label_Top_bc);
  tree->Branch("label_Top_bq", &label_Top_bq);
  tree->Branch("label_W_cq", &label_W_cq);
  tree->Branch("label_W_qq", &label_W_qq);
  tree->Branch("label_Z_bb", &label_Z_bb);
  tree->Branch("label_Z_cc", &label_Z_cc);
  tree->Branch("label_Z_qq", &label_Z_qq);
  tree->Branch("label_H_bb", &label_H_bb);
  tree->Branch("label_H_cc", &label_H_cc);
  tree->Branch("label_H_qqqq", &label_H_qqqq);
  tree->Branch("label_H_tautau", &label_H_tautau);
  tree->Branch("label_H_qq", &label_H_qq);
  tree->Branch("label_QCD_all", &label_QCD_all);
}

ScoutingNanoAOD::~ScoutingNanoAOD() {
}

void ScoutingNanoAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace fastjet;
  using namespace fastjet::contrib;

  Handle<vector<Run3ScoutingParticleParticleNet> > pfcandsParticleNetH;
  iEvent.getByToken(pfcandsParticleNetToken, pfcandsParticleNetH);

  Handle<GenParticleCollection> genpartH;
  iEvent.getByToken(genpartsToken, genpartH);

  // Create AK8 Jet
  vector<PseudoJet> fj_part;
  fj_part.reserve(pfcandsParticleNetH->size());
  int pfcand_i = 0;
  for (auto pfcands_iter = pfcandsParticleNetH->begin(); pfcands_iter != pfcandsParticleNetH->end(); ++pfcands_iter) {
    math::PtEtaPhiMLorentzVector p4(pfcands_iter->pt(), pfcands_iter->eta(), pfcands_iter->phi(), pfcands_iter->m());
    fj_part.emplace_back(p4.px(), p4.py(), p4.pz(), p4.energy());
    fj_part.back().set_user_index(pfcand_i);
    pfcand_i++;
  }

  JetDefinition ak8_def = JetDefinition(antikt_algorithm, 0.8);
  fastjet::GhostedAreaSpec area_spec(5.0,1,0.01);
  fastjet::AreaDefinition area_def(fastjet::active_area, area_spec);

  // Substructure
  double sd_z_cut = 0.10;
  double sd_beta = 0;
  SoftDrop sd_groomer = SoftDrop(sd_z_cut, sd_beta, 1.0);
  EnergyCorrelatorN2 N2 = EnergyCorrelatorN2(1.0);

  ClusterSequenceArea ak8_cs(fj_part, ak8_def, area_def);
  vector<PseudoJet> ak8_jets = sorted_by_pt(ak8_cs.inclusive_jets(170.0));

  cout << "Number of jets: " << ak8_jets.size() << endl;

  for(auto &j: ak8_jets) {

    // Match AK8 jet to truth label
    auto ak8_label = ak8_match.flavorLabel(j, *genpartH, 0.8);
    cout << "Label: " << ak8_label.first << endl;
    if ((ak8_label.first == FatJetMatching::QCD_all && !isQCD) || (ak8_label.first != FatJetMatching::QCD_all && isQCD)) continue;

    float etasign = j.eta() > 0 ? 1 : -1;

    // The following is needed to compute btagEtaRel, btagPtRatio and btagPParRatio
    float jet_px = j.pt() * cos(j.phi());
    float jet_py = j.pt() * sin(j.phi());
    float jet_pz = j.pt() * sinh(j.eta());
    math::XYZVector jet_dir_temp(jet_px, jet_py, jet_pz);
    math::XYZVector jet_dir = jet_dir_temp.Unit();
    TVector3 jet_dir3(jet_px, jet_py, jet_pz);

    // Loop over AK8 jet constituents
    const vector<PseudoJet> constituents = j.constituents();
    for (auto &cand : constituents) {
      // Match PseudoJet constituent to PF candidate
      auto *reco_cand = dynamic_cast<const Run3ScoutingParticleParticleNet*> (&pfcandsParticleNetH->at(cand.user_index()));
      // The following is needed to compute btagEtaRel, btagPtRatio and btagPParRatio
      float trk_px = reco_cand->trk_pt() * cos(reco_cand->trk_phi());
      float trk_py = reco_cand->trk_pt() * sin(reco_cand->trk_phi());
      float trk_pz = reco_cand->trk_pt() * sinh(reco_cand->trk_eta());
      math::XYZVector track_mom(trk_px, trk_py, trk_pz);
      TVector3 track_mom3(trk_px, trk_py, trk_pz);
      double track_mag = sqrt(trk_px * trk_px + trk_py * trk_py + trk_pz * trk_pz);

      float reco_cand_p = reco_cand->pt() * cosh(reco_cand->eta());
      pfcand_pt_log_nopuppi.push_back(log(reco_cand->pt()));
      pfcand_e_log_nopuppi.push_back(log(sqrt(reco_cand_p*reco_cand_p + reco_cand->m()*reco_cand->m())));
      pfcand_etarel.push_back(etasign * (reco_cand->eta() - j.eta()));
      pfcand_phirel.push_back(deltaPhi(reco_cand->phi(), j.phi()));
      pfcand_abseta.push_back(abs(reco_cand->eta()));
      pfcand_charge.push_back(reco_cand->charge());
      pfcand_isEl.push_back(abs(reco_cand->pdgId()) == 11);
      pfcand_isMu.push_back(abs(reco_cand->pdgId()) == 13);
      pfcand_isGamma.push_back(abs(reco_cand->pdgId()) == 22);
      pfcand_isChargedHad.push_back(abs(reco_cand->pdgId()) == 211);
      pfcand_isNeutralHad.push_back(abs(reco_cand->pdgId()) == 130);
      pfcand_lostInnerHits.push_back(reco_cand->lostInnerHits());
      pfcand_normchi2.push_back(reco_cand->normchi2());
      pfcand_quality.push_back(reco_cand->quality());
      pfcand_dz.push_back(reco_cand->dz());
      pfcand_dzsig.push_back(reco_cand->dzsig());
      pfcand_dxy.push_back(reco_cand->dxy());
      pfcand_dxysig.push_back(reco_cand->dxysig());
      pfcand_btagEtaRel.push_back(reco::btau::etaRel(jet_dir, track_mom));
      pfcand_btagPtRatio.push_back(track_mom3.Perp(jet_dir3) / track_mag);
      pfcand_btagPParRatio.push_back(jet_dir.Dot(track_mom) / track_mag);
    }

    fj_pt = j.pt();
    fj_eta = j.eta();
    fj_phi = j.phi();
    fj_mass = j.m();

    PseudoJet sd_ak8 = sd_groomer(j);
    fj_msd = sd_ak8.m();
    fj_n2b1 = N2(sd_ak8);

    label_Top_bcq = (ak8_label.first == FatJetMatching::Top_bcq);
    label_Top_bqq = (ak8_label.first == FatJetMatching::Top_bqq);
    label_Top_bc = (ak8_label.first == FatJetMatching::Top_bc);
    label_Top_bq = (ak8_label.first == FatJetMatching::Top_bq);
    label_W_cq = (ak8_label.first == FatJetMatching::W_cq);
    label_W_qq = (ak8_label.first == FatJetMatching::W_qq);
    label_Z_bb = (ak8_label.first == FatJetMatching::Z_bb);
    label_Z_cc = (ak8_label.first == FatJetMatching::Z_cc);
    label_Z_qq = (ak8_label.first == FatJetMatching::Z_qq);
    label_H_bb = (ak8_label.first == FatJetMatching::H_bb);
    label_H_cc = (ak8_label.first == FatJetMatching::H_cc);
    label_H_qqqq = (ak8_label.first == FatJetMatching::H_qqqq);
    label_H_tautau = (ak8_label.first == FatJetMatching::H_tautau);
    label_H_qq = (ak8_label.first == FatJetMatching::H_qq);
    label_QCD_all = (ak8_label.first == FatJetMatching::QCD_all);

    tree->Fill();	
    clearVars();
  }
}

void ScoutingNanoAOD::clearVars(){
  pfcand_pt_log_nopuppi.clear();
  pfcand_e_log_nopuppi.clear();
  pfcand_etarel.clear();
  pfcand_phirel.clear();
  pfcand_abseta.clear();
  pfcand_charge.clear();
  pfcand_isEl.clear();
  pfcand_isMu.clear();
  pfcand_isGamma.clear();
  pfcand_isChargedHad.clear();
  pfcand_isNeutralHad.clear();
  pfcand_lostInnerHits.clear();
  pfcand_normchi2.clear();
  pfcand_quality.clear();
  pfcand_dz.clear();
  pfcand_dzsig.clear();
  pfcand_dxy.clear();
  pfcand_dxysig.clear();
  pfcand_btagEtaRel.clear();
  pfcand_btagPtRatio.clear();
  pfcand_btagPParRatio.clear();
}

void ScoutingNanoAOD::beginJob() {
  
}

void ScoutingNanoAOD::endJob() {
}

void ScoutingNanoAOD::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
  // HLT paths

  triggerPathsVector.push_back("DST_DoubleMu1_noVtx_CaloScouting_v*");
  triggerPathsVector.push_back("DST_DoubleMu3_noVtx_CaloScouting_v*");
  triggerPathsVector.push_back("DST_DoubleMu3_noVtx_Mass10_PFScouting_v*");
  triggerPathsVector.push_back("DST_L1HTT_CaloScouting_PFScouting_v*");
  triggerPathsVector.push_back("DST_CaloJet40_CaloScouting_PFScouting_v*");
  triggerPathsVector.push_back("DST_HT250_CaloScouting_v*");
  triggerPathsVector.push_back("DST_HT410_PFScouting_v*");
  triggerPathsVector.push_back("DST_HT450_PFScouting_v*");

  HLTConfigProvider hltConfig;
  bool changedConfig = false;
  hltConfig.init(iRun, iSetup, triggerResultsTag.process(), changedConfig);

  for (size_t i = 0; i < triggerPathsVector.size(); i++) {
    triggerPathsMap[triggerPathsVector[i]] = -1;
  }

  for(size_t i = 0; i < triggerPathsVector.size(); i++){
    TPRegexp pattern(triggerPathsVector[i]);
    for(size_t j = 0; j < hltConfig.triggerNames().size(); j++){
      std::string pathName = hltConfig.triggerNames()[j];
      if(TString(pathName).Contains(pattern)){
	triggerPathsMap[triggerPathsVector[i]] = j;
      }
    }
  }
}

void ScoutingNanoAOD::endRun(edm::Run const&, edm::EventSetup const&) {
}

void ScoutingNanoAOD::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {
}

void ScoutingNanoAOD::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void ScoutingNanoAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScoutingNanoAOD);
