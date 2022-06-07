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
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "Run3ScoutingAnalysisTools/Analysis/interface/FatJetMatching.h"

#include "DataFormats/JetMatching/interface/JetFlavourInfo.h"
#include "DataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Common/interface/View.h"

using namespace std;
using namespace deepntuples;

FatJetMatching ak4_match;

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
  bool isNeutralPdg(int);
  //const edm::EDGetTokenT<std::vector<Run3ScoutingParticle> >  	pfcandsParticleNetToken;
  const edm::EDGetTokenT<edm::View<pat::Jet>> jetToken_;
  const edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> jetFlavourInfosToken_;

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

  float jet_pt;
  float jet_eta;
  float jet_phi;
  float jet_mass;
  
  float jet_nCHadrons;
  float jet_nBHadrons;
  float jet_partonFlavour;
  float jet_hadronFlavour;

  int event_no;

  bool debug_match = false;
  int debug_match_numJets = 5;
};

ScoutingNanoAOD::ScoutingNanoAOD(const edm::ParameterSet& iConfig):
  //pfcandsParticleNetToken  (consumes<std::vector<Run3ScoutingParticle> > (iConfig.getParameter<edm::InputTag>("pfcandsParticleNet"))),
  jetToken_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"))),
  jetFlavourInfosToken_    (consumes<reco::JetFlavourInfoMatchingCollection> (iConfig.getParameter<edm::InputTag>("jetFlavourInfos")))
{
  usesResource("TFileService");

  edm::Service<TFileService> fs;

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

  tree->Branch("jet_pt", &jet_pt);
  tree->Branch("jet_eta", &jet_eta);
  tree->Branch("jet_phi", &jet_phi);
  tree->Branch("jet_mass", &jet_mass);
  // NEW:  nc, nb, partonFlavour
  tree->Branch("jet_nCHadrons", &jet_nCHadrons);
  tree->Branch("jet_nBHadrons", &jet_nBHadrons);
  tree->Branch("jet_partonFlavour", &jet_partonFlavour);
  tree->Branch("jet_hadronFlavour", &jet_hadronFlavour);

  tree->Branch("event_no", &event_no);
}

ScoutingNanoAOD::~ScoutingNanoAOD() {
}

// https://github.com/cms-sw/cmssw/blob/6d2f66057131baacc2fcbdd203588c41c885b42c/PhysicsTools/JetMCAlgos/plugins/GenHFHadronMatcher.cc#L1022
bool ScoutingNanoAOD::isNeutralPdg(int pdgId) {
   const int neutralPdgs_array[] = {9, 21, 22, 23, 25, 12, 14, 16, 111, 130, 310, 311, 421, 511, 2112}; // gluon, gluon, gamma, Z0, higgs, electron neutrino, muon neutrino, tau neutrino, pi0, K0_L, K0_S; K0, neutron
   const std::vector<int> neutralPdgs(neutralPdgs_array, neutralPdgs_array + sizeof(neutralPdgs_array) / sizeof(int));
   if (std::find(neutralPdgs.begin(), neutralPdgs.end(), std::abs(pdgId)) == neutralPdgs.end())
     return false;
 
   return true;
}

void ScoutingNanoAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace pat;
  using namespace fastjet;
  using namespace fastjet::contrib;

  //Handle<vector<Run3ScoutingParticle> > pfcandsParticleNetH;
  //iEvent.getByToken(pfcandsParticleNetToken, pfcandsParticleNetH);

  Handle<edm::View<pat::Jet>> ak4_jets;
  iEvent.getByToken(jetToken_, ak4_jets);

  Handle<reco::JetFlavourInfoMatchingCollection> theJetFlavourInfos;
  iEvent.getByToken(jetFlavourInfosToken_, theJetFlavourInfos );

  // Create AK4 Jets
  //vector<PseudoJet> jet_part;
  //jet_part.reserve(pfcandsParticleNetH->size());
  //int pfcand_i = 0;
  //for (auto pfcands_iter = pfcandsParticleNetH->begin(); pfcands_iter != pfcandsParticleNetH->end(); ++pfcands_iter) {
  //  math::PtEtaPhiMLorentzVector p4(pfcands_iter->pt(), pfcands_iter->eta(), pfcands_iter->phi(), pfcands_iter->m());
  //  jet_part.emplace_back(p4.px(), p4.py(), p4.pz(), p4.energy());
  //  jet_part.back().set_user_index(pfcand_i);
  //  pfcand_i++;
  //}

  //JetDefinition ak4_def = JetDefinition(antikt_algorithm, 0.4);  // was 0.8
  //fastjet::GhostedAreaSpec area_spec(5.0,1,0.01);
  //fastjet::AreaDefinition area_def(fastjet::active_area, area_spec);

  //ClusterSequenceArea ak4_cs(jet_part, ak4_def, area_def);
  //vector<PseudoJet> ak4_jets = sorted_by_pt(ak4_cs.inclusive_jets(15.0));

  // NEW:  Match slimjets w/ ak4_jets, record pdgID info
  // Note:  must do before main loop to avoid previous bug
  // outline:  (nc, nb, parton flavor, partflav, hadronFlavour)
  // - create pairList (vector of int vectors):  All dRs < 0.4 of all possible (slimjet, ak4) pairs in increasing order
  // - create map of results:  std::map<pat::Jet, PseudoJet>
  // - Loop:  Grab smallest dR, store match in map, remove jet from pairList
  //set_user_index, user_index
  //NEW:  use the user_index defined above to reference each jet!
  std::map<int, reco::JetFlavourInfo> resultMap;
  std::vector<int> unmatchedJets;  //vector of IDs of unmatched jets
  // note:  pairList will have duplicates; this is okay
  std::vector<std::tuple<int, int, float> > pairList;  // jet, ak4_, dR

  int ak4_jet_idx = 0;
  if (debug_match) std::cout << "dR loop:" << std::endl;

  for(unsigned int i=0; i<ak4_jets->size(); i++) {
    if (debug_match && ak4_jet_idx > debug_match_numJets) break;
    bool found_match = false;
    for(unsigned int j=0; j<theJetFlavourInfos->size(); j++) {
      // calc dR
      float dR = reco::deltaR(ak4_jets->at(i).eta(), ak4_jets->at(i).phi(), (*theJetFlavourInfos)[j].first.get()->eta(), (*theJetFlavourInfos)[j].first.get()->phi());
      if (debug_match) std::cout << i << " " << j << " " << dR << " " << ak4_jets->at(i).pt() << std::endl;
      if(dR < 0.4) {
        pairList.push_back(std::make_tuple(i, j, dR));
        found_match = true;
      }
    }
    ak4_jet_idx++;
    if(!found_match) {
      //unmatchedJets.push_back(ak4_jets[i].user_index());
      unmatchedJets.push_back(i);
    }
  }

  if (debug_match) {
    std::cout << "\nunmatchedJets loop:" << std::endl;
    for(auto &j: unmatchedJets) {
      std::cout << j << std::endl;
    }
    std::cout << "\npairList loop:" << std::endl;
  }

  // next, sort in order of increasing dR:
  std::sort(pairList.begin(), pairList.end(), [](std::tuple<int, int, float> t1, std::tuple<int, int, float> t2){ return std::get<2>(t1) < std::get<2>(t2); });
  // go through each jet and match:
  while(pairList.size() > 0) {
    if (debug_match) std::cout << std::get<0>(pairList[0]) << " " << std::get<1>(pairList[0]) << " " << std::get<2>(pairList[0]) << std::endl;

    // grab first element, assign pair
    //auto ak4_assn = ak4_jets[std::get<0>(pairList[0])];
    //int uindex = ak4_assn.user_index();
    reco::JetFlavourInfo jetFlavour_assn = (*theJetFlavourInfos)[std::get<1>(pairList[0])].second;
    //resultMap[uindex] = jetFlavour_assn;
    resultMap[std::get<0>(pairList[0])] = jetFlavour_assn;
    // remove all particles matched to that jet
    for(unsigned int k=1; k<pairList.size(); k++) {
      if(std::get<0>(pairList[k]) == std::get<0>(pairList[0]) ||
         std::get<1>(pairList[k]) == std::get<1>(pairList[0])) {
        pairList.erase(pairList.begin() + k);
      }
    }
    pairList.erase(pairList.begin());
  } // matching finished; results stored in resultMap

  if (debug_match) {
    std::cout << "\nresultMap loop:" << std::endl;
    for(auto &r: resultMap) {
      std::cout << r.first << std::endl;
    }
  }

  int num_unmatched_geq30 = 0;
  int num_matched_geq30 = 0;

  int num_unmatched_leq30 = 0;
  int num_matched_leq30 = 0;
  ak4_jet_idx = 0;
  for(const pat::Jet & j : *ak4_jets) {  // was ak8, etc
    if (debug_match && ak4_jet_idx > debug_match_numJets) break;

    //float etasign = j.eta() > 0 ? 1 : -1;

    //// The following is needed to compute btagEtaRel, btagPtRatio and btagPParRatio
    //float jet_px = j.pt() * cos(j.phi());
    //float jet_py = j.pt() * sin(j.phi());
    //float jet_pz = j.pt() * sinh(j.eta());
    //math::XYZVector jet_dir_temp(jet_px, jet_py, jet_pz);
    //math::XYZVector jet_dir = jet_dir_temp.Unit();
    //TVector3 jet_dir3(jet_px, jet_py, jet_pz);

    //// Loop over AK4 jet constituents
    //const vector<PseudoJet> constituents = j.constituents();
    //for (auto &cand : constituents) {
    //  // Match PseudoJet constituent to PF candidate
    //  auto *reco_cand = dynamic_cast<const Run3ScoutingParticle*> (&pfcandsParticleNetH->at(cand.user_index()));
    //  // The following is needed to compute btagEtaRel, btagPtRatio and btagPParRatio
    //  float trk_px = reco_cand->trk_pt() * cos(reco_cand->trk_phi());
    //  float trk_py = reco_cand->trk_pt() * sin(reco_cand->trk_phi());
    //  float trk_pz = reco_cand->trk_pt() * sinh(reco_cand->trk_eta());
    //  math::XYZVector track_mom(trk_px, trk_py, trk_pz);
    //  TVector3 track_mom3(trk_px, trk_py, trk_pz);
    //  double track_mag = sqrt(trk_px * trk_px + trk_py * trk_py + trk_pz * trk_pz);

    //  float reco_cand_p = reco_cand->pt() * cosh(reco_cand->eta());
    //  pfcand_pt_log_nopuppi.push_back(log(reco_cand->pt()));
    //  pfcand_e_log_nopuppi.push_back(log(sqrt(reco_cand_p*reco_cand_p + reco_cand->m()*reco_cand->m())));
    //  pfcand_etarel.push_back(etasign * (reco_cand->eta() - j.eta()));
    //  pfcand_phirel.push_back(deltaPhi(reco_cand->phi(), j.phi()));
    //  pfcand_abseta.push_back(abs(reco_cand->eta()));
    //  if (isNeutralPdg(reco_cand->pdgId())) {
    //     pfcand_charge.push_back(0);
    //  } else {
    //     pfcand_charge.push_back(abs(reco_cand->pdgId())/reco_cand->pdgId());
    //  }
    //  pfcand_isEl.push_back(abs(reco_cand->pdgId()) == 11);
    //  pfcand_isMu.push_back(abs(reco_cand->pdgId()) == 13);
    //  pfcand_isGamma.push_back(abs(reco_cand->pdgId()) == 22);
    //  pfcand_isChargedHad.push_back(abs(reco_cand->pdgId()) == 211);
    //  pfcand_isNeutralHad.push_back(abs(reco_cand->pdgId()) == 130);
    //  pfcand_lostInnerHits.push_back(reco_cand->lostInnerHits());
    //  pfcand_normchi2.push_back(reco_cand->normchi2());
    //  pfcand_quality.push_back(reco_cand->quality());
    //  pfcand_dz.push_back(reco_cand->dz());
    //  pfcand_dzsig.push_back(reco_cand->dzsig());
    //  pfcand_dxy.push_back(reco_cand->dxy());
    //  pfcand_dxysig.push_back(reco_cand->dxysig());
    //  pfcand_btagEtaRel.push_back(reco::btau::etaRel(jet_dir, track_mom));
    //  pfcand_btagPtRatio.push_back(track_mom3.Perp(jet_dir3) / track_mag);
    //  pfcand_btagPParRatio.push_back(jet_dir.Dot(track_mom) / track_mag);
    //}

    jet_pt = j.pt();
    jet_eta = j.eta();
    jet_phi = j.phi();
    jet_mass = j.mass();
    // check whether jet has been matched first
    //if(std::find(unmatchedJets.begin(), unmatchedJets.end(), j.user_index()) != unmatchedJets.end()) {
    if(std::find(unmatchedJets.begin(), unmatchedJets.end(), ak4_jet_idx) != unmatchedJets.end()) {

      if (debug_match) std::cout << "\nUnmatched!" << std::endl;
      if (j.pt() > 30) {
        num_unmatched_geq30++;
      } else {
        num_unmatched_leq30++;
      }

      //unmatched; assign dummy value
      jet_nCHadrons = -99;
      jet_nBHadrons = -99;
      jet_partonFlavour = -99;
      jet_hadronFlavour = -99;
    } else {
	
      if (debug_match) std::cout << "\nMatched!" << std::endl;
      if (j.pt() > 30) {
        num_matched_geq30++;
      } else {
        num_matched_leq30++;
      }

      reco::JetFlavourInfo flavJet = resultMap[ak4_jet_idx];
      jet_nCHadrons = flavJet.getcHadrons().size();
      jet_nBHadrons = flavJet.getbHadrons().size();
      jet_partonFlavour = flavJet.getPartonFlavour();
      jet_hadronFlavour = flavJet.getHadronFlavour();
    }

    ak4_jet_idx++;

    event_no = iEvent.id().event();

    tree->Fill();	
    clearVars();
  }

  if (!debug_match) {
    std::cout << "\nMatched (pT > 30 GeV): " << num_matched_geq30 << "\n(" << (float)num_matched_geq30/(float)(num_unmatched_geq30 + num_matched_geq30) * 100 << "%)" << std::endl; 
    std::cout << "Matched (pT < 30 GeV): " << num_matched_leq30 << "\n(" << (float)num_matched_leq30/(float)(num_unmatched_leq30 + num_matched_leq30) * 100 << "%)" << std::endl; 
    std::cout << "Unmatched (pT > 30 GeV): " << num_unmatched_geq30 << "\n(" << (float)num_unmatched_geq30/(float)(num_unmatched_geq30 + num_matched_geq30) * 100 << "%)" << std::endl;
    std::cout << "Unmatched (pT < 30 GeV): " << num_unmatched_leq30 << "\n(" << (float)num_unmatched_leq30/(float)(num_unmatched_leq30 + num_matched_leq30) * 100 << "%)\n" << std::endl;
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
