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
  const edm::EDGetTokenT<std::vector<Run3ScoutingParticle> >  	pfcandsParticleNetToken;
  const edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> jetFlavourInfosToken_;

  // TTree carrying the event weight information
  TTree* tree;

  //PFCand ParticleNet
  vector<Float16_t> jet_pfcand_pt_log;
  vector<Float16_t> jet_pfcand_energy_log;
  vector<Float16_t> jet_pfcand_deta;
  vector<Float16_t> jet_pfcand_dphi;
  vector<Float16_t> jet_pfcand_eta;
  vector<Float16_t> jet_pfcand_charge;
  vector<Float16_t> jet_pfcand_nlostinnerhits;
  vector<Float16_t> jet_pfcand_track_chi2;
  vector<Float16_t> jet_pfcand_track_qual;
  vector<Float16_t> jet_pfcand_dz;
  vector<Float16_t> jet_pfcand_dzsig;
  vector<Float16_t> jet_pfcand_dxy;
  vector<Float16_t> jet_pfcand_dxysig;
  vector<Float16_t> jet_pfcand_etarel;
  vector<Float16_t> jet_pfcand_pperp_ratio;
  vector<Float16_t> jet_pfcand_ppara_ratio;

  //Jet kinematics
  float jet_pt;
  float jet_eta;
  float jet_phi;
  float jet_mass;
  // NEW:
  float jet_nCHadrons;
  float jet_nBHadrons;
  float jet_partonFlavour;
  float jet_hadronFlavour;

  int event_no;
};

ScoutingNanoAOD::ScoutingNanoAOD(const edm::ParameterSet& iConfig):
  pfcandsParticleNetToken  (consumes<std::vector<Run3ScoutingParticle> > (iConfig.getParameter<edm::InputTag>("pfcandsParticleNet"))),
  jetFlavourInfosToken_    (consumes<reco::JetFlavourInfoMatchingCollection> (iConfig.getParameter<edm::InputTag>("jetFlavourInfos")))
{
  usesResource("TFileService");

 // Access the TFileService
  edm::Service<TFileService> fs;

  // Create the TTree
  tree = fs->make<TTree>("tree", "tree");

  tree->Branch("jet_pfcand_pt_log", &jet_pfcand_pt_log);
  tree->Branch("jet_pfcand_energy_log", &jet_pfcand_energy_log);
  tree->Branch("jet_pfcand_deta", &jet_pfcand_deta);
  tree->Branch("jet_pfcand_dphi", &jet_pfcand_dphi);
  tree->Branch("jet_pfcand_eta", &jet_pfcand_eta);
  tree->Branch("jet_pfcand_charge", &jet_pfcand_charge);
  tree->Branch("jet_pfcand_nlostinnerhits", &jet_pfcand_nlostinnerhits);
  tree->Branch("jet_pfcand_track_chi2", &jet_pfcand_track_chi2);
  tree->Branch("jet_pfcand_track_qual", &jet_pfcand_track_qual);
  tree->Branch("jet_pfcand_dz", &jet_pfcand_dz);
  tree->Branch("jet_pfcand_dzsig", &jet_pfcand_dzsig);
  tree->Branch("jet_pfcand_dxy", &jet_pfcand_dxy);
  tree->Branch("jet_pfcand_dxysig", &jet_pfcand_dxysig);
  tree->Branch("jet_pfcand_etarel", &jet_pfcand_etarel);
  tree->Branch("jet_pfcand_pperp_ratio", &jet_pfcand_pperp_ratio);
  tree->Branch("jet_pfcand_ppara_ratio", &jet_pfcand_ppara_ratio);

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
  using namespace pat;  // new
  using namespace fastjet;
  using namespace fastjet::contrib;

  Handle<vector<Run3ScoutingParticle> > pfcandsParticleNetH;
  iEvent.getByToken(pfcandsParticleNetToken, pfcandsParticleNetH);

  Handle<reco::JetFlavourInfoMatchingCollection> theJetFlavourInfos;
  iEvent.getByToken(jetFlavourInfosToken_, theJetFlavourInfos );

  // Create AK4 Jets
  vector<PseudoJet> jet_part;
  jet_part.reserve(pfcandsParticleNetH->size());
  int pfcand_i = 0;
  for (auto pfcands_iter = pfcandsParticleNetH->begin(); pfcands_iter != pfcandsParticleNetH->end(); ++pfcands_iter) {
    math::PtEtaPhiMLorentzVector p4(pfcands_iter->pt(), pfcands_iter->eta(), pfcands_iter->phi(), pfcands_iter->m());
    jet_part.emplace_back(p4.px(), p4.py(), p4.pz(), p4.energy());
    jet_part.back().set_user_index(pfcand_i);
    pfcand_i++;
  }

  JetDefinition ak4_def = JetDefinition(antikt_algorithm, 0.4);  // was 0.8
  fastjet::GhostedAreaSpec area_spec(5.0,1,0.01);
  fastjet::AreaDefinition area_def(fastjet::active_area, area_spec);

  ClusterSequenceArea ak4_cs(jet_part, ak4_def, area_def);
  vector<PseudoJet> ak4_jets = sorted_by_pt(ak4_cs.inclusive_jets(15.0));

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
  for(unsigned int i=0; i<ak4_jets.size(); i++) {
    bool found_match = false;
    for(unsigned int j=0; j<theJetFlavourInfos->size(); j++) {
      // calc dR
      float dR = reco::deltaR(ak4_jets[i].eta(), ak4_jets[i].phi(), (*theJetFlavourInfos)[j].first.get()->eta(), (*theJetFlavourInfos)[j].first.get()->phi());
      if(dR < 0.4) {
        //std::cout << "chadrons is " << (*slimjetH)[j].jetFlavourInfo().getcHadrons().size() << std::endl;
        //std::cout << "eta is " << (*slimjetH)[j].eta() << std::endl;
        pairList.push_back(std::make_tuple(i, j, dR));
        found_match = true;
      }
    }
    if(!found_match) {
      unmatchedJets.push_back(ak4_jets[i].user_index());
    }
  }
  // next, sort in order of increasing dR:
  std::sort(pairList.begin(), pairList.end(), [](std::tuple<int, int, float> t1, std::tuple<int, int, float> t2){ return std::get<2>(t1) < std::get<2>(t2); });
  // go through each jet and match:
  while(pairList.size() > 0) {
    // grab first element, assign pair
    PseudoJet ak4_assn = ak4_jets[std::get<0>(pairList[0])];
    int uindex = ak4_assn.user_index();
    reco::JetFlavourInfo jetFlavour_assn = (*theJetFlavourInfos)[std::get<1>(pairList[0])].second;
    resultMap[uindex] = jetFlavour_assn;
    //std::cout << "Adding, nchadrons= " << jetFlavour_assn.jetFlavourInfo().getcHadrons().size() << std::endl;
    // remove all particles matched to that jet
    for(unsigned int k=1; k<pairList.size(); k++) {
      if(std::get<0>(pairList[k]) == std::get<0>(pairList[0]) ||
         std::get<1>(pairList[k]) == std::get<1>(pairList[0])) {
        pairList.erase(pairList.begin() + k);
      }
    }
    pairList.erase(pairList.begin());
  } // matching finished; results stored in resultMap


  for(auto &j: ak4_jets) {  // was ak8, etc


    // The following is needed to compute btagEtaRel, btagPtRatio and btagPParRatio
    float jet_px = j.pt() * cos(j.phi());
    float jet_py = j.pt() * sin(j.phi());
    float jet_pz = j.pt() * sinh(j.eta());
    math::XYZVector jet_dir_temp(jet_px, jet_py, jet_pz);
    math::XYZVector jet_dir = jet_dir_temp.Unit();
    TVector3 jet_direction(jet_px, jet_py, jet_pz);

    // Loop over AK4 jet constituents
    const vector<PseudoJet> constituents = j.constituents();
    for (auto &cand : constituents) {
      // Match PseudoJet constituent to PF candidate
      auto *reco_cand = dynamic_cast<const Run3ScoutingParticle*> (&pfcandsParticleNetH->at(cand.user_index()));
      float pfcand_px = reco_cand->pt() * cos(reco_cand->phi());
      float pfcand_py = reco_cand->pt() * sin(reco_cand->phi());
      float pfcand_pz = reco_cand->pt() * sinh(reco_cand->eta());
      math::XYZVector cand_dir(pfcand_px, pfcand_py, pfcand_pz);
      TVector3 cand_direction(cand_dir.x(), cand_dir.y(), cand_dir.z());

      float reco_cand_p = reco_cand->pt() * cosh(reco_cand->eta());
      jet_pfcand_pt_log.push_back(log(reco_cand->pt()));
      jet_pfcand_energy_log.push_back(log(sqrt(reco_cand_p*reco_cand_p + reco_cand->m()*reco_cand->m())));
      jet_pfcand_deta.push_back(j.eta() - reco_cand->eta());
      jet_pfcand_dphi.push_back(deltaPhi(j.phi(), reco_cand->phi()));
      jet_pfcand_eta.push_back(reco_cand->eta());
      if (isNeutralPdg(reco_cand->pdgId())) {
         jet_pfcand_charge.push_back(0);
      } else {
         jet_pfcand_charge.push_back(abs(reco_cand->pdgId())/reco_cand->pdgId());
      }
      jet_pfcand_nlostinnerhits.push_back(reco_cand->lostInnerHits());
      jet_pfcand_track_chi2.push_back(reco_cand->normchi2());
      jet_pfcand_track_qual.push_back(reco_cand->quality());
      jet_pfcand_dz.push_back(reco_cand->dz());
      jet_pfcand_dzsig.push_back(reco_cand->dzsig());
      jet_pfcand_dxy.push_back(reco_cand->dxy());
      jet_pfcand_dxysig.push_back(reco_cand->dxysig());
      jet_pfcand_etarel.push_back(reco::btau::etaRel(jet_dir, cand_dir));
      jet_pfcand_pperp_ratio.push_back(jet_direction.Perp(cand_direction) / cand_direction.Mag());
      jet_pfcand_ppara_ratio.push_back(jet_direction.Dot(cand_direction) / cand_direction.Mag());
    }

    jet_pt = j.pt();
    jet_eta = j.eta();
    jet_phi = j.phi();
    jet_mass = j.m();
    // NEW
    // check whether jet has been matched first
    if(std::find(unmatchedJets.begin(), unmatchedJets.end(), j.user_index()) == unmatchedJets.end()) {
      //unmatched; assign dummy value
      jet_nCHadrons = -99;
      jet_nBHadrons = -99;
      jet_partonFlavour = -99;
      jet_hadronFlavour = -99;
    } else {
      reco::JetFlavourInfo flavJet = resultMap[j.user_index()];
      jet_nCHadrons = flavJet.getcHadrons().size();
      jet_nBHadrons = flavJet.getbHadrons().size();
      jet_partonFlavour = flavJet.getPartonFlavour();
      jet_hadronFlavour = flavJet.getHadronFlavour();
    }

    event_no = iEvent.id().event();

    tree->Fill();	
    clearVars();
  }
}

void ScoutingNanoAOD::clearVars(){
  jet_pfcand_pt_log.clear();
  jet_pfcand_energy_log.clear();
  jet_pfcand_deta.clear();
  jet_pfcand_dphi.clear();
  jet_pfcand_eta.clear();
  jet_pfcand_charge.clear();
  jet_pfcand_nlostinnerhits.clear();
  jet_pfcand_track_chi2.clear();
  jet_pfcand_track_qual.clear();
  jet_pfcand_dz.clear();
  jet_pfcand_dzsig.clear();
  jet_pfcand_dxy.clear();
  jet_pfcand_dxysig.clear();
  jet_pfcand_etarel.clear();
  jet_pfcand_pperp_ratio.clear();
  jet_pfcand_ppara_ratio.clear();
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
