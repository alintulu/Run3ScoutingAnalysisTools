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
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

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
std::unordered_set<const reco::GenParticle*> processed_;

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
  std::pair<const reco::GenParticle*, double> closest(fastjet::PseudoJet, std::vector<const reco::GenParticle*>&, bool);
  const edm::InputTag triggerResultsTag;
  const edm::EDGetTokenT<std::vector<Run3ScoutingParticle> >  	pfcandsParticleNetToken;
  const edm::EDGetTokenT<reco::GenParticleCollection>      genpartsToken;
  const edm::EDGetTokenT<double>  	pfMetToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingMuon> >      muonsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingPFJet> >  	pfjetsToken;


  // TTree carrying the event weight information
  TTree* tree;

  vector<Float16_t> fatjet_1_dr_T;
  vector<Float16_t> fatjet_1_dr_T_Wq_max;
  vector<int> fatjet_1_T_Wq_max_pdgId;
  vector<Float16_t> fatjet_1_dr_W_daus;
  vector<Float16_t> fatjet_1_dr_T_b;
  vector<Float_t> muon_pt;
  vector<Float_t> muon_eta;
  vector<Float_t> muon_trk_dxy;
  vector<Float_t> muon_trk_dz;
  vector<Float_t> muon_m;
  vector<Float_t> muon_charge;
  vector<Float_t> muon_type;
  vector<Float_t> muon_phi;
  vector<Float16_t> fatjet_area;
  vector<Float16_t> fatjet_eta;
  vector<Float16_t> fatjet_n2b1;
  vector<Float16_t> fatjet_n3b1;
  vector<Float16_t> fatjet_phi;
  vector<Float16_t> fatjet_pt;
  vector<Float16_t> fatjet_mass;
  vector<Float16_t> fatjet_msoftdrop;
  vector<Float16_t> fatjet_mregressed;
  vector<Float16_t> fatjet_doubleBTag;
  vector<Float_t> jet_pt;
  vector<Float_t> jet_eta;
  vector<Float_t> jet_phi;
  vector<Float_t> jet_m;
  vector<Float_t> jet_jetArea;
  vector<Float_t> jet_bTagScore;
  double met;
  
  bool debug = false;
};

ScoutingNanoAOD::ScoutingNanoAOD(const edm::ParameterSet& iConfig):
  pfcandsParticleNetToken  (consumes<std::vector<Run3ScoutingParticle> > (iConfig.getParameter<edm::InputTag>("pfcandsParticleNet"))),
  genpartsToken            (consumes<reco::GenParticleCollection> (iConfig.getParameter<edm::InputTag>("genpart"))),
  pfMetToken             (consumes<double>(iConfig.getParameter<edm::InputTag>("met"))),
  muonsToken               (consumes<std::vector<Run3ScoutingMuon> >             (iConfig.getParameter<edm::InputTag>("muons"))),
  pfjetsToken              (consumes<std::vector<Run3ScoutingPFJet> >            (iConfig.getParameter<edm::InputTag>("pfjets")))
{
  usesResource("TFileService");

 // Access the TFileService
  edm::Service<TFileService> fs;

  // Create the TTree
  tree = fs->make<TTree>("tree", "tree");

  tree->Branch("fatjet_1_dr_T",&fatjet_1_dr_T);
  tree->Branch("fatjet_1_dr_T_Wq_max",&fatjet_1_dr_T_Wq_max);
  tree->Branch("fatjet_1_T_Wq_max_pdgId",&fatjet_1_T_Wq_max_pdgId);
  tree->Branch("fatjet_1_dr_W_daus",&fatjet_1_dr_W_daus);
  tree->Branch("fatjet_1_dr_T_b",&fatjet_1_dr_T_b);
  tree->Branch("muon_pt", &muon_pt);
  tree->Branch("muon_eta", &muon_eta);
  tree->Branch("muon_trk_dxy", &muon_trk_dxy);
  tree->Branch("muon_trk_dz", &muon_trk_dz);
  tree->Branch("muon_m", &muon_m);
  tree->Branch("muon_charge", &muon_charge);
  tree->Branch("muon_type", &muon_type);
  tree->Branch("muon_phi", &muon_phi);
  tree->Branch("fatjet_area", &fatjet_area);
  tree->Branch("fatjet_eta", &fatjet_eta);
  tree->Branch("fatjet_n2b1", &fatjet_n2b1);
  tree->Branch("fatjet_n3b1", &fatjet_n3b1);
  tree->Branch("fatjet_phi", &fatjet_phi);
  tree->Branch("fatjet_pt", &fatjet_pt);
  tree->Branch("fatjet_mass", &fatjet_mass);
  tree->Branch("fatjet_msoftdrop", &fatjet_msoftdrop);
  tree->Branch("fatjet_mregressed", &fatjet_mregressed);
  tree->Branch("fatjet_doubleBTag", &fatjet_doubleBTag);
  tree->Branch("jet_pt", &jet_pt);
  tree->Branch("jet_eta", &jet_eta);
  tree->Branch("jet_phi", &jet_phi);
  tree->Branch("jet_m", &jet_m);
  tree->Branch("jet_jetArea", &jet_jetArea);
  tree->Branch("jet_bTagScore", &jet_bTagScore);
  tree->Branch("met", &met );
}

ScoutingNanoAOD::~ScoutingNanoAOD() {
}

std::pair<const reco::GenParticle*, double> ScoutingNanoAOD::closest(fastjet::PseudoJet jet, std::vector<const reco::GenParticle*>& gps, bool debug){

 const reco::GenParticle *ret = nullptr;
 double drMin = 1e6;
 for (const auto &gp : gps) {
   double dr = reco::deltaR(jet.eta(), jet.phi(), gp->eta(), gp->phi());
   if (dr < drMin) {
     ret = dynamic_cast<const reco::GenParticle*>(&(*gp));
     drMin = dr;
   }
   if (debug and drMin == 1e6) std::cout << "Did not find a smallest distance" << std::endl;
 }
 return std::make_pair(ret, drMin);
}

void ScoutingNanoAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace fastjet;
  using namespace fastjet::contrib;

  Handle<vector<Run3ScoutingParticle> > pfcandsParticleNetH;
  iEvent.getByToken(pfcandsParticleNetToken, pfcandsParticleNetH);

  Handle<GenParticleCollection> genpartH;
  iEvent.getByToken(genpartsToken, genpartH);

  Handle<double> pfMetH;
  iEvent.getByToken(pfMetToken, pfMetH);
  met = *pfMetH;

  Handle<vector<Run3ScoutingPFJet> > pfjetsH;
  iEvent.getByToken(pfjetsToken, pfjetsH);

  Handle<vector<Run3ScoutingMuon> > muonsH;
  iEvent.getByToken(muonsToken, muonsH);

  for (auto iter = muonsH->begin(); iter != muonsH->end(); ++iter) {
    muon_pt.push_back(iter->pt());
    muon_eta.push_back(iter->eta());
    muon_phi.push_back(iter->phi());
    muon_m.push_back(iter->m());
    muon_type.push_back(iter->type());
    muon_charge.push_back(iter->charge());
    muon_trk_dxy.push_back(iter->trk_dxy());
    muon_trk_dz.push_back(iter->trk_dz());
  }

  for (auto iter = pfjetsH->begin(); iter != pfjetsH->end(); ++iter) {
    jet_pt.push_back(iter->pt());
    jet_eta.push_back(iter->eta());
    jet_phi.push_back(iter->phi());
    jet_m.push_back(iter->m());
    jet_jetArea.push_back(iter->jetArea()); 
    jet_bTagScore.push_back(99.0);
  }
  
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
  SoftDrop sd_groomer = SoftDrop(sd_beta, sd_z_cut, 0.8);
  EnergyCorrelatorN2 N2 = EnergyCorrelatorN2(1.0);
  EnergyCorrelatorN3 N3=EnergyCorrelatorN3(1.0);

  ClusterSequenceArea ak8_cs(fj_part, ak8_def, area_def);
  //vector<PseudoJet> ak8_jets = sorted_by_pt(ak8_cs.inclusive_jets(170.0));
  vector<PseudoJet> ak8_jets = sorted_by_pt(ak8_cs.inclusive_jets());

  std::vector<const reco::GenParticle*> hadGenTops;
  std::vector<const reco::GenParticle*> hadGenWs;
 
  if (debug) {
    ak8_match.printGenInfoHeader();
    for(unsigned int j=0; j<genpartH->size(); j++) {
      const auto *gp = &(*genpartH)[j];
      ak8_match.printGenParticleInfo(gp, j);
    }
   }

  // Loop over gen particles
  for(unsigned int j=0; j<genpartH->size(); j++) {
    const auto *gp = &(*genpartH)[j];

    if (processed_.count(gp)) continue;
    processed_.insert(gp);

    // Check if it is top
    if (std::abs(gp->pdgId()) == 6 and gp->numberOfDaughters() != 1) {
      auto top = ak8_match.getFinal(gp);

      // Check the daughters of the top
      const reco::GenParticle *w_from_top = nullptr;
      for (const auto &dau : top->daughterRefVector()) {
         if (std::abs(dau->pdgId()) == 24) w_from_top = ak8_match.getFinal(&(*dau));
      }
      // If the daughter is the W
      if (w_from_top) {
         // Check if the W decays hadronically
         if (ak8_match.isHadronic(w_from_top)) {
           hadGenTops.push_back(gp);
           if (debug) std::cout << j << " Top -> Wx -> qqx" << std::endl;
         }
      }
    // Check if it is the W
    } else if (std::abs(gp->pdgId()) == 24) {

        // Check if the W decays hadronically
        if(ak8_match.isHadronic(gp)) {
          hadGenWs.push_back(gp);
          if (debug) std::cout << j << " W -> qq" << std::endl;
        }
     }
  }

  if (debug) std::cout << "Number of jets: " << ak8_jets.size() << std::endl;

  int i = 0;
  // Loop over jets
  for (auto &jet: ak8_jets) {

    fatjet_area.push_back(jet.area());
    fatjet_eta.push_back(jet.pseudorapidity());
    fatjet_mass.push_back(jet.m());
    fatjet_phi.push_back(jet.phi_std());
    fatjet_pt.push_back(jet.pt());

    PseudoJet sd_ak8 = sd_groomer(jet);
    fatjet_msoftdrop.push_back(sd_ak8.m());
    fatjet_n2b1.push_back(N2(sd_ak8));
    fatjet_n3b1.push_back(N3(sd_ak8));
    
    fatjet_mregressed.push_back(-99.0);
    fatjet_doubleBTag.push_back(-99.0);

    i++;
    if (debug) std::cout << i << " Jet pT: " << jet.pt() << std::endl;
    // Find closest top to jet
    auto closestT = closest(jet, hadGenTops, debug);
    auto genT = closestT.first;
    auto drT = closestT.second;
    if (debug) std::cout << i <<" deltaR(jet, top)   : " << drT << std::endl;

    // Find closest W to jet
    auto closestW = closest(jet, hadGenWs, debug);
    auto genW = closestW.first;
    auto drW = closestW.second;
    if (debug) std::cout << i << " deltaR(jet, W)   : " << drW << std::endl;

    if (genW) {
      auto wdaus = ak8_match.getDaughterQuarks(genW);
      double dr_q1 = reco::deltaR(jet, *wdaus.at(0));
      double dr_q2 = reco::deltaR(jet, *wdaus.at(1));
      // If the 1st daughter is closer to the jet than the 2nd, then swap their places
      if (dr_q1 < dr_q2) {
          std::swap(dr_q1, dr_q2);
      }
      if (debug) {
        using namespace std;
        cout << "Wx -> qqx" << endl;
        cout << "deltaR(jet, q1)    : " << dr_q1 << endl;
        cout << "deltaR(jet, q2)    : " << dr_q2 << endl;
      }
      fatjet_1_dr_W_daus.push_back(dr_q1);
    } else {
      fatjet_1_dr_W_daus.push_back(99.0); 
    }

    fatjet_1_dr_T.push_back(drT);

    // If there is no top
    if (!genT) {
      fatjet_1_dr_T_Wq_max.push_back(99.0);
      fatjet_1_T_Wq_max_pdgId.push_back(0);
      fatjet_1_dr_T_b.push_back(99.0);
      continue;
    }

    // Check the daughters of the top
    for (const auto &dau : genT->daughterRefVector()){
      const reco::GenParticle *w_from_top = nullptr;
      const reco::GenParticle *b_from_top = nullptr;
      if (std::abs(dau->pdgId()) == 24) w_from_top = ak8_match.getFinal(&(*dau));
      if (std::abs(dau->pdgId()) == 5) b_from_top = dynamic_cast<const reco::GenParticle*>(&(*dau));

      if (w_from_top) {
         // Get daughters of top
        auto wdaus = ak8_match.getDaughterQuarks(w_from_top);
        double dr_q1 = reco::deltaR(jet, *wdaus.at(0));
        double dr_q2 = reco::deltaR(jet, *wdaus.at(1));
        // If the 1st daughter is closer to the jet than the 2nd, then swap their places
        if (dr_q1 < dr_q2){
          std::swap(dr_q1, dr_q2);
          std::swap(wdaus.at(0), wdaus.at(1));
        }
        if (debug) {
          using namespace std;
          cout << i << " Top -> Wx -> qqx" << endl;
          cout << "deltaR(jet, q1)    : " << dr_q1 << endl;
          cout << "deltaR(jet, q2)    : " << dr_q2 << endl;
        }
        fatjet_1_dr_T_Wq_max.push_back(dr_q1);
        fatjet_1_T_Wq_max_pdgId.push_back(wdaus.at(0)->pdgId());
      }
      if (b_from_top) {
        double dr_b = reco::deltaR(jet, *b_from_top);
        fatjet_1_dr_T_b.push_back(dr_b);
        if (debug) { 
          cout << i << " Top -> bx" << endl;
          std::cout << "deltaR(jet, b)   : " << dr_b << std::endl;
        }
      } 
    }
  }
 
  tree->Fill();	
  clearVars();
}

void ScoutingNanoAOD::clearVars(){
  fatjet_1_dr_T.clear();
  fatjet_1_dr_T_Wq_max.clear();
  fatjet_1_T_Wq_max_pdgId.clear();
  fatjet_1_dr_W_daus.clear();
  fatjet_1_dr_T_b.clear();
  fatjet_area.clear();
  fatjet_eta.clear();
  fatjet_mass.clear();
  fatjet_phi.clear();
  fatjet_pt.clear();
  fatjet_msoftdrop.clear();
  fatjet_n2b1.clear();
  fatjet_n3b1.clear();
  fatjet_doubleBTag.clear();
  fatjet_mregressed.clear();
  jet_pt.clear();
  jet_eta.clear();
  jet_phi.clear();
  jet_m.clear();
  jet_jetArea.clear();
  jet_bTagScore.clear();
  muon_pt.clear();
  muon_eta.clear();
  muon_phi.clear();
  muon_m.clear();
  muon_type.clear();
  muon_charge.clear();
  muon_trk_dxy.clear();
  muon_trk_dz.clear();
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
