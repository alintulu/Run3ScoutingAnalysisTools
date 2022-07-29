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


  // TTree carrying the event weight information
  TTree* tree;

  vector<Float16_t> fj_1_dr_T;
  vector<Float16_t> fj_1_dr_T_Wq_max;
  vector<int> fj_1_T_Wq_max_pdgId;
  vector<Float16_t> fj_1_dr_W_daus;
  vector<Float16_t> fj_1_dr_T_b;
  
  bool debug = true;
};

ScoutingNanoAOD::ScoutingNanoAOD(const edm::ParameterSet& iConfig):
  pfcandsParticleNetToken  (consumes<std::vector<Run3ScoutingParticle> > (iConfig.getParameter<edm::InputTag>("pfcandsParticleNet"))),
  genpartsToken            (consumes<reco::GenParticleCollection> (iConfig.getParameter<edm::InputTag>("genpart")))
{
  usesResource("TFileService");

 // Access the TFileService
  edm::Service<TFileService> fs;

  // Create the TTree
  tree = fs->make<TTree>("tree", "tree");

  tree->Branch("fj_1_dr_T",&fj_1_dr_T);
  tree->Branch("fj_1_dr_T_Wq_max",&fj_1_dr_T_Wq_max);
  tree->Branch("fj_1_T_Wq_max_pdgId",&fj_1_T_Wq_max_pdgId);
  tree->Branch("fj_1_dr_W_daus",&fj_1_dr_W_daus);
  tree->Branch("fj_1_dr_T_b",&fj_1_dr_T_b);
}

ScoutingNanoAOD::~ScoutingNanoAOD() {
}

std::pair<const reco::GenParticle*, double> ScoutingNanoAOD::closest(fastjet::PseudoJet jet, std::vector<const reco::GenParticle*>& gps, bool debug){

 const reco::GenParticle *ret = nullptr;
 double drMin = 1e6;
 for (const auto &gp : gps) {
   double dr = reco::deltaR(jet.eta(), jet.phi(), gp->eta(), gp->phi());
   if (debug) std::cout << gp->pdgId() << dr << std::endl;
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
  vector<PseudoJet> ak8_jets = sorted_by_pt(ak8_cs.inclusive_jets(170.0));

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

    // Check if it is top
    if (std::abs(gp->pdgId()) == 6) {
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
           if (debug) std::cout << "Top -> Wx -> qqx" << std::endl;
         }
      }
    // Check if it is the W
    } else if (std::abs(gp->pdgId()) == 24) {

        // Check if the W decays hadronically
        if(ak8_match.isHadronic(gp)) {
          hadGenWs.push_back(gp);
          if (debug) std::cout << "W -> qq" << std::endl;
        }
     }
  }

  // Loop over jets
  for (auto &jet: ak8_jets) {
    // Find closest top to jet
    auto closestT = closest(jet, hadGenTops, debug);
    auto genT = closestT.first;
    auto drT = closestT.second;
    if (debug) std::cout << "deltaR(jet, top)   : " << drT << std::endl;

    // Find closest W to jet
    auto closestW = closest(jet, hadGenWs, debug);
    auto genW = closestW.first;
    auto drW = closestW.second;
    if (debug) std::cout << "deltaR(jet, W)   : " << drW << std::endl;

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
        cout << "deltaR(jet, q1)    : " << dr_q1 << endl;
        cout << "deltaR(jet, q2)    : " << dr_q2 << endl;
      }
      fj_1_dr_W_daus.push_back(dr_q1);
    } else {
      fj_1_dr_W_daus.push_back(99.0); 
    }

    fj_1_dr_T.push_back(drT);

    // If there is no top
    if (!genT) {
      fj_1_dr_T_Wq_max.push_back(99.0);
      fj_1_T_Wq_max_pdgId.push_back(0);
      fj_1_dr_T_b.push_back(99.0);
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
          cout << "deltaR(jet, q1)    : " << dr_q1 << endl;
          cout << "deltaR(jet, q2)    : " << dr_q2 << endl;
        }
        fj_1_dr_T_Wq_max.push_back(dr_q1);
        fj_1_T_Wq_max_pdgId.push_back(wdaus.at(0)->pdgId());
      }
      if (b_from_top) {
        double dr_b = reco::deltaR(jet, *b_from_top);
        fj_1_dr_T_b.push_back(dr_b);
        if (debug) std::cout << "deltaR(jet, b)   : " << dr_b << std::endl;
      } 
    }
  }
 
  tree->Fill();	
  clearVars();
}

void ScoutingNanoAOD::clearVars(){
  fj_1_dr_T.clear();
  fj_1_dr_T_Wq_max.clear();
  fj_1_T_Wq_max_pdgId.clear();
  fj_1_dr_W_daus.clear();
  fj_1_dr_T_b.clear();
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
