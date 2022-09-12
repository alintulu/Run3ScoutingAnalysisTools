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
#include "DataFormats/JetReco/interface/PFJet.h"

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

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

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
  bool isNeutralPdg(int);
  std::pair<std::map<int, int>, std::vector<int>> match_pseudo(std::vector<fastjet::PseudoJet>, const std::vector<reco::GenJet>&);
  std::pair<std::map<int, int>, std::vector<int>> match_reco(const std::vector<reco::PFJet>&, const std::vector<reco::GenJet>&);
  const edm::InputTag triggerResultsTag;
  const edm::EDGetTokenT<std::vector<Run3ScoutingParticle> >  	scoutpartToken;
  const edm::EDGetTokenT<reco::GenJetCollection> genjetToken;
  const edm::EDGetTokenT<double>  	scoutrhoToken;
  const edm::EDGetTokenT<std::vector<reco::PFJet> >  	hltjetToken;
  const edm::EDGetTokenT<std::vector<reco::PFJet> >  	hltjetcorrToken;
  const edm::EDGetTokenT<double>  	hltrhoToken;
  const edm::ESGetToken<HepPDT::ParticleDataTable, edm::DefaultRecord> particleTableToken;

  // TTree carrying the event weight information
  TTree* tree;

  vector<Float16_t> HLTJet_pt;
  vector<Float16_t> HLTJet_mass;
  vector<Float16_t> HLTJet_eta;
  vector<Float16_t> HLTJet_phi;
  vector<Float16_t> HLTJet_ptRaw;
  vector<Float16_t> HLTJet_massRaw;
  vector<Float16_t> HLTJet_area;
  vector<int> HLTJet_genJetAK8Idx;
  double HLTJet_rho;
  vector<Float16_t> Jet_pt;
  vector<Float16_t> Jet_mass;
  vector<Float16_t> Jet_eta;
  vector<Float16_t> Jet_phi;
  vector<Float16_t> Jet_area;
  vector<int> Jet_genJetAK8Idx;
  double Jet_rho;
  vector<Float16_t> GenJet_pt;
  vector<Float16_t> GenJet_mass;
  vector<Float16_t> GenJet_eta;
  vector<Float16_t> GenJet_phi;

};

ScoutingNanoAOD::ScoutingNanoAOD(const edm::ParameterSet& iConfig):
  scoutpartToken  (consumes<std::vector<Run3ScoutingParticle> > (iConfig.getParameter<edm::InputTag>("scoutpart"))),
  genjetToken                (consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("ak8genjet"))),
  scoutrhoToken             (consumes<double>(iConfig.getParameter<edm::InputTag>("scoutrho"))),
  hltjetToken  (consumes<std::vector<reco::PFJet> > (iConfig.getParameter<edm::InputTag>("ak8hltjet"))),
  hltjetcorrToken  (consumes<std::vector<reco::PFJet> > (iConfig.getParameter<edm::InputTag>("ak8hltcorrjet"))),
  hltrhoToken             (consumes<double>(iConfig.getParameter<edm::InputTag>("hltrho"))),
  particleTableToken       (esConsumes<HepPDT::ParticleDataTable, edm::DefaultRecord>())
{
  usesResource("TFileService");

 // Access the TFileService
  edm::Service<TFileService> fs;

  // Create the TTree
  tree = fs->make<TTree>("tree", "tree");

  tree->Branch("HLTJet_pt", &HLTJet_pt);
  tree->Branch("HLTJet_mass", &HLTJet_mass);
  tree->Branch("HLTJet_eta", &HLTJet_eta);
  tree->Branch("HLTJet_phi", &HLTJet_phi);
  tree->Branch("HLTJet_ptRaw", &HLTJet_ptRaw);
  tree->Branch("HLTJet_massRaw", &HLTJet_massRaw);
  tree->Branch("HLTJet_rho", &HLTJet_rho);
  tree->Branch("HLTJet_area", &HLTJet_area);
  tree->Branch("HLTJet_genJetAK8Idx", &HLTJet_genJetAK8Idx);
  tree->Branch("Jet_pt", &Jet_pt);
  tree->Branch("Jet_mass", &Jet_mass);
  tree->Branch("Jet_eta", &Jet_eta);
  tree->Branch("Jet_phi", &Jet_phi);
  tree->Branch("Jet_rho", &Jet_rho);
  tree->Branch("Jet_area", &Jet_area);
  tree->Branch("Jet_genJetAK8Idx", &Jet_genJetAK8Idx);
  tree->Branch("GenJet_pt", &GenJet_pt);
  tree->Branch("GenJet_mass", &GenJet_mass);
  tree->Branch("GenJet_eta", &GenJet_eta);
  tree->Branch("GenJet_phi", &GenJet_phi);

}

ScoutingNanoAOD::~ScoutingNanoAOD() {
}

std::pair<std::map<int, int>, std::vector<int>> ScoutingNanoAOD::match_pseudo(std::vector<fastjet::PseudoJet> jets, const std::vector<reco::GenJet>& genjets) {

  std::map<int, int> resultMap;
  std::vector<int> unmatchedJets;
  std::vector<std::tuple<int, int, float> > pairList;

  for(unsigned int i=0; i < jets.size(); i++) {
    bool found_match = false;
    for(unsigned int j=0; j < genjets.size(); j++) {
      float dR = reco::deltaR(jets[i].eta(), jets[i].phi(), genjets[j].eta(), genjets[j].phi());
      if (dR < 0.4) {
        pairList.push_back(std::make_tuple(i, j, dR));
        found_match = true;
      }
    }
    if (!found_match) {
      unmatchedJets.push_back(i);
    }
  }

  std::sort(pairList.begin(), pairList.end(), [](std::tuple<int, int, float> t1, std::tuple<int, int, float> t2){ return std::get<2>(t1) < std::get<2>(t2); });

  while(pairList.size() > 0) {

    resultMap[std::get<0>(pairList[0])] = std::get<1>(pairList[0]);
    
    for(unsigned int k=1; k < pairList.size(); k++) {
      if (std::get<0>(pairList[k]) == std::get<0>(pairList[0]) ||
         std::get<1>(pairList[k]) == std::get<1>(pairList[0])) {
        pairList.erase(pairList.begin() + k);
      }
    }
    pairList.erase(pairList.begin());
  }
  return std::make_pair(resultMap, unmatchedJets);
}

std::pair<std::map<int, int>, std::vector<int>> ScoutingNanoAOD::match_reco(const std::vector<reco::PFJet>& jets, const std::vector<reco::GenJet>& genjets) {

  std::map<int, int> resultMap;
  std::vector<int> unmatchedJets;
  std::vector<std::tuple<int, int, float> > pairList;

  for(unsigned int i=0; i < jets.size(); i++) {
    bool found_match = false;
    for(unsigned int j=0; j < genjets.size(); j++) {
      float dR = reco::deltaR(jets[i].eta(), jets[i].phi(), genjets[j].eta(), genjets[j].phi());
      if (dR < 0.4) {
        pairList.push_back(std::make_tuple(i, j, dR));
        found_match = true;
      }
    }
    if (!found_match) {
      unmatchedJets.push_back(i);
    }
  }

  std::sort(pairList.begin(), pairList.end(), [](std::tuple<int, int, float> t1, std::tuple<int, int, float> t2){ return std::get<2>(t1) < std::get<2>(t2); });

  while(pairList.size() > 0) {

    resultMap[std::get<0>(pairList[0])] = std::get<1>(pairList[0]);
    
    for(unsigned int k=1; k < pairList.size(); k++) {
      if (std::get<0>(pairList[k]) == std::get<0>(pairList[0]) ||
         std::get<1>(pairList[k]) == std::get<1>(pairList[0])) {
        pairList.erase(pairList.begin() + k);
      }
    }
    pairList.erase(pairList.begin());
  }
  return std::make_pair(resultMap, unmatchedJets);
}
void ScoutingNanoAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace fastjet;
  using namespace fastjet::contrib;

  Handle<vector<Run3ScoutingParticle> > scoutpartH;
  iEvent.getByToken(scoutpartToken, scoutpartH);

  Handle<GenJetCollection> genjetH;
  iEvent.getByToken(genjetToken, genjetH);

  Handle<double> scoutrhoH;
  iEvent.getByToken(scoutrhoToken, scoutrhoH);
  Jet_rho = *scoutrhoH;

  Handle<vector<reco::PFJet> > hltjetH;
  iEvent.getByToken(hltjetToken, hltjetH);

  Handle<vector<reco::PFJet> > hltjetcorrH;
  iEvent.getByToken(hltjetcorrToken, hltjetcorrH);

  Handle<double> hltrhoH;
  iEvent.getByToken(hltrhoToken, hltrhoH);
  HLTJet_rho = *hltrhoH;

  auto pdt = iSetup.getHandle(particleTableToken);
  const HepPDT::ParticleDataTable* pdTable = pdt.product();

  //std::cout << "HLT jets: " << hltjetH->size() << std::endl;
  //std::cout << "HLT corr jets: " << hltjetcorrH->size() << std::endl;
  
  auto hltjet_match = match_reco(*hltjetH, *genjetH);
  auto hltjet_resultMap = hltjet_match.first;
  auto hltjet_unmatchedJets = hltjet_match.second;

  for(unsigned int j=0; j<hltjetcorrH->size(); j++) {
    HLTJet_pt.push_back((*hltjetcorrH)[j].pt());
    HLTJet_mass.push_back((*hltjetcorrH)[j].mass());
  }

  for(unsigned int j=0; j<hltjetH->size(); j++) {
    HLTJet_eta.push_back((*hltjetH)[j].eta());
    HLTJet_phi.push_back((*hltjetH)[j].phi());
    HLTJet_ptRaw.push_back((*hltjetH)[j].pt());
    HLTJet_massRaw.push_back((*hltjetH)[j].mass());
    HLTJet_area.push_back((*hltjetH)[j].jetArea());

    if (std::find(hltjet_unmatchedJets.begin(), hltjet_unmatchedJets.end(), j) != hltjet_unmatchedJets.end()) {
      HLTJet_genJetAK8Idx.push_back(-99);
    } else {
      HLTJet_genJetAK8Idx.push_back(hltjet_resultMap[j]);
    }
  }

  for(unsigned int j=0; j<genjetH->size(); j++) {
    GenJet_eta.push_back((*genjetH)[j].eta());
    GenJet_phi.push_back((*genjetH)[j].phi());
    GenJet_pt.push_back((*genjetH)[j].pt());
    GenJet_mass.push_back((*genjetH)[j].mass()); 
  }
 
  // Create Scouting AK8 Jet
  vector<PseudoJet> fj_part;
  fj_part.reserve(scoutpartH->size());
  int pfcand_i = 0;
  for (auto pfcands_iter = scoutpartH->begin(); pfcands_iter != scoutpartH->end(); ++pfcands_iter) {

    auto m = pdTable->particle(HepPDT::ParticleID(pfcands_iter->pdgId())) != nullptr
                        ? pdTable->particle(HepPDT::ParticleID(pfcands_iter->pdgId()))->mass()
                        : -99.f;

    if (m < 0) continue;

    math::PtEtaPhiMLorentzVector p4(pfcands_iter->pt(), pfcands_iter->eta(), pfcands_iter->phi(), m);
    fj_part.emplace_back(p4.px(), p4.py(), p4.pz(), p4.energy());
    fj_part.back().set_user_index(pfcand_i);
    pfcand_i++;
  }

  JetDefinition ak8_def = JetDefinition(antikt_algorithm, 0.8);
  fastjet::GhostedAreaSpec area_spec(5.0,1,0.01);
  fastjet::AreaDefinition area_def(fastjet::active_area, area_spec);

  ClusterSequenceArea ak8_cs(fj_part, ak8_def, area_def);
  vector<PseudoJet> ak8_jets = sorted_by_pt(ak8_cs.inclusive_jets());

  auto jet_match = match_pseudo(ak8_jets, *genjetH);
  auto jet_resultMap = jet_match.first;
  auto jet_unmatchedJets = jet_match.second;

  //std::cout << "Scouting jets: " << ak8_jets.size() << std::endl;

  for(unsigned int i=0; i < ak8_jets.size(); i++) {

    Jet_pt.push_back(ak8_jets[i].pt());
    Jet_mass.push_back(ak8_jets[i].m());
    Jet_eta.push_back(ak8_jets[i].eta());
    Jet_phi.push_back(ak8_jets[i].phi());
    Jet_area.push_back(ak8_jets[i].area());

    if (std::find(jet_unmatchedJets.begin(), jet_unmatchedJets.end(), i) != jet_unmatchedJets.end()) {
      Jet_genJetAK8Idx.push_back(-99);
    } else {
      Jet_genJetAK8Idx.push_back(jet_resultMap[i]);
    }
  }

  tree->Fill();	
  clearVars();
}

void ScoutingNanoAOD::clearVars(){
  HLTJet_pt.clear();
  HLTJet_mass.clear();
  HLTJet_eta.clear();
  HLTJet_phi.clear();
  HLTJet_ptRaw.clear();
  HLTJet_massRaw.clear();
  HLTJet_area.clear();
  HLTJet_genJetAK8Idx.clear();
  Jet_pt.clear();
  Jet_mass.clear();
  Jet_eta.clear();
  Jet_phi.clear();
  Jet_area.clear();
  Jet_genJetAK8Idx.clear();
  GenJet_pt.clear();
  GenJet_mass.clear();
  GenJet_eta.clear();
  GenJet_phi.clear();
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
