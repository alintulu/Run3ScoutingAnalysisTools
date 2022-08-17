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

  vector<Float16_t> JetHLT_pt;
  vector<Float16_t> JetHLT_mass;
  vector<Float16_t> JetHLT_eta;
  vector<Float16_t> JetHLT_ptGen;
  vector<Float16_t> JetHLT_ptRaw;
  vector<Float16_t> JetHLT_massRaw;
  double JetHLT_rho;
  vector<Float16_t> Jet_pt;
  vector<Float16_t> Jet_mass;
  vector<Float16_t> Jet_eta;
  vector<Float16_t> Jet_ptGen;
  vector<Float16_t> Jet_ptRaw;
  vector<Float16_t> Jet_massRaw;
  double Jet_rho;

  bool debug = false;
  unsigned int debug_match_numJets_HLT = 1000;
  unsigned int debug_match_numJets_scout = 1000;

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

  tree->Branch("JetHLT_pt", &JetHLT_pt);
  tree->Branch("JetHLT_mass", &JetHLT_mass);
  tree->Branch("JetHLT_eta", &JetHLT_eta);
  tree->Branch("JetHLT_ptGen", &JetHLT_ptGen);
  tree->Branch("JetHLT_ptRaw", &JetHLT_ptRaw);
  tree->Branch("JetHLT_massRaw", &JetHLT_massRaw);
  tree->Branch("JetHLT_rho", &JetHLT_rho);
  tree->Branch("Jet_pt", &Jet_pt);
  tree->Branch("Jet_mass", &Jet_mass);
  tree->Branch("Jet_eta", &Jet_eta);
  tree->Branch("Jet_ptGen", &Jet_ptGen);
  tree->Branch("Jet_ptRaw", &Jet_ptRaw);
  tree->Branch("Jet_massRaw", &Jet_massRaw);
  tree->Branch("Jet_rho", &Jet_rho);

}

ScoutingNanoAOD::~ScoutingNanoAOD() {
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
  JetHLT_rho = *hltrhoH;
  //JetHLT_rho = 10;

  auto pdt = iSetup.getHandle(particleTableToken);
  const HepPDT::ParticleDataTable* pdTable = pdt.product();

  std::cout << "HLT jets: " << hltjetH->size() << std::endl;
  std::cout << "HLT corr jets: " << hltjetcorrH->size() << std::endl;

  for(unsigned int j=0; j<hltjetcorrH->size(); j++) {
    JetHLT_pt.push_back((*hltjetcorrH)[j].pt());
    JetHLT_mass.push_back((*hltjetcorrH)[j].mass());
  }

  // Match HLT AK8 jets to GEN jets
  std::map<int, reco::GenJet> resultMapHLT;
  std::vector<int> unmatchedJetsHLT;
  std::vector<std::tuple<int, int, float> > pairListHLT;

  if (debug) std::cout << "dR loop:" << std::endl;

  for(unsigned int i=0; i<hltjetH->size(); i++) {
    if (debug && i > debug_match_numJets_HLT) break;
    bool found_match = false;
    for(unsigned int j=0; j<genjetH->size(); j++) {
      float dR = reco::deltaR((*hltjetH)[i].eta(), (*hltjetH)[i].phi(), (*genjetH)[j].eta(), (*genjetH)[j].phi());
      if (debug) std::cout << i << " " << j << " " << dR << " " << (*hltjetH)[i].pt() << " " << (*genjetH)[j].pt() << std::endl;
      if(dR < 0.8) {
        pairListHLT.push_back(std::make_tuple(i, j, dR));
        found_match = true;
      }
    }
    if(!found_match) {
       unmatchedJetsHLT.push_back(i);
    }
  }

  if (debug) {
    std::cout << "\nunmatchedJets HLT loop:" << std::endl;
    for(auto &j: unmatchedJetsHLT) {
      std::cout << j << std::endl;
    }
    std::cout << "\npairList HLT loop:" << std::endl;
  }

  std::sort(pairListHLT.begin(), pairListHLT.end(), [](std::tuple<int, int, float> t1, std::tuple<int, int, float> t2){ return std::get<2>(t1) < std::get<2>(t2); });

  while(pairListHLT.size() > 0) {
    if (debug) std::cout << std::get<0>(pairListHLT[0]) << " " << std::get<1>(pairListHLT[0]) << " " << std::get<2>(pairListHLT[0]) << std::endl;

    reco::GenJet genjet_assn = (*genjetH)[std::get<1>(pairListHLT[0])];
    resultMapHLT[std::get<0>(pairListHLT[0])] = genjet_assn;
    for(unsigned int k=1; k<pairListHLT.size(); k++) {
      if(std::get<0>(pairListHLT[k]) == std::get<0>(pairListHLT[0]) ||
         std::get<1>(pairListHLT[k]) == std::get<1>(pairListHLT[0])) {
        pairListHLT.erase(pairListHLT.begin() + k);
      }
    }
    pairListHLT.erase(pairListHLT.begin());
  }

  if (debug) {
    std::cout << "\nresultMap HLT loop:" << std::endl;
    for(auto &r: resultMapHLT) {
      std::cout << r.first << std::endl;
    }
  }

  for(unsigned int j=0; j<hltjetH->size(); j++) {
    JetHLT_eta.push_back((*hltjetH)[j].eta());
    JetHLT_ptRaw.push_back((*hltjetH)[j].pt());
    JetHLT_massRaw.push_back((*hltjetH)[j].mass());
   
    if (debug && j > debug_match_numJets_HLT) continue;
 
    if (std::find(unmatchedJetsHLT.begin(), unmatchedJetsHLT.end(), j) != unmatchedJetsHLT.end()) {

      JetHLT_ptGen.push_back(-99);

    } else {
      
      JetHLT_ptGen.push_back(resultMapHLT[j].pt());

    }
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

  // Substructure
  double sd_z_cut = 0.10;
  double sd_beta = 0;
  SoftDrop sd_groomer = SoftDrop(sd_beta, sd_z_cut, 0.8);
  EnergyCorrelatorN2 N2 = EnergyCorrelatorN2(1.0);

  ClusterSequenceArea ak8_cs(fj_part, ak8_def, area_def);
  vector<PseudoJet> ak8_jets = sorted_by_pt(ak8_cs.inclusive_jets());

  std::cout << "Scouting jets: " << ak8_jets.size() << std::endl;

  // Match Scouting AK8 jets to GEN jets
  std::map<int, reco::GenJet> resultMap;
  std::vector<int> unmatchedJets;
  std::vector<std::tuple<int, int, float> > pairList;

  if (debug) std::cout << "dR loop:" << std::endl;

  for(unsigned int i=0; i<ak8_jets.size(); i++) {
    if (debug && i > debug_match_numJets_scout) break;
    bool found_match = false;
    for(unsigned int j=0; j<genjetH->size(); j++) {
      float dR = reco::deltaR(ak8_jets[i].eta(), ak8_jets[i].phi(), (*genjetH)[j].eta(), (*genjetH)[j].phi());
      if (debug) std::cout << i << " " << j << " " << dR << " " << ak8_jets[i].pt() << (*genjetH)[j].pt() << std::endl;
      if(dR < 0.8) {
        pairList.push_back(std::make_tuple(i, j, dR));
        found_match = true;
      }
    }
    if(!found_match) {
       unmatchedJets.push_back(i);
    }
  }

  if (debug) {
    std::cout << "\nunmatchedJets loop:" << std::endl;
    for(auto &j: unmatchedJets) {
      std::cout << j << std::endl;
    }
    std::cout << "\npairList loop:" << std::endl;
  }

  std::sort(pairList.begin(), pairList.end(), [](std::tuple<int, int, float> t1, std::tuple<int, int, float> t2){ return std::get<2>(t1) < std::get<2>(t2); });

  while(pairList.size() > 0) {
    if (debug) std::cout << std::get<0>(pairList[0]) << " " << std::get<1>(pairList[0]) << " " << std::get<2>(pairList[0]) << std::endl;

    reco::GenJet genjet_assn = (*genjetH)[std::get<1>(pairList[0])];
    resultMap[std::get<0>(pairList[0])] = genjet_assn;
    for(unsigned int k=1; k<pairList.size(); k++) {
      if(std::get<0>(pairList[k]) == std::get<0>(pairList[0]) ||
         std::get<1>(pairList[k]) == std::get<1>(pairList[0])) {
        pairList.erase(pairList.begin() + k);
      }
    }
    pairList.erase(pairList.begin());
  }

  if (debug) {
    std::cout << "\nresultMap loop:" << std::endl;
    for(auto &r: resultMap) {
      std::cout << r.first << std::endl;
    }
  }

  for(unsigned int i=0; i<ak8_jets.size(); i++) {

    Jet_pt.push_back(ak8_jets[i].pt());
    Jet_mass.push_back(ak8_jets[i].m());
    Jet_eta.push_back(ak8_jets[i].eta());
    Jet_ptRaw.push_back(ak8_jets[i].pt());
    Jet_massRaw.push_back(ak8_jets[i].m());

    if (debug && i > debug_match_numJets_scout) continue;
    
    if (std::find(unmatchedJets.begin(), unmatchedJets.end(), i) != unmatchedJets.end()) {

      Jet_ptGen.push_back(-99.0);

    } else {
      
      Jet_ptGen.push_back(resultMap[i].pt());

    }
  }

  tree->Fill();	
  clearVars();
}

void ScoutingNanoAOD::clearVars(){
  JetHLT_pt.clear();
  JetHLT_mass.clear();
  JetHLT_eta.clear();
  JetHLT_ptGen.clear();
  JetHLT_ptRaw.clear();
  JetHLT_massRaw.clear();
  Jet_pt.clear();
  Jet_mass.clear();
  Jet_eta.clear();
  Jet_ptGen.clear();
  Jet_ptRaw.clear();
  Jet_massRaw.clear();
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
