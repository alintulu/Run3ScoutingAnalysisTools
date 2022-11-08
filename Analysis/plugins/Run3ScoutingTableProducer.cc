#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include <string> 

#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPhoton.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingElectron.h"
#include "DataFormats/Scouting/interface/Run3ScoutingTrack.h"

class Run3ScoutingTableProducer : public edm::stream::EDProducer<> {
public:
  explicit Run3ScoutingTableProducer(const edm::ParameterSet &);
  ~Run3ScoutingTableProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void produce(edm::Event &, const edm::EventSetup &) override;
   
  const edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> input_primaryvertex_token_;
  const edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> input_displacedvertex_token_;
  const edm::EDGetTokenT<std::vector<Run3ScoutingPhoton>> input_photon_token_;
  const edm::EDGetTokenT<std::vector<Run3ScoutingMuon>> input_muon_token_;
  const edm::EDGetTokenT<std::vector<Run3ScoutingElectron>> input_electron_token_;
  const edm::EDGetTokenT<std::vector<Run3ScoutingTrack>> input_track_token_;
  const edm::EDGetTokenT<double> input_metpt_token_;
  const edm::EDGetTokenT<double> input_metphi_token_;
  const edm::EDGetTokenT<double> input_rho_token_;
  
};

//
// constructors and destructor
//
Run3ScoutingTableProducer::Run3ScoutingTableProducer(const edm::ParameterSet &iConfig)
    : input_primaryvertex_token_(consumes(iConfig.getParameter<edm::InputTag>("primaryvertex"))), 
      input_displacedvertex_token_(consumes(iConfig.getParameter<edm::InputTag>("displacedvertex"))), 
      input_photon_token_(consumes(iConfig.getParameter<edm::InputTag>("photon"))),
      input_muon_token_(consumes(iConfig.getParameter<edm::InputTag>("muon"))),
      input_electron_token_(consumes(iConfig.getParameter<edm::InputTag>("electron"))),
      input_track_token_(consumes(iConfig.getParameter<edm::InputTag>("track"))),
      input_metpt_token_(consumes<double>(iConfig.getParameter<edm::InputTag>("metpt"))),
      input_metphi_token_ (consumes<double>(iConfig.getParameter<edm::InputTag>("metphi"))),
      input_rho_token_ (consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))) {
  produces<nanoaod::FlatTable>("ScoutingPrimaryVertex");
  produces<nanoaod::FlatTable>("ScoutingDisplacedVertex");
  produces<nanoaod::FlatTable>("ScoutingPhoton");
  produces<nanoaod::FlatTable>("ScoutingMuon");
  produces<nanoaod::FlatTable>("ScoutingElectron");
  produces<nanoaod::FlatTable>("ScoutingTrack");
  produces<nanoaod::FlatTable>("ScoutingMET");
  produces<nanoaod::FlatTable>("ScoutingRho");
}

Run3ScoutingTableProducer::~Run3ScoutingTableProducer() {}

void Run3ScoutingTableProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  
  auto vertex = iEvent.getHandle(input_primaryvertex_token_);
  std::vector<float> vertex_x;
  std::vector<float> vertex_y;
  std::vector<float> vertex_z;
  std::vector<float> vertex_xError;
  std::vector<float> vertex_yError;
  std::vector<float> vertex_zError;
  std::vector<int> vertex_tracksSize;
  std::vector<float> vertex_chi2;
  std::vector<bool> vertex_isValidVtx;

  for (unsigned i = 0; i < vertex->size(); ++i) {
     vertex_x.push_back((*vertex)[i].x());
     vertex_y.push_back((*vertex)[i].y());
     vertex_z.push_back((*vertex)[i].z());
     vertex_xError.push_back((*vertex)[i].xError());
     vertex_yError.push_back((*vertex)[i].yError());
     vertex_zError.push_back((*vertex)[i].zError());
     vertex_tracksSize.push_back((*vertex)[i].tracksSize());
     vertex_chi2.push_back((*vertex)[i].chi2());
     vertex_isValidVtx.push_back((*vertex)[i].isValidVtx());
  }
      
  auto vertexTable = std::make_unique<nanoaod::FlatTable>(vertex->size(), "ScoutingPrimaryVertex", false, false);
  vertexTable->setDoc("PFScouting primary vertex, i.e. hltPixelVertices");
  vertexTable->addColumn<float>("x", vertex_x, "vertex x");
  vertexTable->addColumn<float>("y", vertex_y, "vertex y");
  vertexTable->addColumn<float>("z", vertex_z, "vertex z");
  vertexTable->addColumn<float>("xError", vertex_xError, "vertex x error");
  vertexTable->addColumn<float>("yError", vertex_yError, "vertex y error");
  vertexTable->addColumn<float>("zError", vertex_zError, "vertex z error");
  vertexTable->addColumn<float>("tracksSize", vertex_tracksSize, "track size");
  vertexTable->addColumn<float>("chi2", vertex_chi2, "vertex chi2");
  vertexTable->addColumn<float>("isValidVertex", vertex_isValidVtx, "boolean if vertex is valid");

  auto dispvertex = iEvent.getHandle(input_displacedvertex_token_);
  std::vector<float> dispvertex_x;
  std::vector<float> dispvertex_y;
  std::vector<float> dispvertex_z;
  std::vector<float> dispvertex_xError;
  std::vector<float> dispvertex_yError;
  std::vector<float> dispvertex_zError;
  std::vector<int> dispvertex_tracksSize;
  std::vector<float> dispvertex_chi2;
  std::vector<bool> dispvertex_isValidVtx;

  for (unsigned i = 0; i < dispvertex->size(); ++i) {
     dispvertex_x.push_back((*dispvertex)[i].x());
     dispvertex_y.push_back((*dispvertex)[i].y());
     dispvertex_z.push_back((*dispvertex)[i].z());
     dispvertex_xError.push_back((*dispvertex)[i].xError());
     dispvertex_yError.push_back((*dispvertex)[i].yError());
     dispvertex_zError.push_back((*dispvertex)[i].zError());
     dispvertex_tracksSize.push_back((*dispvertex)[i].tracksSize());
     dispvertex_chi2.push_back((*dispvertex)[i].chi2());
     dispvertex_isValidVtx.push_back((*dispvertex)[i].isValidVtx());
  }
      
  auto dispvertexTable = std::make_unique<nanoaod::FlatTable>(dispvertex->size(), "ScoutingDisplacedVertex", false, false);
  dispvertexTable->setDoc("PFScouting displaced vertex, i.e. hltScoutingMuonPacker");
  dispvertexTable->addColumn<float>("x", dispvertex_x, "vertex x");
  dispvertexTable->addColumn<float>("y", dispvertex_y, "vertex y");
  dispvertexTable->addColumn<float>("z", dispvertex_z, "vertex z");
  dispvertexTable->addColumn<float>("xError", dispvertex_xError, "vertex x error");
  dispvertexTable->addColumn<float>("yError", dispvertex_yError, "vertex y error");
  dispvertexTable->addColumn<float>("zError", dispvertex_zError, "vertex z error");
  dispvertexTable->addColumn<float>("tracksSize", dispvertex_tracksSize, "track size");
  dispvertexTable->addColumn<float>("chi2", dispvertex_chi2, "vertex chi2");
  dispvertexTable->addColumn<float>("isValidVertex", dispvertex_isValidVtx, "boolean if vertex is valid");

  auto photon = iEvent.getHandle(input_photon_token_);
  std::vector<float> photon_pt;
  std::vector<float> photon_eta;
  std::vector<float> photon_phi;
  std::vector<float> photon_m;
  std::vector<float> photon_sigmaIetaIeta;
  std::vector<float> photon_hOverE;
  std::vector<float> photon_ecalIso;
  std::vector<float> photon_hcalIso;
  std::vector<float> photon_trkIso;
  std::vector<float> photon_r9;
  std::vector<float> photon_sMin;
  std::vector<float> photon_sMaj;
  std::vector<uint32_t> photon_seedId;
  std::vector<std::vector<float>> photon_energyMatrix;
  std::vector<std::vector<uint32_t>> photon_detIds;
  std::vector<std::vector<float>> photon_timingMatrix;
  std::vector<bool> photon_rechitZeroSuppression;
    
  for (unsigned i = 0; i < photon->size(); ++i) {
     photon_pt.push_back((*photon)[i].pt());
     photon_eta.push_back((*photon)[i].eta());
     photon_phi.push_back((*photon)[i].phi());
     photon_m.push_back((*photon)[i].m());
     photon_sigmaIetaIeta.push_back((*photon)[i].sigmaIetaIeta());
     photon_hOverE.push_back((*photon)[i].hOverE());
     photon_ecalIso.push_back((*photon)[i].ecalIso());
     photon_hcalIso.push_back((*photon)[i].hcalIso());
     photon_trkIso.push_back((*photon)[i].trkIso());
     photon_r9.push_back((*photon)[i].r9());
     photon_sMin.push_back((*photon)[i].sMin());
     photon_sMaj.push_back((*photon)[i].sMaj());
     photon_seedId.push_back((*photon)[i].seedId());
     photon_energyMatrix.push_back((*photon)[i].energyMatrix());
     photon_detIds.push_back((*photon)[i].detIds());
     photon_timingMatrix.push_back((*photon)[i].timingMatrix());
     photon_rechitZeroSuppression.push_back((*photon)[i].rechitZeroSuppression());
  }
      
  auto photonTable = std::make_unique<nanoaod::FlatTable>(photon->size(), "ScoutingPhoton", false, false);
  photonTable->setDoc("PFScouting photon, i.e. hltScoutingEgammaPacker");
  photonTable->addColumn<float>("pt", photon_pt, "photon pt");
  photonTable->addColumn<float>("eta", photon_eta, "photon eta");
  photonTable->addColumn<float>("phi", photon_phi, "photon phi");
  photonTable->addColumn<float>("m", photon_m, "photon x error");
  photonTable->addColumn<float>("sigmaIetaIeta", photon_sigmaIetaIeta, "photon sigmaIetaIeta");
  photonTable->addColumn<float>("hOverE", photon_hOverE, "photon hOverE");
  photonTable->addColumn<float>("ecalIso", photon_ecalIso, "photon ecalIso");
  photonTable->addColumn<float>("hcalIso", photon_hcalIso, "photon hcalIso");
  photonTable->addColumn<float>("trkiso", photon_trkIso, "photon trkIso");
  photonTable->addColumn<float>("r9", photon_r9, "photon ecalIso");
  photonTable->addColumn<float>("sMin", photon_sMin, "photon sMaj");
  photonTable->addColumn<float>("sMaj", photon_sMaj, "photon sMin");
  photonTable->addColumn<uint32_t>("seedId", photon_seedId, "photon seedId");
  //photonTable->addColumn<std::vector<float>>("energyMatrix", photon_energyMatrix, "photon energyMatrix");
  //photonTable->addColumn<std::vector<uint32_t>>("detIds", photon_detIds, "photon detIds");
  //photonTable->addColumn<std::vector<float>>("timingMatrix", photon_timingMatrix, "photon timingMatrix");
  //photonTable->addColumn<bool>("rechitZeroSuppression", photon_rechitZeroSuppression, "photon rechitZeroSuppression");

  auto muon = iEvent.getHandle(input_muon_token_);
  std::vector<float> muon_pt;
  std::vector<float> muon_eta;
  std::vector<float> muon_phi;
  std::vector<float> muon_m;
  std::vector<unsigned int> muon_type;
  std::vector<int> muon_charge;
  std::vector<float> muon_normchi2;
  std::vector<float> muon_ecalIso;
  std::vector<float> muon_hcalIso;
  std::vector<float> muon_trkIso;
  std::vector<int> muon_nValidStandAloneMuonHits;
  std::vector<int> muon_nStandAloneMuonMatchedStations;
  std::vector<int> muon_nValidRecoMuonHits;
  std::vector<int> muon_nRecoMuonChambers;
  std::vector<int> muon_nRecoMuonChambersCSCorDT;
  std::vector<int> muon_nRecoMuonMatches;
  std::vector<int> muon_nRecoMuonMatchedStations;
  std::vector<unsigned int> muon_nRecoMuonExpectedMatchedStations;
  std::vector<unsigned int> muon_recoMuonStationMask;
  std::vector<int> muon_nRecoMuonMatchedRPCLayers;
  std::vector<unsigned int> muon_recoMuonRPClayerMask;
  std::vector<int> muon_nValidPixelHits;
  std::vector<int> muon_nValidStripHits;
  std::vector<int> muon_nPixelLayersWithMeasurement;
  std::vector<int> muon_nTrackerLayersWithMeasurement;
  std::vector<float> muon_trk_chi2;
  std::vector<float> muon_trk_ndof;
  std::vector<float> muon_trk_dxy;
  std::vector<float> muon_trk_dz;
  std::vector<float> muon_trk_qoverp;
  std::vector<float> muon_trk_lambda;
  std::vector<float> muon_trk_pt;
  std::vector<float> muon_trk_phi;
  std::vector<float> muon_trk_eta;
  std::vector<float> muon_trk_dxyError;
  std::vector<float> muon_trk_dzError;
  std::vector<float> muon_trk_qoverpError;
  std::vector<float> muon_trk_lambdaError;
  std::vector<float> muon_trk_phiError;
  std::vector<float> muon_trk_dsz;
  std::vector<float> muon_trk_dszError;
  std::vector<float> muon_trk_qoverp_lambda_cov;
  std::vector<float> muon_trk_qoverp_phi_cov;
  std::vector<float> muon_trk_qoverp_dxy_cov;
  std::vector<float> muon_trk_qoverp_dsz_cov;
  std::vector<float> muon_trk_lambda_phi_cov;
  std::vector<float> muon_trk_lambda_dxy_cov;
  std::vector<float> muon_trk_lambda_dsz_cov;
  std::vector<float> muon_trk_phi_dxy_cov;
  std::vector<float> muon_trk_phi_dsz_cov;
  std::vector<float> muon_trk_dxy_dsz_cov;
  std::vector<float> muon_trk_vx;
  std::vector<float> muon_trk_vy;
  std::vector<float> muon_trk_vz;
    
  for (unsigned i = 0; i < muon->size(); ++i) {
     muon_pt.push_back((*muon)[i].pt());
     muon_eta.push_back((*muon)[i].eta());
     muon_phi.push_back((*muon)[i].phi());
     muon_m.push_back((*muon)[i].m());
     muon_type.push_back((*muon)[i].type());
     muon_charge.push_back((*muon)[i].charge());
     muon_normchi2.push_back((*muon)[i].normalizedChi2());
     muon_ecalIso.push_back((*muon)[i].ecalIso());
     muon_hcalIso.push_back((*muon)[i].hcalIso());
     muon_trkIso.push_back((*muon)[i].trackIso());
     muon_nValidStandAloneMuonHits.push_back((*muon)[i].nValidStandAloneMuonHits());
     muon_nStandAloneMuonMatchedStations.push_back((*muon)[i].nStandAloneMuonMatchedStations());
     muon_nValidRecoMuonHits.push_back((*muon)[i].nValidRecoMuonHits());
     muon_nRecoMuonChambers.push_back((*muon)[i].nRecoMuonChambers());
     muon_nRecoMuonChambersCSCorDT.push_back((*muon)[i].nRecoMuonChambersCSCorDT());
     muon_nRecoMuonMatches.push_back((*muon)[i].nRecoMuonMatches());
     muon_nRecoMuonMatchedStations.push_back((*muon)[i].nRecoMuonMatchedStations());
     muon_nRecoMuonExpectedMatchedStations.push_back((*muon)[i].nRecoMuonExpectedMatchedStations());
     muon_recoMuonStationMask.push_back((*muon)[i].recoMuonStationMask());
     muon_nRecoMuonMatchedRPCLayers.push_back((*muon)[i].nRecoMuonMatchedRPCLayers());
     muon_recoMuonRPClayerMask.push_back((*muon)[i].recoMuonRPClayerMask());
     muon_nValidPixelHits.push_back((*muon)[i].nValidPixelHits());
     muon_nValidStripHits.push_back((*muon)[i].nValidStripHits());
     muon_nPixelLayersWithMeasurement.push_back((*muon)[i].nPixelLayersWithMeasurement());
     muon_nTrackerLayersWithMeasurement.push_back((*muon)[i].nTrackerLayersWithMeasurement());
     muon_trk_chi2.push_back((*muon)[i].trk_chi2());
     muon_trk_ndof.push_back((*muon)[i].trk_ndof());
     muon_trk_dxy.push_back((*muon)[i].trk_dxy());
     muon_trk_dz.push_back((*muon)[i].trk_dz());
     muon_trk_qoverp.push_back((*muon)[i].trk_qoverp());
     muon_trk_lambda.push_back((*muon)[i].trk_lambda());
     muon_trk_pt.push_back((*muon)[i].trk_pt());
     muon_trk_phi.push_back((*muon)[i].trk_phi());
     muon_trk_eta.push_back((*muon)[i].trk_eta());
     muon_trk_dxyError.push_back((*muon)[i].trk_dxyError());
     muon_trk_dzError.push_back((*muon)[i].trk_dzError());
     muon_trk_qoverpError.push_back((*muon)[i].trk_qoverpError());
     muon_trk_lambdaError.push_back((*muon)[i].trk_lambdaError());
     muon_trk_phiError.push_back((*muon)[i].trk_phiError());
     muon_trk_dsz.push_back((*muon)[i].trk_dsz());
     muon_trk_dszError.push_back((*muon)[i].trk_dszError());
     muon_trk_qoverp_lambda_cov.push_back((*muon)[i].trk_qoverp_lambda_cov());
     muon_trk_qoverp_phi_cov.push_back((*muon)[i].trk_qoverp_phi_cov());
     muon_trk_qoverp_dxy_cov.push_back((*muon)[i].trk_qoverp_dxy_cov());
     muon_trk_qoverp_dsz_cov.push_back((*muon)[i].trk_qoverp_dsz_cov());
     muon_trk_lambda_phi_cov.push_back((*muon)[i].trk_lambda_phi_cov());
     muon_trk_lambda_dxy_cov.push_back((*muon)[i].trk_lambda_dxy_cov());
     muon_trk_lambda_dsz_cov.push_back((*muon)[i].trk_lambda_dsz_cov());
     muon_trk_vx.push_back((*muon)[i].trk_vx());
     muon_trk_vy.push_back((*muon)[i].trk_vy());
     muon_trk_vz.push_back((*muon)[i].trk_vz());
  }
      
  auto muonTable = std::make_unique<nanoaod::FlatTable>(muon->size(), "ScoutingMuon", false, false);
  muonTable->setDoc("PFScouting muon, i.e. hltScoutingMuonPacker");
  muonTable->addColumn<float>("pt", muon_pt, "muon pt");
  muonTable->addColumn<float>("eta", muon_eta, "muon eta");
  muonTable->addColumn<float>("phi", muon_phi, "muon phi");
  muonTable->addColumn<float>("m", muon_m, "muon x error");
  muonTable->addColumn<unsigned int>("type", muon_type, "muon type");
  muonTable->addColumn<int>("charge", muon_charge, "muon charge");
  muonTable->addColumn<float>("normchi2", muon_normchi2, "muon normchi2");
  muonTable->addColumn<float>("ecalIso", muon_ecalIso, "muon ecalIso");
  muonTable->addColumn<float>("hcalIso", muon_hcalIso, "muon hcalIso");
  muonTable->addColumn<float>("trkiso", muon_trkIso, "muon trkIso");
  muonTable->addColumn<int>("nValidStandAloneMuonHits", muon_nValidStandAloneMuonHits, "muon nValidStandAloneMuonHits");
  muonTable->addColumn<int>("nStandAloneMuonMatchedStations", muon_nStandAloneMuonMatchedStations, "muon nStandAloneMuonMatchedStations");
  muonTable->addColumn<int>("nValidRecoMuonHits", muon_nValidRecoMuonHits, "muon nValidRecoMuonHits");
  muonTable->addColumn<int>("nRecoMuonChambers", muon_nRecoMuonChambers, "muon nRecoMuonChambers");
  muonTable->addColumn<int>("nRecoMuonChambersCSCorDT", muon_nRecoMuonChambersCSCorDT, "muon nRecoMuonChambersCSCorDT");
  muonTable->addColumn<int>("nRecoMuonMatches", muon_nRecoMuonMatches, "muon nRecoMuonMatches");
  muonTable->addColumn<int>("nRecoMuonMatchesStations", muon_nRecoMuonMatchedStations, "muon nRecoMuonMatchesStations");
  muonTable->addColumn<unsigned int>("nRecoMuonExpectedMatchedStations", muon_nRecoMuonExpectedMatchedStations, "muon nRecoMuonExpectedMatchedStations");
  muonTable->addColumn<unsigned int>("recoMuonStationMask", muon_recoMuonStationMask, "muon recoMuonStationMask");
  muonTable->addColumn<int>("nRecoMuonMatchedRPCLayers", muon_nRecoMuonMatchedRPCLayers, "muon nRecoMuonMatchedRPCLayers");
  muonTable->addColumn<unsigned int>("recoMuonRPClayerMask", muon_recoMuonRPClayerMask, "muon recoMuonRPClayerMask");
  muonTable->addColumn<int>("nValidPixelHits", muon_nValidPixelHits, "muon nValidPixelHits");
  muonTable->addColumn<int>("nValidStripHits", muon_nValidStripHits, "muon nValidStripHits");
  muonTable->addColumn<int>("nPixelLayersWithMeasurement", muon_nPixelLayersWithMeasurement, "muon nPixelLayersWithMeasurement");
  muonTable->addColumn<float>("nTrackerLayersWithMeasurement", muon_nTrackerLayersWithMeasurement, "muon nTrackerLayersWithMeasurement");
  muonTable->addColumn<float>("trk_chi2", muon_trk_chi2, "muon trk chi2");
  muonTable->addColumn<float>("trk_ndof", muon_trk_ndof, "muon trk ndof");
  muonTable->addColumn<float>("trk_dxy", muon_trk_dxy, "muon trk dxy");
  muonTable->addColumn<float>("trk_dz", muon_trk_dz, "muon trk dz");
  muonTable->addColumn<float>("trk_qoverp", muon_trk_qoverp, "muon trk qoverp");
  muonTable->addColumn<float>("trk_lambda", muon_trk_lambda, "muon trk lambda");
  muonTable->addColumn<float>("trk_pt", muon_trk_pt, "muon trk pt");
  muonTable->addColumn<float>("trk_phi", muon_trk_phi, "muon trk phi");
  muonTable->addColumn<float>("trk_dxyError", muon_trk_dxyError, "muon trk dxyError");
  muonTable->addColumn<float>("trk_dzError", muon_trk_dzError, "muon trk dzError");
  muonTable->addColumn<float>("trk_qoverpError", muon_trk_qoverpError, "muon trk qoverpError");
  muonTable->addColumn<float>("trk_lambdaError", muon_trk_lambdaError, "muon trk lambdaError");
  muonTable->addColumn<float>("trk_phiError", muon_trk_phiError, "muon trk phiError");
  muonTable->addColumn<float>("trk_dsz", muon_trk_dsz, "muon trk dsz");
  muonTable->addColumn<float>("trk_dszError", muon_trk_dszError, "muon trk dszError");
  muonTable->addColumn<float>("trk_qoverp_lambda_cov", muon_trk_qoverp_lambda_cov, "muon trk qoverp lambda cov");
  muonTable->addColumn<float>("trk_qoverp_phi_cov", muon_trk_qoverp_phi_cov, "muon trk qoverp phi cov");
  muonTable->addColumn<float>("trk_qoverp_dxy_cov", muon_trk_qoverp_dxy_cov, "muon trk qoverp dxy cov");
  muonTable->addColumn<float>("trk_qoverp_dsz_cov", muon_trk_qoverp_dsz_cov, "muon trk qoverp dsz cov");
  muonTable->addColumn<float>("trk_lambda_phi_cov", muon_trk_lambda_phi_cov, "muon trk lambda phi cov");
  muonTable->addColumn<float>("trk_lambda_dxy_cov", muon_trk_lambda_dxy_cov, "muon trk lambda dxy cov");
  muonTable->addColumn<float>("trk_lambda_dsz_cov", muon_trk_lambda_dsz_cov, "muon trk lambda dsz cov");
  muonTable->addColumn<float>("trk_vx", muon_trk_vx, "muon trk vx");
  muonTable->addColumn<float>("trk_vy", muon_trk_vy, "muon trk vy");
  muonTable->addColumn<float>("trk_vz", muon_trk_vz, "muon trk vz");

  auto electron = iEvent.getHandle(input_electron_token_);
  std::vector<float> electron_pt;
  std::vector<float> electron_eta;
  std::vector<float> electron_phi;
  std::vector<float> electron_m;
  std::vector<float> electron_d0;
  std::vector<float> electron_dz;
  std::vector<float> electron_dEtaIn;
  std::vector<float> electron_dPhiIn;
  std::vector<float> electron_sigmaIetaIeta;
  std::vector<float> electron_hOverE;
  std::vector<float> electron_ooEMOop;
  std::vector<float> electron_missingHits;
  std::vector<float> electron_charge;
  std::vector<float> electron_ecalIso;
  std::vector<float> electron_hcalIso;
  std::vector<float> electron_trkIso;
  std::vector<float> electron_r9;
  std::vector<float> electron_sMin;
  std::vector<float> electron_sMaj;
  std::vector<uint32_t> electron_seedId;
  std::vector<std::vector<float>> electron_energyMatrix;
  std::vector<std::vector<uint32_t>> electron_detIds;
  std::vector<std::vector<float>> electron_timingMatrix;
  std::vector<bool> electron_rechitZeroSuppression;
    
  for (unsigned i = 0; i < electron->size(); ++i) {
     electron_pt.push_back((*electron)[i].pt());
     electron_eta.push_back((*electron)[i].eta());
     electron_phi.push_back((*electron)[i].phi());
     electron_m.push_back((*electron)[i].m());
     electron_d0.push_back((*electron)[i].d0());
     electron_dz.push_back((*electron)[i].dz());
     electron_dEtaIn.push_back((*electron)[i].dEtaIn());
     electron_dPhiIn.push_back((*electron)[i].dPhiIn());
     electron_sigmaIetaIeta.push_back((*electron)[i].sigmaIetaIeta());
     electron_hOverE.push_back((*electron)[i].hOverE());
     electron_ooEMOop.push_back((*electron)[i].ooEMOop());
     electron_missingHits.push_back((*electron)[i].missingHits());
     electron_charge.push_back((*electron)[i].charge());
     electron_ecalIso.push_back((*electron)[i].ecalIso());
     electron_hcalIso.push_back((*electron)[i].hcalIso());
     electron_trkIso.push_back((*electron)[i].trackIso());
     electron_r9.push_back((*electron)[i].r9());
     electron_sMin.push_back((*electron)[i].sMin());
     electron_sMaj.push_back((*electron)[i].sMaj());
     electron_seedId.push_back((*electron)[i].seedId());
     electron_energyMatrix.push_back((*electron)[i].energyMatrix());
     electron_detIds.push_back((*electron)[i].detIds());
     electron_timingMatrix.push_back((*electron)[i].timingMatrix());
     electron_rechitZeroSuppression.push_back((*electron)[i].rechitZeroSuppression());
  }
      
  auto electronTable = std::make_unique<nanoaod::FlatTable>(electron->size(), "ScoutingElectron", false, false);
  electronTable->setDoc("PFScouting electron, i.e. hltScoutingEgammaPacker");
  electronTable->addColumn<float>("pt", electron_pt, "electron pt");
  electronTable->addColumn<float>("eta", electron_eta, "electron eta");
  electronTable->addColumn<float>("phi", electron_phi, "electron phi");
  electronTable->addColumn<float>("m", electron_m, "electron x error");
  electronTable->addColumn<float>("d0", electron_d0, "electron d0");
  electronTable->addColumn<float>("dz", electron_dz, "electron dz");
  electronTable->addColumn<float>("dEtaIn", electron_dEtaIn, "electron dEtaIn");
  electronTable->addColumn<float>("dPhiIn", electron_dPhiIn, "electron x error");
  electronTable->addColumn<float>("sigmaIetaIeta", electron_sigmaIetaIeta, "electron sigmaIetaIeta");
  electronTable->addColumn<float>("hOverE", electron_hOverE, "electron hOverE");
  electronTable->addColumn<float>("ooEMOop", electron_ooEMOop, "electron ooEMOop");
  electronTable->addColumn<float>("missingHits", electron_missingHits, "electron missingHits");
  electronTable->addColumn<float>("charge", electron_charge, "electron charge");
  electronTable->addColumn<float>("ecalIso", electron_ecalIso, "electron ecalIso");
  electronTable->addColumn<float>("hcalIso", electron_hcalIso, "electron hcalIso");
  electronTable->addColumn<float>("trkiso", electron_trkIso, "electron trkIso");
  electronTable->addColumn<float>("r9", electron_r9, "electron ecalIso");
  electronTable->addColumn<float>("sMin", electron_sMin, "electron sMaj");
  electronTable->addColumn<float>("sMaj", electron_sMaj, "electron sMin");
  electronTable->addColumn<uint32_t>("seedId", electron_seedId, "electron seedId");

  auto track = iEvent.getHandle(input_track_token_);
  std::vector<float> track_pt;
  std::vector<float> track_eta;
  std::vector<float> track_phi;
  std::vector<float> track_chi2;
  std::vector<float> track_ndof;
  std::vector<int> track_charge;
  std::vector<float> track_dxy;
  std::vector<float> track_dz;
  std::vector<int> track_nValidPixelHits;
  std::vector<int> track_nTrackerLayersWithMeasurement;
  std::vector<int> track_nValidStripHits;
  std::vector<float> track_qoverp;
  std::vector<float> track_lambda;
  std::vector<float> track_dxyError;
  std::vector<float> track_dzError;
  std::vector<float> track_qoverpError;
  std::vector<float> track_lambdaError;
  std::vector<float> track_phiError;
  std::vector<float> track_dsz;
  std::vector<float> track_dszError;
  std::vector<float> track_qoverp_lambda_cov;
  std::vector<float> track_qoverp_phi_cov;
  std::vector<float> track_qoverp_dxy_cov;
  std::vector<float> track_qoverp_dsz_cov;
  std::vector<float> track_lambda_phi_cov;
  std::vector<float> track_lambda_dxy_cov;
  std::vector<float> track_lambda_dsz_cov;
  std::vector<float> track_phi_dxy_cov;
  std::vector<float> track_phi_dsz_cov;
  std::vector<float> track_dxy_dsz_cov;
  std::vector<int> track_vtxInd;
  std::vector<float> track_vx;
  std::vector<float> track_vy;
  std::vector<float> track_vz;
    
  for (unsigned i = 0; i < track->size(); ++i) {
     track_pt.push_back((*track)[i].tk_pt());
     track_eta.push_back((*track)[i].tk_eta());
     track_phi.push_back((*track)[i].tk_phi());
     track_chi2.push_back((*track)[i].tk_chi2());
     track_ndof.push_back((*track)[i].tk_ndof());
     track_charge.push_back((*track)[i].tk_charge());
     track_dxy.push_back((*track)[i].tk_dxy());
     track_dz.push_back((*track)[i].tk_dz());
     track_nValidPixelHits.push_back((*track)[i].tk_nValidPixelHits());
     track_nTrackerLayersWithMeasurement.push_back((*track)[i].tk_nTrackerLayersWithMeasurement());
     track_nValidStripHits.push_back((*track)[i].tk_nValidStripHits());
     track_qoverp.push_back((*track)[i].tk_qoverp());
     track_lambda.push_back((*track)[i].tk_lambda());
     track_dxyError.push_back((*track)[i].tk_dxy_Error());
     track_dzError.push_back((*track)[i].tk_dz_Error());
     track_qoverpError.push_back((*track)[i].tk_qoverp_Error());
     track_lambdaError.push_back((*track)[i].tk_lambda_Error());
     track_phiError.push_back((*track)[i].tk_phi_Error());
     track_dsz.push_back((*track)[i].tk_dsz());
     track_dszError.push_back((*track)[i].tk_dsz_Error());
     track_qoverp_lambda_cov.push_back((*track)[i].tk_qoverp_lambda_cov());
     track_qoverp_phi_cov.push_back((*track)[i].tk_qoverp_phi_cov());
     track_qoverp_dxy_cov.push_back((*track)[i].tk_qoverp_dxy_cov());
     track_qoverp_dsz_cov.push_back((*track)[i].tk_qoverp_dsz_cov());
     track_lambda_phi_cov.push_back((*track)[i].tk_lambda_phi_cov());
     track_lambda_dxy_cov.push_back((*track)[i].tk_lambda_dxy_cov());
     track_lambda_dsz_cov.push_back((*track)[i].tk_lambda_dsz_cov());
     track_vtxInd.push_back((*track)[i].tk_vtxInd());
     track_vx.push_back((*track)[i].tk_vx());
     track_vy.push_back((*track)[i].tk_vy());
     track_vz.push_back((*track)[i].tk_vz());
  }
      
  auto trackTable = std::make_unique<nanoaod::FlatTable>(track->size(), "ScoutingTrack", false, false);
  trackTable->setDoc("PFScouting track, i.e. hltPixelOnlyPFMuonMerging");
  trackTable->addColumn<float>("pt", track_pt, "track pt");
  trackTable->addColumn<float>("eta", track_eta, "track eta");
  trackTable->addColumn<float>("phi", track_phi, "track phi");
  trackTable->addColumn<float>("chi2", track_chi2, "track trk chi2");
  trackTable->addColumn<float>("ndof", track_ndof, "track trk ndof");
  trackTable->addColumn<int>("charge", track_charge, "track charge");
  trackTable->addColumn<float>("dxy", track_dxy, "track trk dxy");
  trackTable->addColumn<float>("dz", track_dz, "track trk dz");
  trackTable->addColumn<int>("nValidPixelHits", track_nValidPixelHits, "track nValidPixelHits");
  trackTable->addColumn<int>("nTrackerLayersWithMeasurement", track_nTrackerLayersWithMeasurement, "track nTrackerLayersWithMeasurement");
  trackTable->addColumn<int>("nValidStripHits", track_nValidStripHits, "track nValidStripHits");
  trackTable->addColumn<float>("qoverp", track_qoverp, "track trk qoverp");
  trackTable->addColumn<float>("lambda", track_lambda, "track trk lambda");
  trackTable->addColumn<float>("dxyError", track_dxyError, "track trk dxyError");
  trackTable->addColumn<float>("dzError", track_dzError, "track trk dzError");
  trackTable->addColumn<float>("qoverpError", track_qoverpError, "track trk qoverpError");
  trackTable->addColumn<float>("lambdaError", track_lambdaError, "track trk lambdaError");
  trackTable->addColumn<float>("phiError", track_phiError, "track trk phiError");
  trackTable->addColumn<float>("dsz", track_dsz, "track trk dsz");
  trackTable->addColumn<float>("dszError", track_dszError, "track trk dszError");
  trackTable->addColumn<float>("qoverp_lambda_cov", track_qoverp_lambda_cov, "track trk qoverp lambda cov");
  trackTable->addColumn<float>("qoverp_phi_cov", track_qoverp_phi_cov, "track trk qoverp phi cov");
  trackTable->addColumn<float>("qoverp_dxy_cov", track_qoverp_dxy_cov, "track trk qoverp dxy cov");
  trackTable->addColumn<float>("qoverp_dsz_cov", track_qoverp_dsz_cov, "track trk qoverp dsz cov");
  trackTable->addColumn<float>("lambda_phi_cov", track_lambda_phi_cov, "track trk lambda phi cov");
  trackTable->addColumn<float>("lambda_dxy_cov", track_lambda_dxy_cov, "track trk lambda dxy cov");
  trackTable->addColumn<float>("lambda_dsz_cov", track_lambda_dsz_cov, "track trk lambda dsz cov");
  trackTable->addColumn<float>("vtxInd", track_vtxInd, "track trk vtxInd");
  trackTable->addColumn<float>("vx", track_vx, "track trk vx");
  trackTable->addColumn<float>("vy", track_vy, "track trk vy");
  trackTable->addColumn<float>("vz", track_vz, "track trk vz");
 
  auto metTable = std::make_unique<nanoaod::FlatTable>(2, "ScoutingMET", true, false);
  metTable->setDoc("PFScouting MET, i.e. hltPixelOnlyPFMETProducer");
  metTable->addColumnValue<double>("pt", *iEvent.getHandle(input_metpt_token_), "met pt");
  metTable->addColumnValue<double>("phi", *iEvent.getHandle(input_metphi_token_), "met phi");

  auto rhoTable = std::make_unique<nanoaod::FlatTable>(1, "ScoutingRho", true, false);
  rhoTable->setDoc("PFScouting rho, i.e. hltFixedGridRhoFastjetPixelOnlyAll");
  rhoTable->addColumnValue<double>("", *iEvent.getHandle(input_rho_token_), "rho");
  
  iEvent.put(std::move(vertexTable), "ScoutingPrimaryVertex");
  iEvent.put(std::move(dispvertexTable), "ScoutingDisplacedVertex");
  iEvent.put(std::move(photonTable), "ScoutingPhoton");
  iEvent.put(std::move(muonTable), "ScoutingMuon");
  iEvent.put(std::move(electronTable), "ScoutingElectron");
  iEvent.put(std::move(trackTable), "ScoutingTrack");
  iEvent.put(std::move(metTable), "ScoutingMET");
  iEvent.put(std::move(rhoTable), "ScoutingRho");
}

void Run3ScoutingTableProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("primaryvertex", edm::InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx"));
  desc.add<edm::InputTag>("displacedvertex", edm::InputTag("hltScoutingMuonPacker","displacedVtx"));
  desc.add<edm::InputTag>("photon", edm::InputTag("hltScoutingEgammaPacker"));
  desc.add<edm::InputTag>("muon", edm::InputTag("hltScoutingMuonPacker"));
  desc.add<edm::InputTag>("electron", edm::InputTag("hltScoutingEgammaPacker"));
  desc.add<edm::InputTag>("track", edm::InputTag("hltScoutingTrackPacker"));
  desc.add<edm::InputTag>("metpt", edm::InputTag("hltScoutingPFPacker","pfMetPt"));
  desc.add<edm::InputTag>("metphi", edm::InputTag("hltScoutingPFPacker","pfMetPhi"));
  desc.add<edm::InputTag>("rho", edm::InputTag("hltScoutingPFPacker", "rho"));
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(Run3ScoutingTableProducer);
