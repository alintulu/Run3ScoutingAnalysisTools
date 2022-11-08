// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Common/interface/OrphanHandle.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/Math/interface/libminifloat.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

class Run3ScoutingToPFCandidateProducer : public edm::global::EDProducer<> {
public:
  explicit Run3ScoutingToPFCandidateProducer(const edm::ParameterSet &);
  ~Run3ScoutingToPFCandidateProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
  void produce(edm::StreamID sid, edm::Event & iEvent, edm::EventSetup const & setup) const final;
  void print(Run3ScoutingParticle s, reco::PFCandidate r) const;

private:
  const edm::EDGetTokenT<std::vector<Run3ScoutingParticle>> input_scoutingparticle_token_;
  const edm::ESGetToken<HepPDT::ParticleDataTable, edm::DefaultRecord> particletable_token_;
  bool debug_;
};

//
// constructors and destructor
//
Run3ScoutingToPFCandidateProducer::Run3ScoutingToPFCandidateProducer(const edm::ParameterSet &iConfig)
    : input_scoutingparticle_token_(consumes(iConfig.getParameter<edm::InputTag>("scoutingparticle"))),
      particletable_token_(esConsumes<HepPDT::ParticleDataTable, edm::DefaultRecord>()),
      debug_(iConfig.existsAs<bool>("debug") ? iConfig.getParameter<bool>("debug") : false) {
  //register products
  produces<reco::PFCandidateCollection>();
  produces<edm::ValueMap<int>>("vertexIndex");
  produces<edm::ValueMap<float>>("normchi2");
  produces<edm::ValueMap<float>>("dz");
  produces<edm::ValueMap<float>>("dxy");
  produces<edm::ValueMap<float>>("dzsig");
  produces<edm::ValueMap<float>>("dxysig");
  produces<edm::ValueMap<int>>("lostInnerHits");
  produces<edm::ValueMap<int>>("quality");
  produces<edm::ValueMap<float>>("trkPt");
  produces<edm::ValueMap<float>>("trkEta");
  produces<edm::ValueMap<float>>("trkPhi");
}

Run3ScoutingToPFCandidateProducer::~Run3ScoutingToPFCandidateProducer() = default;

void Run3ScoutingToPFCandidateProducer::print(Run3ScoutingParticle s, reco::PFCandidate r) const {

  std::cout << "pdgId" << std::endl; 
  std::cout << s.pdgId() << " " << r.pdgId() << std::endl;  
  std::cout << "pt" << std::endl; 
  std::cout << s.pt() << " " << r.pt() << std::endl;  
  std::cout << "phi" << std::endl; 
  std::cout << s.phi() << " " << r.phi() << std::endl;  
  std::cout << "eta" << std::endl; 
  std::cout << s.eta() << " " << r.eta() << std::endl;  

}

// ------------ method called to produce the data  ------------
void Run3ScoutingToPFCandidateProducer::produce(edm::StreamID sid, edm::Event & iEvent, edm::EventSetup const & setup) const {
  using namespace edm;

  auto pdt = setup.getHandle(particletable_token_);
  const HepPDT::ParticleDataTable* pdTable = pdt.product();

  Handle<std::vector<Run3ScoutingParticle>> scoutingparticleHandle;
  iEvent.getByToken(input_scoutingparticle_token_, scoutingparticleHandle);

  std::vector<int8_t> vertexIndex(scoutingparticleHandle->size());
  std::vector<float> normchi2(scoutingparticleHandle->size());
  std::vector<float> dz(scoutingparticleHandle->size());
  std::vector<float> dxy(scoutingparticleHandle->size());
  std::vector<float> dzsig(scoutingparticleHandle->size());
  std::vector<float> dxysig(scoutingparticleHandle->size());
  std::vector<int> lostInnerHits(scoutingparticleHandle->size());
  std::vector<int> quality(scoutingparticleHandle->size());
  std::vector<float> trkPt(scoutingparticleHandle->size());
  std::vector<float> trkEta(scoutingparticleHandle->size());
  std::vector<float> trkPhi(scoutingparticleHandle->size());

  auto pfcands = std::make_unique<reco::PFCandidateCollection>(scoutingparticleHandle->size());
  for (unsigned int icand = 0; icand < scoutingparticleHandle->size(); ++icand) {

      auto& pfcand = dynamic_cast<reco::PFCandidate&>((*pfcands)[icand]);
      auto& scoutingparticle = (*scoutingparticleHandle)[icand];

      auto m = pdTable->particle(HepPDT::ParticleID(scoutingparticle.pdgId())) != nullptr
                        ? pdTable->particle(HepPDT::ParticleID(scoutingparticle.pdgId()))->mass()
                        : -99.f;
      auto q = pdTable->particle(HepPDT::ParticleID(scoutingparticle.pdgId())) != nullptr
                        ? pdTable->particle(HepPDT::ParticleID(scoutingparticle.pdgId()))->charge()
                        : -99.f;
      if (m < -90 ||  q < -90) continue;
      
      float px = scoutingparticle.pt() * cos(scoutingparticle.phi());
      float py = scoutingparticle.pt() * sin(scoutingparticle.phi());
      float pz = scoutingparticle.pt() * sinh(scoutingparticle.eta());
      float p = scoutingparticle.pt() * cosh(scoutingparticle.eta());
      float energy = std::sqrt(p*p + m*m);
      reco::Particle::LorentzVector p4(px, py, pz, energy); 
 
      pfcand = reco::PFCandidate(q, p4, pfcand.translatePdgIdToType(scoutingparticle.pdgId()));

      bool relativeTrackVars = scoutingparticle.relative_trk_vars();
      vertexIndex[icand] = scoutingparticle.vertex();
      normchi2[icand] = scoutingparticle.normchi2();
      dz[icand] = scoutingparticle.dz();
      dxy[icand] = scoutingparticle.dxy();
      dzsig[icand] = scoutingparticle.dzsig();
      dxysig[icand] = scoutingparticle.dxysig();
      lostInnerHits[icand] = scoutingparticle.lostInnerHits();
      quality[icand] = scoutingparticle.quality();
      trkPt[icand] = relativeTrackVars ? scoutingparticle.trk_pt() + scoutingparticle.pt() : scoutingparticle.trk_pt();
      trkEta[icand] = relativeTrackVars ? scoutingparticle.trk_eta() + scoutingparticle.eta() : scoutingparticle.trk_eta();
      trkPhi[icand] = relativeTrackVars ? scoutingparticle.trk_phi() + scoutingparticle.phi() : scoutingparticle.trk_phi();

      if (debug_) print(scoutingparticle, pfcand);
  }

  //ProductID const pidK0(1, 0);
  //reco::PFCandidateCollection const *const_pfcands = &(*pfcands);
  edm::OrphanHandle<reco::PFCandidateCollection> oh = iEvent.put(std::move(pfcands));
  //edm::OrphanHandle<reco::PFCandidateCollection> oh = iEvent.emplace(ptokenPFCandidates_, pfcands);
  
  std::unique_ptr<edm::ValueMap<int>> vertexIndex_VM(new edm::ValueMap<int>());
  edm::ValueMap<int>::Filler filler_vertexIndex(*vertexIndex_VM);
  filler_vertexIndex.insert(oh, vertexIndex.begin(), vertexIndex.end());
  filler_vertexIndex.fill();
  iEvent.put(std::move(vertexIndex_VM), "vertexIndex");
  
  std::unique_ptr<edm::ValueMap<float>> normchi2_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_normchi2(*normchi2_VM);
  filler_normchi2.insert(oh, normchi2.begin(), normchi2.end());
  filler_normchi2.fill();
  iEvent.put(std::move(normchi2_VM), "normchi2");

  std::unique_ptr<edm::ValueMap<float>> dz_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_dz(*dz_VM);
  filler_dz.insert(oh, dz.begin(), dz.end());
  filler_dz.fill();
  iEvent.put(std::move(dz_VM), "dz");

  std::unique_ptr<edm::ValueMap<float>> dxy_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_dxy(*dxy_VM);
  filler_dxy.insert(oh, dxy.begin(), dxy.end());
  filler_dxy.fill();
  iEvent.put(std::move(dxy_VM), "dxy");

  std::unique_ptr<edm::ValueMap<float>> dzsig_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_dzsig(*dzsig_VM);
  filler_dzsig.insert(oh, dzsig.begin(), dzsig.end());
  filler_dzsig.fill();
  iEvent.put(std::move(dzsig_VM), "dzsig");

  std::unique_ptr<edm::ValueMap<float>> dxysig_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_dxysig(*dxysig_VM);
  filler_dxysig.insert(oh, dxysig.begin(), dxysig.end());
  filler_dxysig.fill();
  iEvent.put(std::move(dxysig_VM), "dxysig");

  std::unique_ptr<edm::ValueMap<int>> lostInnerHits_VM(new edm::ValueMap<int>());
  edm::ValueMap<int>::Filler filler_lostInnerHits(*lostInnerHits_VM);
  filler_lostInnerHits.insert(oh, lostInnerHits.begin(), lostInnerHits.end());
  filler_lostInnerHits.fill();
  iEvent.put(std::move(lostInnerHits_VM), "lostInnerHits");

  std::unique_ptr<edm::ValueMap<int>> quality_VM(new edm::ValueMap<int>());
  edm::ValueMap<int>::Filler filler_quality(*quality_VM);
  filler_quality.insert(oh, quality.begin(), quality.end());
  filler_quality.fill();
  iEvent.put(std::move(quality_VM), "quality");

  std::unique_ptr<edm::ValueMap<float>> trkPt_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_trkPt(*trkPt_VM);
  filler_trkPt.insert(oh, trkPt.begin(), trkPt.end());
  filler_trkPt.fill();
  iEvent.put(std::move(trkPt_VM), "trkPt");

  std::unique_ptr<edm::ValueMap<float>> trkEta_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_trkEta(*trkEta_VM);
  filler_trkEta.insert(oh, trkEta.begin(), trkEta.end());
  filler_trkEta.fill();
  iEvent.put(std::move(trkEta_VM), "trkEta");

  std::unique_ptr<edm::ValueMap<float>> trkPhi_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_trkPhi(*trkPhi_VM);
  filler_trkPhi.insert(oh, trkPhi.begin(), trkPhi.end());
  filler_trkPhi.fill();
  iEvent.put(std::move(trkPhi_VM), "trkPhi");
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Run3ScoutingToPFCandidateProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("scoutingparticle", edm::InputTag("hltScoutingPFPacker"));
  descriptions.addWithDefaultLabel(desc);
}

// declare this class as a framework plugin
DEFINE_FWK_MODULE(Run3ScoutingToPFCandidateProducer);
