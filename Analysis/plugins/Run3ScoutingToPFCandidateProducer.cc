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
  produces<std::vector<reco::PFCandidate>>();
  produces<edm::ValueMap<float>>("normchi2");
  produces<edm::ValueMap<float>>("dz");
  produces<edm::ValueMap<float>>("dxy");
  produces<edm::ValueMap<float>>("dzsig");
  produces<edm::ValueMap<float>>("dxysig");
  produces<edm::ValueMap<uint8_t>>("lostInnerHits");
  produces<edm::ValueMap<uint8_t>>("quality");
  produces<edm::ValueMap<float>>("trk_pt");
  produces<edm::ValueMap<float>>("trk_eta");
  produces<edm::ValueMap<float>>("trk_phi");
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

  std::vector<float> normchi2(scoutingparticleHandle->size());
  std::vector<float> dz(scoutingparticleHandle->size());
  std::vector<float> dxy(scoutingparticleHandle->size());
  std::vector<float> dzsig(scoutingparticleHandle->size());
  std::vector<float> dxysig(scoutingparticleHandle->size());
  std::vector<uint8_t> lostInnerHits(scoutingparticleHandle->size());
  std::vector<uint8_t> quality(scoutingparticleHandle->size());
  std::vector<float> trk_pt(scoutingparticleHandle->size());
  std::vector<float> trk_eta(scoutingparticleHandle->size());
  std::vector<float> trk_phi(scoutingparticleHandle->size());

  auto pfcands = std::make_unique<std::vector<reco::PFCandidate>>(scoutingparticleHandle->size());
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

      normchi2[icand] = scoutingparticle.normchi2();
      dz[icand] = scoutingparticle.dz();
      dxy[icand] = scoutingparticle.dxy();
      dzsig[icand] = scoutingparticle.dzsig();
      dxysig[icand] = scoutingparticle.dxysig();
      lostInnerHits[icand] = scoutingparticle.lostInnerHits();
      quality[icand] = scoutingparticle.quality();
      trk_pt[icand] = scoutingparticle.trk_pt();
      trk_eta[icand] = scoutingparticle.trk_eta();
      trk_phi[icand] = scoutingparticle.trk_phi();

      if (debug_) print(scoutingparticle, pfcand);
  }

  ProductID const pidK0(1, 0);
  std::vector<reco::PFCandidate> const *const_pfcands = &(*pfcands);
  iEvent.put(std::move(pfcands)); 
  
  std::unique_ptr<edm::ValueMap<float>> normchi2V(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_normchi2(*normchi2V);
  filler_normchi2.insert(OrphanHandle<std::vector<reco::PFCandidate>>(const_pfcands, pidK0), normchi2.begin(), normchi2.end());
  filler_normchi2.fill();
  iEvent.put(std::move(normchi2V), "normchi2");

  std::unique_ptr<edm::ValueMap<float>> dzV(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_dz(*dzV);
  filler_dz.insert(OrphanHandle<std::vector<reco::PFCandidate>>(const_pfcands, pidK0), dz.begin(), dz.end());
  filler_dz.fill();
  iEvent.put(std::move(dzV), "dz");

  std::unique_ptr<edm::ValueMap<float>> dxyV(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_dxy(*dxyV);
  filler_dxy.insert(OrphanHandle<std::vector<reco::PFCandidate>>(const_pfcands, pidK0), dxy.begin(), dxy.end());
  filler_dxy.fill();
  iEvent.put(std::move(dxyV), "dxy");

  std::unique_ptr<edm::ValueMap<float>> dzsigV(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_dzsig(*dzsigV);
  filler_dzsig.insert(OrphanHandle<std::vector<reco::PFCandidate>>(const_pfcands, pidK0), dzsig.begin(), dzsig.end());
  filler_dzsig.fill();
  iEvent.put(std::move(dzsigV), "dzsig");

  std::unique_ptr<edm::ValueMap<float>> dxysigV(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_dxysig(*dxysigV);
  filler_dxysig.insert(OrphanHandle<std::vector<reco::PFCandidate>>(const_pfcands, pidK0), dxysig.begin(), dxysig.end());
  filler_dxysig.fill();
  iEvent.put(std::move(dxysigV), "dxysig");

  std::unique_ptr<edm::ValueMap<uint8_t>> lostInnerHitsV(new edm::ValueMap<uint8_t>());
  edm::ValueMap<uint8_t>::Filler filler_lostInnerHits(*lostInnerHitsV);
  filler_lostInnerHits.insert(OrphanHandle<std::vector<reco::PFCandidate>>(const_pfcands, pidK0), lostInnerHits.begin(), lostInnerHits.end());
  filler_lostInnerHits.fill();
  iEvent.put(std::move(lostInnerHitsV), "lostInnerHits");

  std::unique_ptr<edm::ValueMap<uint8_t>> qualityV(new edm::ValueMap<uint8_t>());
  edm::ValueMap<uint8_t>::Filler filler_quality(*qualityV);
  filler_quality.insert(OrphanHandle<std::vector<reco::PFCandidate>>(const_pfcands, pidK0), quality.begin(), quality.end());
  filler_quality.fill();
  iEvent.put(std::move(qualityV), "quality");

  std::unique_ptr<edm::ValueMap<float>> trk_ptV(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_trk_pt(*trk_ptV);
  filler_trk_pt.insert(OrphanHandle<std::vector<reco::PFCandidate>>(const_pfcands, pidK0), trk_pt.begin(), trk_pt.end());
  filler_trk_pt.fill();
  iEvent.put(std::move(trk_ptV), "trk_pt");

  std::unique_ptr<edm::ValueMap<float>> trk_etaV(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_trk_eta(*trk_etaV);
  filler_trk_eta.insert(OrphanHandle<std::vector<reco::PFCandidate>>(const_pfcands, pidK0), trk_eta.begin(), trk_eta.end());
  filler_trk_eta.fill();
  iEvent.put(std::move(trk_etaV), "trk_eta");

  std::unique_ptr<edm::ValueMap<float>> trk_phiV(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_trk_phi(*trk_phiV);
  filler_trk_phi.insert(OrphanHandle<std::vector<reco::PFCandidate>>(const_pfcands, pidK0), trk_phi.begin(), trk_phi.end());
  filler_trk_phi.fill();
  iEvent.put(std::move(trk_phiV), "trk_phi");
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Run3ScoutingToPFCandidateProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("scoutingparticle", edm::InputTag("hltScoutingPFPacker"));
  descriptions.addWithDefaultLabel(desc);
}

// declare this class as a framework plugin
DEFINE_FWK_MODULE(Run3ScoutingToPFCandidateProducer);
