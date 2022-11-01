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

      if (debug_) print(scoutingparticle, pfcand);
  }

  //put output
  iEvent.put(std::move(pfcands));
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Run3ScoutingToPFCandidateProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("scoutingparticle", edm::InputTag("hltScoutingPFPacker"));
  descriptions.addWithDefaultLabel(desc);
}

// declare this class as a framework plugin
DEFINE_FWK_MODULE(Run3ScoutingToPFCandidateProducer);
