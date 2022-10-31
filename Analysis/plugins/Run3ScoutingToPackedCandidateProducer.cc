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
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/Math/interface/libminifloat.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

class Run3ScoutingToPackedCandidateProducer : public edm::global::EDProducer<> {
public:
  explicit Run3ScoutingToPackedCandidateProducer(const edm::ParameterSet &);
  ~Run3ScoutingToPackedCandidateProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
  void produce(edm::StreamID sid, edm::Event & iEvent, edm::EventSetup const & setup) const final;

private:

  const edm::EDGetTokenT<std::vector<Run3ScoutingParticle>> input_scoutingparticle_token_;
  const edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> input_scoutingvertex_input_;
  const edm::ESGetToken<HepPDT::ParticleDataTable, edm::DefaultRecord> particleTableToken;
};

//
// constructors and destructor
//
Run3ScoutingToPackedCandidateProducer::Run3ScoutingToPackedCandidateProducer(const edm::ParameterSet &iConfig)
    : input_scoutingparticle_token_(consumes(iConfig.getParameter<edm::InputTag>("scoutingparticle"))),
      input_scoutingvertex_input_(consumes(iConfig.getParameter<edm::InputTag>("scoutingvertex"))),
      particleTableToken       (esConsumes<HepPDT::ParticleDataTable, edm::DefaultRecord>()) {
  //register products
  produces<std::vector<pat::PackedCandidate>>();
}

Run3ScoutingToPackedCandidateProducer::~Run3ScoutingToPackedCandidateProducer() = default;

// ------------ method called to produce the data  ------------
void Run3ScoutingToPackedCandidateProducer::produce(edm::StreamID sid, edm::Event & iEvent, edm::EventSetup const & setup) const {
  using namespace edm;

  auto pdt = setup.getHandle(particleTableToken);
  const HepPDT::ParticleDataTable* pdTable = pdt.product();

  Handle<std::vector<Run3ScoutingParticle>> scoutingparticleHandle;
  iEvent.getByToken(input_scoutingparticle_token_, scoutingparticleHandle);

  Handle<std::vector<Run3ScoutingVertex>> scoutingvertexHandle;
  iEvent.getByToken(input_scoutingvertex_input_, scoutingvertexHandle);

  auto packs = std::make_unique<std::vector<pat::PackedCandidate>>(scoutingparticleHandle->size());
  for (unsigned int icand = 0; icand < scoutingparticleHandle->size(); ++icand) {

      auto& pack = dynamic_cast<pat::PackedCandidate&>((*packs)[icand]);
      auto& scoutingparticle = (*scoutingparticleHandle)[icand];
      int ivertex = scoutingparticle.vertex();

      auto m = pdTable->particle(HepPDT::ParticleID(scoutingparticle.pdgId())) != nullptr
                        ? pdTable->particle(HepPDT::ParticleID(scoutingparticle.pdgId()))->mass()
                        : -99.f;
      if (m < -90) continue;
      
      pat::PackedCandidate::LorentzVector p4(scoutingparticle.pt(), scoutingparticle.eta(), scoutingparticle.phi(), m);
      Run3ScoutingVertex vertex = (*scoutingvertexHandle)[ivertex];
      pat::PackedCandidate::Point v(vertex.x(), vertex.y(), vertex.z());
 
      pack = pat::PackedCandidate(p4, v, scoutingparticle.pt(), scoutingparticle.eta(), scoutingparticle.phi(), scoutingparticle.pdgId(), reco::VertexRefProd(), reco::VertexRef().key());
  }

  //put output
  iEvent.put(std::move(packs));
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Run3ScoutingToPackedCandidateProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("scoutingparticle", edm::InputTag("hltScoutingPFPacker"));
  desc.add<edm::InputTag>("scoutingvertex", edm::InputTag("hltScoutingPrimaryVertexPacker"));
  descriptions.addWithDefaultLabel(desc);
}

// declare this class as a framework plugin
DEFINE_FWK_MODULE(Run3ScoutingToPackedCandidateProducer);
