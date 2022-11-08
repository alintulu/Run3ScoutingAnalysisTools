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
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/OrphanHandle.h"

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
  produces<edm::ValueMap<float>>("normchi2");
}

Run3ScoutingToPackedCandidateProducer::~Run3ScoutingToPackedCandidateProducer() = default;

// ------------ method called to produce the data  ------------
void Run3ScoutingToPackedCandidateProducer::produce(edm::StreamID sid, edm::Event & iEvent, edm::EventSetup const & setup) const {
  using namespace edm;

  auto pdt = setup.getHandle(particleTableToken);
  const HepPDT::ParticleDataTable* pdTable = pdt.product();

  auto scoutingparticleHandle = iEvent.getHandle(input_scoutingparticle_token_);
  auto scoutingvertexHandle = iEvent.getHandle(input_scoutingvertex_input_);

  std::vector<float> normchi2(scoutingparticleHandle->size(), -1);
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
 
      pack = pat::PackedCandidate(p4, v, scoutingparticle.trk_pt(), scoutingparticle.trk_eta(), scoutingparticle.trk_phi(), scoutingparticle.pdgId(), reco::VertexRefProd(), reco::VertexRef().key());

      normchi2[icand] = scoutingparticle.normchi2();
  }

  std::vector<pat::PackedCandidate> const *const_packs = &(*packs);

  //put output
  iEvent.put(std::move(packs));
  std::unique_ptr<edm::ValueMap<float>> normchi2V(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_normchi2(*normchi2V);
  ProductID const pidK(1, 3);
  filler_normchi2.insert(OrphanHandle<std::vector<pat::PackedCandidate>>(const_packs, pidK), normchi2.begin(), normchi2.end());
  filler_normchi2.fill();
  iEvent.put(std::move(normchi2V), "normchi2");
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Run3ScoutingToPackedCandidateProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("scoutingparticle", edm::InputTag("hltScoutingPFPacker"));
  desc.add<edm::InputTag>("scoutingvertex", edm::InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx"));
  descriptions.addWithDefaultLabel(desc);
}

// declare this class as a framework plugin
DEFINE_FWK_MODULE(Run3ScoutingToPackedCandidateProducer);
