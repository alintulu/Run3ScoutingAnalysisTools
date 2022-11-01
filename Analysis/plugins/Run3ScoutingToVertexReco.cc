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
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/Math/interface/libminifloat.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

class Run3ScoutingToVertexReco : public edm::global::EDProducer<> {
public:
  explicit Run3ScoutingToVertexReco(const edm::ParameterSet &);
  ~Run3ScoutingToVertexReco() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
  void produce(edm::StreamID sid, edm::Event & iEvent, edm::EventSetup const & setup) const final;

private:

  const edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> input_scoutingvertex_token_;
};

//
// constructors and destructor
//
Run3ScoutingToVertexReco::Run3ScoutingToVertexReco(const edm::ParameterSet &iConfig)
    : input_scoutingvertex_token_(consumes(iConfig.getParameter<edm::InputTag>("scoutingvertex"))) {
  //register products
  produces<std::vector<reco::Vertex>>();
}

Run3ScoutingToVertexReco::~Run3ScoutingToVertexReco() = default;

// ------------ method called to produce the data  ------------
void Run3ScoutingToVertexReco::produce(edm::StreamID sid, edm::Event & iEvent, edm::EventSetup const & setup) const {
  using namespace edm;

  Handle<std::vector<Run3ScoutingVertex>> scoutingvertexHandle;
  iEvent.getByToken(input_scoutingvertex_token_, scoutingvertexHandle);

  auto vertexs = std::make_unique<std::vector<reco::Vertex>>(scoutingvertexHandle->size());
  for (unsigned int icand = 0; icand < scoutingvertexHandle->size(); ++icand) {

      auto& vertex = dynamic_cast<reco::Vertex&>((*vertexs)[icand]);
      auto& scoutingvertex = (*scoutingvertexHandle)[icand];
      double x = scoutingvertex.x(), y = scoutingvertex.y(), z = scoutingvertex.z();
      reco::Vertex::Error err;
      err(0, 0) = scoutingvertex.xError();
      err(1, 1) = scoutingvertex.yError();
      err(2, 2) = scoutingvertex.zError();
      double chi2 = scoutingvertex.chi2();
      double ndof = scoutingvertex.ndof();
      size_t size = scoutingvertex.tracksSize(); 

      vertex = reco::Vertex(reco::Vertex::Point(x, y, z), err, chi2, ndof, size);
  }

  //put output
  iEvent.put(std::move(vertexs));
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Run3ScoutingToVertexReco::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("scoutingvertex", edm::InputTag("hltScoutingPrimaryVertexPacker"));
  descriptions.addWithDefaultLabel(desc);
}

// declare this class as a framework plugin
DEFINE_FWK_MODULE(Run3ScoutingToVertexReco);
