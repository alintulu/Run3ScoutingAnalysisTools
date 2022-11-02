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

class Run3ScoutingProducer : public edm::stream::EDProducer<> {
public:
  explicit Run3ScoutingProducer(const edm::ParameterSet &);
  ~Run3ScoutingProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void produce(edm::Event &, const edm::EventSetup &) override;
    
  const edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> input_vertex_token_;
  
};

//
// constructors and destructor
//
Run3ScoutingProducer::Run3ScoutingProducer(const edm::ParameterSet &iConfig)
    : input_vertex_token_(consumes(iConfig.getParameter<edm::InputTag>("vertex"))) {
  produces<nanoaod::FlatTable>("Scouting");
}

Run3ScoutingProducer::~Run3ScoutingProducer() {}

void Run3ScoutingProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  
  auto vertex = iEvent.getHandle(input_vertex_token_);
  std::vector<float> x;
  std::vector<float> y;
  std::vector<float> z;
  std::vector<float> xError;
  std::vector<float> yError;
  std::vector<float> zError;
  std::vector<int> tracksSize;
  std::vector<float> chi2;
  std::vector<bool> isValidVtx;

  for (unsigned i = 0; i < vertex->size(); ++i) {
     x.push_back((*vertex)[i].x());
     y.push_back((*vertex)[i].y());
     z.push_back((*vertex)[i].z());
     xError.push_back((*vertex)[i].xError());
     yError.push_back((*vertex)[i].yError());
     zError.push_back((*vertex)[i].zError());
     tracksSize.push_back((*vertex)[i].tracksSize());
     chi2.push_back((*vertex)[i].chi2());
     isValidVtx.push_back((*vertex)[i].isValidVtx());
  }
    
  auto table = std::make_unique<nanoaod::FlatTable>(vertex->size(), "Scouting", false, true);
    
  table->addColumn<float>("ScoutingVertex_x", x, "vertex x");
  table->addColumn<float>("ScoutingVertex_y", y, "vertex y");
  table->addColumn<float>("ScoutingVertex_z", z, "vertex z");
  table->addColumn<float>("ScoutingVertex_x", xError, "vertex x error");
  table->addColumn<float>("ScoutingVertex_y", yError, "vertex y error");
  table->addColumn<float>("ScoutingVertex_z", zError, "vertex z error");
  table->addColumn<float>("ScoutingVertex_tracksSize", tracksSize, "track size");
  table->addColumn<float>("ScoutingVertex_chi2", chi2, "vertex chi2");
  table->addColumn<float>("ScoutingVertex_isValidVertex", isValidVtx, "boolean if vertex is valid");
    
  iEvent.put(std::move(table), "Scouting");
}

void Run3ScoutingProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("vertex", edm::InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx"));
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(Run3ScoutingProducer);
