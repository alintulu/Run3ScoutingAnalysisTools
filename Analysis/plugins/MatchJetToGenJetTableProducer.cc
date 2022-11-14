#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include <string> 

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Math/interface/deltaR.h"

class MatchJetToGenJetTableProducer : public edm::stream::EDProducer<> {
public:
  explicit MatchJetToGenJetTableProducer(const edm::ParameterSet &);
  ~MatchJetToGenJetTableProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void produce(edm::Event &, const edm::EventSetup &) override;
    
  const std::string nameTable_;
  const edm::EDGetTokenT<edm::View<reco::Jet>> input_scoutingjet_token_;
  const edm::EDGetTokenT<std::vector<reco::GenJet>> input_genjet_token_;
  
};

//
// constructors and destructor
//
MatchJetToGenJetTableProducer::MatchJetToGenJetTableProducer(const edm::ParameterSet &iConfig)
    : nameTable_(iConfig.getParameter<std::string>("nameTable")),
      input_scoutingjet_token_(consumes(iConfig.getParameter<edm::InputTag>("jets"))),
      input_genjet_token_(consumes(iConfig.getParameter<edm::InputTag>("genjets"))) {
  produces<nanoaod::FlatTable>(nameTable_);
}

MatchJetToGenJetTableProducer::~MatchJetToGenJetTableProducer() {}

void MatchJetToGenJetTableProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  auto jets = iEvent.getHandle(input_scoutingjet_token_);
  auto genjets = iEvent.getHandle(input_genjet_token_);

  // match jet to genjet  
  std::map<int, int> resultMap; // jet, genjet
  std::vector<int> unmatchedJets;
  std::vector<std::tuple<int, int, float> > pairList;  // jet, genjet, dR

  for(unsigned int i = 0; i < jets->size(); i++) {
    bool found_match = false;
    for(unsigned int j = 0; j < genjets->size(); j++) {
      float dR = reco::deltaR((*jets)[i].eta(), (*jets)[i].phi(), (*genjets)[j].eta(), (*genjets)[j].phi());
    if(dR < 0.4) {
        pairList.push_back(std::make_tuple(i, j, dR));
        found_match = true;
      }
    }
    if (!found_match) unmatchedJets.push_back(i);
  }

  std::sort(pairList.begin(), pairList.end(), [](std::tuple<int, int, float> t1, std::tuple<int, int, float> t2){ return std::get<2>(t1) < std::get<2>(t2); });

  while(pairList.size() > 0) {
    resultMap[std::get<0>(pairList[0])] = std::get<1>(pairList[0]);

    for(unsigned int k = 1; k < pairList.size(); k++) {
      if(std::get<0>(pairList[k]) == std::get<0>(pairList[0]) ||
         std::get<1>(pairList[k]) == std::get<1>(pairList[0])) {
        pairList.erase(pairList.begin() + k);
      }
    }
    pairList.erase(pairList.begin());
  } // matching finished; results stored in resultMap

  std::vector<int> jet_genJetIdx(jets->size());  

  for (unsigned i_jet = 0; i_jet < jets->size(); ++i_jet) {
      // new version of the jet loop which reads tag info instead of constituent info

      if (std::find(unmatchedJets.begin(), unmatchedJets.end(), i_jet) != unmatchedJets.end()) {             
          jet_genJetIdx[i_jet] = -1;
      } else {
          jet_genJetIdx[i_jet] = resultMap[i_jet];
      }
  }
    
  // DeepJetInputs table
  auto matchTable = std::make_unique<nanoaod::FlatTable>(jet_genJetIdx.size(), nameTable_, false, true);
    
  matchTable->addColumn<int>("genJetIdx",
                          jet_genJetIdx,
                          "index of matched gen jet");
    
  iEvent.put(std::move(matchTable), nameTable_);  
  
}

void MatchJetToGenJetTableProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("genjets", edm::InputTag("slimmedGenJets"));
  desc.add<edm::InputTag>("jets", edm::InputTag("ak4Jets"));
  desc.add<std::string>("nameTable", "ScoutingJet");
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(MatchJetToGenJetTableProducer);
