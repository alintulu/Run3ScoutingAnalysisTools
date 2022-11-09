#include <memory>

// user include files
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/Math/interface/libminifloat.h"

class Run3ScoutingToRecoJetProducer : public edm::stream::EDProducer<> {
//
//  construction/destruction
// 
public:
  explicit Run3ScoutingToRecoJetProducer(const edm::ParameterSet &);
  ~Run3ScoutingToRecoJetProducer() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

//
//  member functions
//
public:
  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;
  void print(Run3ScoutingPFJet s, reco::PFJet r) const;

//
//  member data
//
protected:
  std::vector<edm::Ptr<reco::Candidate>> pfcands_;

private:
  const edm::EDGetTokenT<std::vector<Run3ScoutingPFJet>> input_scoutingjet_token_;
  edm::EDGetTokenT<reco::CandidateView> input_candidateview_token_;
  bool debug_;
};

//
// constructors and destructor
//
Run3ScoutingToRecoJetProducer::Run3ScoutingToRecoJetProducer(const edm::ParameterSet &iConfig)
    : input_scoutingjet_token_(consumes(iConfig.getParameter<edm::InputTag>("scoutingjet"))),
      input_candidateview_token_(consumes(iConfig.getParameter<edm::InputTag>("recopfcand"))),
      debug_(iConfig.existsAs<bool>("debug") ? iConfig.getParameter<bool>("debug") : false) {
  //register products
  produces<reco::PFJetCollection>();
  produces<edm::ValueMap<float>>("jetArea");
  produces<edm::ValueMap<float>>("mass");
}

Run3ScoutingToRecoJetProducer::~Run3ScoutingToRecoJetProducer() = default;

void Run3ScoutingToRecoJetProducer::print(Run3ScoutingPFJet s, reco::PFJet r) const {

  std::cout << "pt: " << s.pt() << " " << r.pt() << std::endl;
  std::cout << "phi: " << s.phi() << " " << r.phi() << std::endl;
  std::cout << "eta: " << s.eta() << " " << r.eta() << std::endl;
  std::cout << "che: " << s.chargedHadronEnergy() << " " << r.chargedHadronEnergy() << std::endl;
  std::cout << "nhe: " << s.neutralHadronEnergy() << " " << r.neutralHadronEnergy() << std::endl;
  std::cout << "mue: " << s.muonEnergy() << " " << r.muonEnergy() << std::endl;
  std::cout << "phoe: " << s.photonEnergy() << " " << r.photonEnergy() << std::endl;

}

// ------------ method called to produce the data  ------------
void Run3ScoutingToRecoJetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  Handle<std::vector<Run3ScoutingPFJet>> scoutingjetHandle;
  iEvent.getByToken(input_scoutingjet_token_, scoutingjetHandle);

  Handle<reco::CandidateView> candidateviewHandle;
  bool isView = iEvent.getByToken(input_candidateview_token_, candidateviewHandle);

  std::vector<float> jetArea(scoutingjetHandle->size());
  std::vector<float> mass(scoutingjetHandle->size());

  auto pfjets = std::make_unique<reco::PFJetCollection>(scoutingjetHandle->size());
  for (unsigned int ijet = 0; ijet < scoutingjetHandle->size(); ++ijet) {

      auto& pfjet = dynamic_cast<reco::PFJet&>((*pfjets)[ijet]);
      auto& scoutingjet = (*scoutingjetHandle)[ijet];
      std::vector<int> constituents = scoutingjet.constituents();

      if (isView) {
         for (size_t iconst = 0; iconst < constituents.size(); ++iconst) {
           pfcands_.push_back(candidateviewHandle->ptrAt(iconst));
         }   
      }

      float px = scoutingjet.pt() * cos(scoutingjet.phi());
      float py = scoutingjet.pt() * sin(scoutingjet.phi());
      float pz = scoutingjet.pt() * sinh(scoutingjet.eta());
      float energy = scoutingjet.chargedHadronEnergy() + scoutingjet.neutralHadronEnergy() + scoutingjet.photonEnergy() + scoutingjet.electronEnergy() + scoutingjet.muonEnergy();

      reco::PFJet::Specific specific;
      //TO DO
      //specific.mJetArea = scoutingjet.jetArea();
      //specific.mMass = scoutingjet.m();
      specific.mChargedHadronEnergy = scoutingjet.chargedHadronEnergy();
      specific.mNeutralHadronEnergy = scoutingjet.neutralHadronEnergy();
      specific.mPhotonEnergy = scoutingjet.photonEnergy();
      specific.mElectronEnergy = scoutingjet.electronEnergy();
      specific.mMuonEnergy = scoutingjet.muonEnergy();
      specific.mHFHadronEnergy = scoutingjet.HFHadronEnergy();
      specific.mHFEMEnergy = scoutingjet.HFEMEnergy();
      specific.mHOEnergy = scoutingjet.HOEnergy();
      specific.mChargedHadronMultiplicity = scoutingjet.chargedHadronMultiplicity();
      specific.mNeutralHadronMultiplicity = scoutingjet.neutralHadronMultiplicity();
      specific.mPhotonMultiplicity = scoutingjet.photonMultiplicity();
      specific.mElectronMultiplicity = scoutingjet.electronMultiplicity();
      specific.mMuonMultiplicity = scoutingjet.muonMultiplicity();
      specific.mHFHadronMultiplicity = scoutingjet.HFHadronMultiplicity();
      specific.mHFEMMultiplicity = scoutingjet.HFEMMultiplicity();

      reco::Particle::LorentzVector p4(px, py, pz, energy);
      reco::Jet::Point vertex(0, 0, 0);
   
      pfjet = reco::PFJet(p4, vertex, specific, pfcands_);

      jetArea[ijet] = scoutingjet.jetArea();
      mass[ijet] = scoutingjet.m();

      if (debug_) print(scoutingjet, pfjet);
  }

  //put output
  edm::OrphanHandle<reco::PFJetCollection> oh = iEvent.put(std::move(pfjets));

  std::unique_ptr<edm::ValueMap<float>> jetArea_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_jetArea(*jetArea_VM);
  filler_jetArea.insert(oh, jetArea.begin(), jetArea.end());
  filler_jetArea.fill();
  iEvent.put(std::move(jetArea_VM), "jetArea");

  std::unique_ptr<edm::ValueMap<float>> mass_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_mass(*mass_VM);
  filler_mass.insert(oh, mass.begin(), mass.end());
  filler_mass.fill();
  iEvent.put(std::move(mass_VM), "mass");
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Run3ScoutingToRecoJetProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("scoutingjet", edm::InputTag("hltScoutingPFPacker"));
  desc.add<edm::InputTag>("recopfcand", edm::InputTag("pfcans"));
  descriptions.addWithDefaultLabel(desc);
}

// declare this class as a framework plugin
DEFINE_FWK_MODULE(Run3ScoutingToRecoJetProducer);
