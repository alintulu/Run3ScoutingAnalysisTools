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
#include "Run3ScoutingAnalysisTools/Analysis/interface/Run3ScoutingAK8PFJet.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/XConePlugin.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/contrib/RecursiveSoftDrop.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"
#include "fastjet/JadePlugin.hh"
#include "fastjet/contrib/SoftKiller.hh"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/Math/interface/libminifloat.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

class Run3ScoutingAK8PFJetProducer : public edm::global::EDProducer<> {
public:
  explicit Run3ScoutingAK8PFJetProducer(const edm::ParameterSet &);
  ~Run3ScoutingAK8PFJetProducer() override;

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
Run3ScoutingAK8PFJetProducer::Run3ScoutingAK8PFJetProducer(const edm::ParameterSet &iConfig)
    : input_scoutingparticle_token_(consumes(iConfig.getParameter<edm::InputTag>("scoutingparticle"))),
      particletable_token_(esConsumes<HepPDT::ParticleDataTable, edm::DefaultRecord>()),
      debug_(iConfig.existsAs<bool>("debug") ? iConfig.getParameter<bool>("debug") : false) {
  //register products
  produces<std::vector<Run3ScoutingAK8PFJet>>();
}

Run3ScoutingAK8PFJetProducer::~Run3ScoutingAK8PFJetProducer() = default;

void Run3ScoutingAK8PFJetProducer::print(Run3ScoutingParticle s, reco::PFCandidate r) const {

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
void Run3ScoutingAK8PFJetProducer::produce(edm::StreamID sid, edm::Event & iEvent, edm::EventSetup const & setup) const {
  using namespace edm;

  auto pdt = setup.getHandle(particletable_token_);
  const HepPDT::ParticleDataTable* pdTable = pdt.product();

  auto particles = iEvent.getHandle(input_scoutingparticle_token_);

  std::vector<fastjet::PseudoJet> fj_part;
  fj_part.reserve(particles->size());
  int idx = 0;
  for (auto part_iter = particles->begin(); part_iter != particles->end(); ++part_iter) {

    auto m = pdTable->particle(HepPDT::ParticleID(part_iter->pdgId())) != nullptr
                        ? pdTable->particle(HepPDT::ParticleID(part_iter->pdgId()))->mass()
                        : -99.f;
    if (m < -90) continue;

    math::PtEtaPhiMLorentzVector p4(part_iter->pt(), part_iter->eta(), part_iter->phi(), m);
    fj_part.emplace_back(p4.px(), p4.py(), p4.pz(), p4.energy());
    fj_part.back().set_user_index(idx);
    idx++;
  }

  fastjet::JetDefinition ak8_def = fastjet::JetDefinition(fastjet::antikt_algorithm, 0.8);
  fastjet::GhostedAreaSpec area_spec(5.0,1,0.01);
  fastjet::AreaDefinition area_def(fastjet::active_area, area_spec);

  double sd_z_cut = 0.10;
  double sd_beta = 0;
  fastjet::contrib::SoftDrop sd_groomer = fastjet::contrib::SoftDrop(sd_beta, sd_z_cut, 0.8);

  fastjet::ClusterSequenceArea ak8_cs(fj_part, ak8_def, area_def);
  std::vector<fastjet::PseudoJet> ak8_jets = fastjet::sorted_by_pt(ak8_cs.inclusive_jets(170.0));

  auto jets = std::make_unique<std::vector<Run3ScoutingAK8PFJet>>(ak8_jets.size());
  for (unsigned int i_jet = 0; i_jet < ak8_jets.size(); ++i_jet) {

      auto& jet = dynamic_cast<Run3ScoutingAK8PFJet&>((*jets)[i_jet]);
      auto& ak8_jet = ak8_jets[i_jet];
      fastjet::PseudoJet sd_ak8 = sd_groomer(ak8_jet);

      std::vector<fastjet::PseudoJet> constituents = ak8_jet.constituents();
      std::vector<int> constituents_idx(ak8_jets.size());
      for (unsigned int i_const = 0; i_const < constituents.size(); ++i_const) {
         constituents_idx[i_const] = constituents[i_const].user_index();
      }

      jet = Run3ScoutingAK8PFJet(ak8_jet.pt(), ak8_jet.eta(), ak8_jet.phi(), ak8_jet.m(), sd_ak8.m(), ak8_jet.area(), constituents_idx);

  }

  //put output
  iEvent.put(std::move(jets));
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Run3ScoutingAK8PFJetProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("scoutingparticle", edm::InputTag("hltScoutingPFPacker"));
  descriptions.addWithDefaultLabel(desc);
}

// declare this class as a framework plugin
DEFINE_FWK_MODULE(Run3ScoutingAK8PFJetProducer);
