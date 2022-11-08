#include <TVector3.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "RecoBTag/FeatureTools/interface/TrackInfoBuilder.h"
#include "RecoBTag/FeatureTools/interface/deep_helpers.h"
#include "RecoBTag/FeatureTools/interface/sorting_modules.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/BTauReco/interface/DeepBoostedJetTagInfo.h"

using namespace btagbtvdeep;

class DeepBoostedJetTagInfoScoutingProducer : public edm::stream::EDProducer<> {
public:
  explicit DeepBoostedJetTagInfoScoutingProducer(const edm::ParameterSet &);
  ~DeepBoostedJetTagInfoScoutingProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  typedef std::vector<reco::DeepBoostedJetTagInfo> DeepBoostedJetTagInfoCollection;
  typedef reco::VertexCompositePtrCandidateCollection SVCollection;
  typedef reco::VertexCollection VertexCollection;
  typedef edm::View<reco::Candidate> CandidateView;

  void beginStream(edm::StreamID) override {}
  void produce(edm::Event &, const edm::EventSetup &) override;
  void endStream() override {}

  void fillParticleFeatures(DeepBoostedJetFeatures &fts, const reco::Jet &jet);

  const double jet_radius_;
  const double min_jet_pt_;
  const double max_jet_eta_;
  const double min_pt_for_track_properties_;
  const double min_pt_for_pfcandidates_;
  const bool use_puppiP4_;
  const double min_puppi_wgt_;
  const bool include_neutrals_;
  const bool sort_by_sip2dsig_;
  const bool flip_ip_sign_;
  const double max_sip3dsig_;
  const bool use_hlt_features_;

  edm::EDGetTokenT<edm::View<reco::Jet>> jet_token_;
  edm::EDGetTokenT<CandidateView> pfcand_token_;

  bool use_puppi_value_map_;
  bool use_pvasq_value_map_;
  bool is_packed_pf_candidate_collection_;

  edm::EDGetTokenT<edm::ValueMap<float>> puppi_value_map_token_;
  edm::EDGetTokenT<edm::ValueMap<int>> pvasq_value_map_token_;
  edm::EDGetTokenT<edm::Association<VertexCollection>> pvas_token_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> track_builder_token_;

  edm::Handle<VertexCollection> vtxs_;
  edm::Handle<SVCollection> svs_;
  edm::Handle<CandidateView> pfcands_;
  edm::ESHandle<TransientTrackBuilder> track_builder_;
  edm::Handle<edm::ValueMap<float>> puppi_value_map_;
  edm::Handle<edm::ValueMap<int>> pvasq_value_map_;
  edm::Handle<edm::Association<VertexCollection>> pvas_;

  const static std::vector<std::string> particle_features_;
  const static std::vector<std::string> sv_features_;
  const static std::vector<std::string> particle_features_hlt_;
  const static std::vector<std::string> sv_features_hlt_;
  const reco::Vertex *pv_ = nullptr;
  const static float min_track_pt_property_;
  const static int min_valid_pixel_hits_;

  std::map<reco::CandidatePtr::key_type, float> puppi_wgt_cache;

  edm::EDGetTokenT<edm::ValueMap<float>> normchi2_value_map_token_;
  edm::EDGetTokenT<edm::ValueMap<float>> dz_value_map_token_;
  edm::EDGetTokenT<edm::ValueMap<float>> dxy_value_map_token_;
  edm::EDGetTokenT<edm::ValueMap<float>> dzsig_value_map_token_;
  edm::EDGetTokenT<edm::ValueMap<float>> dxysig_value_map_token_;
  edm::EDGetTokenT<edm::ValueMap<int>> lostInnerHits_value_map_token_;
  edm::EDGetTokenT<edm::ValueMap<int>> quality_value_map_token_;
  edm::EDGetTokenT<edm::ValueMap<float>> trkPt_value_map_token_;
  edm::EDGetTokenT<edm::ValueMap<float>> trkEta_value_map_token_;
  edm::EDGetTokenT<edm::ValueMap<float>> trkPhi_value_map_token_;

  edm::Handle<edm::ValueMap<float>> normchi2_value_map_;
  edm::Handle<edm::ValueMap<float>> dz_value_map_;
  edm::Handle<edm::ValueMap<float>> dxy_value_map_;
  edm::Handle<edm::ValueMap<float>> dzsig_value_map_;
  edm::Handle<edm::ValueMap<float>> dxysig_value_map_;
  edm::Handle<edm::ValueMap<int>> lostInnerHits_value_map_;
  edm::Handle<edm::ValueMap<int>> quality_value_map_;
  edm::Handle<edm::ValueMap<float>> trkPt_value_map_;
  edm::Handle<edm::ValueMap<float>> trkEta_value_map_;
  edm::Handle<edm::ValueMap<float>> trkPhi_value_map_;
};

const std::vector<std::string> DeepBoostedJetTagInfoScoutingProducer::particle_features_{
    "pfcand_quality",       "pfcand_charge",   "pfcand_isEl",           "pfcand_isMu",
    "pfcand_isChargedHad",  "pfcand_isGamma",  "pfcand_isNeutralHad",   "pfcand_phirel",
    "pfcand_etarel",        "pfcand_deltaR",   "pfcand_abseta",         "pfcand_ptrel_log",
    "pfcand_erel_log",      "pfcand_pt_log",   "pfcand_normchi2",       "pfcand_dz", 
    "pfcand_dxy",           "pfcand_dxysig",   "pfcand_btagEtaRel",     "pfcand_btagPtRatio",
    "pfcand_btagPParRatio", "pfcand_mask",     "pfcand_pt_log_nopuppi", "pfcand_dzsig",
    "pfcand_e_log_nopuppi", "pfcand_ptrel",    "pfcand_erel",           "pfcand_lostInnerHits"};

const std::vector<std::string> DeepBoostedJetTagInfoScoutingProducer::particle_features_hlt_{"jet_pfcand_pt_log",
                                                                                     "jet_pfcand_energy_log",
                                                                                     "jet_pfcand_deta",
                                                                                     "jet_pfcand_dphi",
                                                                                     "jet_pfcand_eta",
                                                                                     "jet_pfcand_charge",
                                                                                     "jet_pfcand_frompv",
                                                                                     "jet_pfcand_nlostinnerhits",
                                                                                     "jet_pfcand_track_chi2",
                                                                                     "jet_pfcand_track_qual",
                                                                                     "jet_pfcand_dz",
                                                                                     "jet_pfcand_dzsig",
                                                                                     "jet_pfcand_dxy",
                                                                                     "jet_pfcand_dxysig",
                                                                                     "jet_pfcand_etarel",
                                                                                     "jet_pfcand_pperp_ratio",
                                                                                     "jet_pfcand_ppara_ratio",
                                                                                     "jet_pfcand_trackjet_d3d",
                                                                                     "jet_pfcand_trackjet_d3dsig",
                                                                                     "jet_pfcand_trackjet_dist",
                                                                                     "jet_pfcand_nhits",
                                                                                     "jet_pfcand_npixhits",
                                                                                     "jet_pfcand_nstriphits",
                                                                                     "jet_pfcand_trackjet_decayL",
                                                                                     "jet_pfcand_puppiw",
                                                                                     "pfcand_mask"};

const std::vector<std::string> DeepBoostedJetTagInfoScoutingProducer::sv_features_{"sv_mask",
                                                                           "sv_ptrel",
                                                                           "sv_erel",
                                                                           "sv_phirel",
                                                                           "sv_etarel",
                                                                           "sv_deltaR",
                                                                           "sv_abseta",
                                                                           "sv_mass",
                                                                           "sv_ptrel_log",
                                                                           "sv_erel_log",
                                                                           "sv_pt_log",
                                                                           "sv_pt",
                                                                           "sv_ntracks",
                                                                           "sv_normchi2",
                                                                           "sv_dxy",
                                                                           "sv_dxysig",
                                                                           "sv_d3d",
                                                                           "sv_d3dsig",
                                                                           "sv_costhetasvpv"};

const std::vector<std::string> DeepBoostedJetTagInfoScoutingProducer::sv_features_hlt_{"jet_sv_pt_log",
                                                                               "jet_sv_mass",
                                                                               "jet_sv_deta",
                                                                               "jet_sv_dphi",
                                                                               "jet_sv_eta",
                                                                               "jet_sv_ntrack",
                                                                               "jet_sv_chi2",
                                                                               "jet_sv_dxy",
                                                                               "jet_sv_dxysig",
                                                                               "jet_sv_d3d",
                                                                               "jet_sv_d3dsig",
                                                                               "sv_mask"};

const float DeepBoostedJetTagInfoScoutingProducer::min_track_pt_property_ = 0.5;
const int DeepBoostedJetTagInfoScoutingProducer::min_valid_pixel_hits_ = 0;

DeepBoostedJetTagInfoScoutingProducer::DeepBoostedJetTagInfoScoutingProducer(const edm::ParameterSet &iConfig)
    : jet_radius_(iConfig.getParameter<double>("jet_radius")),
      min_jet_pt_(iConfig.getParameter<double>("min_jet_pt")),
      max_jet_eta_(iConfig.getParameter<double>("max_jet_eta")),
      min_pt_for_track_properties_(iConfig.getParameter<double>("min_pt_for_track_properties")),
      min_pt_for_pfcandidates_(iConfig.getParameter<double>("min_pt_for_pfcandidates")),
      use_puppiP4_(iConfig.getParameter<bool>("use_puppiP4")),
      min_puppi_wgt_(iConfig.getParameter<double>("min_puppi_wgt")),
      include_neutrals_(iConfig.getParameter<bool>("include_neutrals")),
      sort_by_sip2dsig_(iConfig.getParameter<bool>("sort_by_sip2dsig")),
      flip_ip_sign_(iConfig.getParameter<bool>("flip_ip_sign")),
      max_sip3dsig_(iConfig.getParameter<double>("sip3dSigMax")),
      use_hlt_features_(iConfig.getParameter<bool>("use_hlt_features")),
      jet_token_(consumes<edm::View<reco::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
      pfcand_token_(consumes<CandidateView>(iConfig.getParameter<edm::InputTag>("pf_candidates"))),
      use_puppi_value_map_(false),
      use_pvasq_value_map_(false),
      track_builder_token_(
          esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))) {
  const auto &puppi_value_map_tag = iConfig.getParameter<edm::InputTag>("puppi_value_map");
  if (!puppi_value_map_tag.label().empty()) {
    puppi_value_map_token_ = consumes<edm::ValueMap<float>>(puppi_value_map_tag);
    use_puppi_value_map_ = true;
  }

  normchi2_value_map_token_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("normchi2_value_map"));
  dz_value_map_token_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("dz_value_map"));
  dxy_value_map_token_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("dxy_value_map"));
  dzsig_value_map_token_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("dzsig_value_map"));
  dxysig_value_map_token_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("dxysig_value_map"));
  lostInnerHits_value_map_token_ = consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("lostInnerHits_value_map"));
  quality_value_map_token_ = consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("quality_value_map"));
  trkPt_value_map_token_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trkPt_value_map"));
  trkEta_value_map_token_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trkEta_value_map"));
  trkPhi_value_map_token_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trkPhi_value_map"));

  produces<DeepBoostedJetTagInfoCollection>();
}

DeepBoostedJetTagInfoScoutingProducer::~DeepBoostedJetTagInfoScoutingProducer() {}

void DeepBoostedJetTagInfoScoutingProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  // pfDeepBoostedJetTagInfos
  edm::ParameterSetDescription desc;
  desc.add<double>("jet_radius", 0.8);
  desc.add<double>("min_jet_pt", 150);
  desc.add<double>("max_jet_eta", 99);
  desc.add<double>("min_pt_for_track_properties", -1);
  desc.add<double>("min_pt_for_pfcandidates", -1);
  desc.add<bool>("use_puppiP4", true);
  desc.add<bool>("include_neutrals", true);
  desc.add<bool>("sort_by_sip2dsig", false);
  desc.add<double>("min_puppi_wgt", 0.01);
  desc.add<bool>("flip_ip_sign", false);
  desc.add<double>("sip3dSigMax", -1);
  desc.add<bool>("use_hlt_features", false);
  desc.add<edm::InputTag>("pf_candidates", edm::InputTag("particleFlow"));
  desc.add<edm::InputTag>("jets", edm::InputTag("ak8PFJetsPuppi"));
  desc.add<edm::InputTag>("puppi_value_map", edm::InputTag("puppi"));
  descriptions.add("pfDeepBoostedJetTagInfos", desc);
}

void DeepBoostedJetTagInfoScoutingProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  // output collection
  auto output_tag_infos = std::make_unique<DeepBoostedJetTagInfoCollection>();
  // Input jets
  auto jets = iEvent.getHandle(jet_token_);
  // Input pfcands
  iEvent.getByToken(pfcand_token_, pfcands_);

  for (std::size_t jet_n = 0; jet_n < jets->size(); jet_n++) {
    const auto &jet = (*jets)[jet_n];
    edm::RefToBase<reco::Jet> jet_ref(jets, jet_n);

    // create jet features
    DeepBoostedJetFeatures features;
    for (const auto &name : particle_features_) {
      features.add(name);
    }

    // fill values only if above pt threshold and has daughters, otherwise left
    bool fill_vars = true;
    if (jet.pt() < min_jet_pt_ or std::abs(jet.eta()) > max_jet_eta_)
      fill_vars = false;
    if (jet.numberOfDaughters() == 0)
      fill_vars = false;

    // fill features
    if (fill_vars) {
      fillParticleFeatures(features, jet);
      features.check_consistency(particle_features_);
    }
    // this should always be done even if features are not filled
    output_tag_infos->emplace_back(features, jet_ref);
  }
  // move output collection
  iEvent.put(std::move(output_tag_infos));
}

void DeepBoostedJetTagInfoScoutingProducer::fillParticleFeatures(DeepBoostedJetFeatures &fts, const reco::Jet &jet) {
  // some jet properties
  math::XYZVector jet_dir = jet.momentum().Unit();
  TVector3 jet_direction(jet.momentum().Unit().x(), jet.momentum().Unit().y(), jet.momentum().Unit().z());
  //GlobalVector jet_ref_track_dir(jet.px(), jet.py(), jet.pz());
  const float etasign = jet.eta() > 0 ? 1 : -1;

  // make list of pf-candidates to be considered
  std::vector<reco::CandidatePtr> daughters;
  for (const auto &dau : jet.daughterPtrVector()) {
    auto cand = pfcands_->ptrAt(dau.key());
    daughters.push_back(cand);
  }

  // sort daughters according to pt
  std::sort(daughters.begin(), daughters.end(), [](const auto &a, const auto &b) { return a->pt() > b->pt(); });
  
  for (const auto &name : particle_features_)
    fts.reserve(name, daughters.size());

  // Build observables
  size_t icand = 0;
  for (const auto &cand : daughters) {
    const auto *packed_cand = dynamic_cast<const pat::PackedCandidate *>(&(*cand));
    const auto *reco_cand = dynamic_cast<const reco::PFCandidate *>(&(*cand));

    if (not packed_cand and not reco_cand)
      throw edm::Exception(edm::errors::InvalidReference)
          << "Cannot convert to either reco::PFCandidate or pat::PackedCandidate";

    const float ip_sign = flip_ip_sign_ ? -1 : 1;

    // input particle is a packed PF candidate
    auto candP4 = cand->p4();

    // Building offline features
    fts.fill("pfcand_lostInnerHits", (*lostInnerHits_value_map_)[cand]);
    fts.fill("pfcand_quality", (*quality_value_map_)[cand]);
    fts.fill("pfcand_charge", reco_cand->charge());
    fts.fill("pfcand_isEl", std::abs(reco_cand->pdgId()) == 11);
    fts.fill("pfcand_isMu", std::abs(reco_cand->pdgId()) == 13);
    fts.fill("pfcand_isChargedHad", std::abs(reco_cand->pdgId()) == 211);
    fts.fill("pfcand_isGamma", std::abs(reco_cand->pdgId()) == 22);
    fts.fill("pfcand_isNeutralHad", std::abs(reco_cand->pdgId()) == 130);
    fts.fill("pfcand_dz", (*dz_value_map_)[cand]);
    fts.fill("pfcand_dzsig", (*dzsig_value_map_)[cand]);
    fts.fill("pfcand_dxy", (*dxy_value_map_)[cand]);
    fts.fill("pfcand_dxysig", (*dxysig_value_map_)[cand]);

    // generic candidate observables
    fts.fill("pfcand_phirel", reco::deltaPhi(candP4, jet));
    fts.fill("pfcand_etarel", etasign * (candP4.eta() - jet.eta()));
    fts.fill("pfcand_deltaR", reco::deltaR(candP4, jet));
    fts.fill("pfcand_abseta", std::abs(candP4.eta()));
    fts.fill("pfcand_ptrel_log", std::log(candP4.pt() / jet.pt()));
    fts.fill("pfcand_ptrel", candP4.pt() / jet.pt());
    fts.fill("pfcand_erel_log", std::log(candP4.energy() / jet.energy()));
    fts.fill("pfcand_erel", candP4.energy() / jet.energy());
    fts.fill("pfcand_pt_log", std::log(candP4.pt()));
    fts.fill("pfcand_mask", 1);
    fts.fill("pfcand_pt_log_nopuppi", std::log(cand->pt()));
    fts.fill("pfcand_e_log_nopuppi", std::log(cand->energy()));
    fts.fill("pfcand_normchi2", (*normchi2_value_map_)[cand]);
 
    float trk_px = (*trkPt_value_map_)[cand] * std::cos((*trkPhi_value_map_)[cand]);
    float trk_py = (*trkPt_value_map_)[cand] * std::sin((*trkPhi_value_map_)[cand]);
    float trk_pz = (*trkPt_value_map_)[cand] * std::sinh((*trkEta_value_map_)[cand]);
    math::XYZVector track_mom(trk_px, trk_py, trk_pz);
    TVector3 track_direction(trk_px, trk_py, trk_pz);
    double track_mag = sqrt(trk_px * trk_px + trk_py * trk_py + trk_pz * trk_pz);
    fts.fill("pfcand_btagEtaRel", reco::btau::etaRel(jet_dir, track_mom));
    fts.fill("pfcand_btagPtRatio", track_direction.Perp(jet_direction) / track_mag);
    fts.fill("pfcand_btagPParRatio", jet_dir.Dot(track_mom) / track_mag);

    icand++;
  }
}

// define this as a plug-in
DEFINE_FWK_MODULE(DeepBoostedJetTagInfoScoutingProducer);
