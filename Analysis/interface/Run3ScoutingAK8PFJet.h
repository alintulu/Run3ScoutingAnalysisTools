#ifndef DataFormats_Run3ScoutingAK8PFJet_h
#define DataFormats_Run3ScoutingAK8PFJet_h

#include <vector>

//class for holding PF jet information, for use in data scouting
//IMPORTANT: the content of this class should be changed only in backwards compatible ways!
class Run3ScoutingAK8PFJet {
public:
  //constructor with values for all data fields
  Run3ScoutingAK8PFJet(float pt,
                    float eta,
                    float phi,
                    float m,
                    float msoftdrop,
                    float jetArea,
                    std::vector<int> constituents)
      : pt_(pt),
        eta_(eta),
        phi_(phi),
        m_(m),
        msoftdrop_(m),
        jetArea_(jetArea),
        constituents_(std::move(constituents)) {}

  //default constructor
  Run3ScoutingAK8PFJet()
      : pt_(0),
        eta_(0),
        phi_(0),
        m_(0),
        msoftdrop_(0),
        jetArea_(0) {}

  //accessor functions
  float pt() const { return pt_; }
  float eta() const { return eta_; }
  float phi() const { return phi_; }
  float m() const { return m_; }
  float msoftdrop() const { return msoftdrop_; }
  float jetArea() const { return jetArea_; }
  std::vector<int> const& constituents() const { return constituents_; }

private:
  float pt_;
  float eta_;
  float phi_;
  float m_;
  float msoftdrop_;
  float jetArea_;
  std::vector<int> constituents_;
};

typedef std::vector<Run3ScoutingAK8PFJet> Run3ScoutingAK8PFJetCollection;

#endif
