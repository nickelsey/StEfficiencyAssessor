
#ifndef STEFFICIENCYASSESSOR__HH
#define STEFFICIENCYASSESSOR__HH

#include "centrality_def.hh"

#include <string>
#include <vector>
#include <set>

#include "StMaker.h"
#include "StMiniMcEvent/StMiniMcEvent.h"
#include "TChain.h"
#include "TFile.h"
#include "TH3F.h"
#include "TH2D.h"

#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"

struct axisDef {
  unsigned nBins;
  double low;
  double high;
  
  axisDef() : nBins(1), low(0), high(1) {}
  
  axisDef(unsigned n, double l, double h)
  : nBins(n), low(l), high(h) {}
  
  axisDef(const axisDef& rhs)
  : nBins(rhs.nBins), low(rhs.low), high(rhs.high) {};
  
  double width() const {return (high - low) / nBins;}
  
  bool valid() const {return nBins > 0 && width() > 0.0;}
  
  int bin(double val) const {
    for (unsigned i = 0; i < nBins; ++i) {
      if (val > low + i * width() &&
          val <= low + (i + 1) * width()) {
        return i;
      }
    }
    return -1;
  }
};

class StEfficiencyAssessor : public StMaker {
public:
  StEfficiencyAssessor(TChain* chain, std::string outputFile = "StEfficiencyAssessor.root");
  
  ~StEfficiencyAssessor();
  
  // loads a new chain
  bool LoadTree(TChain* chain);
  
  // set axis bounds
  void SetDefaultAxes();
  void SetLuminosityAxis(unsigned n, double low, double high);
  void SetCentralityAxis(unsigned n, double low, double high);
  void SetVzAxis(unsigned n, double low, double high);
  void SetPtAxis(unsigned n, double low, double high);
  void SetEtaAxis(unsigned n, double low, double high);
  void SetPhiAxis(unsigned n, double low, double high);
  
  // allows you to modify the centrality and StRefMultCorr definitions
  CentralityDef& CentralityDefinition() {return cent_def_;}
  
  // set track cuts for matched tracks
  void SetDCAMax(double dca) {maxDCA_ = dca;}
  double DCAMax() const      {return maxDCA_;}
  
  void SetMinFitPoints(unsigned fit) {minFit_ = fit;}
  unsigned MinFitPoints() const      {return minFit_;}
  
  void AddGeantId(int id)   {geant_ids_.insert(id);}
  std::set<int>& GeantIds() {return geant_ids_;}
  
  // (re)creates histograms from current axisDefs
  Int_t Init();
  
  // process event
  Int_t Make();
  
  // save result histograms to disk
  Int_t Finish();
  
private:
  
  CentralityDef cent_def_;
  
  int InitInput();
  int InitOutput();
  bool LoadEvent();
  
  bool CheckAxes();
  
  TChain* chain_;
  TFile* out_;
  
  unsigned current_;
  
  StMuDstMaker* muDstMaker_;
  StMuDst* muDst_;
  StMuEvent* muInputEvent_;
  
  StMiniMcEvent* event_;
  
  axisDef lumi_axis_;
  axisDef cent_axis_;
  axisDef vz_axis_;
  axisDef pt_axis_;
  axisDef eta_axis_;
  axisDef phi_axis_;
  
  std::vector<std::vector<std::vector<std::vector<TH3F*>>>> matched_dca_;
  std::vector<std::vector<std::vector<std::vector<TH3F*>>>> data_dca_;
  std::vector<std::vector<std::vector<std::vector<TH3F*>>>> data_pion_dca_;
  
  std::vector<std::vector<std::vector<std::vector<TH3F*>>>> matched_nhit_;
  std::vector<std::vector<std::vector<std::vector<TH3F*>>>> data_nhit_;
  std::vector<std::vector<std::vector<std::vector<TH3F*>>>> data_pion_nhit_;
  
  std::vector<std::vector<std::vector<std::vector<TH3F*>>>> matched_nhit_pos_;
  std::vector<std::vector<std::vector<std::vector<TH3F*>>>> data_nhit_pos_;
  std::vector<std::vector<std::vector<std::vector<TH3F*>>>> data_pion_nhit_pos_;
  
  std::vector<int> runids;
  std::vector<int> eventids;
  
  TH1D* nPrimaries_;
  TH1D* refMult_;
  
  TH1D* nMC_;
  TH1D* nMCPrim_;
  TH1D* nMCRefMult_;
  TH1D* nMCNoParent_;
  TH1D* nMatched_;
  
  TH2D* etaMatch;
  TH2D* phiMatch;
  
  TH2D* etaPrim;
  TH2D* phiPrim;
  
  TH2D* etaPrimPion;
  TH2D* phiPrimPion;
  
  TH3F* weighted_dca_match;
  TH3F* weighted_dca_data;
  TH3F* weighted_dca_data_pion;
  
  int minFit_;
  double maxDCA_;
  std::set<int> geant_ids_;
  
  ClassDef(StEfficiencyAssessor,1)
};

#endif // STEFFICIENCYASSESSOR__HH
