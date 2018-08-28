
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
        
        void SetMinFitFrac(double frac) {minFitFrac_ = frac;}
        double MinFitFrac() const       {return minFitFrac_;}

        void AddGeantId(int id)   {geant_ids_.insert(id);}
        std::set<int>& GeantIds() {return geant_ids_;}

        // (re)creates histograms from current axisDefs
        Int_t Init();

        // process event
        Int_t Make();

        // save result histograms to disk
        Int_t Finish();

    private:

        CentralityDef* p17id_cent_def_;
        StRefMultCorr* p16id_cent_def_;

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
        
        TH3D* mc_eta_;
        TH3D* mc_phi_;

        TH3D* reco_nhit_;
        TH3D* reco_dca_;
        TH3D* reco_nhitposs_;
        TH3D* reco_eta_;
        TH3D* reco_phi_;

        TH3D* reco_cut_nhit_;
        TH3D* reco_cut_dca_;
        TH3D* reco_cut_nhitposs_;
        TH3D* reco_eta_;
        TH3D* reco_phi_;

        TH3D* data_nhit_;
        TH3D* data_dca_;
        TH3D* data_nhitposs_;
        TH3D* data_eta_;
        TH3D* data_phi_;

        TH1D* vz_;
        TH1D* refmult_;
        TH1D* grefmult_;
        TH1D* centrality_;

        TH3D* mc_reco_tracks_;

        TH2D* mc_tracks_;
        TH2D* reco_tracks_;

        int minFit_;
        double minFitFrac_;
        double maxDCA_;
        std::set<int> geant_ids_;

        ClassDef(StEfficiencyAssessor,1)
};

#endif // STEFFICIENCYASSESSOR__HH
