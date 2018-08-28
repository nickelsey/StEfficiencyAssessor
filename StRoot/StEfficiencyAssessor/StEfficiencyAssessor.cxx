#include "StEfficiencyAssessor.hh"

#include "St_base/StMessMgr.h"
#include "StMiniMcEvent/StMiniMcPair.h"
#include "StMiniMcEvent/StTinyMcTrack.h"
#include "StMiniMcEvent/StContamPair.h"

#include "StMuDSTMaker/COMMON/StMuTrack.h"

#include "TMath.h"
#include "TTree.h"

ClassImp(StEfficiencyAssessor);

StEfficiencyAssessor::StEfficiencyAssessor(TChain* mcTree, std::string outputFile) {
    if (!LoadTree(mcTree)) {
        LOG_ERROR << "load chain failed" << endm;
    }

    SetDefaultAxes();

    current_ = 0;

    muDstMaker_ = nullptr;
    muDst_ = nullptr;
    muInputEvent_ = nullptr;

    minFit_ = 20;
    minFitFrac_ = 0.52;
    maxDCA_ = 3.0;

    out_ = new TFile(outputFile.c_str(), "RECREATE");
}

StEfficiencyAssessor::~StEfficiencyAssessor() {

}

int StEfficiencyAssessor::Init() {
    if (InitInput() != kStOK)
        return kStFatal;
    if (InitOutput() != kStOK)
        return kStFatal;
    return kStOK;
}

void StEfficiencyAssessor::SetDefaultAxes() {
    lumi_axis_ = axisDef(3, 0.0, 1e5);
    cent_axis_ = axisDef(9, -0.5, 8.5);
    vz_axis_ = axisDef(5, -30, 30);
    pt_axis_   = axisDef(10, 0.0, 5.0);
    eta_axis_  = axisDef(5, -1.0, 1.0);
    phi_axis_  = axisDef(6, -TMath::Pi(), TMath::Pi());
}

void StEfficiencyAssessor::SetLuminosityAxis(unsigned n, double low, double high) {
    lumi_axis_ = axisDef(n, low, high);
}

void StEfficiencyAssessor::SetCentralityAxis(unsigned n, double low, double high) {
    cent_axis_ = axisDef(n, low, high);
}

void StEfficiencyAssessor::SetVzAxis(unsigned n, double low, double high) {
    vz_axis_ = axisDef(n, low, high);
}

void StEfficiencyAssessor::SetPtAxis(unsigned n, double low, double high) {
    pt_axis_ = axisDef(n, low, high);
}

void StEfficiencyAssessor::SetEtaAxis(unsigned n, double low, double high) {
    eta_axis_ = axisDef(n, low, high);
}

void StEfficiencyAssessor::SetPhiAxis(unsigned n, double low, double high) {
    phi_axis_ = axisDef(n, low, high);
}


bool StEfficiencyAssessor::LoadTree(TChain* chain) {
    if (chain == nullptr) {
        LOG_INFO << "chain does not exist" << endm;
        return false;
    }
    if (chain->GetBranch("StMiniMcEvent") == nullptr) {
        LOG_ERROR << "chain does not contain StMiniMcEvent branch" << endm;
        return false;
    }

    chain_ = chain;
    event_ = new StMiniMcEvent;

    chain_->SetBranchAddress("StMiniMcEvent", &event_);
    chain_->GetEntry(0);
    return true;
}

bool StEfficiencyAssessor::CheckAxes() {
    return lumi_axis_.valid() && cent_axis_.valid() && pt_axis_.valid()
        && eta_axis_.valid() && phi_axis_.valid();
}

Int_t StEfficiencyAssessor::Make() {

    if (event_ == nullptr) {
        LOG_ERROR << "StMiniMcEvent Branch not loaded properly: exiting run loop" << endm;
        return kStFatal;
    }


    // load the matching miniMC event
    if (LoadEvent() == false) {
        LOG_ERROR << "Could not find miniMC event matching muDST event" << endm;
        return kStErr;
    }

    int centrality = 0;

    if (p17id_cent_def_) {
        p17id_cent_def_->setEvent(muInputEvent_->runId(), muInputEvent_->refMult(), muInputEvent_->runInfo().zdcCoincidenceRate(), event_->vertexZ());
        centrality = p17id_cent_def_->centrality9();
    }
    else if (p16id_cent_def_) {
        p16id_cent_def_->init(muInputEvent_->runId());
        p16id_cent_def_->initEvent(muInputEvent_->grefmult(), muInputEvent_->primaryVertexPosition().z(), muInputEvent_->runInfo().zdcCoincidenceRate());
        centrality = p16id_cent_def_->getCentralityBin9();
    }
    else {
        LOG_ERROR << "error: no centrality definition created" << endm;
        return kStFatal;
    }
    if (centrality < 0 || centrality > 8)
        continue;

    if (fabs(muInputEvent_->primaryVertexPosition().z()) > 30)
        return kStOK;

    vz_->Fill(muInputEvent_->primaryVertexPosition().z());
    refmult_->Fill(muInputEvent_->refMult());
    grefmult_->Fill(muInputEvent_->grefmult());
    centrality_->Fill(centrality);

    TClonesArray* mc_array = event_->tracks(MC);
    TIter next_mc(mc_array);
    StTinyMcTrack* track = nullptr;
    unsigned count_mc = 0;
    while ((track = (StTinyMcTrack*) next_mc())) {

        if (geant_ids_.size() && geant_ids_.find(track->geantId()) == geant_ids_.end())
            continue;

        if (track->parentGeantId() != 0)
            continue;

        count_mc++;
        mc_tracks_->Fill(centrality, track->ptMc());
        mc_eta_->Fill(centrality, track->ptMc(), track->etaMc());
        mc_phi_->Fill(centrality, track->ptMc(), track->phiMc());
    }


    TClonesArray* match_array = event_->tracks(MATCHED);
    TIter next_match(match_array);
    StMiniMcPair* pair = nullptr;
    unsigned count_pair = 0;
    while ((pair = (StMiniMcPair*) next_match())) {

        if (geant_ids_.size() && geant_ids_.find(pair->geantId()) == geant_ids_.end())
            continue;

        if (pair->parentGeantId() != 0)
            continue;
        reco_tracks_->Fill(centrality, pair->ptPr());
        reco_nhit_->Fill(centrality, pair->ptPr(), pair->fitPts());
        reco_dca_->Fill(centrality, pair->ptPr(), pair->dcaGl());
        reco_eta_->Fill(centrality, pair->ptPr(), pair->etaPr());
        reco_phi_->Fill(centrality, pair->ptPr(), pair->phiPr());
        reco_nhitposs_->Fill(centrality, pair->ptPr(), pair->nPossiblePts()+1);
        reco_fitfrac_->Fill(centrality, pair->ptPr(), (double)pair->fitPts()/pair->nPossiblePts());

        if (pair->dcaGl() > maxDCA_ || pair->fitPts() < minFit_)
            continue;

        if (fabs(pair->etaPr()) > 1.0)
            continue;

        if ((double) pair->fitPts() / (pair->nPossiblePts()+1) < minFitFrac_)
            continue;
        count_pair++;
        reco_cut_nhit_->Fill(centrality, pair->ptPr(), pair->fitPts());
        reco_cut_dca_->Fill(centrality, pair->ptPr(), pair->dcaGl());
        reco_cut_eta_->Fill(centrality, pair->ptPr(), pair->etaPr());
        reco_cut_phi_->Fill(centrality, pair->ptPr(), pair->phiPr());
        reco_cut_nhitposs_->Fill(centrality, pair->ptPr(), pair->nPossiblePts()+1);
        reco_cut_fitfrac_->Fill(centrality, pair->ptPr(), (double)pair->fitPts()/pair->nPossiblePts());
    }
    mc_reco_track_->Fill(centrality, count_mc, count_pair);

    for (int i = 0; i < muDst_->primaryTracks()->GetEntries(); ++i) {
        StMuTrack* muTrack = (StMuTrack*) muDst_->primaryTracks(i);
        if (muTrack->flag() < 0)
            continue;
        if (muTrack->dcaGlobal().mag() > maxDCA_)
            continue;
        if (muTrack->nHitsFit() < minFit_)
            continue;
        if ((double) muTrack->nHitsFit() / (muTrack->nHitsPoss(kTpcId) + 1) < minFitFrac_)
            continue;
        if (fabs(muTrack->eta()) > 1.0)
            continue;
        
        data_nhit_->Fill(centrality, muTrack->pt(), pair->fitPts());
        data_dca_->Fill(centrality, muTrack->pt(), pair->dcaGl());
        data_eta_->Fill(centrality, muTrack->pt(), pair->etaPr());
        data_phi_->Fill(centrality, muTrack->pt(), pair->phiPr());
        data_nhitposs_->Fill(centrality, muTrack->pt(), pair->nPossiblePts()+1);
        data_fitfrac_->Fill(centrality, muTrack->pt(), (double)pair->fitPts()/pair->nPossiblePts());

    }


    return kStOK;
}


Int_t StEfficiencyAssessor::Finish() {
    if (out_ == nullptr) {
        out_ = new TFile("stefficiencyassessor.root", "RECREATE");
    }

    out_->cd();



    out_->Close();
    return kStOk;
}


int StEfficiencyAssessor::InitInput() {
    muDstMaker_ = (StMuDstMaker*) GetMakerInheritsFrom("StMuDstMaker");
    if (muDstMaker_ == nullptr) {
        LOG_ERROR << "No muDstMaker found in chain: StEfficiencyAssessor init failed" << endm;
        return kStFatal;
    }

    if (TString(muDstMaker->GetFile()).contains("SL17d")) {
        p17id_cent_def_ = new CentralityDef();
    }
    else if (TString(muDstMaker->GetFile()).contains("SL16d")) {
        p16id_cent_def_ = CentralityMaker::instance()->getgRefMultCorr_P16id();
        p16id_cent_def->setVzForWeight(6, -6.0, 6.0);
        p16id_cent_def->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_vpd30_vpd5_Run14_P16id.txt");
    }
    else {
        LOG_ERROR << "Library could not be discovered: exiting" << endm;
        return kSFatal;
    }   


    return kStOK;
}

int StEfficiencyAssessor::InitOutput() {
    if (!CheckAxes()) {
        LOG_ERROR << "axes not valid: could not initialize histograms";
        return kStFatal;
    }

    mc_eta_ = new TH3D("mceta", ";cent;pt;#eta", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, -1, 1);
    mc_phi_ = new TH3D("mcphi", ";cent;pt;#phi", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, -1, 1);

    reco_nhit_ = new TH3D("reconhit", ";cent;pt;nhit", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, 0, 50);
    reco_dca_ = new TH3D("recodca", ";cent;pt;DCA[cm]", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, 0, 3.0);
    reco_nhitposs_ = new TH3D("reconhitposs", ";cent;pt;nhitposs", cent_axis_.nBins, cent_axis_.low, cent_axis_.high pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, 0, 50);
    reco_eta_ = new TH3D("recoeta", ";cent;pt;#eta", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, -1, 1);
    reco_phi_ = new TH3D("recophi", ";cent;pt;#phi", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, -TMath::Pi(), TMath::Pi());

    reco_cut_nhit_ = new TH3D("reconhitcut", ";cent;pt;nhit", cent_axis_.nBins, cent_axis_.low, cent_axis_.high pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, 0, 50);
    reco_cut_dca_ = new TH3D("recodcacut", ";cent;pt;DCA[cm]", cent_axis_.nBins, cent_axis_.low, cent_axis_.high pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, 0, 3.0);
    reco_cut_nhitposs_ = new TH3D("reconhitposscut", ";cent;pt;nhitposs", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, 0, 50);
    reco_cut_eta_ = new TH3D("recoetacut", ";cent;pt;#eta", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, -1, 1);
    reco_cut_phi_ = new TH3D("recophicut", ";cent;pt;#phi", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, -TMath::Pi(), TMath::Pi());

    data_nhit_ = new TH3D("datanhit", ";cent;pt;nhit", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, 0, 50);
    data_dca_ = new TH3D("datadca", ";cent;pt;DCA[cm]", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, 0, 3.0);
    data_nhitposs_ = new TH3D("datanhitposs", ";cent;pt;nhitposs", cent_axis_.nBins, cent_axis_.low, cent_axis_.high pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, 0, 50);
    data_eta_ = new TH3D("recoetacut", ";cent;pt;#eta", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, -1, 1);
    data_phi_ = new TH3D("recophicut", ";cent;pt;#phi", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, -TMath::Pi(), TMath::Pi());

    vz_ = new TH1D("vz", ";v_{z}[cm]", 60, -30, 30);
    refmult_ = new TH1D("refmult", ";refmult", 800, 0, 800);
    grefmult_ = new TH1D("grefmult", ";grefmult", 800, 0, 800);
    centrality_ = new TH1D("centrality", ";centrality", cent_axis_.nBins, cent_axis_.low, cent_axis_.high);

    mc_reco_tracks_ = new TH3D("mcrecotracks", ";cent;mc tracks;reco tracks", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, 50, 0, 50, 50, 0, 50);

    mc_tracks_ = new TH2D("mctracks", ";cent;pt", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high);
    reco_tracks_ = new TH2D("recotracks", ";cent;pt", cent_axis.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high);



    return kStOK;
}

bool StEfficiencyAssessor::LoadEvent() {
    muDst_ = muDstMaker_->muDst();
    if (muDst_ == nullptr) {
        LOG_ERROR << "Could not load MuDst" << endm;
        return kStErr;
    }
    muInputEvent_ = muDst_->event();
    if (muInputEvent_ == nullptr) {
        LOG_ERROR << "Could not load MuDstEvent" << endm;
        return kStErr;
    }

    int eventID = muInputEvent_->eventId();
    int runID = muInputEvent_->runId();

    // now try to match the event to a miniMC event in the chain

    int nTries = chain_->GetEntries();

    if (event_->eventId() == eventID &&
            event_->runId() == runID)
        return true;

    while (nTries >= 0) {

        current_++;
        if (current_ >= chain_->GetEntries())
            current_ = 0;
        nTries--;

        chain_->GetEntry(current_);

        if (event_->eventId() == eventID &&
                event_->runId() == runID)
            return true;
    }

    if (nTries < 0) {
        LOG_ERROR << "could not match event to miniMC" << endm;
        return false;
    }

    return true;
}


