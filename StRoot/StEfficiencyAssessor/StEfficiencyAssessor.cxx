#include "StEfficiencyAssessor.hh"

#include "St_base/StMessMgr.h"
#include "StMiniMcEvent/StMiniMcPair.h"
#include "StMiniMcEvent/StTinyMcTrack.h"
#include "StMiniMcEvent/StContamPair.h"

#include "StMuDSTMaker/COMMON/StMuTrack.h"

#include "TMath.h"
#include "TTree.h"

#include "StRefMultCorr/CentralityMaker.h"

#include <iostream>

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
    pt_axis_   = axisDef(20, 0.0, 5.0);
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

    // check event cuts 
    if (!cuts_.AcceptEvent(muInputEvent_))
        return kStOk;

    int centrality = 0;
    if (p18ih_cent_def_ != nullptr) {
        p18ih_cent_def_->setEvent(muInputEvent_->runId(), muInputEvent_->refMult(), muInputEvent_->runInfo().zdcCoincidenceRate(), event_->vertexZ());
        centrality = p18ih_cent_def_->centrality9();
    }
    else if (p16id_cent_def_ != nullptr) {
        p16id_cent_def_->init(muInputEvent_->runId());
        p16id_cent_def_->initEvent(muInputEvent_->grefmult(), muInputEvent_->primaryVertexPosition().z(), muInputEvent_->runInfo().zdcCoincidenceRate());
        centrality = p16id_cent_def_->getCentralityBin9();
    }
    else {
        LOG_ERROR << "error: no centrality definition created" << endm;
        return kStFatal;
    }
    if (centrality < 0 || centrality > 8)
        return kStOK;
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
    
        reco_nhit_->Fill(centrality, pair->ptPr(), pair->fitPts()+1);
        reco_dca_->Fill(centrality, pair->ptPr(), pair->dcaGl());
        reco_eta_->Fill(centrality, pair->ptPr(), pair->etaPr());
        reco_phi_->Fill(centrality, pair->ptPr(), pair->phiPr());
        reco_nhitposs_->Fill(centrality, pair->ptPr(), pair->nPossiblePts()+1);
        reco_fitfrac_->Fill(centrality, pair->ptPr(), (double)(pair->fitPts()+1)/(pair->nPossiblePts()+1));

        if (fabs(pair->etaPr()) > 1.0)
            continue;

        if ((double) (pair->fitPts()+1) / (pair->nPossiblePts()+1) < minFitFrac_)
            continue;
      
        reco_dca_scale_->Fill(centrality, pair->ptPr(), pair->dcaGl());
      
        if (pair->dcaGl() > maxDCA_ || pair->fitPts() < minFit_)
          continue;
      
        count_pair++;
        reco_tracks_->Fill(centrality, pair->ptPr());
        reco_cut_nhit_->Fill(centrality, pair->ptPr(), pair->fitPts()+1);
        reco_cut_dca_->Fill(centrality, pair->ptPr(), pair->dcaGl());
        reco_cut_eta_->Fill(centrality, pair->ptPr(), pair->etaPr());
        reco_cut_phi_->Fill(centrality, pair->ptPr(), pair->phiPr());
        reco_cut_nhitposs_->Fill(centrality, pair->ptPr(), pair->nPossiblePts()+1);
        reco_cut_fitfrac_->Fill(centrality, pair->ptPr(), (double)(pair->fitPts()+1)/(pair->nPossiblePts()+1));
        dca_reco_cut_ext_->Fill(centrality, pair->ptPr(), pair->dcaGl());
    }
    mc_reco_tracks_->Fill(centrality, count_mc, count_pair);

    for (int i = 0; i < muDst_->primaryTracks()->GetEntries(); ++i) {
        StMuTrack* muTrack = (StMuTrack*) muDst_->primaryTracks(i);
        if (muTrack->flag() < 0)
            continue;
        if (muTrack->nHitsFit() < minFit_)
            continue;
        if ((double) muTrack->nHitsFit() / (muTrack->nHitsPoss(kTpcId) + 1) < minFitFrac_)
            continue;
        if (fabs(muTrack->eta()) > 1.0)
            continue;
      
        data_dca_scale_->Fill(centrality, muTrack->pt(), muTrack->dcaGlobal().mag());
      
        if (muTrack->dcaGlobal().mag() > maxDCA_)
          continue;
        
        data_nhit_->Fill(centrality, muTrack->pt(), muTrack->nHitsFit());
        data_dca_->Fill(centrality, muTrack->pt(), muTrack->dcaGlobal().mag());
        data_eta_->Fill(centrality, muTrack->pt(), muTrack->eta());
        data_phi_->Fill(centrality, muTrack->pt(), muTrack->phi());
        data_nhitposs_->Fill(centrality, muTrack->pt(), muTrack->nHitsPoss(kTpcId)+1);
        data_fitfrac_->Fill(centrality, muTrack->pt(), (double)(muTrack->nHitsFit())/(muTrack->nHitsPoss(kTpcId)+1));
        dca_data_cut_ext_->Fill(centrality, muTrack->pt(), muTrack->dcaGlobal().mag());
    }


    return kStOK;
}


Int_t StEfficiencyAssessor::Finish() {
    if (out_ == nullptr) {
        out_ = new TFile("stefficiencyassessor.root", "RECREATE");
    }

    out_->cd();

    vz_->Write();
    refmult_->Write();
    grefmult_->Write();
    centrality_->Write();
    mc_reco_tracks_->Write();
    mc_tracks_->Write();
    reco_tracks_->Write();

    mc_eta_->Write();
    mc_phi_->Write();

    reco_nhit_->Write();
    reco_dca_->Write();
    reco_nhitposs_->Write();
    reco_eta_->Write();
    reco_phi_->Write();
    reco_fitfrac_->Write();
    reco_dca_scale_->Write();
  
    reco_cut_nhit_->Write();
    reco_cut_dca_->Write();
    reco_cut_nhitposs_->Write();
    reco_cut_eta_->Write();
    reco_cut_phi_->Write();
    reco_cut_fitfrac_->Write();

    data_nhit_->Write();
    data_dca_->Write();
    data_nhitposs_->Write();
    data_eta_->Write();
    data_phi_->Write();
    data_fitfrac_->Write();
    data_dca_scale_->Write();

    dca_reco_cut_ext_->Write();
    dca_data_cut_ext_->Write();

    out_->Close();
    return kStOk;
}


int StEfficiencyAssessor::InitInput() {
    muDstMaker_ = (StMuDstMaker*) GetMakerInheritsFrom("StMuDstMaker");
    if (muDstMaker_ == nullptr) {
        LOG_ERROR << "No muDstMaker found in chain: StEfficiencyAssessor init failed" << endm;
        return kStFatal;
    }
    if (TString(muDstMaker_->GetFile()).Contains("SL17d") ||
        TString(muDstMaker_->GetFile()).Contains("SL18f") ||
        TString(muDstMaker_->GetFile()).Contains("SL18h")) {
        p18ih_cent_def_ = new CentralityDef();
        p16id_cent_def_ = nullptr;
    }
    else if (TString(muDstMaker_->GetFile()).Contains("SL16d")) {
        p16id_cent_def_ = CentralityMaker::instance()->getgRefMultCorr_P16id();
        p16id_cent_def_->setVzForWeight(6, -6.0, 6.0);
        p16id_cent_def_->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_vpd30_vpd5_Run14_P16id.txt");
        p18ih_cent_def_ = nullptr;
    }
    else {
        LOG_ERROR << "Library could not be discovered: exiting" << endm;
        return kStFatal;
    }   
    return kStOK;
}

int StEfficiencyAssessor::InitOutput() {
    if (!CheckAxes()) {
        LOG_ERROR << "axes not valid: could not initialize histograms";
        return kStFatal;
    }
    
    mc_eta_ = new TH3D("mceta", ";cent;pt;#eta", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, -1, 1);
    mc_phi_ = new TH3D("mcphi", ";cent;pt;#phi", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, -TMath::Pi(), TMath::Pi());

    reco_nhit_ = new TH3D("reconhit", ";cent;pt;nhit", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, 0, 50);
    reco_dca_ = new TH3D("recodca", ";cent;pt;DCA[cm]", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, 0, 3.0);
    reco_nhitposs_ = new TH3D("reconhitposs", ";cent;pt;nhitposs", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, 0, 50);
    reco_eta_ = new TH3D("recoeta", ";cent;pt;#eta", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, -1, 1);
    reco_phi_ = new TH3D("recophi", ";cent;pt;#phi", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, -TMath::Pi(), TMath::Pi());
    reco_fitfrac_ = new TH3D("recofitfrac", ";cent;pt;fitfrac", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, 0, 1);
  
    reco_dca_scale_ = new TH3D("recodcascale", ";cent;pt;DCA[cm]", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, 0, 3.0);
  
    reco_cut_nhit_ = new TH3D("reconhitcut", ";cent;pt;nhit", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, 0, 50);
    reco_cut_dca_ = new TH3D("recodcacut", ";cent;pt;DCA[cm]", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, 0, 3.0);
    reco_cut_nhitposs_ = new TH3D("reconhitposscut", ";cent;pt;nhitposs", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, 0, 50);
    reco_cut_eta_ = new TH3D("recoetacut", ";cent;pt;#eta", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, -1, 1);
    reco_cut_phi_ = new TH3D("recophicut", ";cent;pt;#phi", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, -TMath::Pi(), TMath::Pi());
    reco_cut_fitfrac_ = new TH3D("recocutfitfrac", ";cent;pt;fitfrac", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, 0, 1);
    
    data_nhit_ = new TH3D("datanhit", ";cent;pt;nhit", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, 0, 50);
    data_dca_ = new TH3D("datadca", ";cent;pt;DCA[cm]", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, 0, 3.0);
    data_nhitposs_ = new TH3D("datanhitposs", ";cent;pt;nhitposs", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, 0, 50);
    data_eta_ = new TH3D("dataeta", ";cent;pt;#eta", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, -1, 1);
    data_phi_ = new TH3D("dataphi", ";cent;pt;#phi", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, -TMath::Pi(), TMath::Pi());
    data_fitfrac_ = new TH3D("datafitfrac", ";cent;pt;fitfrac", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, 0, 1);
  
    data_dca_scale_ = new TH3D("datadcascale", ";cent;pt;DCA[cm]", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high, 50, 0, 3.0);
  
    vz_ = new TH1D("vz", ";v_{z}[cm]", 60, -30, 30);
    refmult_ = new TH1D("refmult", ";refmult", 800, 0, 800);
    grefmult_ = new TH1D("grefmult", ";grefmult", 800, 0, 800);
    centrality_ = new TH1D("centrality", ";centrality", cent_axis_.nBins, cent_axis_.low, cent_axis_.high);

    mc_reco_tracks_ = new TH3D("mcrecotracks", ";cent;mc tracks;reco tracks", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, 50, 0, 50, 50, 0, 50);

    mc_tracks_ = new TH2D("mctracks", ";cent;pt", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high);
    reco_tracks_ = new TH2D("recotracks", ";cent;pt", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, pt_axis_.nBins, pt_axis_.low, pt_axis_.high);

    dca_reco_cut_ext_ = new TH3D("recocutdcaext", ";cent;pt;DCA[cm]", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, 100, pt_axis_.low, pt_axis_.high, 50, 0, 3.0);
    dca_data_cut_ext_ = new TH3D("datadcaext", ";cent;pt;DCA[cm]", cent_axis_.nBins, cent_axis_.low, cent_axis_.high, 100, pt_axis_.low, pt_axis_.high, 50, 0, 3.0);

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


