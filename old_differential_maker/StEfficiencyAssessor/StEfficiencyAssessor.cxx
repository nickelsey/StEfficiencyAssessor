#include "StEfficiencyAssessor.hh"

#include "St_base/StMessMgr.h"
#include "StMiniMcEvent/StMiniMcPair.h"
#include "StMiniMcEvent/StTinyMcTrack.h"
#include "StMiniMcEvent/StContamPair.h"

#include "StMuDSTMaker/COMMON/StMuTrack.h"

#include "TMath.h"
#include "TTree.h"

ClassImp(StEfficiencyAssessor)

StEfficiencyAssessor::StEfficiencyAssessor(TChain* mcTree, std::string outputFile) {
  if (!LoadTree(mcTree)) {
    LOG_ERROR << "load chain failed" << endm;
  }
  
  SetDefaultAxes();
  
  current_ = 0;
  
  muDstMaker_ = nullptr;
  muDst_ = nullptr;
  muInputEvent_ = nullptr;
  
  minFit_ = 10;
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
  cent_axis_ = axisDef(16, 0, 1);
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
  
  // get luminosity bin
  double zdcAnd = muInputEvent_->runInfo().zdcCoincidenceRate();
  int zdcBin = lumi_axis_.bin(zdcAnd);
  
  // get centrality
  cent_def_.setEvent(muInputEvent_->runId(), muInputEvent_->refMult(), zdcAnd, event_->vertexZ());
  int centBin = cent_def_.centrality16();
  
  // finally get vz bin
  double vz = event_->vertexZ();
  int vzBin = vz_axis_.bin(vz);
  
  if (zdcBin < 0 || centBin < 0 || vzBin < 0)
    return kStOK;
  
  // add event to our list of saved events
  runids.push_back(muInputEvent_->runId());
  eventids.push_back(muInputEvent_->eventId());
  
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

    if (pair->dcaGl() > maxDCA_ || pair->fitPts() < minFit_)
      continue;
    
    if (fabs(pair->etaPr()) > 1.0)
      continue;
    
    int ptBin = pt_axis_.bin(pair->ptPr());
    
    if (ptBin < 0)
      continue;
    
    etaMatch->Fill(pair->etaPr(), pair->ptPr());
    phiMatch->Fill(pair->phiPr(), pair->ptPr());
    
    double jacobian_inverse = 1.0 / (pow(pair->dcaGl(), 2.0) * asin(pair->dcaXYGl() / pair->dcaGl()));
    weighted_dca_match->Fill(centBin, pair->ptPr(), pair->dcaGl(), fabs(jacobian_inverse));
    
    matched_dca_[zdcBin][centBin][vzBin][ptBin]->Fill(pair->phiPr(), pair->etaPr(), pair->dcaGl());
    matched_nhit_[zdcBin][centBin][vzBin][ptBin]->Fill(pair->phiPr(), pair->etaPr(), pair->fitPts());
    matched_nhit_pos_[zdcBin][centBin][vzBin][ptBin]->Fill(pair->phiPr(), pair->etaPr(), pair->nPossiblePts());

    count_pair++;

  }
  
  StMuTrack* muTrack;
  for (int i = 0; i < muDst_->primaryTracks()->GetEntries(); ++i) {
    muTrack = (StMuTrack*) muDst_->primaryTracks(i);
    if (muTrack->flag() < 0 || muTrack->dcaGlobal().mag() > maxDCA_ || muTrack->nHitsFit() < minFit_)
      continue;
    
    if (fabs(muTrack->eta()) > 1.0)
      continue;
    
    int ptBin = pt_axis_.bin(muTrack->pt());
    
    if (ptBin < 0)
      continue;
    
    etaPrim->Fill(muTrack->eta(), muTrack->pt());
    phiPrim->Fill(muTrack->phi(), muTrack->pt());
    
    data_dca_[zdcBin][centBin][vzBin][ptBin]->Fill(muTrack->phi(), muTrack->eta(), muTrack->dcaGlobal().mag());
    data_nhit_[zdcBin][centBin][vzBin][ptBin]->Fill(muTrack->phi(), muTrack->eta(), muTrack->nHitsFit());
    data_nhit_pos_[zdcBin][centBin][vzBin][ptBin]->Fill(muTrack->phi(), muTrack->eta(), muTrack->nHitsPoss());
    
    double jacobian_inverse = 1.0 / (pow(muTrack->dcaGlobal().mag(), 2.0) * asin(muTrack->dcaGlobal().perp() / muTrack->dcaGlobal().mag()));
    weighted_dca_data->Fill(centBin, muTrack->pt(), muTrack->dcaGlobal().mag(), fabs(jacobian_inverse));
    
    if (fabs(muTrack->nSigmaPion()) < 2.0) {
      
      etaPrimPion->Fill(muTrack->eta(), muTrack->pt());
      phiPrimPion->Fill(muTrack->phi(), muTrack->pt());
      
      weighted_dca_data_pion->Fill(centBin, muTrack->pt(), muTrack->dcaGlobal().mag(), fabs(jacobian_inverse));
      
      data_pion_dca_[zdcBin][centBin][vzBin][ptBin]->Fill(muTrack->phi(), muTrack->eta(), muTrack->dcaGlobal().mag());
      data_pion_nhit_[zdcBin][centBin][vzBin][ptBin]->Fill(muTrack->phi(), muTrack->eta(), muTrack->nHitsFit());
      data_pion_nhit_pos_[zdcBin][centBin][vzBin][ptBin]->Fill(muTrack->phi(), muTrack->eta(), muTrack->nHitsPoss());
    }
  }
  
  
  nPrimaries_->Fill(muDst_->primaryTracks()->GetEntries());
  refMult_->Fill(muDst_->event()->refMult());
  
  nMC_->Fill(mc_array->GetEntries());
  nMCPrim_->Fill(event_->mcMult());
  nMCRefMult_->Fill(event_->nMcNch());
  nMCNoParent_->Fill(count_mc);
  nMatched_->Fill(match_array->GetEntries());
  
  return kStOK;
}

Int_t StEfficiencyAssessor::Finish() {
  if (out_ == nullptr) {
    out_ = new TFile("stefficiencyassessor.root", "RECREATE");
  }
  
  out_->cd();
  
  // write runids and eventids to file
  TTree* tree = new TTree("eventid","event id tree");
  tree->Branch("runids", &runids);
  tree->Branch("eventids", &eventids);
  tree->Fill();
  
  tree->Write();
  
  nPrimaries_->Write();
  refMult_->Write();
  
  nMC_->Write();
  nMCPrim_->Write();
  nMCRefMult_->Write();
  nMCNoParent_->Write();
  nMatched_->Write();
  
  etaMatch->Write();
  phiMatch->Write();
  
  etaPrim->Write();
  phiPrim->Write();
  
  etaPrimPion->Write();
  phiPrimPion->Write();
  
  weighted_dca_match->Write();
  weighted_dca_data->Write();
  weighted_dca_data_pion->Write();
  
  for (unsigned i = 0; i < matched_dca_.size(); ++i) {
    for (unsigned j = 0; j < matched_dca_[i].size(); ++j) {
      for (unsigned k = 0; k < matched_dca_[i][j].size(); ++k) {
        for (unsigned l = 0; l < matched_dca_[i][j][k].size(); ++l) {
          matched_dca_[i][j][k][l]->Write();
          matched_nhit_[i][j][k][l]->Write();
          matched_nhit_pos_[i][j][k][l]->Write();
          
          data_dca_[i][j][k][l]->Write();
          data_nhit_[i][j][k][l]->Write();
          data_nhit_pos_[i][j][k][l]->Write();

          data_pion_dca_[i][j][k][l]->Write();
          data_pion_nhit_[i][j][k][l]->Write();
          data_pion_nhit_pos_[i][j][k][l]->Write();
        }
      }
    }
  }
  
  out_->Close();
  return kStOk;
}


int StEfficiencyAssessor::InitInput() {
  muDstMaker_ = (StMuDstMaker*) GetMakerInheritsFrom("StMuDstMaker");
  if (muDstMaker_ == nullptr) {
    LOG_ERROR << "No muDstMaker found in chain: StEfficiencyAssessor init failed" << endm;
    return kStFatal;
  }
  return kStOK;
}

int StEfficiencyAssessor::InitOutput() {
  if (!CheckAxes()) {
    LOG_ERROR << "axes not valid: could not initialize histograms";
    return kStFatal;
  }
  
  nPrimaries_ = new TH1D("nprimaries", ";primaries", 400, 0, 2000);
  refMult_ = new TH1D("refmult", ";refmult", 800, 0, 800);
  
  nMC_ = new TH1D("nmctracks", ";mc tracks", 200, 0, 200);
  nMCPrim_ = new TH1D("nmcprimtracks", ";primary mc", 100, 0, 100);
  nMCRefMult_ = new TH1D("nmcrefmult", ";mc refmult", 100, 0, 100);
  nMCNoParent_ = new TH1D("nembed", ";embedded tracks", 100, 0, 100);
  nMatched_ = new TH1D("nmatched", ";matched tracks", 100, 0, 100);
  
  etaMatch = new TH2D("etamatch", ";#eta", 100, -1 , 1, 50, 0, 5.0);
  phiMatch = new TH2D("phimatch", ";#phi", 100, -TMath::Pi(), TMath::Pi(), 50, 0, 5.0);
  
  etaPrim = new TH2D("etaprim", ";#eta", 100, -1 , 1, 50, 0, 5.0);
  phiPrim = new TH2D("phiprim", ";#phi", 100, -TMath::Pi(), TMath::Pi(), 50, 0, 5.0);
  
  etaPrimPion = new TH2D("etaprimpion", ";#eta", 100, -1 , 1, 50, 0, 5.0);
  phiPrimPion = new TH2D("phiprimpion", ";#phi", 100, -TMath::Pi(), TMath::Pi(), 50, 0, 5.0);
  
  weighted_dca_match = new TH3F("weighteddcamatch", ";cent;pt;dca", 16, 0, 16, 50, 0, 5, 30, 0, 3);
  weighted_dca_data = new TH3F("weighteddcadata", ";cent;pt;dca", 16, 0, 16, 50, 0, 5, 30, 0, 3);
  weighted_dca_data_pion = new TH3F("weighteddcadatapion", ";cent;pt;dca", 16, 0, 16, 50, 0, 5, 30, 0, 3);
  
  
  matched_dca_ = std::vector<std::vector<std::vector<std::vector<TH3F*>>>>(lumi_axis_.nBins, std::vector<std::vector<std::vector<TH3F*>>>(cent_axis_.nBins, std::vector<std::vector<TH3F*>>(vz_axis_.nBins, std::vector<TH3F*>(pt_axis_.nBins, nullptr))));
  data_dca_ = std::vector<std::vector<std::vector<std::vector<TH3F*>>>>(lumi_axis_.nBins, std::vector<std::vector<std::vector<TH3F*>>>(cent_axis_.nBins, std::vector<std::vector<TH3F*>>(vz_axis_.nBins, std::vector<TH3F*>(pt_axis_.nBins, nullptr))));
  data_pion_dca_ = std::vector<std::vector<std::vector<std::vector<TH3F*>>>>(lumi_axis_.nBins, std::vector<std::vector<std::vector<TH3F*>>>(cent_axis_.nBins, std::vector<std::vector<TH3F*>>(vz_axis_.nBins, std::vector<TH3F*>(pt_axis_.nBins, nullptr))));
  
  matched_nhit_ = std::vector<std::vector<std::vector<std::vector<TH3F*>>>>(lumi_axis_.nBins, std::vector<std::vector<std::vector<TH3F*>>>(cent_axis_.nBins, std::vector<std::vector<TH3F*>>(vz_axis_.nBins, std::vector<TH3F*>(pt_axis_.nBins, nullptr))));
  data_nhit_ = std::vector<std::vector<std::vector<std::vector<TH3F*>>>>(lumi_axis_.nBins, std::vector<std::vector<std::vector<TH3F*>>>(cent_axis_.nBins, std::vector<std::vector<TH3F*>>(vz_axis_.nBins, std::vector<TH3F*>(pt_axis_.nBins, nullptr))));
  data_pion_nhit_ = std::vector<std::vector<std::vector<std::vector<TH3F*>>>>(lumi_axis_.nBins, std::vector<std::vector<std::vector<TH3F*>>>(cent_axis_.nBins, std::vector<std::vector<TH3F*>>(vz_axis_.nBins, std::vector<TH3F*>(pt_axis_.nBins, nullptr))));
  
  matched_nhit_pos_ = std::vector<std::vector<std::vector<std::vector<TH3F*>>>>(lumi_axis_.nBins, std::vector<std::vector<std::vector<TH3F*>>>(cent_axis_.nBins, std::vector<std::vector<TH3F*>>(vz_axis_.nBins, std::vector<TH3F*>(pt_axis_.nBins, nullptr))));
  data_nhit_pos_ = std::vector<std::vector<std::vector<std::vector<TH3F*>>>>(lumi_axis_.nBins, std::vector<std::vector<std::vector<TH3F*>>>(cent_axis_.nBins, std::vector<std::vector<TH3F*>>(vz_axis_.nBins, std::vector<TH3F*>(pt_axis_.nBins, nullptr))));
  data_pion_nhit_pos_ = std::vector<std::vector<std::vector<std::vector<TH3F*>>>>(lumi_axis_.nBins, std::vector<std::vector<std::vector<TH3F*>>>(cent_axis_.nBins, std::vector<std::vector<TH3F*>>(vz_axis_.nBins, std::vector<TH3F*>(pt_axis_.nBins, nullptr))));
  
  for (unsigned lumi = 0; lumi < lumi_axis_.nBins; ++lumi) {
    for (unsigned cent = 0; cent < cent_axis_.nBins; ++cent) {
      for (unsigned vz = 0; vz < vz_axis_.nBins; ++vz) {
        for (unsigned pt = 0; pt < pt_axis_.nBins; ++pt) {
          std::string matched_dca_name = "match_lumi_" + std::to_string(lumi) +
          "_cent_" + std::to_string(cent) + "_vz_" + std::to_string(vz) + "_pt_" + std::to_string(pt) + "_dca";
          std::string matched_nhit_name = "match_lumi_" + std::to_string(lumi) +
          "_cent_" + std::to_string(cent) + "_vz_" + std::to_string(vz) + "_pt_" + std::to_string(pt) + "_nhit";
          std::string matched_nhit_pos_name = "match_lumi_" + std::to_string(lumi) +
          "_cent_" + std::to_string(cent) + "_vz_" + std::to_string(vz) + "_pt_" + std::to_string(pt) + "_nhit_pos";

          matched_dca_[lumi][cent][vz][pt] = new TH3F(matched_dca_name.c_str(), ";#phi;#eta;dca",
                                                       phi_axis_.nBins, phi_axis_.low, phi_axis_.high,
                                                       eta_axis_.nBins, eta_axis_.low, eta_axis_.high,
                                                       30, 0, 3);
          matched_dca_[lumi][cent][vz][pt]->Sumw2();
          matched_nhit_[lumi][cent][vz][pt] = new TH3F(matched_nhit_name.c_str(), ";#phi;#eta;dca",
                                               phi_axis_.nBins, phi_axis_.low, phi_axis_.high,
                                               eta_axis_.nBins, eta_axis_.low, eta_axis_.high,
                                               40, 10, 50);
          matched_nhit_[lumi][cent][vz][pt]->Sumw2();

          matched_nhit_pos_[lumi][cent][vz][pt] = new TH3F(matched_nhit_pos_name.c_str(), ";#phi;#eta;dca",
                                                            phi_axis_.nBins, phi_axis_.low, phi_axis_.high,
                                                            eta_axis_.nBins, eta_axis_.low, eta_axis_.high,
                                                            40, 10, 50);
          matched_nhit_pos_[lumi][cent][vz][pt]->Sumw2();

          std::string data_dca_name = "data_lumi_" + std::to_string(lumi) +
          "_cent_" + std::to_string(cent) + "_vz_" + std::to_string(vz) + "_pt_" + std::to_string(pt) + "_dca";
          std::string data_nhit_name = "data_lumi_" + std::to_string(lumi) +
          "_cent_" + std::to_string(cent) + "_vz_" + std::to_string(vz) + "_pt_" + std::to_string(pt) + "_nhit";
          std::string data_nhit_pos_name = "data_lumi_" + std::to_string(lumi) +
          "_cent_" + std::to_string(cent) + "_vz_" + std::to_string(vz) + "_pt_" + std::to_string(pt) + "_nhit_pos";

          data_dca_[lumi][cent][vz][pt] = new TH3F(data_dca_name.c_str(), ";phi;#eta;dca",
                                                    phi_axis_.nBins, phi_axis_.low, phi_axis_.high,
                                                    eta_axis_.nBins, eta_axis_.low, eta_axis_.high,
                                                    30, 0, 3);
          data_dca_[lumi][cent][vz][pt]->Sumw2();
          data_nhit_[lumi][cent][vz][pt] = new TH3F(data_nhit_name.c_str(), ";#phi;#eta;dca",
                                                     phi_axis_.nBins, phi_axis_.low, phi_axis_.high,
                                                     eta_axis_.nBins, eta_axis_.low, eta_axis_.high,
                                                     40, 10, 50);
          data_nhit_[lumi][cent][vz][pt]->Sumw2();

          data_nhit_pos_[lumi][cent][vz][pt] = new TH3F(data_nhit_pos_name.c_str(), ";#phi;#eta;dca",
                                                         phi_axis_.nBins, phi_axis_.low, phi_axis_.high,
                                                         eta_axis_.nBins, eta_axis_.low, eta_axis_.high,
                                                         40, 10, 50);
          data_nhit_pos_[lumi][cent][vz][pt]->Sumw2();

          std::string data_pion_dca_name = "data_pion_lumi_" + std::to_string(lumi) +
          "_cent_" + std::to_string(cent) + "_vz_" + std::to_string(vz) + "_pt_" + std::to_string(pt) + "_dca";
          std::string data_pion_nhit_name = "data_pion_lumi_" + std::to_string(lumi) +
          "_cent_" + std::to_string(cent) + "_vz_" + std::to_string(vz) + "_pt_" + std::to_string(pt) + "_nhit";
          std::string data_pion_nhit_pos_name = "data_pion_lumi_" + std::to_string(lumi) +
          "_cent_" + std::to_string(cent) + "_vz_" + std::to_string(vz) + "_pt_" + std::to_string(pt) + "_nhit_pos";

          data_pion_dca_[lumi][cent][vz][pt] = new TH3F(data_pion_dca_name.c_str(), ";phi;#eta;dca",
                                                    phi_axis_.nBins, phi_axis_.low, phi_axis_.high,
                                                    eta_axis_.nBins, eta_axis_.low, eta_axis_.high,
                                                    30, 0, 3);
          data_pion_dca_[lumi][cent][vz][pt]->Sumw2();
          data_pion_nhit_[lumi][cent][vz][pt] = new TH3F(data_pion_nhit_name.c_str(), ";#phi;#eta;dca",
                                                     phi_axis_.nBins, phi_axis_.low, phi_axis_.high,
                                                     eta_axis_.nBins, eta_axis_.low, eta_axis_.high,
                                                     40, 10, 50);
          data_pion_nhit_[lumi][cent][vz][pt]->Sumw2();

          data_pion_nhit_pos_[lumi][cent][vz][pt] = new TH3F(data_pion_nhit_pos_name.c_str(), ";#phi;#eta;dca",
                                                              phi_axis_.nBins, phi_axis_.low, phi_axis_.high,
                                                              eta_axis_.nBins, eta_axis_.low, eta_axis_.high,
                                                              40, 10, 50);
          data_pion_nhit_pos_[lumi][cent][vz][pt]->Sumw2();

        }
      }
    }
  }
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


