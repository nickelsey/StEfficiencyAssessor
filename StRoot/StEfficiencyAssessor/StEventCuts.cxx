#include "StEventCuts.hh"
#include "St_base/StMessMgr.h"

#include <algorithm>
#include <string>
#include <sstream>
#include <istream>

template <typename T> string tostr(const T& t) {
  std::ostringstream os;
  os<<t;
  return os.str();
}

ClassImp(StEventCuts)

StEventCuts::StEventCuts()
  : mNEvents(0), mEventsFailed(0), mEventsFailedVx(0),
    mEventsFailedVy(0),  mEventsFailedVz(0),
    mEventsFailedVr(0), mEventsFailedRef(0),
    mEventsFailedTrigger(), mEventsFailedTriggerTotal(0),
    mCheckVx(kFALSE), mCheckVy(kFALSE), mCheckVz(kFALSE),
    mCheckVr(kFALSE), mCheckRefMult(kFALSE),
    mCheckTrigger(kFALSE), mUseGrefMult(kFALSE), mMinVx(0),
    mMaxVx(0), mMinVy(0), mMaxVy(0), mMinVz(0), mMaxVz(0),
    mMinVr(0), mMaxVr(0), mMinRef(0), mMaxRef(0),
    mTriggers(), mRuns() {}

StEventCuts::~StEventCuts() {
  
}

bool StEventCuts::AcceptEvent(StMuEvent* event) {
  /* first see if the run is masked out, if so, reject */
  if (!AcceptRunId(event))
    return kFALSE;
  
  /* checks all cuts for all events, event if it 
     has already failed one - gives a better idea of 
     relative acceptance rates 
   */
  mNEvents++;
  bool accept_event = kTRUE;
  if (mCheckVx && !AcceptVx(event))           accept_event = kFALSE;
  if (mCheckVy && !AcceptVy(event))           accept_event = kFALSE;
  if (mCheckVz && !AcceptVz(event))           accept_event = kFALSE;
  if (mCheckVr && !AcceptVr(event))           accept_event = kFALSE;
  if (mCheckRefMult && !AcceptRefMult(event)) accept_event = kFALSE;
  if (mCheckTrigger && !AcceptTrigger(event)) accept_event = kFALSE;
  if (accept_event != kTRUE) {mEventsFailed++;
                              LOG_DEBUG << "Event Failed" << endm;}
  else {LOG_DEBUG << "Event Accepted" << endm;}
  return accept_event;
}

Bool_t StEventCuts::AcceptVx(StMuEvent* event) {
  Double_t vx = event->primaryVertexPosition().x();
  if (vx > mMaxVx || vx < mMinVx) {
    LOG_DEBUG << "Fail: Vx=" << vx << " minVx=" << mMinVx << " maxVx=" << mMaxVx << endm;
    mEventsFailedVx++;
    return kFALSE;
  }
  LOG_DEBUG << "Pass: Vx=" << vx << " minVx=" << mMinVx << " maxVx=" << mMaxVx << endm;
  return kTRUE;
}

Bool_t StEventCuts::AcceptVy(StMuEvent* event) {
  Double_t vy = event->primaryVertexPosition().y();
  if (vy > mMaxVy || vy < mMinVy) {
    LOG_DEBUG << "Fail: Vy=" << vy << " minVy=" << mMinVy << " maxVy=" << mMaxVy << endm;
    mEventsFailedVy++;
    return kFALSE;
  }
  LOG_DEBUG << "Pass: Vy=" << vy << " minVy=" << mMinVy << " maxVy=" << mMaxVy << endm;
  return kTRUE;
}

Bool_t StEventCuts::AcceptVz(StMuEvent* event) {
  Double_t vz = event->primaryVertexPosition().z();
  if (vz > mMaxVz || vz < mMinVz) {
    LOG_DEBUG << "Fail: Vz=" << vz << " minVz=" << mMinVz << " maxVz=" << mMaxVz << endm;
    mEventsFailedVz++;
    return kFALSE;
  }
  LOG_DEBUG << "Pass: Vz=" << vz << " minVz=" << mMinVz << " maxVz=" << mMaxVz << endm;
  return kTRUE;
}

Bool_t StEventCuts::AcceptVr(StMuEvent* event) {
  Double_t vr = sqrt(pow(event->primaryVertexPosition().x(), 2) +
                     pow(event->primaryVertexPosition().y(), 2));
  if (vr > mMaxVr || vr < mMinVr) {
    LOG_DEBUG << "Fail: Vr=" << vr << " minVr=" << mMinVr << " maxVr=" << mMaxVr << endm;
    mEventsFailedVr++;
    return kFALSE;
  }
  LOG_DEBUG << "Pass: Vr=" << vr << " minVr=" << mMinVr << " maxVr=" << mMaxVr << endm;
  return kTRUE;
}

Bool_t StEventCuts::AcceptTrigger(StMuEvent* event) {
  Bool_t accept = kFALSE;
  for (unsigned int i = 0; i < mTriggers.size(); ++i) {
    if (event->triggerIdCollection().nominal().isTrigger(mTriggers[i]) == kTRUE) {
      accept = kTRUE;
      LOG_DEBUG << "Pass: trigger accepted=" << mTriggers[i] << endm;
    } else {mEventsFailedTrigger[i]++;}
  }
  if (accept == kFALSE) {mEventsFailedTriggerTotal++; LOG_DEBUG << "Fail: no trigger accepted" << endm;}
  return accept;
}

Bool_t StEventCuts::AcceptRefMult(StMuEvent* event) {
  UInt_t refmult = mUseGrefMult ? event->grefmult() : event->refMult();
  if (refmult > mMaxRef || refmult < mMinRef) {
    LOG_DEBUG << "Fail: refmult=" << refmult << " min refmult=" << mMinRef << " max refmult=" << mMaxRef << endm;
    mEventsFailedRef++;
    return kFALSE;
  }
  LOG_DEBUG << "Pass: refmult=" << refmult << " min refmult=" << mMinRef << " max refmult=" << mMaxRef << endm;
  return kTRUE;
}

Bool_t StEventCuts::AcceptRunId(StMuEvent* event) {
  UInt_t runid = event->runId();
  if (mRuns.find(runid) != mRuns.end()) {
    LOG_DEBUG <<"Fail: runid=" << runid << " is masked out" << endm;
    return kFALSE;
  }
  
  LOG_DEBUG << "Pass: runid=" << runid << endm;
  return kTRUE;
}

void  StEventCuts::SetVxRange(Double_t min, Double_t max) {
  mCheckVx = kTRUE;
  mMinVx = min;
  mMaxVx = max;
}

void   StEventCuts::SetVyRange(Double_t min, Double_t max) {
  mCheckVy = kTRUE;
  mMinVy = min;
  mMaxVy = max;
}

void   StEventCuts::SetVzRange(Double_t min, Double_t max) {
  mCheckVz = kTRUE;
  mMinVz = min;
  mMaxVz = max;
}

void   StEventCuts::SetVrRange(Double_t min, Double_t max) {
  mCheckVr = kTRUE;
  mMinVr = min;
  mMaxVr = max;
}

void   StEventCuts::SetRefMultRange(UInt_t min, UInt_t max) {
  mCheckRefMult = kTRUE;
  mMinRef = min;
  mMaxRef = max;
}

void   StEventCuts::UsegRefMult(bool use_gref) {
  mUseGrefMult = use_gref;
}

void StEventCuts::AddTrigger(unsigned int trigger) {
  if(std::find(mTriggers.begin(), mTriggers.end(), trigger) == mTriggers.end()) {
    mCheckTrigger = kTRUE;
    mTriggers.push_back(trigger);
    mEventsFailedTrigger.push_back(0);
  }
}

void StEventCuts::AddTrigger(std::vector<unsigned int> triggers) {
  for (unsigned id = 0; id < triggers.size(); ++id) {
    UInt_t trigger = triggers[id];
    if(std::find(mTriggers.begin(), mTriggers.end(), trigger) == mTriggers.end()) {
      mCheckTrigger = kTRUE;
      mTriggers.push_back(trigger);
      mEventsFailedTrigger.push_back(0);
    }
  }
}

bool StEventCuts::MaskRuns(std::string filename) {
  LOG_DEBUG << "Mask Run file: " << filename << endm;
  std::ifstream file(filename.c_str());
  
  if (!file.good()) {
    LOG_ERROR << "Can't read file: " << filename << endm;
    return false;
  }
  
  std::string line;
  while (std::getline (file, line)){
    if (line.size()==0) continue; // skip empty lines
    if (line[0] == '#') continue; // skip comments
    
    std::istringstream ss(line);
    while(ss){
      std::string entry;
      std::getline(ss, entry, ',');
      int id = atoi(entry.c_str());
      if (id) {
        mRuns.insert(id);
        LOG_DEBUG << "Added masked Run: " << id << endm;
      }
    }
  }
  
  return true;
}

void StEventCuts::MaskRuns(UInt_t run) {
  mRuns.insert(run);
}

void StEventCuts::PrintCuts() {
  LOG_INFO << "// ------------------ StEventCuts Cuts ------------------ //" << endm;
  LOG_INFO << endm;
  if (mCheckVx || mCheckVy || mCheckVz) {LOG_INFO << "Vertex cuts: " << endm;}
  if (mCheckVx)                         {LOG_INFO << mMinVx << " < Vx < " << mMaxVx << endm;}
  if (mCheckVy)                         {LOG_INFO << mMinVy << " < Vy < " << mMaxVy << endm;}
  if (mCheckVz)                         {LOG_INFO << mMinVz << " < Vz < " << mMaxVz << endm;}
  if (mCheckVr)                         {LOG_INFO << mMinVr << " < Vr < " << mMaxVr << endm;}
  if (mCheckRefMult)                    {LOG_INFO << " Refmult Cuts: " << endm;}
  if (mCheckRefMult && mUseGrefMult)    {LOG_INFO << "use grefmult" << endm;
                                         LOG_INFO << mMinRef << " < gRefMult < " << mMaxRef << endm;}
  if (mCheckRefMult && !mUseGrefMult)   {LOG_INFO << "use refmult" << endm;
                                         LOG_INFO << mMinRef << " < RefMult < " << mMaxRef << endm;}
  
  if (mCheckTrigger && mTriggers.size() > 0) {
    std::string trigger_string = "[ " + tostr(mTriggers[0]);
    for (unsigned i = 1; i < mTriggers.size(); ++i) {
      trigger_string += ", " + tostr(mTriggers[i]);
    }
    LOG_INFO << "Using triggers: " << trigger_string << endm;
  }
  
  if (mRuns.size() > 0) {
    std::string run_string = "[ ";
    UInt_t count = 0;
    for (std::set<UInt_t>::iterator run = mRuns.begin(); run != mRuns.end(); ++run) {
      count++;
      run_string += tostr(*run) + (count == mRuns.size() ? " ]" : ", ");
    }
    LOG_INFO << "Masking runs: " << run_string << endm;
  }
  
  LOG_INFO << " // ------------------        End Cuts        ------------------ //" << endm;
}

void  StEventCuts::PrintStats() {
  LOG_INFO << "// ------------------ StEventCuts Stats ------------------ //" << endm;
  LOG_INFO << "after removing masked runs" << endm;
  LOG_INFO << "number of events:   " << mNEvents << endm;
  
  if (mCheckVx) {LOG_INFO << "events rejected by Vx cut: " << mEventsFailedVx << endm;
                 LOG_INFO << "\t percent loss: " << (Double_t) mEventsFailedVx / mNEvents << endm;}
  if (mCheckVy) {LOG_INFO << "events rejected by Vy cut: " << mEventsFailedVy << endm;
                 LOG_INFO << "\t percent loss: " << (Double_t) mEventsFailedVy / mNEvents << endm;}
  if (mCheckVz) {LOG_INFO << "events rejected by Vz cut: " << mEventsFailedVz << endm;
                 LOG_INFO << "\t percent loss: " << (Double_t) mEventsFailedVz / mNEvents << endm;}
  if (mCheckVz) {LOG_INFO << "events rejected by Vr cut: " << mEventsFailedVr << endm;
                 LOG_INFO << "\t percent loss: " << (Double_t) mEventsFailedVr / mNEvents << endm;}
  
  if (mCheckRefMult) {
    LOG_INFO << "events rejected by RefMult cut: " << mEventsFailedVx << endm;
    LOG_INFO << "\t using grefmult: "; if (mUseGrefMult) {LOG_INFO << "true" << endm;} else {LOG_INFO << "false" << endm;}
    LOG_INFO << "\t percent loss: " << (Double_t) mEventsFailedRef / mNEvents << endm;
 }
  
  if (mCheckTrigger && mTriggers.size() > 0) {
    std::string trigger_string = "[ " + tostr(mTriggers[0]);
    std::string loss_string = "[ " + tostr(mEventsFailedTrigger[0]);
    for (unsigned i = 1; i < mTriggers.size(); ++i) {
      trigger_string += ", " + tostr(mTriggers[i]);
      loss_string += ", " + tostr(mEventsFailedTrigger[i]);
    }
    trigger_string += " ]"; loss_string += " ]";
    LOG_INFO << "events rejected by triggers: " << endm;
    LOG_INFO << "using triggers: " << trigger_string << endm;
    LOG_INFO << "rejected: " << loss_string << endm;
    LOG_INFO << "total events rejected by trigger: " << mEventsFailedTriggerTotal << endm;
  }
  
  LOG_INFO << "total number of events rejected: " << mEventsFailed << endm;
  LOG_INFO << "\t percent loss: " << (Double_t) mEventsFailed / mNEvents << endm;
  LOG_INFO << " // ------------------        End Stats        ------------------ //" << endm;
  
}

