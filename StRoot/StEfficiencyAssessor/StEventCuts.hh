/* internal class for StEfficiencyAssessor
   implements cuts on vertex, trigger selection,
   etc.
 */


#ifndef STEVENTCUTS_HH
#define STEVENTCUTS_HH

#include "TObject.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"

#include <vector>
#include <set>

class StEventCuts : public TObject {
  
public:
  
  StEventCuts();
  ~StEventCuts();
  
  /* function called by StEfficiencyAssessor, returns true if
     event should be accepted based on the defined cuts
   */
  Bool_t AcceptEvent(StMuEvent* event);
  
  /* there is no default cut value - cuts are turned off 
     until set by the user. Each must be set individually
   */
  void SetVxRange(Double_t min, Double_t max);
  void SetVyRange(Double_t min, Double_t max);
  void SetVzRange(Double_t min, Double_t max);
  void SetVrRange(Double_t min, Double_t max);
  void SetRefMultRange(UInt_t min, UInt_t max);
  void UsegRefMult(Bool_t use_grefmult);
  
  // set runs to be masked out
  bool MaskRuns(std::string filename);
  void MaskRuns(UInt_t run);
  
  std::set<UInt_t> MaskRuns() const {return mRuns;}
  
  /* by default all triggers are accepted. Once a trigger is
     added, only those triggers are used
   */
  void AddTrigger(UInt_t trigger);
  void AddTrigger(std::vector<UInt_t> triggers);
  
  std::vector<UInt_t> Triggers() const {return mTriggers;}
  
  /* print a list of the cuts & triggers used */
  void PrintCuts();
  
  /* print statistics on the number of events rejected. 
     Called in Finish() of TStarJetPicoMaker */
  void PrintStats();
  
  // access to cuts
  inline Double_t MinVz() const {return mMinVz;}
  inline Double_t MaxVz() const {return mMaxVz;}
  inline Double_t MinVy() const {return mMinVy;}
  inline Double_t MaxVy() const {return mMaxVy;}
  inline Double_t MinVx() const {return mMinVx;}
  inline Double_t MaxVx() const {return mMaxVx;}
  inline Double_t MinVr() const {return mMinVr;}
  inline Double_t MaxVr() const {return mMaxVr;}
  inline UInt_t  MinRef() const {return mMinRef;}
  inline UInt_t  MaxRef() const {return mMaxRef;}
  
private:
  
  /* functions used by AcceptEvent() to determine which
     criteria the event passes or fails. Each give their
     own output & iterate their respective counters
   */
  Bool_t AcceptVx(StMuEvent* event);
  Bool_t AcceptVy(StMuEvent* event);
  Bool_t AcceptVz(StMuEvent* event);
  Bool_t AcceptVr(StMuEvent* event);
  Bool_t AcceptRefMult(StMuEvent* event);
  Bool_t AcceptTrigger(StMuEvent* event);
  Bool_t AcceptRunId(StMuEvent* event);

  UInt_t         mNEvents;
  UInt_t         mEventsFailed;
  UInt_t         mEventsFailedVx;
  UInt_t         mEventsFailedVy;
  UInt_t         mEventsFailedVz;
  UInt_t         mEventsFailedVr;
  UInt_t         mEventsFailedRef;
  std::vector<UInt_t> mEventsFailedTrigger;
  UInt_t         mEventsFailedTriggerTotal;
  
  Bool_t    mCheckVx;
  Bool_t    mCheckVy;
  Bool_t    mCheckVz;
  Bool_t    mCheckVr;
  Bool_t    mCheckRefMult;
  Bool_t    mCheckTrigger;
  Bool_t    mUseGrefMult;
  
  Double_t  mMinVx,  mMaxVx;
  Double_t  mMinVy,  mMaxVy;
  Double_t  mMinVz,  mMaxVz;
  Double_t  mMinVr,  mMaxVr;
  UInt_t    mMinRef, mMaxRef;
  
  std::vector<UInt_t> mTriggers;
  std::set<UInt_t>    mRuns;
  
  
  ClassDef(StEventCuts, 1)
};

#endif /* STEVENTCUTS_HH */
