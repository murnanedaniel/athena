// Dear emacs, this is -*- c++ -*-

/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// -------------------------------------------------------------
//  author: Tim Martin <Tim.Martin@cern.ch>
// -------------------------------------------------------------
#ifndef TrigCostRootAnalysis_CounterROI_H
#define TrigCostRootAnalysis_CounterROI_H

// STL include(s):
#include <map>
#include <string>

// Local include(s):
#include "CounterBase.h"

// ROOT include(s):
#include <Rtypes.h>

namespace TrigCostRootAnalysis {

  //Forward declaration
  class TrigCostData;

  /**
   * @class CounterROI
   * One counter ROI is created per ROI type (None,Muon,EM/Tau,Jet,JetEt,Energy) along with one 'All' counter.
   * These are processed the relevent ROIs, 'All' is processed on all ROIs
   * @see CounterBase
   * @see MonitorROI
   */
  class CounterROI : public CounterBase {

   public:

    CounterROI( const TrigCostData* _costData, const std::string& _name, Int_t _ID, UInt_t _detailLevel = 10, MonitorBase* _parent = 0 );
    ~CounterROI();
    void startEvent();
    void processEventCounter(UInt_t _e, UInt_t _f, Float_t _weight = 1.);
    void endEvent(Float_t _weight = 1.);
    void debug(UInt_t _e);

   private:

    Double_t getPrescaleFactor(UInt_t _e = INT_MAX);

  }; //class CounterROI

} // namespace TrigCostRootAnalysis

#endif //TrigCostRootAnalysis_CounterROI_H
