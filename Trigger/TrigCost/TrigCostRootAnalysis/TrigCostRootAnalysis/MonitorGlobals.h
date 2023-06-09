// Dear emacs, this is -*- c++ -*-

/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// -------------------------------------------------------------
//  author: Tim Martin <Tim.Martin@cern.ch>
// -------------------------------------------------------------
#ifndef TrigCostRootAnalysis_MonitorGlobals_H
#define TrigCostRootAnalysis_MonitorGlobals_H

// STL include(s):
#include <map>
#include <string>
#include <vector>

// Local include(s):
#include "MonitorBase.h"

// ROOT include(s):
#include <TH1.h>
#include <TCanvas.h>

namespace TrigCostRootAnalysis {

  
  /**
   * @class MonitorGlobals
   * Keep track of global quantities for an entire event
   */
  class MonitorGlobals : public MonitorBase {
  
   public:
   
    MonitorGlobals(const TrigCostData* _costData);
    void newEvent(Float_t _weight = 1.);
    CounterBase* newCounter( const std::string &_name, Int_t _ID );
    Bool_t getIfActive(ConfKey_t _mode);
    void saveOutput();
    
  }; //class MonitorGlobals
  
} // namespace TrigCostRootAnalysis

#endif //TrigCostRootAnalysis_MonitorGlobals_H
