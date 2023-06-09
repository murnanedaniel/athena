// Dear emacs, this is -*- c++ -*-

/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// -------------------------------------------------------------
//  author: Tim Martin <Tim.Martin@cern.ch>
// -------------------------------------------------------------
#ifndef TrigCostRootAnalysis_MonitorChain_H
#define TrigCostRootAnalysis_MonitorChain_H

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
   * @class MonitorChain
   * Chain monitor implimentation
   */
  class MonitorChain : public MonitorBase {
  
   public:
   
    MonitorChain(const TrigCostData* _costData);
    void newEvent(Float_t _weight = 1.);
    CounterBase* newCounter( const std::string &_name, Int_t _ID );
    Bool_t getIfActive(ConfKey_t _mode);
    void saveOutput();
    
  }; //class MonitorChain
  
} // namespace TrigCostRootAnalysis

#endif //TrigCostRootAnalysis_MonitorChain_H
