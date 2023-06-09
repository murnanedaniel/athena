// Dear emacs, this is -*- c++ -*-

/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// -------------------------------------------------------------
//  author: Tim Martin <Tim.Martin@cern.ch>
// -------------------------------------------------------------
#ifndef TrigCostRootAnalysis_CounterAlgorithm_H
#define TrigCostRootAnalysis_CounterAlgorithm_H

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
   * @class CounterAlgorithm
   * A Counter object, responsible for all the low-level monitoring performed on individual HLT algorithms.
   * Used by MonitorAlgorithm and by MonitorFullEvent, where one instance of CounterAlgorithm is used for each alg
   * in each monitored full event.
   * @see CounterBase
   */
  class CounterAlgorithm : public CounterBase {

   public:

    CounterAlgorithm( const TrigCostData* _costData, const std::string& _name, Int_t _ID, UInt_t _detailLevel = 10, MonitorBase* _parent = 0 );
    ~CounterAlgorithm();
    void startEvent();
    void processEventCounter(UInt_t _e, UInt_t _f,  Float_t _weight = 1.);
    void endEvent(Float_t _weight = 1);
    void debug(UInt_t _e, UInt_t _a);

   private:

    Double_t getPrescaleFactor(UInt_t _e = INT_MAX);

    void fullExecutionInformation(UInt_t _e, UInt_t _f, Float_t _weight);
    Float_t getElapsedTime(UInt_t _e, UInt_t _f);

    static Float_t s_eventTimeExecute; //!< Static tabulation of execute time for all algorithms in an event.
    Float_t m_firstAlgStartTime; //!< Used to keep track of the earliest (first) execution of this counter in any one event
    Float_t m_firstAlgTime; //<! Used to buffer the execution time for the first call in an event
    //Float_t m_firstAlgTimeWeight; //<! DEFUNCT Used to buffer the weight of m_firstAlgTime

    //Float_t m_eventWeight; //<! DEFUNCT What weight to use for the active events variable
    //std::set< std::string > m_chainsSeen; //<! DEFUNCT List of chains which have used this in a given event

  }; //class CounterAlgorithm

} // namespace TrigCostRootAnalysis

#endif //TrigCostRootAnalysis_CounterAlgorithm_H
