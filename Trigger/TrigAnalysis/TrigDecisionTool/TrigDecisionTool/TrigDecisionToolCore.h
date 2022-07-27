/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TrigDecision_TrigDecisionToolCore_h
#define TrigDecision_TrigDecisionToolCore_h
/**********************************************************************************
 * @Project: TrigDecisionTool
 * @Package: TrigDecisionTool
 * @class  : TrigDecisionTool
 *
 * @brief main tool
 *    This class defines basic functionality of TDT irrespectively 
 *      of the environment in which it is used (ARA, Athena).
 * 
 *
 * @author Michael Begel  <michael.begel@cern.ch>  - Brookhaven National Laboratory
 * @author Tomasz Bold    <Tomasz.Bold@cern.ch>    - UC Irvine - AGH Krakow
 * @author Joerg Stelzer  <Joerg.Stelzer@cern.ch>  - DESY
 *
 ***********************************************************************************/
#include "AsgMessaging/StatusCode.h"
#include "AsgTools/SlotSpecificObj.h"
#include "TrigDecisionTool/ChainGroupFunctions.h"
#include "TrigDecisionTool/Conditions.h"
#include "TrigDecisionTool/ChainGroup.h"
#include "TrigDecisionTool/CacheGlobalMemory.h"
#include "TrigDecisionTool/ConfigurationAccess.h"
#include "TrigDecisionTool/DecisionAccess.h"
#include "TrigDecisionTool/ExpertMethods.h"

namespace Trig {
  
  class TrigDecisionToolCore : public ChainGroupFunctions, 
                               public ConfigurationAccess, 
                               public DecisionAccess, 
                               public virtual Logger
  {
  public:
    // constructors, destructor
    TrigDecisionToolCore();
    virtual ~TrigDecisionToolCore();
    
    // initialize routine as required for an Algorithm
    virtual StatusCode initialize();
    virtual StatusCode finalize();

    using ChainGroupFunctions::getChainGroup;
    using ConfigurationAccess::getListOfTriggers;
    using ConfigurationAccess::getListOfStreams;
    using ConfigurationAccess::getListOfGroups;
    using ConfigurationAccess::getListOfTriggerElements;
    using ConfigurationAccess::getPrescale;
    using DecisionAccess::isPassed;
    using DecisionAccess::isPassedBits;
    using DecisionAccess::features;
    using DecisionAccess::ancestor;

    const Trig::ExpertMethods& ExperimentalAndExpertMethods() const { return m_expertMethods; }

  protected:
    virtual Trig::CacheGlobalMemory* cgm();
    virtual const Trig::CacheGlobalMemory* cgm() const;

    
  private:

    SG::SlotSpecificObj<Trig::CacheGlobalMemory> m_cacheGlobalMemory;

    Trig::ExpertMethods m_expertMethods;
    TrigDecisionToolCore (const TrigDecisionToolCore&);
    TrigDecisionToolCore& operator= (const TrigDecisionToolCore&);

  protected:
    // Set by TrigDecisionTool.
    HLT::TrigNavStructure* m_navigation = nullptr;
  };
  
} // End of namespace

#endif
