/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TRIG_IDECISIONUNPACKER_H
#define TRIG_IDECISIONUNPACKER_H

#include "AsgTools/StatusCode.h"
#include <map>
#include <string>

namespace HLT {
  class TrigNavStructure;
  class Chain;
}

namespace LVL1CTP{
  class Lvl1Item;
}

namespace Trig{
  class IDecisionUnpacker{
  public:
    typedef unsigned CTPID;
    typedef unsigned CHAIN_COUNTER;
    virtual ~IDecisionUnpacker();
    virtual StatusCode unpackDecision(std::map<std::string, const LVL1CTP::Lvl1Item*>&,
				      std::map<CTPID, LVL1CTP::Lvl1Item*>& itemsCache,
				      std::map<std::string, const HLT::Chain*>&,
				      std::map<CHAIN_COUNTER, HLT::Chain*>&,
				      std::map<std::string, const HLT::Chain*>&,
				      std::map<CHAIN_COUNTER, HLT::Chain*>&,
				      char&,
				      bool
				      ) = 0;
    virtual StatusCode unpackNavigation(HLT::TrigNavStructure*) = 0;
    virtual bool assert_handle() = 0;
    virtual void validate_handle() = 0;
    virtual void invalidate_handle() = 0;
  };
} //end of Trig namespace

#endif
