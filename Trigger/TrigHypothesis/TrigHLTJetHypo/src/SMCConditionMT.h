/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TRIGHLTJETHYPO_SMCCONDITIONMT_H
#define TRIGHLTJETHYPO_SMCCONDITIONMT_H

/********************************************************************
 *
 * NAME:     SMCConditionMT.h (mass of jet consition
 * PACKAGE:  Trigger/TrigHypothesis/TrigHLTJetHypo
 *
 * AUTHOR:   P. Sherwood
 *********************************************************************/

#include <string>
#include "./IConditionMT.h"

namespace HypoJet{
  class IJet;
}

class ITrigJetHypoInfoCollector;

class SMCConditionMT: public IConditionMT{
 public:
  SMCConditionMT(double massMin,
		 double massMax);

  bool isSatisfied(const HypoJetVector&,
                   const std::unique_ptr<ITrigJetHypoInfoCollector>&) const override;

  virtual unsigned int capacity() const override{return s_capacity;}
  std::string toString() const override;
  
 private:

  double m_min;
  double m_max;
  bool isSatisfied(const pHypoJet&,
                   const std::unique_ptr<ITrigJetHypoInfoCollector>&) const;
  
  const static  unsigned int s_capacity{1};
  
};

#endif
