/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TRIGHLTJETHYPO_JETGROUPPRODUCT_H
#define TRIGHLTJETHYPO_JETGROUPPRODUCT_H

#include "./IJetGroupProduct.h"
#include "./ProductGen.h"
#include "./JetStreamer.h"
#include "./DebugInfoCollector.h"
#include "TrigHLTJetHypo/TrigHLTJetHypoUtils/HypoJetDefs.h"
#include <vector>

using CondInd2JetGroupsInds = std::map<int, std::vector<std::size_t>>;
using JetGroupInd2ElemInds = std::map<int, std::vector<std::size_t>>;

typedef std::unique_ptr<ITrigJetHypoInfoCollector> Collector;

class JetGroupProduct: public IJetGroupProduct{
  /*
   * Iterate through the combinations of jet groups.
   * The jet groups are those that satisfied a set up siblings.
   * Their parent is tested against all the combioinations of jet groups
   * that satisfy the siblings.
   * 
   * eg for a vector of siblings [s1, s2]
   * satisfiedBy[s1] is a  vector<vector<size_t>> say <[jg1, jg2], [jg3, jg4]>
   * satisfiedBy[s2] is a  vector<vector<size_t>> say <[jg5, jg6]>
   * the products are then [jg1, jg2, jg5, jg6], [jg3, jg4, jg5, jg6]
   * jg1 is an key in a map that has a value of an input jet group (typically
   * containing a single jet.
   */
 public:
  JetGroupProduct(const std::vector<std::size_t>& siblings,
		  const CondInd2JetGroupsInds& satisfiedBy,
		  const std::vector<std::size_t>& condMult,
		  const JetGroupInd2ElemInds&);

  virtual std::vector<std::size_t> next(const Collector&) override;
  virtual bool valid() const override;
  
 private:

  JetGroupInd2ElemInds m_jg2elemjgs;

  std::vector<bool>  m_jetMask;
  std::size_t  m_jetEnd{0};
  // ProductGen m_productGen;
  std::vector<std::vector<std::size_t>> m_seenIndices;
  std::unique_ptr<JetStreamer> m_jetstreamer{nullptr};
  bool m_valid{false};

  void init(const std::vector<std::size_t>& siblings,
	    const CondInd2JetGroupsInds& satisfiedBy,
	    const std::vector<std::size_t>& condMult
	    );

};

#endif
