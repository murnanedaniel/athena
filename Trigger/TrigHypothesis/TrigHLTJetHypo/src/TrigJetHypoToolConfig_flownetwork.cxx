/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/


#include "TrigJetHypoToolConfig_flownetwork.h"

#include "GaudiKernel/StatusCode.h"

#include "./conditionsFactoryMT.h"

#include "TrigHLTJetHypo/TrigHLTJetHypoUtils/SingleJetGrouper.h"
#include "TrigHLTJetHypo/TrigHLTJetHypoUtils/xAODJetAsIJetFactory.h"
#include "TrigHLTJetHypo/TrigHLTJetHypoUtils/groupsMatcherFactory.h"
#include "TrigHLTJetHypo/TrigHLTJetHypoUtils/CleanerFactory.h"
#include "TrigHLTJetHypo/TrigHLTJetHypoUtils/TrigHLTJetHypoHelper2.h"
#include "./groupsMatcherFactoryMT.h"

#include "DecisionHandling/TrigCompositeUtils.h"

using TrigCompositeUtils::DecisionID;
using TrigCompositeUtils::Decision;
using TrigCompositeUtils::DecisionContainer;

TrigJetHypoToolConfig_flownetwork::TrigJetHypoToolConfig_flownetwork(const std::string& type,
                                                 const std::string& name,
                                                 const IInterface* parent) :
  base_class(type, name, parent){

}


TrigJetHypoToolConfig_flownetwork::~TrigJetHypoToolConfig_flownetwork(){
}

StatusCode TrigJetHypoToolConfig_flownetwork::initialize() {

  // gymanastics as cannot pass vecor<vecotr<int>> as a Gaudi::Property
  if(m_sharedNodesVec.empty()){return StatusCode::FAILURE;}

  if(m_sharedNodesVec.size() == 1){
    if(m_sharedNodesVec[0] >= 0){
      m_sharedNodes.push_back(m_sharedNodesVec);
      return StatusCode::SUCCESS;
    } else  {
      return StatusCode::FAILURE;
    }
  }
  
  auto be = m_sharedNodesVec.begin();
  auto en = m_sharedNodesVec.end();
  auto it = std::find(be, en, -1);
  while(true){
    m_sharedNodes.push_back(std::vector<int>(be, it));
    if(it >= en - 1){break;}
    be = it + 1;
    it = std::find(be, it, -1);
  };
    
  return StatusCode::SUCCESS;
}




std::optional<ConditionsMT>
TrigJetHypoToolConfig_flownetwork::getConditions() const {

  ConditionsMT conditions;

  // collect the Conditions objects from the various sources
  // return an invalid optional if any src signals a problem
  for(const auto& cm : m_conditionMakers){
    conditions.push_back(cm->getCondition());
  }
    
  return std::make_optional<ConditionsMT>(std::move(conditions));
}


// following function not used for treeless hypos
std::size_t
TrigJetHypoToolConfig_flownetwork::requiresNJets() const {
  return 0;
}

 
std::unique_ptr<IJetGrouper>
TrigJetHypoToolConfig_flownetwork::getJetGrouper() const {
  return std::make_unique<SingleJetGrouper>();
}

std::unique_ptr<IGroupsMatcherMT>
TrigJetHypoToolConfig_flownetwork::getMatcher () const {

  auto opt_conds = getConditions();

  if(!opt_conds.has_value()){
    return std::unique_ptr<IGroupsMatcherMT>(nullptr);
  }

  return groupsMatcherFactoryMT_Unified(std::move(*opt_conds),
					m_treeVec,
					m_sharedNodes);
}

StatusCode TrigJetHypoToolConfig_flownetwork::checkVals() const {
  return StatusCode::SUCCESS;
}

std::vector<std::shared_ptr<ICleaner>> 
TrigJetHypoToolConfig_flownetwork::getCleaners() const {
  std::vector<std::shared_ptr<ICleaner>> v;
  return v;
}


