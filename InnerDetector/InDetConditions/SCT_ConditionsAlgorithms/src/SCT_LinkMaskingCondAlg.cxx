/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "SCT_LinkMaskingCondAlg.h"

#include <memory>

SCT_LinkMaskingCondAlg::SCT_LinkMaskingCondAlg(const std::string& name, ISvcLocator* pSvcLocator)
  : ::AthReentrantAlgorithm(name, pSvcLocator)
{
}

StatusCode SCT_LinkMaskingCondAlg::initialize() {
  ATH_MSG_DEBUG("initialize " << name());

  // Read Cond Handle
  ATH_CHECK(m_readKey.initialize());
  // Write Cond Handles
  ATH_CHECK(m_writeKey.initialize());

  return StatusCode::SUCCESS;
}

StatusCode SCT_LinkMaskingCondAlg::execute(const EventContext& ctx) const {
  ATH_MSG_DEBUG("execute " << name());

  // Write Cond Handle
  SG::WriteCondHandle<SCT_ModuleVetoCondData> writeHandle{m_writeKey, ctx};
  // Do we have a valid Write Cond Handle for current time?
  if (writeHandle.isValid()) {
    ATH_MSG_DEBUG("CondHandle " << writeHandle.fullKey() << " is already valid."
                  << ". In theory this should not be called, but may happen"
                  << " if multiple concurrent events are being processed out of order.");
    return StatusCode::SUCCESS; 
  }

  // Read Cond Handle
  SG::ReadCondHandle<CondAttrListCollection> readHandle{m_readKey, ctx};
  const CondAttrListCollection* readCdo{*readHandle}; 
  if (readCdo==nullptr) {
    ATH_MSG_FATAL("Null pointer to the read conditions object");
    return StatusCode::FAILURE;
  }
  // Add dependency
  writeHandle.addDependency(readHandle);
  ATH_MSG_INFO("Size of CondAttrListCollection " << readHandle.fullKey() << " readCdo->size()= " << readCdo->size());
  ATH_MSG_INFO("Range of input is " << readHandle.getRange());
  
  // Construct the output Cond Object and fill it in
  std::unique_ptr<SCT_ModuleVetoCondData> writeCdo{std::make_unique<SCT_ModuleVetoCondData>()};

  // Read bad wafer info
  CondAttrListCollection::const_iterator linkItr{readCdo->begin()};
  CondAttrListCollection::const_iterator linkEnd{readCdo->end()};
  for (;linkItr != linkEnd; ++linkItr) {
    //A CondAttrListCollection is a map of ChanNum and AttributeList
    Identifier waferId{(*linkItr).first};
    const CondAttrListCollection::AttributeList &payload{(*linkItr).second};
    bool lastProbedState{payload[0].data<bool>()};
    if (not lastProbedState) writeCdo->setBadWaferId(waferId);
    ATH_MSG_INFO("LINK " << waferId << " (" << waferId.get_identifier32().get_compact() << " in 32 bit): " << lastProbedState);
  }
  if (writeCdo->size()>0) writeCdo->setFilled();

  // Record the output cond object
  if (writeHandle.record(std::move(writeCdo)).isFailure()) {
    ATH_MSG_FATAL("Could not record SCT_ModuleVetoCondData " << writeHandle.key()
                  << " with EventRange " << writeHandle.getRange()
                  << " into Conditions Store");
    return StatusCode::FAILURE;
  }
  ATH_MSG_INFO("recorded new CDO " << writeHandle.key() << " with range " << writeHandle.getRange() << " into Conditions Store");

  return StatusCode::SUCCESS;
}

StatusCode SCT_LinkMaskingCondAlg::finalize() {
  ATH_MSG_DEBUG("finalize " << name());
  return StatusCode::SUCCESS;
}
