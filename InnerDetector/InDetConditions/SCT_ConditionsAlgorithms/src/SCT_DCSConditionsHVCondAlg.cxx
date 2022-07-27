/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "SCT_DCSConditionsHVCondAlg.h"

#include "Identifier/IdentifierHash.h"
#include "AthenaKernel/IOVInfiniteRange.h"

#include <memory>

SCT_DCSConditionsHVCondAlg::SCT_DCSConditionsHVCondAlg(const std::string& name, ISvcLocator* pSvcLocator)
  : ::AthReentrantAlgorithm(name, pSvcLocator)
{
}

StatusCode SCT_DCSConditionsHVCondAlg::initialize() {
  ATH_MSG_DEBUG("initialize " << name());

  // Read Cond Handle
  ATH_CHECK(m_readKey.initialize(m_returnHVTemp));
  // Write Cond Handle
  ATH_CHECK(m_writeKey.initialize(m_returnHVTemp));

  return StatusCode::SUCCESS;
}

StatusCode SCT_DCSConditionsHVCondAlg::execute(const EventContext& ctx) const {
  ATH_MSG_DEBUG("execute " << name());

  if (not m_returnHVTemp) {
    return StatusCode::SUCCESS;
  }

  // Write Cond Handle
  SG::WriteCondHandle<SCT_DCSFloatCondData> writeHandle{m_writeKey, ctx};
  // Do we have a valid Write Cond Handle for current time?
  if (writeHandle.isValid()) {
    ATH_MSG_DEBUG("CondHandle " << writeHandle.fullKey() << " is already valid."
                  << ". In theory this should not be called, but may happen"
                  << " if multiple concurrent events are being processed out of order.");
    return StatusCode::SUCCESS; 
  }

  writeHandle.addDependency(IOVInfiniteRange::infiniteMixed());

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
  std::unique_ptr<SCT_DCSFloatCondData> writeCdo{std::make_unique<SCT_DCSFloatCondData>()};

  // Read temperature info
  std::string param{"HVCHVOLT_RECV"};
  CondAttrListCollection::const_iterator attrList{readCdo->begin()};
  CondAttrListCollection::const_iterator end{readCdo->end()};
  // CondAttrListCollection doesn't support C++11 type loops, no generic 'begin'
  for (; attrList!=end; ++attrList) {
    // A CondAttrListCollection is a map of ChanNum and AttributeList
    CondAttrListCollection::ChanNum channelNumber{attrList->first};
    const CondAttrListCollection::AttributeList &payload{attrList->second};
    if (payload.exists(param) and not payload[param].isNull()) {
      float val{payload[param].data<float>()};
      writeCdo->setValue(channelNumber, val);
    } else {
      ATH_MSG_WARNING(param << " does not exist for ChanNum " << channelNumber);
    }
  }

  // Record the output cond object
  if (writeHandle.record(std::move(writeCdo)).isFailure()) {
    ATH_MSG_FATAL("Could not record SCT_DCSFloatCondData " << writeHandle.key() 
                  << " with EventRange " << writeHandle.getRange()
                  << " into Conditions Store");
    return StatusCode::FAILURE;
  }
  ATH_MSG_INFO("recorded new CDO " << writeHandle.key() << " with range " << writeHandle.getRange() << " into Conditions Store");

  return StatusCode::SUCCESS;
}

StatusCode SCT_DCSConditionsHVCondAlg::finalize()
{
  ATH_MSG_DEBUG("finalize " << name());
  return StatusCode::SUCCESS;
}
