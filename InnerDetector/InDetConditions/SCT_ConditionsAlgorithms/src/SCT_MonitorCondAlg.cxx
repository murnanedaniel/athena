/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "SCT_MonitorCondAlg.h"

#include "InDetIdentifier/SCT_ID.h"

#include <memory>

SCT_MonitorCondAlg::SCT_MonitorCondAlg(const std::string& name, ISvcLocator* pSvcLocator)
  : ::AthReentrantAlgorithm(name, pSvcLocator)
{
}

StatusCode SCT_MonitorCondAlg::initialize()
{
  ATH_MSG_DEBUG("initialize " << name());

  ATH_CHECK(detStore()->retrieve(m_helper, "SCT_ID"));

  // Read Cond Handle
  ATH_CHECK(m_readKey.initialize());

  // Write Cond Handle
  ATH_CHECK(m_writeKey.initialize());

  return StatusCode::SUCCESS;
}

StatusCode SCT_MonitorCondAlg::execute(const EventContext& ctx) const
{
  ATH_MSG_DEBUG("execute " << name());

  // Write Cond Handle
  SG::WriteCondHandle<SCT_MonitorCondData> writeHandle{m_writeKey, ctx};

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
  ATH_MSG_INFO("Size of CondAttrListCollection readCdo->size()= " << readCdo->size());

  // Add dependency
  writeHandle.addDependency(readHandle);

  // Construct the output Cond Object and fill it in
  std::unique_ptr<SCT_MonitorCondData> writeCdo{std::make_unique<SCT_MonitorCondData>()};

  // Fill Write Cond Handle
  static const unsigned int defectListIndex{7};
  CondAttrListCollection::const_iterator iter{readCdo->begin()};
  CondAttrListCollection::const_iterator last{readCdo->end()};
  for (; iter!=last; ++iter) {
    const coral::AttributeList& list{iter->second};
    if (list.size()>defectListIndex) {
      const Identifier moduleId{m_helper->module_id(Identifier{iter->first})};
      const IdentifierHash moduleHash{m_helper->wafer_hash(moduleId)};
      writeCdo->insert(moduleHash, list[defectListIndex].data<std::string>());
    }
  }

  // Record validity of the output cond obbject
  if (writeHandle.record(std::move(writeCdo)).isFailure()) {
    ATH_MSG_FATAL("Could not record SCT_MonitorCondData " << writeHandle.key()
                  << " with EventRange " << writeHandle.getRange()
                  << " into Conditions Store");
    return StatusCode::FAILURE;
  }
  ATH_MSG_INFO("recorded new CDO " << writeHandle.key() << " with range " << writeHandle.getRange() << " into Conditions Store");

  return StatusCode::SUCCESS;
}

StatusCode SCT_MonitorCondAlg::finalize()
{
  ATH_MSG_DEBUG("finalize " << name());
  return StatusCode::SUCCESS;
}
