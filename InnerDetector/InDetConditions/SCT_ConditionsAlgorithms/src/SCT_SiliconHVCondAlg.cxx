/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "SCT_SiliconHVCondAlg.h"

#include "Identifier/IdentifierHash.h"
#include "InDetIdentifier/SCT_ID.h"

#include <memory>

SCT_SiliconHVCondAlg::SCT_SiliconHVCondAlg(const std::string& name, ISvcLocator* pSvcLocator)
  : ::AthReentrantAlgorithm(name, pSvcLocator)
{
}

StatusCode SCT_SiliconHVCondAlg::initialize() {
  ATH_MSG_DEBUG("initialize " << name());

  // SCT DCS tool
  ATH_CHECK(m_sctDCSTool.retrieve());
  // SCT ID helper
  ATH_CHECK(detStore()->retrieve(m_pHelper, "SCT_ID"));

  // Read Cond Handles
  ATH_CHECK(m_readKeyState.initialize(m_useState));
  ATH_CHECK(m_readKeyHV.initialize());
  // Write Cond Handle
  ATH_CHECK(m_writeKey.initialize());

  return StatusCode::SUCCESS;
}

StatusCode SCT_SiliconHVCondAlg::execute(const EventContext& ctx) const {
  ATH_MSG_DEBUG("execute " << name());

  // Write Cond Handle
  SG::WriteCondHandle<SCT_DCSFloatCondData> writeHandle{m_writeKey, ctx};
  // Do we have a valid Write Cond Handle for current time?
  if (writeHandle.isValid()) {
    ATH_MSG_DEBUG("CondHandle " << writeHandle.fullKey() << " is already valid."
                  << " In theory this should not be called, but may happen"
                  << " if multiple concurrent events are being processed out of order.");
    return StatusCode::SUCCESS; 
  }

  // Read Cond Handle (HV)
  SG::ReadCondHandle<SCT_DCSFloatCondData> readHandleHV{m_readKeyHV, ctx};
  const SCT_DCSFloatCondData* readCdoHV{*readHandleHV};
  if (readCdoHV==nullptr) {
    ATH_MSG_FATAL("Null pointer to the read conditions object");
    return StatusCode::FAILURE;
  }
  writeHandle.addDependency(readHandleHV);
  ATH_MSG_INFO("Input is " << readHandleHV.fullKey() << " with the range of " << readHandleHV.getRange());

  if (m_useState) {
    // Read Cond Handle (state)
    SG::ReadCondHandle<SCT_DCSStatCondData> readHandleState{m_readKeyState, ctx};
    const SCT_DCSStatCondData* readCdoState{*readHandleState};
    if (readCdoState==nullptr) {
      ATH_MSG_FATAL("Null pointer to the read conditions object");
      return StatusCode::FAILURE;
    }
    writeHandle.addDependency(readHandleState);
    ATH_MSG_INFO("Input is " << readHandleState.fullKey() << " with the range of " << readHandleState.getRange());
  }
  
  // Construct the output Cond Object and fill it in
  std::unique_ptr<SCT_DCSFloatCondData> writeCdo{std::make_unique<SCT_DCSFloatCondData>()};
  const SCT_ID::size_type wafer_hash_max{m_pHelper->wafer_hash_max()};
  for (SCT_ID::size_type hash{0}; hash<wafer_hash_max; hash++) {
    writeCdo->setValue(hash, m_sctDCSTool->modHV(IdentifierHash(hash), ctx));
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

StatusCode SCT_SiliconHVCondAlg::finalize()
{
  ATH_MSG_DEBUG("finalize " << name());
  return StatusCode::SUCCESS;
}
