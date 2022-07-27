/*
 * Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration.
 */
/**
 * @file CaloRec/ToolConstantsCondAlg.h
 * @author scott snyder <snyder@bnl.gov>
 * @date Jun, 2020
 * @brief Convert from COOL inline data / POOL data to ToolConstants.
 */


#include "ToolConstantsCondAlg.h"
#include "AthenaKernel/IOVInfiniteRange.h"

/**
 * @brief Algorithm initialize method.
 */
StatusCode ToolConstantsCondAlg::initialize()
{
  if (!m_coolFolderKey.empty()) {
    if (!m_detStoreKey.empty()) {
      ATH_MSG_ERROR( "Configuration error: both COOL folder and det store key specified." );
      return StatusCode::FAILURE;
    }
    ATH_CHECK( m_blobTool.retrieve() );
    ATH_CHECK( m_coolFolderKey.initialize() );
  }
  else if (m_detStoreKey.empty()) {
    ATH_MSG_ERROR( "Configuration error: neither COOL folder nor det store key specified." );
    return StatusCode::FAILURE;
  }

  ATH_CHECK( m_toolConstantsKey.initialize() );

  return StatusCode::SUCCESS;
}


/**
 * @brief Algorithm execute method.
 * @param ctx Event context.
 */
StatusCode ToolConstantsCondAlg::execute (const EventContext& ctx) const
{
  SG::WriteCondHandle<CaloRec::ToolConstants> toolConstants
    (m_toolConstantsKey, ctx);

  auto tc = std::make_unique<CaloRec::ToolConstants>();

  if (!m_coolFolderKey.empty()) {
    SG::ReadCondHandle<CondAttrListCollection> coolFolder(m_coolFolderKey, ctx);

    const std::string key = m_toolConstantsKey.key();

    const unsigned chNbr = m_blobTool->nameToChannelNumber (key);
    // Check that this channel actually exits
    const std::string& chanName = coolFolder->chanName (chNbr);
    if (!chanName.empty() && key!=chanName) {
      ATH_MSG_ERROR( "Channel name does not match! Expected " << key << " found " << chanName );
      return StatusCode::FAILURE;
    }
    else {
      ATH_MSG_DEBUG( "Found channel number " << chNbr << " named " << key );
    }

    const coral::AttributeList& attrList = coolFolder->attributeList (chNbr);
    ATH_CHECK( m_blobTool->AttrListToToolConstants (attrList, *tc) );

    toolConstants.addDependency (coolFolder);
  }
  else if (!m_detStoreKey.empty()) {
    const CaloRec::ToolConstants* tc_in = nullptr;
    ATH_CHECK( detStore()->retrieve (tc_in, m_detStoreKey) );
    *tc = *tc_in;

    const EventIDRange fullRange=IOVInfiniteRange::infiniteRunLB();
    toolConstants.addDependency (fullRange);
  }
  else {
    ATH_MSG_ERROR( "Bad configuration." );
    return StatusCode::FAILURE;
  }

  ATH_CHECK( toolConstants.record (std::move (tc)) );
  return StatusCode::SUCCESS;
}
