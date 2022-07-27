/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "SCT_MajorityCondAlg.h"

#include "SCT_ConditionsData/SCT_ConditionsParameters.h"

#include <memory>

SCT_MajorityCondAlg::SCT_MajorityCondAlg(const std::string& name, ISvcLocator* pSvcLocator)
  : ::AthReentrantAlgorithm(name, pSvcLocator)
{
}

StatusCode SCT_MajorityCondAlg::initialize()
{
  ATH_MSG_DEBUG("initialize " << name());

  // Read Cond Handle
  ATH_CHECK(m_readKey.initialize());

  // Write Cond Handle
  ATH_CHECK(m_writeKey.initialize());

  return StatusCode::SUCCESS;
}

StatusCode SCT_MajorityCondAlg::execute(const EventContext& ctx) const
{
  ATH_MSG_DEBUG("execute " << name());

  // Write Cond Handle
  SG::WriteCondHandle<SCT_MajorityCondData> writeHandle{m_writeKey, ctx};

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

  int numFilled{0};

  // Construct the output Cond Object and fill it in
  std::unique_ptr<SCT_MajorityCondData> writeCdo{std::make_unique<SCT_MajorityCondData>()};

  CondAttrListCollection::const_iterator majItr{readCdo->begin()};
  CondAttrListCollection::const_iterator majEnd{readCdo->end()};
  for (;majItr != majEnd; ++majItr) {
    // A CondAttrListCollection is a map of ChanNum and AttributeList
    CondAttrListCollection::ChanNum channelNumber{(*majItr).first};
    const CondAttrListCollection::AttributeList &payload{(*majItr).second};
    // Possible components
    if ((channelNumber == SCT_ConditionsData::OVERALL) or 
        (channelNumber == SCT_ConditionsData::BARREL) or 
        (channelNumber == SCT_ConditionsData::ECA) or 
        (channelNumber == SCT_ConditionsData::ECC)) {
      // Reset default to true at callback
      bool majorityState{true};

      // Majority state
      if (not payload[SCT_ConditionsData::INDEX_MajorityState].isNull()) {
	ATH_MSG_DEBUG("Majority state for " << channelNumber << " = " <<
                      payload[SCT_ConditionsData::INDEX_MajorityState].data<int>());
	majorityState = (payload[SCT_ConditionsData::INDEX_MajorityState].data<int>()
                         == SCT_ConditionsData::HighAndLowVoltageOK);
      }
      writeCdo->setMajorityState(channelNumber, majorityState);

      // HV fraction in majority state (>50% by definition) IF majority state is HV and LV on
      float hvFraction{1.};
      if (majorityState and (not payload[SCT_ConditionsData::INDEX_HVfraction].isNull())) {
	ATH_MSG_DEBUG("Majority HV fraction for " << channelNumber << " = " <<
                      payload[SCT_ConditionsData::INDEX_HVfraction].data<float>());
	hvFraction = payload[SCT_ConditionsData::INDEX_HVfraction].data<float>();
	numFilled++;
      }
      writeCdo->setHVFraction(channelNumber, hvFraction);

    } else {
      ATH_MSG_WARNING("Unexpected channel number " << channelNumber);
    }
  }

  // Has data been filled?
  // Four regions (OVERALL, BARREL, ECA, ECC) are needed.
  writeCdo->setFilled(numFilled==SCT_ConditionsData::N_REGIONS);

  // Record the out output Cond Object
  if (writeHandle.record(std::move(writeCdo)).isFailure()) {
    ATH_MSG_FATAL("Could not record SCT_MajorityCondData " << writeHandle.key() 
		  << " with EventRange " << writeHandle.getRange()
		  << " into Conditions Store");
    return StatusCode::FAILURE;
  }
  ATH_MSG_INFO("recorded new CDO " << writeHandle.key() << " with range " << writeHandle.getRange() << " into Conditions Store");

  return StatusCode::SUCCESS;
}

StatusCode SCT_MajorityCondAlg::finalize()
{
  ATH_MSG_DEBUG("finalize " << name());
  return StatusCode::SUCCESS;
}
