/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#include "SCT_MajorityCondAlg.h"

#include "Identifier/IdentifierHash.h"
#include "SCT_Cabling/SCT_OnlineId.h"
#include "EventInfo/EventID.h"

#include "GaudiKernel/EventIDRange.h"

#include "SCT_ConditionsServices/SCT_ConditionsParameters.h"
using namespace SCT_ConditionsServices;

SCT_MajorityCondAlg::SCT_MajorityCondAlg(const std::string& name, ISvcLocator* pSvcLocator)
  : ::AthAlgorithm(name, pSvcLocator)
  , m_readKey{"/SCT/DCS/MAJ"}
  , m_writeKey{"SCT_MajorityCondData"}
  , m_condSvc{"CondSvc", name}
  , m_cablingSvc{"SCT_CablingSvc", name}
{
  declareProperty("ReadKey", m_readKey, "Key of input (raw) conditions folder");
  declareProperty("WriteKey", m_writeKey, "Key of output (derived) conditions folder");
}

SCT_MajorityCondAlg::~SCT_MajorityCondAlg()
{
}

StatusCode SCT_MajorityCondAlg::initialize()
{
  ATH_MSG_DEBUG("initialize " << name());

  // CondSvc
  ATH_CHECK( m_condSvc.retrieve() );
  // SCT cabling service
  ATH_CHECK( m_cablingSvc.retrieve() );

  // Read Cond Handle
  ATH_CHECK( m_readKey.initialize() );

  // Write Cond Handle
  ATH_CHECK( m_writeKey.initialize() );
  // Register write handle
  if(m_condSvc->regHandle(this, m_writeKey, m_writeKey.dbKey()).isFailure()) {
    ATH_MSG_FATAL("unable to register WriteCondHandle " << m_writeKey.fullKey() << " with CondSvc");
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

StatusCode SCT_MajorityCondAlg::execute()
{
  ATH_MSG_DEBUG("execute " << name());

  // Write Cond Handle
  SG::WriteCondHandle<SCT_MajorityCondData> writeHandle{m_writeKey};

  // Do we have a valid Write Cond Handle for current time?
  if(writeHandle.isValid()) {
    // in theory this should never be called in MT
    writeHandle.updateStore();
    ATH_MSG_DEBUG("CondHandle " << writeHandle.fullKey() << " is already valid."
                  << ". In theory this should not be called, but may happen"
                  << " if multiple concurrent events are being processed out of order."
                  << " Forcing update of Store contents");
    return StatusCode::SUCCESS; 
  }

  // Read Cond Handle 
  SG::ReadCondHandle<CondAttrListCollection> readHandle{m_readKey};
  const CondAttrListCollection* readCdo{*readHandle}; 
  if(readCdo==nullptr) {
    ATH_MSG_FATAL("Null pointer to the read conditions object");
    return StatusCode::FAILURE;
  }
  ATH_MSG_INFO("Size of CondAttrListCollection readCdo->size()= " << readCdo->size());

  // Define validity of the output cond object
  EventIDRange rangeW;
  if(!readHandle.range(rangeW)) {
    ATH_MSG_FATAL("Failed to retrieve validity range for " << readHandle.key());
    return StatusCode::FAILURE;
  }

  int numFilled{0};

  // Construct the output Cond Object and fill it in
  SCT_MajorityCondData* writeCdo{new SCT_MajorityCondData()};

  CondAttrListCollection::const_iterator majItr{m_dataMajority->begin()};
  CondAttrListCollection::const_iterator majEnd{m_dataMajority->end()};
  for (;majItr != majEnd; ++majItr) {
    // A CondAttrListCollection is a map of ChanNum and AttributeList
    CondAttrListCollection::ChanNum channelNumber{(*majItr).first};
    CondAttrListCollection::AttributeList payload{(*majItr).second};
    // Possible components
    if ((channelNumber == OVERALL) or (channelNumber == BARREL) or (channelNumber == ECA) or (channelNumber == ECC)) {
      // Reset default to true at callback
      bool majorityState{true};

      // Majority state
      if (not payload[INDEX_MajorityState].isNull()) {
	ATH_MSG_DEBUG("Majority state for " << channelNumber << " = " << payload[INDEX_MajorityState].data<int>());
	majorityState = (payload[INDEX_MajorityState].data<int>() == HighAndLowVoltageOK);
      }
      writeCdo->setMajorityState(channelNumber, majorityState);

      // HV fraction in majority state (>50% by definition) IF majority state is HV and LV on
      float hvFraction{1.};
      if (majorityState and (not payload[INDEX_HVfraction].isNull())) {
	ATH_MSG_DEBUG("Majority HV fraction for " << channelNumber << " = " << payload[INDEX_HVfraction].data<float>());
	hvFraction = payload[INDEX_HVfraction].data<float>();
	numFilled++;
      }
      writeCdo->setHVFraction(channelNumber, hvFraction);

    } else {
      ATH_MSG_WARNING("Unexpected channel number " << channelNumber);
    }
  }

  // Has data been filled?
  // Four regions (OVERALL, BARREL, ECA, ECC) are needed.
  writeCdo->setFilled(numFilled==N_REGIONS);

  // Record the out output Cond Object
  if(writeHandle.record(rangeW, writeCdo).isFailure()) {
    ATH_MSG_FATAL("Could not record SCT_MajorityCondData " << writeHandle.key() 
		  << " with EventRange " << rangeW
		  << " into Conditions Store");
    return StatusCode::FAILURE;
  }
  ATH_MSG_INFO("recorded new CDO " << writeHandle.key() << " with range " << rangeW << " into Conditions Store");

  return StatusCode::SUCCESS;
}

StatusCode SCT_MajorityCondAlg::finalize()
{
  ATH_MSG_DEBUG("finalize " << name());
  return StatusCode::SUCCESS;
}
