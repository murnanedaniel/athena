/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "SCT_SensorsCondAlg.h"

#include <memory>

SCT_SensorsCondAlg::SCT_SensorsCondAlg(const std::string& name, ISvcLocator* pSvcLocator)
  : ::AthReentrantAlgorithm(name, pSvcLocator)
{
}

StatusCode SCT_SensorsCondAlg::initialize()
{
  ATH_MSG_DEBUG("initialize " << name());

  // Read Cond Handle
  ATH_CHECK(m_readKey.initialize());

  // Write Cond Handle
  ATH_CHECK(m_writeKey.initialize());

  return StatusCode::SUCCESS;
}

StatusCode SCT_SensorsCondAlg::execute(const EventContext& ctx) const
{
  ATH_MSG_DEBUG("execute " << name());

  // Write Cond Handle
  SG::WriteCondHandle<SCT_SensorsCondData> writeHandle{m_writeKey, ctx};

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
  std::unique_ptr<SCT_SensorsCondData> writeCdo{std::make_unique<SCT_SensorsCondData>()};

  // Fill write conditions data object
  CondAttrListCollection::const_iterator attrList{readCdo->begin()};
  CondAttrListCollection::const_iterator end{readCdo->end()};
  // CondAttrListCollection doesnt support C++11 type loops, no generic 'begin'
  for (; attrList!=end; ++attrList) {
    CondAttrListCollection::ChanNum truncatedSerialNumber{attrList->first};
    SCT_SensorCondData data;
    bool isOK{false};
    isOK |= data.setTruncatedSerialNumber(truncatedSerialNumber);
    isOK |= data.setManufacturer(attrList->second[0].data<std::string>());
    isOK |= data.setDepletionVoltage(SCT_SensorCondData::SENSOR1, attrList->second[1].data<float>());
    isOK |= data.setDepletionVoltage(SCT_SensorCondData::SENSOR2, attrList->second[2].data<float>());
    isOK |= data.setDepletionVoltage(SCT_SensorCondData::SENSOR3, attrList->second[3].data<float>());
    isOK |= data.setDepletionVoltage(SCT_SensorCondData::SENSOR4, attrList->second[4].data<float>());
    isOK |= data.setCrystalOrientation(SCT_SensorCondData::SENSOR1, static_cast<int>(attrList->second[5].data<long long>()));
    isOK |= data.setCrystalOrientation(SCT_SensorCondData::SENSOR2, static_cast<int>(attrList->second[6].data<long long>()));
    isOK |= data.setCrystalOrientation(SCT_SensorCondData::SENSOR3, static_cast<int>(attrList->second[7].data<long long>()));
    isOK |= data.setCrystalOrientation(SCT_SensorCondData::SENSOR4, static_cast<int>(attrList->second[8].data<long long>()));
    if (not isOK) {
      ATH_MSG_WARNING("At least one element of SCT_SensorCondData for truncatedSerialNumber " << truncatedSerialNumber << " was not correctly stored.");
    }
    (*writeCdo)[truncatedSerialNumber] = data;
  }

  // Record write conditions data object
  if (writeHandle.record(std::move(writeCdo)).isFailure()) {
    ATH_MSG_FATAL("Could not record SCT_SensorsCondData " << writeHandle.key() 
                  << " with EventRange " << writeHandle.getRange()
                  << " into Conditions Store");
    return StatusCode::FAILURE;
  }
  ATH_MSG_INFO("recorded new CDO " << writeHandle.key() << " with range " << writeHandle.getRange() << " into Conditions Store");

  return StatusCode::SUCCESS;
}

StatusCode SCT_SensorsCondAlg::finalize()
{
  ATH_MSG_DEBUG("finalize " << name());
  return StatusCode::SUCCESS;
}
