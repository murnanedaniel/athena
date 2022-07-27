/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "SCT_ModuleVetoCondAlg.h"

#include <memory>

namespace {
  template <class T>
  std::vector<T>
  string2Vector(const std::string& s) {
    std::vector<T> v;
    std::istringstream inputStream{s};
    std::istream_iterator<T> vecRead{inputStream};
    std::istream_iterator<T> endOfString; //relies on default constructor to produce eof
    std::copy(vecRead,endOfString, std::back_inserter(v)); // DOESN'T ALLOW NON-WHITESPACE DELIMITER !
    return v;
  }
}

SCT_ModuleVetoCondAlg::SCT_ModuleVetoCondAlg(const std::string& name, ISvcLocator* pSvcLocator)
  : ::AthReentrantAlgorithm(name, pSvcLocator)
{
}

StatusCode SCT_ModuleVetoCondAlg::initialize() {
  ATH_MSG_DEBUG("initialize " << name());

  // Read Cond Handle
  ATH_CHECK(m_readKey.initialize());
  // Write Cond Handles
  ATH_CHECK(m_writeKey.initialize());

  return StatusCode::SUCCESS;
}

StatusCode SCT_ModuleVetoCondAlg::execute(const EventContext& ctx) const {
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
  SG::ReadCondHandle<AthenaAttributeList> readHandle{m_readKey, ctx};
  const AthenaAttributeList* readCdo{*readHandle}; 
  if (readCdo==nullptr) {
    ATH_MSG_FATAL("Null pointer to the read conditions object");
    return StatusCode::FAILURE;
  }
  // Add dependency
  writeHandle.addDependency(readHandle);
  ATH_MSG_INFO("Size of AthenaAttributeList " << readHandle.fullKey() << " readCdo->size()= " << readCdo->size());
  ATH_MSG_INFO("Range of input is " << readHandle.getRange());
  
  // Construct the output Cond Object and fill it in
  std::unique_ptr<SCT_ModuleVetoCondData> writeCdo{std::make_unique<SCT_ModuleVetoCondData>()};

  // Read bad wafer info
  const std::string &badModuleString{(*readCdo)["ModuleList"].data<std::string>()};
  std::vector<int> v{string2Vector<int>(badModuleString)};
  int numberInDb{static_cast<int>(v.size())};
  ATH_MSG_INFO(numberInDb << " elements were declared bad in the database.");
  for (const int badWaferId: v) {
    writeCdo->setBadWaferId(Identifier{badWaferId});
  }

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

StatusCode SCT_ModuleVetoCondAlg::finalize() {
  ATH_MSG_DEBUG("finalize " << name());
  return StatusCode::SUCCESS;
}
