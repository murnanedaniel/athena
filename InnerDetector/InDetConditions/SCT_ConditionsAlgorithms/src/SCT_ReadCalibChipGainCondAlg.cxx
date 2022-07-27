/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "SCT_ReadCalibChipGainCondAlg.h"

#include "FillFromStringUtility.h"

#include "Identifier/IdentifierHash.h"
#include "InDetIdentifier/SCT_ID.h"
#include "SCT_ConditionsData/SCT_ConditionsParameters.h"
#include "SCT_ConditionsTools/SCT_ReadCalibChipDefs.h"

#include <memory>

using namespace SCT_ConditionsData;
using namespace SCT_ReadCalibChipDefs;

SCT_ReadCalibChipGainCondAlg::SCT_ReadCalibChipGainCondAlg(const std::string& name, ISvcLocator* pSvcLocator)
  : ::AthReentrantAlgorithm(name, pSvcLocator)
{
}

StatusCode SCT_ReadCalibChipGainCondAlg::initialize() {
  ATH_MSG_DEBUG("initialize " << name());

  // Get SCT helper
  ATH_CHECK(detStore()->retrieve(m_id_sct, "SCT_ID"));

  // Read Cond Handle
  ATH_CHECK(m_readKey.initialize());
  // Write Cond Handle
  ATH_CHECK(m_writeKey.initialize());

  return StatusCode::SUCCESS;
}

StatusCode SCT_ReadCalibChipGainCondAlg::execute(const EventContext& ctx) const {
  ATH_MSG_DEBUG("execute " << name());

  // Write Cond Handle
  SG::WriteCondHandle<SCT_GainCalibData> writeHandle{m_writeKey, ctx};
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
  std::unique_ptr<SCT_GainCalibData> writeCdo{std::make_unique<SCT_GainCalibData>()};

  // Initialization
  const float errVal{std::numeric_limits<float>::quiet_NaN()};
  for (int m{0}; m!=NUMBER_OF_MODULES; ++m) {
    for (int p{0}; p!=N_NPTGAIN; ++p) {
      for (int c{0}; c!=CHIPS_PER_MODULE; ++c) {
        (*writeCdo)[m][p][c]=errVal;
      }
    }
  }

  // loop over collection
  CondAttrListCollection::const_iterator itLoop{readCdo->begin()};
  CondAttrListCollection::const_iterator itLoop_end{readCdo->end()};
  for (; itLoop!=itLoop_end; ++itLoop) {
    CondAttrListCollection::ChanNum chanNum{itLoop->first};
    const coral::AttributeList &anAttrList{itLoop->second};
    // Convert chanNum=offlineID into identifier
    Identifier32 moduleId{chanNum};
    //find the corresponding hash
    const IdentifierHash hashId{m_id_sct->wafer_hash(moduleId)};
    //find the index to the module (hash is for each side), to use as index into array
    const unsigned int moduleIdx{hashId/SIDES_PER_MODULE};
    SCT_ModuleGainCalibData& theseCalibData{(*writeCdo)[moduleIdx]};
    insertNptGainFolderData(theseCalibData, anAttrList);
  }

  // Record the output cond object
  if (writeHandle.record(std::move(writeCdo)).isFailure()) {
    ATH_MSG_FATAL("Could not record SCT_GainCalibData " << writeHandle.key() 
                  << " with EventRange " << writeHandle.getRange()
                  << " into Conditions Store");
    return StatusCode::FAILURE;
  }
  ATH_MSG_INFO("recorded new CDO " << writeHandle.key() << " with range " << writeHandle.getRange() << " into Conditions Store");

  return StatusCode::SUCCESS;
}

StatusCode SCT_ReadCalibChipGainCondAlg::finalize() {
  ATH_MSG_DEBUG("finalize " << name());
  return StatusCode::SUCCESS;
}

void 
SCT_ReadCalibChipGainCondAlg::insertNptGainFolderData(SCT_ModuleGainCalibData& theseCalibData, const coral::AttributeList& folderData) {
  for (int i{0}; i!=N_NPTGAIN; ++i) {
    SCT_ModuleCalibParameter& datavec{theseCalibData[i]};
    const std::string &dbData{((folderData)[nPtGainDbParameterNames[i]]).data<std::string>()};
    fillArrayFromString(dbData, datavec);
  }
}
