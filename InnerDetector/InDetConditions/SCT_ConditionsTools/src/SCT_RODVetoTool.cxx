/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

/**
 * @file SCT_RODVetoTool.cxx
 * implementation file for tool allowing one to declare modules as bad
 * @author daiki.hayakawa@cern.ch
 **/

#include "SCT_RODVetoTool.h"
#include "InDetReadoutGeometry/SiDetectorElementCollection.h"
#include "SCT_DetectorElementStatus.h"

//Athena includes
#include "Identifier/IdentifierHash.h"
#include "InDetIdentifier/SCT_ID.h"
#include "StoreGate/DataHandle.h"
#include "StoreGate/ReadHandle.h"

//STL includes
#include <algorithm>
#include <iterator>
#include <list>
#include <sstream>
#include <vector>

// Constructor
SCT_RODVetoTool::SCT_RODVetoTool(const std::string& type, const std::string& name, const IInterface* parent) :
  base_class(type, name, parent)
{
}

//Initialize
StatusCode 
SCT_RODVetoTool::initialize() {
  ATH_CHECK(detStore()->retrieve(m_pHelper, "SCT_ID"));
  ATH_CHECK(m_badModuleIds.initialize());
 
  return  StatusCode::SUCCESS;
}

//Finalize
StatusCode
SCT_RODVetoTool::finalize() {
  return StatusCode::SUCCESS;
}

bool 
SCT_RODVetoTool::canReportAbout(InDetConditions::Hierarchy h) const {
  return ((h==InDetConditions::DEFAULT) or
          (h==InDetConditions::SCT_SIDE) or
          (h==InDetConditions::SCT_MODULE));
}

bool 
SCT_RODVetoTool::isGood(const Identifier& elementId, const EventContext& ctx, InDetConditions::Hierarchy h) const {
  if (not canReportAbout(h)) return true;
  const IdentifierSet* badIds{getCondData(ctx)};
  if (badIds==nullptr) {
    ATH_MSG_ERROR("IdentifierSet cannot be retrieved in isGood. true is returned.");
    return true;
  }
  const Identifier waferId{m_pHelper->wafer_id(elementId)};
  const Identifier moduleId{m_pHelper->module_id(waferId)};
  bool result{badIds->find(moduleId) == badIds->end()};
  return result;
}

bool 
SCT_RODVetoTool::isGood(const Identifier& elementId, InDetConditions::Hierarchy h) const {
  const EventContext& ctx{Gaudi::Hive::currentContext()};

  return isGood(elementId, ctx, h);
}

bool 
SCT_RODVetoTool::isGood(const IdentifierHash& hashId, const EventContext& ctx) const {
  Identifier elementId{m_pHelper->wafer_id(hashId)};
  Identifier moduleId{m_pHelper->module_id(elementId)};
  return isGood(moduleId, ctx, InDetConditions::SCT_MODULE);
}

bool 
SCT_RODVetoTool::isGood(const IdentifierHash& hashId) const {
  const EventContext& ctx{Gaudi::Hive::currentContext()};

  return isGood(hashId, ctx);
}

void
SCT_RODVetoTool::getDetectorElementStatus(const EventContext& ctx, InDet::SiDetectorElementStatus &element_status, 
                                          SG::WriteCondHandle<InDet::SiDetectorElementStatus>* whandle) const  {
  SG::ReadCondHandle<IdentifierSet> condDataHandle{m_badModuleIds, ctx};
  if (not condDataHandle.isValid() || !condDataHandle.cptr()) {
    ATH_MSG_ERROR("Failed to get " << m_badModuleIds.key());
    return;
  }
  if (whandle) {
    whandle->addDependency (condDataHandle);
  }
  const IdentifierSet* badIds{ condDataHandle.cptr() };
  if (badIds==nullptr) {
    ATH_MSG_ERROR("IdentifierSet cannot be retrieved in isGood. true is returned.");
    return;
  }
  std::vector<bool> &status = element_status.getElementStatus();
  if (status.empty()) {
     status.resize(m_pHelper->wafer_hash_max(),true);
  }
  for (const Identifier &module_id: *badIds) {
     for (int side_i=0; side_i<2; ++side_i) {
        Identifier wafer_id( m_pHelper->wafer_id(  m_pHelper->barrel_ec(module_id),
                                                   m_pHelper->layer_disk(module_id),
                                                   m_pHelper->phi_module(module_id),
                                                   m_pHelper->eta_module(module_id),
                                                   side_i));
        status.at( m_pHelper->wafer_hash(wafer_id) ) = false;
     }
  }
}

const IdentifierSet*
SCT_RODVetoTool::getCondData(const EventContext& ctx) const {
  SG::ReadCondHandle<IdentifierSet> condData{m_badModuleIds, ctx};
  if (not condData.isValid()) {
    ATH_MSG_ERROR("Failed to get " << m_badModuleIds.key());
    return nullptr;
  }
  return condData.cptr();
}
