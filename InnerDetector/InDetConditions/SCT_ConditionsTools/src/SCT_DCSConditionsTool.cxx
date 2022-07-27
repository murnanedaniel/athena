/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

// New SCT_DCSConditions Tool, based on existing tool in SCT_ConditionsAlgs
// A. R-Veronneau 26/02/08

#include "SCT_DCSConditionsTool.h"
#include "InDetIdentifier/SCT_ID.h"
#include "SCT_DetectorElementStatus.h"
#include "InDetReadoutGeometry/SiDetectorElementCollection.h"

#include "GaudiKernel/ThreadLocalContext.h"

#include "SCT_IdConverter.h"
using SCT_ConditionsTools::castId;

const Identifier SCT_DCSConditionsTool::s_invalidId;
const float SCT_DCSConditionsTool::s_defaultHV{-30.}; 
const float SCT_DCSConditionsTool::s_defaultTemperature{-40.};

SCT_DCSConditionsTool::SCT_DCSConditionsTool(const std::string& type, const std::string& name, const IInterface* parent) :
  base_class(type, name, parent)
{ 
}

StatusCode SCT_DCSConditionsTool::initialize() {
  if (detStore()->retrieve(m_pHelper,"SCT_ID").isFailure()) {
    ATH_MSG_ERROR("SCT helper failed to retrieve ");
    return StatusCode::FAILURE;
  }
  
  // Read Cond Handle Keys
  ATH_CHECK(m_condKeyHV.initialize(m_returnHVTemp));
  ATH_CHECK(m_condKeyTemp0.initialize(m_returnHVTemp));

  const bool useState = (m_readAllDBFolders and m_returnHVTemp) or (not m_readAllDBFolders and not m_returnHVTemp);
  ATH_CHECK(m_condKeyState.initialize(useState));

  return StatusCode::SUCCESS;
}

StatusCode SCT_DCSConditionsTool::finalize() { 
  return StatusCode::SUCCESS;
}

//Can report about the module as a whole or the wafer
bool SCT_DCSConditionsTool::canReportAbout(InDetConditions::Hierarchy h) const {
  return (h==InDetConditions::SCT_MODULE or h==InDetConditions::SCT_SIDE or h==InDetConditions::SCT_STRIP);
}   

//returns the module ID (int), or returns -1 if not able to report
Identifier SCT_DCSConditionsTool::getModuleID(const Identifier& elementId, InDetConditions::Hierarchy h) const {
  if (not canReportAbout(h)) return s_invalidId;  

  Identifier moduleId;

  if (h==InDetConditions::SCT_MODULE) {
    moduleId = elementId;
  } else if (h==InDetConditions::SCT_SIDE) {
    moduleId = m_pHelper->module_id(elementId); 
  } else if (h==InDetConditions::SCT_STRIP) {
    Identifier waferId = m_pHelper->wafer_id(elementId);
    moduleId = m_pHelper->module_id(waferId); 
  }
  return moduleId;
}

//Returns if element Id is good or bad
bool SCT_DCSConditionsTool::isGood(const Identifier& elementId, const EventContext& ctx, InDetConditions::Hierarchy h) const {
  Identifier moduleId{getModuleID(elementId, h)};
  if (not moduleId.is_valid()) return true; // not canreportabout

  if ((m_readAllDBFolders and m_returnHVTemp) or (not m_readAllDBFolders and not m_returnHVTemp)) {
    const SCT_DCSStatCondData* condDataState{getCondDataState(ctx)};
    if (condDataState==nullptr) return false; // no cond data
    else if (condDataState->output(castId(moduleId))==0) return true; //No params are listed as bad
    else return false;
  } else {
    return true;
  }
}

bool SCT_DCSConditionsTool::isGood(const Identifier& elementId, InDetConditions::Hierarchy h) const {
  const EventContext& ctx{Gaudi::Hive::currentContext()};
  return isGood(elementId, ctx, h);
}

//Does the same for hashIds
bool SCT_DCSConditionsTool::isGood(const IdentifierHash& hashId, const EventContext& ctx) const {
  Identifier waferId{m_pHelper->wafer_id(hashId)};
  Identifier moduleId{m_pHelper->module_id(waferId)};
  return isGood(moduleId, ctx, InDetConditions::SCT_MODULE);
}

void SCT_DCSConditionsTool::getDetectorElementStatus(const EventContext& ctx, InDet::SiDetectorElementStatus &element_status, 
                                                     SG::WriteCondHandle<InDet::SiDetectorElementStatus>* whandle) const {
   if ((m_readAllDBFolders and m_returnHVTemp) or (not m_readAllDBFolders and not m_returnHVTemp)) {
      SG::ReadCondHandle<SCT_DCSStatCondData> condDataHandle{m_condKeyState, ctx};
      if (not condDataHandle.isValid()) {
         ATH_MSG_ERROR("Failed to get " << m_condKeyState.key());
         return;
      }
      if (whandle) {
        whandle->addDependency (condDataHandle);
      }
      const SCT_DCSStatCondData* condDataState{ condDataHandle.cptr() };
      if (condDataState==nullptr) return;
      std::vector<bool> &status = element_status.getElementStatus();
      if (status.empty()) {
         status.resize(m_pHelper->wafer_hash_max(),true);
      }
      std::vector<InDet::ChipFlags_t> &chip_status = element_status.getElementChipStatus();
      constexpr unsigned int N_CHIPS_PER_SIDE = 6;
      constexpr unsigned int N_SIDES = 2;
      constexpr InDet::ChipFlags_t all_flags_set = static_cast<InDet::ChipFlags_t>((1ul<<(N_CHIPS_PER_SIDE*N_SIDES)) - 1ul);
      if (chip_status.empty()) {
         chip_status.resize(status.size(), all_flags_set);
      }
      for ( const std::pair< const CondAttrListCollection::ChanNum, std::vector<std::string> > & dcs_state : condDataState->badChannelsMap() ) {
         Identifier moduleId{ Identifier32{dcs_state.first} };
         //find the corresponding hash
         for (unsigned int side_i=0; side_i<2; ++side_i) {
            Identifier waferId(m_pHelper->wafer_id(m_pHelper->barrel_ec (moduleId),
                                                   m_pHelper->layer_disk(moduleId),
                                                   m_pHelper->phi_module(moduleId),
                                                   m_pHelper->eta_module(moduleId),
                                                   side_i) );
            const IdentifierHash hashId{m_pHelper->wafer_hash(waferId)};
            status.at(hashId.value() ) = status.at(hashId.value() ) && dcs_state.second.empty();
            chip_status.at(hashId.value())  = chip_status.at(hashId.value()) & ( dcs_state.second.empty() ? all_flags_set : 0 );
         }
      }
   }
}

bool SCT_DCSConditionsTool::isGood(const IdentifierHash& hashId) const {
  const EventContext& ctx{Gaudi::Hive::currentContext()};
  return isGood(hashId, ctx);
}

/////////////////////////////////// 

// some lame helper methods: 
// returns HV (s_defaultHV(-30) if there is no information)
float SCT_DCSConditionsTool::modHV(const Identifier& elementId, const EventContext& ctx, InDetConditions::Hierarchy h) const {
  Identifier moduleId{getModuleID(elementId, h)};
  if (not moduleId.is_valid()) return s_defaultHV; // not canreportabout, return s_defaultHV(-30)

  const SCT_DCSFloatCondData* condDataHV{getCondDataHV(ctx)};
  if (condDataHV==nullptr) return s_defaultHV; // no cond data

  float hvval{s_defaultHV};
  if (condDataHV->getValue(castId(moduleId), hvval) and isGood(elementId, h)) {
    return hvval;
  }
  return s_defaultHV; //didn't find the module, return s_defaultHV(-30)
}

float SCT_DCSConditionsTool::modHV(const Identifier& elementId, InDetConditions::Hierarchy h) const {
  const EventContext& ctx{Gaudi::Hive::currentContext()};
  return modHV(elementId, ctx, h);
}

//Does the same for hashIds
float SCT_DCSConditionsTool::modHV(const IdentifierHash& hashId, const EventContext& ctx) const {
  Identifier waferId{m_pHelper->wafer_id(hashId)};
  Identifier moduleId{m_pHelper->module_id(waferId)};
  return modHV(moduleId, ctx, InDetConditions::SCT_MODULE);
}

float SCT_DCSConditionsTool::modHV(const IdentifierHash& hashId) const {
  const EventContext& ctx{Gaudi::Hive::currentContext()};
  return modHV(hashId, ctx);
}

//Returns temp0 (s_defaultTemperature(-40) if there is no information)
float SCT_DCSConditionsTool::hybridTemperature(const Identifier& elementId, const EventContext& ctx, InDetConditions::Hierarchy h) const {
  Identifier moduleId{getModuleID(elementId, h)};
  if (not moduleId.is_valid()) return s_defaultTemperature; // not canreportabout

  const SCT_DCSFloatCondData* condDataTemp0{getCondDataTemp0(ctx)};
  if (condDataTemp0==nullptr) return s_defaultTemperature; // no cond data

  float temperature{s_defaultTemperature};
  if (condDataTemp0->getValue(castId(moduleId), temperature) and isGood(elementId, h)) {
    return temperature;
  }
  return s_defaultTemperature;//didn't find the module, return -40. 
}

float SCT_DCSConditionsTool::hybridTemperature(const Identifier& elementId, InDetConditions::Hierarchy h) const {
  const EventContext& ctx{Gaudi::Hive::currentContext()};
  return hybridTemperature(elementId, ctx, h);
}

//Does the same for hashIds
float SCT_DCSConditionsTool::hybridTemperature(const IdentifierHash& hashId, const EventContext& ctx) const {
  Identifier waferId{m_pHelper->wafer_id(hashId)};
  Identifier moduleId{m_pHelper->module_id(waferId)};
  return hybridTemperature(moduleId, ctx, InDetConditions::SCT_MODULE);
}

float SCT_DCSConditionsTool::hybridTemperature(const IdentifierHash& hashId) const {
  const EventContext& ctx{Gaudi::Hive::currentContext()};
  return hybridTemperature(hashId, ctx);
}

//Returns temp0 + correction for Lorentz angle calculation (s_defaultTemperature(-40) if there is no information)
float SCT_DCSConditionsTool::sensorTemperature(const Identifier& elementId, const EventContext& ctx, InDetConditions::Hierarchy h) const {
  Identifier moduleId{getModuleID(elementId, h)};
  if (not moduleId.is_valid()) return s_defaultTemperature; // not canreportabout

  const SCT_DCSFloatCondData* condDataTemp0{getCondDataTemp0(ctx)};
  if (condDataTemp0==nullptr) return s_defaultTemperature; // no cond data

  float temperature{s_defaultTemperature};
  if (condDataTemp0->getValue(castId(moduleId), temperature) and isGood(elementId, h)) {
    int bec{m_pHelper->barrel_ec(moduleId)};
    if (bec==0) { // Barrel
      return ( temperature + m_barrel_correction);  //return the temp+correction
    } else { // Endcaps
      int modeta{m_pHelper->eta_module(moduleId)};
      if (modeta==2) {
        return (temperature + m_ecInner_correction);  //return the temp+correction
      } else {
        return (temperature + m_ecOuter_correction);  //return the temp+correction
      }
    }
  }
  return s_defaultTemperature;  //didn't find the module, return s_defaultTemperature(-40).
} 

float SCT_DCSConditionsTool::sensorTemperature(const Identifier& elementId, InDetConditions::Hierarchy h) const {
  const EventContext& ctx{Gaudi::Hive::currentContext()};
  return sensorTemperature(elementId, ctx, h);
}

//Does the same for hashIds
float SCT_DCSConditionsTool::sensorTemperature(const IdentifierHash& hashId, const EventContext& ctx) const {
  Identifier waferId{m_pHelper->wafer_id(hashId)};
  Identifier moduleId{m_pHelper->module_id(waferId)};
  return sensorTemperature(moduleId, ctx, InDetConditions::SCT_MODULE);
}

float SCT_DCSConditionsTool::sensorTemperature(const IdentifierHash& hashId) const {
  const EventContext& ctx{Gaudi::Hive::currentContext()};
  return sensorTemperature(hashId, ctx);
}

///////////////////////////////////

const SCT_DCSStatCondData*
SCT_DCSConditionsTool::getCondDataState(const EventContext& ctx) const {
  SG::ReadCondHandle<SCT_DCSStatCondData> condData{m_condKeyState, ctx};
  return condData.retrieve();
}

const SCT_DCSFloatCondData*
SCT_DCSConditionsTool::getCondDataHV(const EventContext& ctx) const {
  SG::ReadCondHandle<SCT_DCSFloatCondData> condData{m_condKeyHV, ctx};
  return condData.retrieve();
}

const SCT_DCSFloatCondData*
SCT_DCSConditionsTool::getCondDataTemp0(const EventContext& ctx) const {
  SG::ReadCondHandle<SCT_DCSFloatCondData> condData{m_condKeyTemp0, ctx};
  return condData.retrieve();
}
