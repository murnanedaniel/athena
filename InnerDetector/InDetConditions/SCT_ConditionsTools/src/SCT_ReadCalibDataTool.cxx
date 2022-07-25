/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

/** @file SCT_ReadCalibDataTool.cxx Implementation file for SCT_ReadCalibDataTool.
    @author Per Johansson, 17/05/08, based on SCT_ReadCalibDataTool.
*/

#include "SCT_ReadCalibDataTool.h"

// Include Athena stuff
#include "InDetIdentifier/SCT_ID.h"
#include "InDetReadoutGeometry/SiDetectorElement.h"
#include "SCT_ReadoutGeometry/SCT_ChipUtils.h"
#include "SCT_DetectorElementStatus.h"
#include "InDetIdentifier/SCT_ID.h"

// Include STL
#include <cstdint>

//----------------------------------------------------------------------
SCT_ReadCalibDataTool::SCT_ReadCalibDataTool(const std::string& type, const std::string& name, const IInterface* parent) :
  base_class(type, name, parent)
{
}

//----------------------------------------------------------------------
StatusCode SCT_ReadCalibDataTool::initialize() {
  // Print where you are
  ATH_MSG_DEBUG("in initialize()");

  // Get SCT helper
  ATH_CHECK(detStore()->retrieve(m_id_sct, "SCT_ID"));

  // Retrieve SCT Cabling tool
  ATH_CHECK(m_cabling.retrieve());

  // Read Cond Handle Key
  ATH_CHECK(m_condKeyGain.initialize());
  ATH_CHECK(m_condKeyNoise.initialize());
  ATH_CHECK(m_condKeyInfo.initialize());
  ATH_CHECK(m_SCTDetEleCollKey.initialize());

  return StatusCode::SUCCESS;
} // SCT_ReadCalibDataTool::initialize()

//----------------------------------------------------------------------
StatusCode SCT_ReadCalibDataTool::finalize() {
  return StatusCode::SUCCESS;
} // SCT_ReadCalibDataTool::finalize()

//----------------------------------------------------------------------
//Can only report good/bad at strip level
bool SCT_ReadCalibDataTool::canReportAbout(InDetConditions::Hierarchy h) const {
  return (h==InDetConditions::SCT_STRIP);
}

//----------------------------------------------------------------------
// Returns a bool summary of the data
bool SCT_ReadCalibDataTool::isGood(const Identifier& elementId, const EventContext& ctx, InDetConditions::Hierarchy h) const {
  // Status of the compId
  bool status{true};
  switch (h) {
  case InDetConditions::SCT_STRIP:
    {
      // Retrieve isGood Wafer data
      const SCT_AllGoodStripInfo* condDataInfo{getCondDataInfo(ctx)};
      if (condDataInfo==nullptr) {
        ATH_MSG_ERROR("In isGood, SCT_AllGoodStripInfo cannot be retrieved");
        return false;
      }
      // Extract the wafer identifier from the strip identifier
      Identifier waferId{m_id_sct->wafer_id(elementId)};
      // Get hashId
      IdentifierHash waferHash{m_id_sct->wafer_hash(waferId)};
      // Get strip on wafer to check
      int strip{m_id_sct->strip(elementId)};
      // Set value
      status = (*condDataInfo)[waferHash][strip];
      break;
    }
  case InDetConditions::SCT_MODULE:
    {
      // Not applicable for Calibration data
      ATH_MSG_WARNING("summary(): Module good/bad is not applicable for Calibration data");
      break;
    }
  case InDetConditions::SCT_SIDE:
    {
      // Not applicable for Calibration data
      ATH_MSG_WARNING("summary(): Wafer good/bad is not applicable for Calibration data");
      break;
    }
  case InDetConditions::SCT_CHIP:
    {
      // Not applicable for Calibration data
      ATH_MSG_WARNING("summary(): Chip good/bad is not applicable for Calibration data");
      break;
    }
  default:
    {
      ATH_MSG_INFO("Unknown component has been asked for, should be Module/Wafer/Chip or Strip; returning 'good' and continuing");
    }    
  } //end of switch structure
  
  // Print status  
  return status;
} //SCT_ReadCalibDataTool::summary()

bool SCT_ReadCalibDataTool::isGood(const Identifier& elementId, InDetConditions::Hierarchy h) const {
  const EventContext& ctx{Gaudi::Hive::currentContext()};

  return isGood(elementId, ctx, h);
}

void SCT_ReadCalibDataTool::getDetectorElementStatus(const EventContext& ctx, InDet::SiDetectorElementStatus &element_status, 
                                                     SG::WriteCondHandle<InDet::SiDetectorElementStatus>* whandle) const {
   SG::ReadCondHandle<SCT_AllGoodStripInfo> condDataHandle{m_condKeyInfo, ctx};
   if (not condDataHandle.isValid()) {
      ATH_MSG_ERROR("Invalid cond data handle " << m_condKeyInfo.key() );
      return;
   }
   if (whandle) {
     whandle->addDependency (condDataHandle);
   }
   const SCT_AllGoodStripInfo* condDataInfo{condDataHandle.cptr()};

   const std::vector<bool> &status = element_status.getElementStatus();
   const std::vector<InDet::ChipFlags_t> &chip_status = element_status.getElementChipStatus();

   std::vector<std::vector<unsigned short> >  &bad_strips = element_status.getBadCells();
   if (bad_strips.empty()) {
      bad_strips.resize(condDataInfo->size());
   }
   unsigned int element_i=0;
   for( const std::array<bool, SCT_ConditionsData::STRIPS_PER_WAFER> &good_strips : *condDataInfo) {
      IdentifierHash moduleHash(element_i);
      if (status.empty() || status.at(element_i)) {
         std::vector<unsigned short>  &bad_module_strips = bad_strips[element_i];
         unsigned int last_geoemtrical_chip_id=SCT::N_CHIPS_PER_SIDE;
         for (unsigned int strip_i=0; strip_i<good_strips.size(); ++strip_i) {
            unsigned int geoemtrical_chip_id = SCT::getGeometricalChipID(strip_i);
            if (geoemtrical_chip_id != last_geoemtrical_chip_id) {
               last_geoemtrical_chip_id=geoemtrical_chip_id;
               if (!chip_status.empty() && !(chip_status.at(element_i) & static_cast<InDet::ChipFlags_t>(1ul<<geoemtrical_chip_id))) {
                  strip_i += (SCT::N_STRIPS_PER_CHIP-1);
                  continue;
               }
            }
            if (!good_strips[strip_i]) {
               std::vector<unsigned short>::const_iterator iter = std::lower_bound(bad_module_strips.begin(),bad_module_strips.end(),strip_i);
               if (iter == bad_module_strips.end() || *iter != strip_i) {
                  bad_module_strips.insert( iter, strip_i);
               }
            }
         }
      }
      ++element_i;
   }
}

//----------------------------------------------------------------------
// Returns a defect summary of a defect strip, scan, type and value
ISCT_ReadCalibDataTool::CalibDefectType SCT_ReadCalibDataTool::defectType(const Identifier& stripId, const EventContext& ctx, InDetConditions::Hierarchy h) const {
  // Print where you are
  ATH_MSG_DEBUG("in defectType()");

  // Create the calibDefectSummary
  CalibDefectType theseSummaryDefects;

  // Retrieve defect data
  const SCT_CalibDefectData* condDataGain{getCondDataGain(ctx)};
  if (condDataGain==nullptr) {
    ATH_MSG_ERROR("In defectType, SCT_CalibDefectData (gain) cannot be retrieved.");
    return theseSummaryDefects;
  }
  const SCT_CalibDefectData* condDataNoise{getCondDataNoise(ctx)};
  if (condDataNoise==nullptr) {
    ATH_MSG_ERROR("In defectType, SCT_CalibDefectData (noise) cannot be retrieved.");
    return theseSummaryDefects;
  }

  // Extract the moduleId from the comp identifier
  Identifier moduleId{m_id_sct->module_id(stripId)};
  ATH_MSG_DEBUG("Summary wanted for component: " << stripId << " on module: " << moduleId);

  // Create the CalibDataDefect objects
  SCT_CalibDefectData::CalibModuleDefects wantedNPGDefects{condDataGain->findModule(moduleId)};
  SCT_CalibDefectData::CalibModuleDefects wantedNODefects{condDataNoise->findModule(moduleId)};

  switch (h) {
  case InDetConditions::SCT_MODULE:
    {
      // Not applicable for Calibration data
      ATH_MSG_WARNING("summary(): Module defect summary is not applicable for Calibration data");
      break;
    }

  case InDetConditions::SCT_SIDE:
    {
      // Not applicable for Calibration data
      ATH_MSG_WARNING("summary(): Wafer defect summary is not applicable for Calibration data");
      break;
    }

  case InDetConditions::SCT_CHIP:
    {
      // Not applicable for Calibration data
      ATH_MSG_WARNING("summary(): Chip defect summary is not applicable for Calibration data");
      break;
    }
  case InDetConditions::SCT_STRIP:
    {
      // Get the strip/channel number to check
      int side{m_id_sct->side(stripId)};
      int strip{m_id_sct->strip(stripId)};
      const Identifier waferId{m_id_sct->wafer_id(stripId)};
      const IdentifierHash waferHash{m_id_sct->wafer_hash(waferId)};
      unsigned int stripNum;
      const InDetDD::SiDetectorElement* p_element{getDetectorElement(waferHash, ctx)};
      if (p_element->swapPhiReadoutDirection()) {
        if (side == 0) {
          stripNum =   STRIPS_PER_WAFER-1 - strip;
        } else {
          stripNum = 2*STRIPS_PER_WAFER-1 - strip;
        }      
      } else {
        stripNum = side * STRIPS_PER_WAFER + strip;
      }
      
      // Find the bad strip and fill calibDefectSummary
      if (wantedNPGDefects.begDefects.empty()) {
        ATH_MSG_VERBOSE("No NPtGain defects in this module");
      } else {
        for (unsigned int i{0}; i<wantedNPGDefects.begDefects.size(); ++i) {
          if (stripNum>=wantedNPGDefects.begDefects[i] and stripNum<=wantedNPGDefects.endDefects[i]) {
            theseSummaryDefects.scan.emplace_back("NPtGain");
            theseSummaryDefects.defect.push_back(wantedNPGDefects.typeOfDefect[i]);
            theseSummaryDefects.value.push_back(wantedNPGDefects.parValue[i]);
            ATH_MSG_VERBOSE("NPtGain defect summary for strip " << stripNum << " filled");
          }
        }
      }

      if (wantedNODefects.begDefects.empty()) {
        ATH_MSG_VERBOSE("No NoiseOccupancy defects in this module");
      } else {
        for (unsigned int i{0}; i != wantedNODefects.begDefects.size(); ++i) {
          if (stripNum>=wantedNODefects.begDefects[i] and stripNum <= wantedNODefects.endDefects[i]) {
            theseSummaryDefects.scan.emplace_back("NoiseOccupancy");
            theseSummaryDefects.defect.push_back(wantedNODefects.typeOfDefect[i]);
            theseSummaryDefects.value.push_back(wantedNODefects.parValue[i]);
            ATH_MSG_VERBOSE("NoiseOccupancy defect summary for strip " << stripNum << "  filled");
          }
        }
      } 
      if (theseSummaryDefects.scan.empty()) {
        ATH_MSG_VERBOSE("defectSummary(): can't retrieve the defects for this strip: " <<  stripNum << " since strip good");
      }     
      break;
    }
  default:
    {
      ATH_MSG_INFO("Unknown component requested, should be one of Module/Side/Chip or Strip");
      return theseSummaryDefects;
    }

  }//end of switch structure

  return theseSummaryDefects;
} //SCT_ReadCalibDataTool::defectType()

ISCT_ReadCalibDataTool::CalibDefectType SCT_ReadCalibDataTool::defectType(const Identifier& stripId, InDetConditions::Hierarchy h) const {
  const EventContext& ctx{Gaudi::Hive::currentContext()};
  return defectType(stripId, ctx, h);
}
//----------------------------------------------------------------------
// Returns a summary of all defects on a module for a given scan
SCT_CalibDefectData::CalibModuleDefects SCT_ReadCalibDataTool::defectsSummary(const Identifier& moduleId, const std::string& scan, const EventContext& ctx) const {
  // Create pointer to the CalibDataDefect object 
  SCT_CalibDefectData::CalibModuleDefects wantedDefects;

  // Retrieve the correct defect map
  if (scan == "NPtGain") {
    const SCT_CalibDefectData* condDataGain{getCondDataGain(ctx)};
    if (condDataGain==nullptr) {
      ATH_MSG_ERROR("In defectType, SCT_CalibDefectData (gain) cannot be retrieved.");
    } else {
      wantedDefects = condDataGain->findModule(moduleId);
    }
  } else if (scan == "NoiseOccupancy") {
    const SCT_CalibDefectData* condDataNoise{getCondDataNoise(ctx)};
    if (condDataNoise==nullptr) {
      ATH_MSG_ERROR("In defectType, SCT_CalibDefectData (noise) cannot be retrieved.");
    } else {
      wantedDefects = condDataNoise->findModule(moduleId);
    }
  } else {
    ATH_MSG_ERROR("defectsSummary(): Module defects for scan" << scan << " does not exist (only NPtGain or NoiseOccupancy).");
  }

  return wantedDefects;
} //SCT_ReadCalibDataTool::defectsSummary()

SCT_CalibDefectData::CalibModuleDefects SCT_ReadCalibDataTool::defectsSummary(const Identifier& moduleId, const std::string& scan) const {
  const EventContext& ctx{Gaudi::Hive::currentContext()};
  return defectsSummary(moduleId, scan, ctx);
}

//---------------------------------------------------------------------- 
//----------------------------------------------------------------------
// Returns a list of all strips with a certain defects
std::list<Identifier> SCT_ReadCalibDataTool::defectList(const std::string& defect, const EventContext& ctx) const {
  std::list<Identifier> defectList;

  // Retrieve defect data
  const SCT_CalibDefectData* condDataGain{getCondDataGain(ctx)};
  if (condDataGain==nullptr) {
    ATH_MSG_ERROR("In defectType, SCT_CalibDefectData (gain) cannot be retrieved.");
    return defectList;
  }
  const SCT_CalibDefectData* condDataNoise{getCondDataNoise(ctx)};
  if (condDataNoise==nullptr) {
    ATH_MSG_ERROR("In defectType, SCT_CalibDefectData (noise) cannot be retrieved.");
    return defectList;
  }

  // Create pointer to the CalibDataDefect object 
  SCT_CalibDefectData::CalibModuleDefects wantedDefects;
  
  //Check which scan the defect belongs
  bool npDefect{false};
  bool noDefect{false};
  if (defect=="NO_HI" or defect=="BAD_OPE" or defect=="DOUBTR_HI") {
    noDefect = true;
  } else {
    npDefect = true;
  }
  
  //Loop over all wafers using hashIds from the cabling service
  std::vector<std::uint32_t> listOfRODs;
  m_cabling->getAllRods(listOfRODs, ctx);
  std::vector<std::uint32_t>::iterator rodIter{listOfRODs.begin()};
  std::vector<std::uint32_t>::iterator rodEnd{listOfRODs.end()};
  for (; rodIter!=rodEnd; ++rodIter) {
    std::vector<IdentifierHash> listOfHashes;
    m_cabling->getHashesForRod(listOfHashes, *rodIter, ctx);
    std::vector<IdentifierHash>::iterator hashIt{listOfHashes.begin()};
    std::vector<IdentifierHash>::iterator hashEnd{listOfHashes.end()};
    for (; hashIt != hashEnd; ++hashIt) {
      Identifier waferId{m_id_sct->wafer_id(*hashIt)};
      int side{m_id_sct->side(waferId)};
      //Only use the hashid for side 0, since data saved per module basis
      if (side!=0) continue;
      Identifier moduleId{m_id_sct->module_id(waferId)};
      if (npDefect) {
        wantedDefects = condDataGain->findModule(moduleId);
      }  else if (noDefect) {
        wantedDefects = condDataNoise->findModule(moduleId);
      }
      if (not wantedDefects.begDefects.empty()) {
        for (unsigned int i{0}; i < wantedDefects.begDefects.size(); ++i) {
          if (wantedDefects.typeOfDefect[i] == defect) {
            // Create identifier for all strips inside begin to end
            int strip_beg{static_cast<int>(wantedDefects.begDefects[i])};
            int strip_end{static_cast<int>(wantedDefects.endDefects[i])};
            // In DB: strip from 0-1535, need to convert to side and 0-767 and take into account phiSwaps
            for (int strip{strip_beg}; strip<strip_end+1; strip++) {
              int nside{(strip<STRIPS_PER_WAFER) ? 0 : 1};
              int strip_cor;
              const InDetDD::SiDetectorElement* p_element;
              if (nside==1) { // if side 1 need waferId of side 1 to get phiswap and correct stripId
                IdentifierHash hashSide1;
                m_id_sct->get_other_side(*hashIt, hashSide1);
                waferId = m_id_sct->wafer_id(hashSide1);
                p_element = (getDetectorElement(hashSide1, ctx));
              } else {
                p_element = (getDetectorElement(*hashIt, ctx));
              }
              if (p_element->swapPhiReadoutDirection()) {
                if (nside == 0) {
                  strip_cor =   STRIPS_PER_WAFER-1 - strip;
                } else {
                  strip_cor = 2*STRIPS_PER_WAFER-1 - strip;
                }      
              } else {
                strip_cor = strip - nside * STRIPS_PER_WAFER;
              }
              Identifier stripId{m_id_sct->strip_id(waferId,strip_cor)};
              defectList.push_back(stripId);
            }
          }
        }
      }
    }
  }
  return defectList;
} //SCT_ReadCalibDataTool::defects()

std::list<Identifier> SCT_ReadCalibDataTool::defectList(const std::string& defect) const {
  const EventContext& ctx{Gaudi::Hive::currentContext()};
  return defectList(defect, ctx);
}
//---------------------------------------------------------------------- 

const SCT_CalibDefectData*
SCT_ReadCalibDataTool::getCondDataGain(const EventContext& ctx) const {
  SG::ReadCondHandle<SCT_CalibDefectData> condData{m_condKeyGain, ctx};
  return condData.retrieve();
}

//----------------------------------------------------------------------

const SCT_CalibDefectData*
SCT_ReadCalibDataTool::getCondDataNoise(const EventContext& ctx) const {
  SG::ReadCondHandle<SCT_CalibDefectData> condData{m_condKeyNoise, ctx};
  return condData.retrieve();
}

//----------------------------------------------------------------------

const SCT_AllGoodStripInfo*
SCT_ReadCalibDataTool::getCondDataInfo(const EventContext& ctx) const {
  SG::ReadCondHandle<SCT_AllGoodStripInfo> condData{m_condKeyInfo, ctx};
  return condData.retrieve();
}

//----------------------------------------------------------------------

const InDetDD::SiDetectorElement*
SCT_ReadCalibDataTool::getDetectorElement(const IdentifierHash& waferHash, const EventContext& ctx) const {
  SG::ReadCondHandle<InDetDD::SiDetectorElementCollection> condData{m_SCTDetEleCollKey, ctx};
  if (not condData.isValid()) return nullptr;
  return condData->getDetectorElement(waferHash);
}

//----------------------------------------------------------------------
