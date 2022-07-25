// -*- C++ -*-

/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

/** @file SCT_ReadCalibDataTool.h Header file for SCT_ReadCalibDataTool.
    @author Per Johansson, 17/05/08, based on ReadSCT_CalibData.
*/

// Multiple inclusion protection
#ifndef SCT_READ_CALIB_DATA_TOOL
#define SCT_READ_CALIB_DATA_TOOL

#include "AthenaBaseComps/AthAlgTool.h"
#include "SCT_ConditionsTools/ISCT_ReadCalibDataTool.h"

#include "InDetReadoutGeometry/SiDetectorElementCollection.h"
#include "SCT_Cabling/ISCT_CablingTool.h"
#include "SCT_ConditionsData/SCT_CalibDefectData.h"
#include "SCT_ConditionsData/SCT_WaferGoodStripInfo.h"
#include "SCT_ConditionsData/SCT_AllGoodStripInfo.h"

// Include Gaudi classes
#include "GaudiKernel/ServiceHandle.h"

// Forward declarations
class SCT_ID;

/** This class contains a Tool that reads SCT calibration data and makes it available to 
    other algorithms. The current implementation reads the data from a COOL database. 
*/

class SCT_ReadCalibDataTool: public extends<AthAlgTool, ISCT_ReadCalibDataTool> {

 public:
  enum {STRIPS_PER_WAFER=768};
  //----------Public Member Functions----------//
  // Structors
  SCT_ReadCalibDataTool(const std::string& type, const std::string& name, const IInterface* parent); //!< Constructor
  virtual ~SCT_ReadCalibDataTool() = default;                               //!< Destructor
  
  // Standard Gaudi functions
  virtual StatusCode initialize() override; //!< Gaudi initialiser
  virtual StatusCode finalize() override; //!< Gaudi finaliser
  
  /// @name Methods to be implemented from virtual baseclass methods, when introduced
  ///Return whether this tool can report on the hierarchy level (e.g. module, chip...)
  virtual bool canReportAbout(InDetConditions::Hierarchy h) const override;
  ///Summarise the result from the tool as good/bad
  virtual bool isGood(const Identifier& elementId, InDetConditions::Hierarchy h=InDetConditions::DEFAULT) const override;
  virtual bool isGood(const Identifier& elementId, const EventContext& ctx, InDetConditions::Hierarchy h=InDetConditions::DEFAULT) const override;
  ///same thing with id hash, introduced by shaun with dummy method for now
  virtual bool isGood(const IdentifierHash& /*hashId*/) const override { return true; }
  virtual bool isGood(const IdentifierHash& hashId, const EventContext& /*ctx*/) const override { return isGood(hashId); }
  virtual void getDetectorElementStatus(const EventContext& ctx, InDet::SiDetectorElementStatus &element_status, 
                                        SG::WriteCondHandle<InDet::SiDetectorElementStatus>* whandle) const override;
  //@}
  
  // Methods to return calibration defect type and summary
  virtual ISCT_ReadCalibDataTool::CalibDefectType defectType(const Identifier& stripId, const EventContext& ctx, InDetConditions::Hierarchy h=InDetConditions::DEFAULT) const override; //!<Return summary of defect type and values for a strip
  virtual ISCT_ReadCalibDataTool::CalibDefectType defectType(const Identifier& stripId, InDetConditions::Hierarchy h=InDetConditions::DEFAULT) const override; //!<Return summary of defect type and values for a strip
  virtual SCT_CalibDefectData::CalibModuleDefects defectsSummary(const Identifier& moduleId, const std::string& scan, const EventContext& ctx) const override; //!<Returns module summary of defect
  virtual SCT_CalibDefectData::CalibModuleDefects defectsSummary(const Identifier& moduleId, const std::string& scan) const override; //!<Returns module summary of defect
  virtual std::list<Identifier> defectList(const std::string& defect, const EventContext& ctx) const override; //!<Returns module summary of defect
  virtual std::list<Identifier> defectList(const std::string& defect) const override; //!<Returns module summary of defect

 private:
  // Read Cond Handles
  SG::ReadCondHandleKey<SCT_CalibDefectData> m_condKeyGain{this, "CondKeyGain", "SCT_CalibDefectNPtGain", "SCT defects due to NPtGain calibration"};
  SG::ReadCondHandleKey<SCT_CalibDefectData> m_condKeyNoise{this, "CondKeyNoise", "SCT_CalibDefectNPtNoise", "SCT defects due to NPtNoise calibration"};
  SG::ReadCondHandleKey<SCT_AllGoodStripInfo> m_condKeyInfo{this, "CondKeyInfo", "SCT_AllGoodStripInfo", "SCT good strips based on calibration defects"};
  SG::ReadCondHandleKey<InDetDD::SiDetectorElementCollection> m_SCTDetEleCollKey{this, "SCTDetEleCollKey", "SCT_DetectorElementCollection", "Key of SiDetectorElementCollection for SCT"};

  const SCT_CalibDefectData* getCondDataGain(const EventContext& ctx) const;
  const SCT_CalibDefectData* getCondDataNoise(const EventContext& ctx) const;
  const SCT_AllGoodStripInfo* getCondDataInfo(const EventContext& ctx) const;
  const InDetDD::SiDetectorElement* getDetectorElement(const IdentifierHash& waferHash, const EventContext& ctx) const;
  
  //----------Private Attributes----------//
  ToolHandle<ISCT_CablingTool> m_cabling{this, "SCT_CablingTool", "SCT_CablingTool", "Tool to retrieve SCT Cabling"};
  const SCT_ID* m_id_sct{nullptr}; //!< Handle to SCT ID helper
};

//---------------------------------------------------------------------- 
#endif // SCT_READ_CALIB_DATA_TOOL
