// -*- C++ -*-

/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/ 

#ifndef SCT_ReadCalibChipGainCondAlg_h
#define SCT_ReadCalibChipGainCondAlg_h

// Include parent class
#include "AthenaBaseComps/AthReentrantAlgorithm.h"

// Include Gaudi classes
#include "Gaudi/Property.h"

// Include Athena classes
#include "StoreGate/ReadCondHandleKey.h"
#include "StoreGate/WriteCondHandleKey.h"
#include "AthenaPoolUtilities/CondAttrListCollection.h"
#include "SCT_ConditionsData/SCT_GainCalibData.h"

// Forward declarations
class SCT_ID;

class SCT_ReadCalibChipGainCondAlg : public AthReentrantAlgorithm 
{  
 public:
  SCT_ReadCalibChipGainCondAlg(const std::string& name, ISvcLocator* pSvcLocator);
  virtual ~SCT_ReadCalibChipGainCondAlg() = default;
  virtual StatusCode initialize() override final;
  virtual StatusCode execute(const EventContext& ctx) const override final;
  virtual StatusCode finalize() override final;
  virtual bool isReEntrant() const override final { return false; }

 private:
  static void insertNptGainFolderData(SCT_ModuleGainCalibData& theseCalibData, const coral::AttributeList& folderData) ;

  SG::ReadCondHandleKey<CondAttrListCollection> m_readKey{this, "ReadKey", "/SCT/DAQ/Calibration/ChipGain", "Key of input (raw) gain conditions folder"};
  SG::WriteCondHandleKey<SCT_GainCalibData> m_writeKey{this, "WriteKey", "SCT_GainCalibData", "Key of output (derived) gain conditions data"};
  const SCT_ID* m_id_sct{nullptr}; //!< Handle to SCT ID helper
};

#endif // SCT_ReadCalibChipGainCondAlg_h
