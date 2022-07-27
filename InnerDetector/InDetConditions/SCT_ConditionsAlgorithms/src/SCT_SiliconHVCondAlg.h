/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/ 

#ifndef SCT_SILICONHVCONDALG
#define SCT_SILICONHVCONDALG

#include "AthenaBaseComps/AthReentrantAlgorithm.h"

#include "StoreGate/ReadCondHandleKey.h"
#include "StoreGate/WriteCondHandleKey.h"
#include "SCT_ConditionsData/SCT_DCSStatCondData.h"
#include "SCT_ConditionsData/SCT_DCSFloatCondData.h"
#include "SCT_ConditionsTools/ISCT_DCSConditionsTool.h"

class SCT_ID;

class SCT_SiliconHVCondAlg : public AthReentrantAlgorithm
{  
 public:
  SCT_SiliconHVCondAlg(const std::string& name, ISvcLocator* pSvcLocator);
  virtual ~SCT_SiliconHVCondAlg() = default;
  virtual StatusCode initialize() override final;
  virtual StatusCode execute(const EventContext& ctx) const override final;
  virtual StatusCode finalize() override final;
  virtual bool isReEntrant() const override final { return false; }

 private:
  BooleanProperty m_useState{this, "UseState", true, "Flag to use state conditions folder"};
  SG::ReadCondHandleKey<SCT_DCSStatCondData> m_readKeyState{this, "ReadKeyState", "SCT_DCSStatCondData", "Key of input state conditions data"};
  SG::ReadCondHandleKey<SCT_DCSFloatCondData> m_readKeyHV{this, "ReadKeyHV", "SCT_DCSHVCondData", "Key of input HV conditions data"};
  SG::WriteCondHandleKey<SCT_DCSFloatCondData> m_writeKey{this, "WriteKey", "SCT_SiliconBiasVoltCondData", "Key of output bias voltage conditions data"};
  ToolHandle<ISCT_DCSConditionsTool> m_sctDCSTool{this, "DCSConditionsTool", "InDetSCT_DCSConditionsTool", "Tool to retrived SCT DCS information"};
  const SCT_ID* m_pHelper{nullptr};//!< ID helper for SCT
};

#endif // SCT_SILICONHVCONDALG
