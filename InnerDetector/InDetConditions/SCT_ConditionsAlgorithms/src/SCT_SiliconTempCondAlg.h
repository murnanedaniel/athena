/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/ 

#ifndef SCT_SILICONTEMPCONDALG
#define SCT_SILICONTEMPCONDALG

#include "AthenaBaseComps/AthReentrantAlgorithm.h"

#include "SCT_ConditionsData/SCT_DCSStatCondData.h"
#include "SCT_ConditionsData/SCT_DCSFloatCondData.h"
#include "SCT_ConditionsTools/ISCT_DCSConditionsTool.h"
#include "StoreGate/ReadCondHandleKey.h"
#include "StoreGate/WriteCondHandleKey.h"

class SCT_ID;

class SCT_SiliconTempCondAlg : public AthReentrantAlgorithm
{  
 public:
  SCT_SiliconTempCondAlg(const std::string& name, ISvcLocator* pSvcLocator);
  virtual ~SCT_SiliconTempCondAlg() = default;
  virtual StatusCode initialize() override final;
  virtual StatusCode execute(const EventContext& ctx) const override final;
  virtual StatusCode finalize() override final;
  virtual bool isReEntrant() const override final { return false; }

 private:
  BooleanProperty m_useState{this, "UseState", true, "Flag to use state conditions folder"};
  SG::ReadCondHandleKey<SCT_DCSStatCondData> m_readKeyState{this, "ReadKeyState", "SCT_DCSStatCondData", "Key of input state conditions data"};
  SG::ReadCondHandleKey<SCT_DCSFloatCondData> m_readKeyTemp0{this, "ReadKeyTemp", "SCT_DCSTemp0CondData", "Key of input (hybrid) temperature conditions data"};
  SG::WriteCondHandleKey<SCT_DCSFloatCondData> m_writeKey{this, "WriteKey", "SCT_SiliconTempCondData", "Key of output (sensor) temperature conditions data"};
  ToolHandle<ISCT_DCSConditionsTool> m_sctDCSTool{this, "DCSConditionsTool", "InDetSCT_DCSConditionsTool", "Tool to retrived SCT DCS information"};
  const SCT_ID* m_pHelper{nullptr}; //!< ID helper for SCT
};

#endif // SCT_SILICONTEMPCONDALG
