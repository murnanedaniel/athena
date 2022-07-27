// -*- C++ -*-

/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/ 

#ifndef SCT_CONDITIONSPARAMETERCONDALG
#define SCT_CONDITIONSPARAMETERCONDALG

#include "AthenaBaseComps/AthReentrantAlgorithm.h"

#include "StoreGate/ReadCondHandleKey.h"
#include "AthenaPoolUtilities/CondAttrListVec.h"
#include "StoreGate/WriteCondHandleKey.h"
#include "SCT_ConditionsData/SCT_CondParameterData.h"
#include "SCT_Cabling/ISCT_CablingTool.h"

#include "Gaudi/Property.h"

class SCT_ConditionsParameterCondAlg : public AthReentrantAlgorithm 
{  
 public:
  SCT_ConditionsParameterCondAlg(const std::string& name, ISvcLocator* pSvcLocator);
  virtual ~SCT_ConditionsParameterCondAlg() = default;
  virtual StatusCode initialize() override final;
  virtual StatusCode execute(const EventContext& ctx) const override final;
  virtual StatusCode finalize() override final;
  virtual bool isReEntrant() const override final { return false; }
 private:
  SG::ReadCondHandleKey<CondAttrListVec> m_readKey{this, "ReadKey", "/SCT/DAQ/Configuration/Chip", "Key of input (raw) chip conditions folder"};
  SG::WriteCondHandleKey<SCT_CondParameterData> m_writeKey{this, "WriteKey", "SCT_CondParameterData", "Key of output (derived) average threshold conditions data"};

  ToolHandle<ISCT_CablingTool> m_cablingTool{this, "SCT_CablingTool", "SCT_CablingTool", "Tool to retrieve SCT Cabling"};
};

#endif // SCT_CONDITIONSPARAMETERCONDALG
