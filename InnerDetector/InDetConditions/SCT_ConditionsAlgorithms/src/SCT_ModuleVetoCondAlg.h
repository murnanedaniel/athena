// -*- C++ -*-

/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/ 

#ifndef SCT_MODULEVETOCONDALG
#define SCT_MODULEVETOCONDALG

#include "AthenaBaseComps/AthReentrantAlgorithm.h"

#include "AthenaPoolUtilities/AthenaAttributeList.h"
#include "SCT_ConditionsData/SCT_ModuleVetoCondData.h"
#include "StoreGate/ReadCondHandleKey.h"
#include "StoreGate/WriteCondHandleKey.h"

#include "Gaudi/Property.h"

class SCT_ModuleVetoCondAlg : public AthReentrantAlgorithm 
{  
 public:
  SCT_ModuleVetoCondAlg(const std::string& name, ISvcLocator* pSvcLocator);
  virtual ~SCT_ModuleVetoCondAlg() = default;
  virtual StatusCode initialize() override final;
  virtual StatusCode execute(const EventContext& ctx) const override final;
  virtual StatusCode finalize() override final;
  virtual bool isReEntrant() const override final { return false; }

 private:
  SG::ReadCondHandleKey<AthenaAttributeList> m_readKey{this, "ReadKey", "/SCT/Manual/BadModules", "Key of input (raw) bad module conditions folder"};
  SG::WriteCondHandleKey<SCT_ModuleVetoCondData> m_writeKey{this, "WriteKey", "SCT_ModuleVetoCondData", "Key of output (derived) bad module conditions data"};
};

#endif // SCT_MODULEVETOCONDALG
