// -*- C++ -*-
/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/ 

#ifndef SCT_MONITORCONDALG
#define SCT_MONITORCONDALG

#include "AthenaBaseComps/AthReentrantAlgorithm.h"
#include "StoreGate/ReadCondHandleKey.h"
#include "AthenaPoolUtilities/CondAttrListCollection.h"
#include "StoreGate/WriteCondHandleKey.h"
#include "SCT_ConditionsData/SCT_MonitorCondData.h"

class SCT_ID;

class SCT_MonitorCondAlg : public AthReentrantAlgorithm 
{
 public:
  SCT_MonitorCondAlg(const std::string& name, ISvcLocator* pSvcLocator);
  virtual ~SCT_MonitorCondAlg() = default;
  virtual StatusCode initialize() override final;
  virtual StatusCode execute(const EventContext& ctx) const override final;
  virtual StatusCode finalize() override final;
  virtual bool isReEntrant() const override final { return false; }

 private:
  SG::ReadCondHandleKey<CondAttrListCollection> m_readKey{this, "ReadKey", "/SCT/Derived/Monitoring", "Key of input (raw) noisy strip conditions folder"};
  SG::WriteCondHandleKey<SCT_MonitorCondData> m_writeKey{this, "WriteKey", "SCT_MonitorCondData", "Key of output (derived) noisy strip conditions data"};
  const SCT_ID* m_helper{nullptr};
};

#endif // SCT_MONITORCONDALG
