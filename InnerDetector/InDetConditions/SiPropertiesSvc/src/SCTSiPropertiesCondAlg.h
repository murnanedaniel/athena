/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
*/ 

#ifndef SCTSIPROPERTIESCONDALG
#define SCTSIPROPERTIESCONDALG

#include "AthenaBaseComps/AthAlgorithm.h"

#include "StoreGate/ReadCondHandleKey.h"
#include "SCT_ConditionsData/SCT_DCSFloatCondData.h"
#include "StoreGate/WriteCondHandleKey.h"
#include "SiPropertiesSvc/SiliconPropertiesVector.h"

#include "GaudiKernel/ICondSvc.h"

class SCT_ID;
class ISiliconConditionsSvc;
namespace InDetDD {
  class SiDetectorManager;
}

class SCTSiPropertiesCondAlg : public AthAlgorithm 
{  
 public:
  SCTSiPropertiesCondAlg(const std::string& name, ISvcLocator* pSvcLocator);
  ~SCTSiPropertiesCondAlg();
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();

 private:
  double m_temperatureMin;
  double m_temperatureMax;
  double m_temperatureDefault;
  SG::ReadCondHandleKey<SCT_DCSFloatCondData> m_readKeyTemp;
  SG::ReadCondHandleKey<SCT_DCSFloatCondData> m_readKeyHV;
  SG::WriteCondHandleKey<InDet::SiliconPropertiesVector> m_writeKey;
  ServiceHandle<ICondSvc> m_condSvc;
  ServiceHandle<ISiliconConditionsSvc> m_siCondSvc;
  const SCT_ID* m_pHelper; //!< ID helper for SCT
  const InDetDD::SiDetectorManager* m_detManager;
};

#endif // SCTSIPROPERTIESCONDALG
