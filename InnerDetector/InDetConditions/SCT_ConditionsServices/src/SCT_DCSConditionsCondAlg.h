/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/ 

#ifndef SCT_DCSCONDITIONSCONDALG
#define SCT_DCSCONDITIONSCONDALG

#include "AthenaBaseComps/AthAlgorithm.h"

#include "StoreGate/ReadCondHandleKey.h"
#include "AthenaPoolUtilities/CondAttrListCollection.h"

#include "StoreGate/WriteCondHandleKey.h"
#include "SCT_ConditionsData/SCT_DCSFloatCondData.h"
#include "SCT_ConditionsData/SCT_DCSStatCondData.h"

#include "GaudiKernel/ICondSvc.h"
#include "SCT_Cabling/ISCT_CablingSvc.h"

#include "GaudiKernel/Property.h"

class SCT_DCSConditionsCondAlg : public AthAlgorithm 
{  
 public:
  SCT_DCSConditionsCondAlg(const std::string& name, ISvcLocator* pSvcLocator);
  ~SCT_DCSConditionsCondAlg();
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();

 private:
  enum FloatKey {
    HV,
    TEMPERATURE0,
    TEMPERATURE1
  };

  StatusCode executeStatus();
  StatusCode recordStatus();
  StatusCode executeFloat(FloatKey key);

  SG::ReadCondHandleKey<CondAttrListCollection> m_readKeyHV;
  SG::WriteCondHandleKey<SCT_DCSFloatCondData> m_writeKeyHV;
  EventIDRange m_rangeHV;
  SG::ReadCondHandleKey<CondAttrListCollection> m_readKeyTemperature;
  SG::WriteCondHandleKey<SCT_DCSFloatCondData> m_writeKeyTemperature0;
  SG::WriteCondHandleKey<SCT_DCSFloatCondData> m_writeKeyTemperature1;
  EventIDRange m_rangeTemperature;
  SG::ReadCondHandleKey<CondAttrListCollection> m_readKeyStatus;
  SG::WriteCondHandleKey<SCT_DCSStatCondData> m_writeKeyStatus;
  SCT_DCSStatCondData* m_writeCdoStatus;
  EventIDRange m_rangeStatus;

  ServiceHandle<ICondSvc> m_condSvc;
  ServiceHandle<ISCT_CablingSvc> m_cablingSvc; //!< Handle on SCT cabling service

  BooleanProperty m_readAllDBFolders;
  BooleanProperty m_returnHVTemp;
  std::string m_chanstatCut;
  float m_hvLowLimit;
  float m_hvUpLimit;
  bool m_useHV;
  float m_useHVLowLimit;
  float  m_useHVUpLimit;
  std::string m_useHVChanCut;
};

#endif // SCT_DCSCONDITIONSCONDALG
