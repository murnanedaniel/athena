/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef CONDALGS_CONDALGX_H
#define CONDALGS_CONDALGX_H 1

#include "AthenaBaseComps/AthAlgorithm.h"
#include "StoreGate/ReadHandle.h"
#include "StoreGate/WriteCondHandleKey.h"

#include "AthExHive/CondDataObj.h"
#include "AthExHive/IASCIICondDbSvc.h"

#include "xAODEventInfo/EventInfo.h"
#include "GaudiKernel/ICondSvc.h"

#include <string>

class CondAlgX  :  public AthAlgorithm {
  
public:
    
  CondAlgX (const std::string& name, ISvcLocator* pSvcLocator);
  virtual ~CondAlgX();
  
  virtual StatusCode initialize() override;
  virtual StatusCode execute() override;

private:
  
  SG::ReadHandleKey<xAOD::EventInfo> m_evt {this,"EvtInfo", "EventInfo", "EventInfo name"};

  SG::WriteCondHandleKey<CondDataObj> m_wchk {this, "Key_CH", "X2", "cond handle key"};

  Gaudi::Property<std::string> m_dbKey {this, "Key_DB", "X2", "explicit dbKey for cond handle"};

  ServiceHandle<IASCIICondDbSvc> m_cds;


};

#endif
