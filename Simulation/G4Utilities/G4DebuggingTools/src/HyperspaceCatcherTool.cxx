/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#include "CxxUtils/make_unique.h"
#include "HyperspaceCatcherTool.h"

namespace G4UA{ 
  
  HyperspaceCatcherTool::HyperspaceCatcherTool(const std::string& type, const std::string& name,const IInterface* parent):
    
    ActionToolBase<HyperspaceCatcher>(type, name, parent),m_config() {
    
    declareProperty("TreatmentLevel",m_config.treatmentLevel);
    declareProperty("KillAfter",m_config.killAfter);
    
  }
  std::unique_ptr<HyperspaceCatcher>  HyperspaceCatcherTool::makeAction(){
    ATH_MSG_DEBUG("makeAction");
    auto action = CxxUtils::make_unique<HyperspaceCatcher>(m_config);
    return std::move(action);
  }
  StatusCode HyperspaceCatcherTool::queryInterface(const InterfaceID& riid, void** ppvIf){
    
    if(riid == IBeginRunActionTool::interfaceID()) {
      *ppvIf = (IBeginRunActionTool*) this;
      addRef();
      return StatusCode::SUCCESS;
    } if(riid == ISteppingActionTool::interfaceID()) {
      *ppvIf = (ISteppingActionTool*) this;
      addRef();
      return StatusCode::SUCCESS;
    } return ActionToolBase<HyperspaceCatcher>::queryInterface(riid, ppvIf);
  }
  
} // namespace G4UA 
