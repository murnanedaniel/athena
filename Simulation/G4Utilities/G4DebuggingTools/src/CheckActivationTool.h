/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef G4DEBUGGINGTOOLS_G4UA__CHECKACTIVATIONTOOL_H 
#define G4DEBUGGINGTOOLS_G4UA__CHECKACTIVATIONTOOL_H 
#include "G4AtlasInterfaces/IBeginEventActionTool.h"
#include "G4AtlasTools/ActionToolBase.h"
#include "CheckActivation.h"

namespace G4UA{ 

/// @class CheckActivationTool
/// @brief Tool which manages the CheckActivation
///
/// @author Andrea Di Simone
///
  

  class CheckActivationTool: 
  public ActionToolBase<CheckActivation>,
    public IBeginEventActionTool
    {
      
    public:
      /// Standard constructor
      CheckActivationTool(const std::string& type, const std::string& name,const IInterface* parent);
      /// Retrieve the BoE action
      virtual IBeginEventAction* getBeginEventAction() override final 
      { return static_cast<IBeginEventAction*>( getAction() ); }
      /// Gaudi interface management
      virtual StatusCode queryInterface(const InterfaceID& riid, void** ppvInterface) override;
    protected:
      /// Create an action for this thread
      virtual std::unique_ptr<CheckActivation> makeAction() override final;
    private:
    }; // class CheckActivationTool
  
  
} // namespace G4UA 
#endif
