/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef G4COSMICFILTER_G4UA__G4COSMICANDFILTERTOOL_H 
#define G4COSMICFILTER_G4UA__G4COSMICANDFILTERTOOL_H 
#include "G4AtlasInterfaces/IEndEventActionTool.h"
#include "G4AtlasTools/ActionToolBase.h"
#include "G4CosmicFilter/G4CosmicAndFilter.h"

namespace G4UA{ 
  
  class G4CosmicAndFilterTool: 
  public ActionToolBaseReport<G4CosmicAndFilter>,
    public IEndEventActionTool
    {
      
    public:

      G4CosmicAndFilterTool(const std::string& type, const std::string& name,const IInterface* parent);
      virtual IEndEventAction* getEndEventAction() override final 

      { return static_cast<IEndEventAction*>( getAction() ); }

      virtual StatusCode queryInterface(const InterfaceID& riid, void** ppvInterface) override;
      virtual StatusCode finalize() override;

    protected:
      virtual std::unique_ptr<G4CosmicAndFilter> makeAction() override final;

    private:
      G4CosmicAndFilter::Config m_config;
    }; // class G4CosmicAndFilterTool
  
  
} // namespace G4UA 
#endif
