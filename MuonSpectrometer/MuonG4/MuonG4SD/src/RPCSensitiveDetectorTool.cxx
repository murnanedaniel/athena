/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#include "RPCSensitiveDetectorTool.h"
#include "RPCSensitiveDetector.h"

RPCSensitiveDetectorTool::RPCSensitiveDetectorTool(const std::string& type, const std::string& name, const IInterface* parent)
  : SensitiveDetectorBase( type , name , parent )
{
  declareInterface<ISensitiveDetector>(this);
}

G4VSensitiveDetector* RPCSensitiveDetectorTool::makeSD()
{
  ATH_MSG_INFO( "Initializing SD" );
  return new RPCSensitiveDetector(name(), m_outputCollectionNames[0]);
}
