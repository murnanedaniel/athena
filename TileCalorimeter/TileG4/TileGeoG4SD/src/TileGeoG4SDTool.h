/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

//************************************************************
//
// Class TileGeoG4SDTool
// AthTool class for holding the Tile sensitive detector
//
//************************************************************

#ifndef TILEGEOG4SD_TILEGEOG4SDTOOL_H
#define TILEGEOG4SD_TILEGEOG4SDTOOL_H

// Base class header
#include "G4AtlasTools/SensitiveDetectorBase.h"

// Package headers
#include "TileGeoG4SD/TileSDOptions.h"

// STL headers
#include <string>

class G4VSensitiveDetector;

class TileGeoG4SDTool : public SensitiveDetectorBase
{
public:
    TileGeoG4SDTool (const std::string& type, const std::string& name, const IInterface *parent);
    // Default destructor is just fine

    /** End of an athena event */
    StatusCode Gather() override final; //FIXME would be good to be able to avoid this.

protected:
    // Make me an SD!
    G4VSensitiveDetector* makeSD() override final;

private:
    // Options for the SD configuration
    TileSDOptions m_options;

};

#endif
