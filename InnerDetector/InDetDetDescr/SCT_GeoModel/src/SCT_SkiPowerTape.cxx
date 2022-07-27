/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

//
// 1st Feb 2005 D.Naito modified.
// 28th Mar 2005 S.Mima modified.
//
#include "SCT_GeoModel/SCT_SkiPowerTape.h"

#include "SCT_GeoModel/SCT_MaterialManager.h"

#include "SCT_GeoModel/SCT_GeometryManager.h"
#include "SCT_GeoModel/SCT_BarrelParameters.h"

#include "SCT_GeoModel/SCT_Ski.h"
#include "SCT_GeoModel/SCT_PowerTape.h"
#include "SCT_GeoModel/SCT_Module.h"

#include "GeoModelKernel/GeoBox.h"
#include "GeoModelKernel/GeoLogVol.h"
#include "GeoModelKernel/GeoPhysVol.h"
#include "GeoModelKernel/GeoNameTag.h"
#include "GeoModelKernel/GeoTransform.h"
#include "GeoModelKernel/GeoMaterial.h"
#include "GeoModelKernel/Units.h"
#include "GeoModelKernel/GeoDefinitions.h"
#include "GaudiKernel/SystemOfUnits.h"
#include <iostream>
#include <sstream>
#include <cmath>


SCT_SkiPowerTape::SCT_SkiPowerTape(const std::string & name,
                                   const SCT_Ski * ski,
                                   double length,
                                   InDetDD::SCT_DetectorManager* detectorManager,
                                   SCT_GeometryManager* geometryManager,
                                   SCT_MaterialManager* materials) :
  SCT_SharedComponentFactory(name, detectorManager, geometryManager, materials),
  m_length(length), m_ski(ski)
{
  getParameters();
  m_physVolume = SCT_SkiPowerTape::build();
}

void 
SCT_SkiPowerTape::getParameters()
{
  const SCT_BarrelParameters * parameters = m_geometryManager->barrelParameters();

  // Width is the same as the width of a single powertape.
  m_powerTapeThickness        = parameters->powerTapeThickness();
  m_width                     = parameters->powerTapeWidth();
  m_powerTapeStartPointOffset = parameters->powerTapeStartPointOffset();
}


GeoVPhysVol * 
SCT_SkiPowerTape::build()
{
  // The power tapes are stacked going from the nearest interlink to the
  // nearset edge of the dogleg.
  // 
  // ---              ---
  // ------        ------
  // ---------  ---------


  int nHalfModules = m_ski->modulesPerSki()/2;

  // Thickness of the stack. 
  m_thickness = nHalfModules * m_powerTapeThickness;


  // This is a volume containing all the power tapes.
  const GeoBox * skiPowerTapeShape = new GeoBox(0.5*m_thickness, 0.5*m_width, 0.5*m_length);
  const GeoLogVol *skiPowerTapeLog = 
    new GeoLogVol(getName(), skiPowerTapeShape, m_materials->gasMaterial());
  GeoPhysVol * skiPowerTape = new GeoPhysVol(skiPowerTapeLog);

  // Loop over modules in ski as we need their z positions.
  for (int iModule = 0; iModule < m_ski->modulesPerSki(); iModule++) {
    
    // Position PowerTapes
    // These run from the nearest interlink to the edge of the dogleg
    double tapeLength, tapeMid, tapeStart, tapeEnd;
    
    // nPos is used to stack the power tapes. Modules closest to interlink are
    // furthest from support. The positive and negative z positions are 
    // syGaudi::Units::mmetric. nPos = 5,4,3,2,1,0,0,1,2,3,4,5 for the 12 modules.
    int nPos;

    // test sign of zpos to determine whether the tape runs to the
    // high z end or low z end.

    if (m_ski->zPos(iModule) > 0) {

      // Tape runs from high z end to edge of dogleg. NB width of dogleg is in z-direction
      tapeStart  = m_ski->zPos(iModule) + m_ski->coolingBlockOffsetZ() + m_powerTapeStartPointOffset;
      tapeEnd    = 0.5 * m_length;
      
      // For 12 module: iModule = 6, 7, 8, ..., 11 --> nPos = 0, 1, 2, ..., 5
      nPos = iModule - nHalfModules;       

    } else {

      // Tape runs from low z end to edge of dogleg. NB width of dogleg is in z-direction
      tapeStart  =  m_ski->zPos(iModule) + m_ski->coolingBlockOffsetZ() - m_powerTapeStartPointOffset;
      tapeEnd    =  -0.5 * m_length;
   
      // For 12 module: iModule = 0, 1, 2, ..., 5 --> nPos = 5, 4, 3, ..., 0
      nPos = nHalfModules - iModule - 1; 

    }

    tapeLength = std::abs(tapeEnd - tapeStart);
    tapeMid = 0.5 * (tapeEnd + tapeStart);


    // Make sure that first half are negative and secand half are positive.
    // Checking that nPos is >= 0 is equivalent to this.
    if (nPos < 0) {  
      std::cout << "SCT_SkiPowerTape: Module position inconsistent with assumption that\n" 
                << "                        first half has z < 0 and second half has z > 0" 
                << std::endl;
    }

    // Create the tape

    // Label tape with M# at end of string
    SCT_PowerTape powerTape(getName()+"PowerTapeM"+intToString(iModule + 1), tapeLength,
                            m_detectorManager, m_geometryManager, m_materials);
    
    // Calculate x position of tape. This will depend on the module number.
    // The reference point is the middle of the stack.
    double xTapePos = - 0.5 * m_thickness + 
                      (nPos + 0.5) * powerTape.thickness();
    double yTapePos = 0;

    
    // Position the tape
    skiPowerTape->add(new GeoTransform(GeoTrf::Translate3D(xTapePos, yTapePos, tapeMid)));
    skiPowerTape->add(powerTape.getVolume());

  }    

  return skiPowerTape;
}
