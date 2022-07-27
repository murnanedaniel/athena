/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#ifndef SCT_GEOMODEL_SCT_SKIPOWERTAPE_H
#define SCT_GEOMODEL_SCT_SKIPOWERTAPE_H

#include "SCT_GeoModel/SCT_ComponentFactory.h"

#include <string>

class SCT_Ski;

class SCT_SkiPowerTape : public  SCT_SharedComponentFactory
{

public:

  SCT_SkiPowerTape(const std::string & name,
		   const SCT_Ski * ski,
		   double length,
                   InDetDD::SCT_DetectorManager* detectorManager,
                   SCT_GeometryManager* geometryManager,
                   SCT_MaterialManager* materials);

  //
  // Methods to return basic and derived parameters. 
  //
  double thickness() const {return m_thickness;}
  double width()     const {return m_width;}
  double length()    const {return m_length;}
  double powerTapeThickness() const {return m_powerTapeThickness;}

private:
  void getParameters();
  virtual GeoVPhysVol * build();

  // Basic and Derived parameters  
  double m_length;
  double m_width = 0.0;
  double m_thickness = 0.0;
  double m_powerTapeThickness = 0.0;

  // Child detector elements
  const SCT_Ski * m_ski;

  // 15th Aug 2005 S.Mima modified.
  // Power Tape Z Offset from Bracket 
  double m_powerTapeStartPointOffset = 0.0;   
};

#endif // SCT_GEOMODEL_SCT_SKIPOWERTAPE_H

