/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "SCT_GeoModel/SCT_Clamp.h"

#include "SCT_GeoModel/SCT_MaterialManager.h"

#include "SCT_GeoModel/SCT_GeometryManager.h"
#include "SCT_GeoModel/SCT_BarrelParameters.h"

#include "GeoModelKernel/GeoTube.h"
#include "GeoModelKernel/GeoLogVol.h"
#include "GeoModelKernel/GeoPhysVol.h"
#include "GeoModelKernel/Units.h"


SCT_Clamp::SCT_Clamp(const std::string & name, int iLayer,
                     double innerRadius,
                     InDetDD::SCT_DetectorManager* detectorManager,
                     SCT_GeometryManager* geometryManager,
                     SCT_MaterialManager* materials)
  : SCT_SharedComponentFactory(name, detectorManager, geometryManager, materials),
    m_iLayer(iLayer), m_innerRadius(innerRadius)
{
  getParameters();
  m_physVolume = SCT_Clamp::build();
}


void 
SCT_Clamp::getParameters()
{
  const SCT_BarrelParameters * parameters = m_geometryManager->barrelParameters();
      
  m_materialName    = parameters->clampMaterial(m_iLayer);
  m_length      = parameters->clampDeltaZ(m_iLayer);
  m_radialWidth = parameters->clampDeltaR(m_iLayer);

  m_outerRadius = m_innerRadius + m_radialWidth;
}

GeoVPhysVol * 
SCT_Clamp::build() 
{
  // Make the clamp. A simple tube.

  const GeoTube * clampShape = new GeoTube(m_innerRadius, m_outerRadius, 0.5 * m_length);
  m_material = m_materials->getMaterialForVolume(m_materialName, clampShape->volume());
  if(!m_material) {m_material = m_materials->getMaterial(m_materialName);}
  const GeoLogVol * clampLog = new GeoLogVol(getName(), clampShape, m_material);

  GeoPhysVol * clamp = new GeoPhysVol(clampLog);

  return clamp;
}
