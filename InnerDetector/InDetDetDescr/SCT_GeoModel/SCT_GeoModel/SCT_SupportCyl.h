/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#ifndef SCT_GEOMODEL_SCT_SUPPORTCYL_H
#define SCT_GEOMODEL_SCT_SUPPORTCYL_H

#include "SCT_GeoModel/SCT_ComponentFactory.h"

#include <string>

class GeoMaterial;

class SCT_SupportCyl : public SCT_SharedComponentFactory
{

public:
  SCT_SupportCyl(const std::string & name, int iLayer, double length,
                 InDetDD::SCT_DetectorManager* detectorManager,
                 SCT_GeometryManager* geometryManager,
                 SCT_MaterialManager* materials);

public:
  const GeoMaterial * material() const {return m_material;}
  double innerRadius() const {return m_innerRadius;} 
  double outerRadius() const {return m_outerRadius;} 
  double length() const {return m_length;} 

 
private:
  void getParameters();
  virtual GeoVPhysVol * build();

  int m_iLayer;

  const GeoMaterial * m_material = nullptr;
  std::string m_materialName;
  double m_innerRadius = 0.0;
  double m_outerRadius = 0.0;
  double m_length;
};

#endif // SCT_GEOMODEL_SCT_SUPPORTCYL_H

