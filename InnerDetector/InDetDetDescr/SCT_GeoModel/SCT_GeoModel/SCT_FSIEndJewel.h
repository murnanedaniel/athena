/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#ifndef SCT_GEOMODEL_SCT_FSIENDJEWEL_H
#define SCT_GEOMODEL_SCT_FSIENDJEWEL_H

#include "SCT_GeoModel/SCT_ComponentFactory.h"

#include <string>

class GeoMaterial;

class SCT_FSIEndJewel : public SCT_SharedComponentFactory
{

public:
  SCT_FSIEndJewel(const std::string & name,
                  InDetDD::SCT_DetectorManager* detectorManager,
                  SCT_GeometryManager* geometryManager,
                  SCT_MaterialManager* materials);

public:
  const GeoMaterial * material() const {return m_material;}
  double radialWidth() const {return m_radialWidth;} 
  double rPhiWidth() const {return m_rPhiWidth;} 
  double length() const {return m_length;} 

 
private:
  void getParameters();
  virtual GeoVPhysVol * build();

  const GeoMaterial * m_material = nullptr;
  std::string m_materialName;
  double m_radialWidth = 0.0;
  double m_rPhiWidth = 0.0;
  double m_length = 0.0;
};

#endif // SCT_GEOMODEL_SCT_FSIENDJEWEL_H

