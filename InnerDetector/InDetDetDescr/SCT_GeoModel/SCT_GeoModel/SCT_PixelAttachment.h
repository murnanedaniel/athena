/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#ifndef SCT_GEOMODEL_SCT_PIXELATTACHMENT_H
#define SCT_GEOMODEL_SCT_PIXELATTACHMENT_H

#include "SCT_GeoModel/SCT_ComponentFactory.h"

#include <string>

class GeoMaterial;

class SCT_PixelAttachment : public SCT_SharedComponentFactory
{

public:
  SCT_PixelAttachment(const std::string & name,
                      InDetDD::SCT_DetectorManager* detectorManager,
                      SCT_GeometryManager* geometryManager,
                      SCT_MaterialManager* materials);

public:
  const GeoMaterial * material() const {return m_material;}
  double innerRadius() const {return m_innerRadius;} 
  double outerRadius() const {return m_outerRadius;} 
  double zMax() const {return m_zMax;}
  double zMin() const {return m_zMin;}

  double length() const {return m_length;} 
  double radialThickness() const {return m_outerRadius - m_innerRadius;}
  double zPosition() const {return 0.5*(m_zMin+m_zMax);}

 
private:

  void getParameters();
  virtual GeoVPhysVol * build();

  const GeoMaterial * m_material = nullptr;

  double m_innerRadius = 0.0;
  double m_outerRadius = 0.0;
  double m_zMin = 0.0;
  double m_zMax = 0.0;

  double m_length = 0.0;
};

#endif // SCT_GEOMODEL_SCT_PIXELATTACHMENT_H

