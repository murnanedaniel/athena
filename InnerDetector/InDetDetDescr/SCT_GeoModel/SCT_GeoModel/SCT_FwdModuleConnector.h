/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#ifndef SCT_GEOMODEL_SCT_FWDMODULECONNECTOR_H
#define SCT_GEOMODEL_SCT_FWDMODULECONNECTOR_H

#include "SCT_GeoModel/SCT_ComponentFactory.h"

class GeoMaterial;

class SCT_FwdModuleConnector : public SCT_SharedComponentFactory
{

public:
  SCT_FwdModuleConnector(const std::string & name, int ringType,
                         InDetDD::SCT_DetectorManager* detectorManager,
                         SCT_GeometryManager* geometryManager,
                         SCT_MaterialManager* materials);

  //
  // Methods to return basic and derived parameters. 
  //
  const GeoMaterial * material() const {return m_material;}
  double deltaR() const {return m_deltaR;} 
  double rphi() const {return m_rphi;} 
  double thickness()   const {return m_thickness;}
 
 
private:
  void getParameters();
  virtual GeoVPhysVol * build();

  int m_ringType;

  double m_deltaR = 0.0;
  double m_rphi = 0.0;
  double m_thickness = 0.0;
  const GeoMaterial * m_material = nullptr;
  std::string m_materialName;

};

#endif // SCT_GEOMODEL_SCT_FWDMODULECONNECTOR_H

