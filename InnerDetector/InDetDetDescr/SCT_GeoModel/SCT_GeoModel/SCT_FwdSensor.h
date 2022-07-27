/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef SCT_GEOMODEL_SCT_FWDSENSOR_H
#define SCT_GEOMODEL_SCT_FWDSENSOR_H


#include "SCT_GeoModel/SCT_ComponentFactory.h"

#include <atomic>
#include <string>

// sensor types
/*
#define SENSOR_NT 0
#define SENSOR_FT 1
#define SENSOR_NB 2
#define SENSOR_FB 3
*/

#define SENSOR_BOTTOM 0
#define SENSOR_TOP  1


class GeoMaterial;
class GeoVPhysVol;
class GeoPhysVol;

namespace InDetDD{class SiDetectorDesign;}

class SCT_FwdSensor : public SCT_UniqueComponentFactory
{
public:
  SCT_FwdSensor(const std::string & name, int ringType,
                InDetDD::SCT_DetectorManager* detectorManager,
                SCT_GeometryManager* geometryManager,
                SCT_MaterialManager* materials);
  ~SCT_FwdSensor();

  // Ring type
  int ringType() const {return m_ringType;}


  // Inner radius of physical sensor
  double innerRadius() const {return m_innerRadius;}

  // Outer radius of physical sensor
  double outerRadius() const {return m_outerRadius;}

  // Width of silicon at inner radius
  double innerWidth() const {return m_innerWidth;}
  
  // Width of silicon at outer radius
  double outerWidth() const {return m_outerWidth;}

  // Thickness of sensor
  double thickness() const {return  m_thickness;}

  // Length of physical sensor
  double length() const {return m_outerRadius - m_innerRadius;} 

  // Radius of center of sensor.
  double centerRadius() const {return m_sensorRadius;}
  double sensorOffset() const {return m_sensorOffset;}

  
  double outerWidthN() const {return m_outerWidthN;}
  double innerWidthN() const {return m_innerWidthN;}
  double lengthN()     const {return m_lengthN;}
  double thicknessN()  const {return m_thicknessN;}
  
  double outerWidthF() const {return m_outerWidthF;}
  double innerWidthF() const {return m_innerWidthF;}
  double lengthF()     const {return m_lengthF;}
  double thicknessF()  const {return m_thicknessF;}

  
  virtual GeoVPhysVol * build(SCT_Identifier id);

  GeoPhysVol * getInactive() {return m_inactive;}

private:
  void getParameters();
  virtual const GeoLogVol * preBuild();
  void makeDesign(); 

  int m_ringType{};

  double m_thickness = 0.0;     // Thickness of sensor
  double m_outerRadius = 0.0;   // Outer radius of physical sensor
  double m_innerRadius = 0.0;   // Outer radius of physical sensor
  double m_innerWidth = 0.0;    // Width of sensor at inner radius
  double m_outerWidth = 0.0;    // Width of sensor at outer radius


  double m_outerWidthN = 0.0; // Outer physical width of near sensor.
  double m_innerWidthN = 0.0; // Inner width
  double m_lengthN = 0.0;     // Length (radially).
  double m_thicknessN = 0.0;
  double m_radiusN = 0.0;
  double m_activeHalfLengthNear = 0.0; // Active half length
	  
  double m_outerWidthF = 0.0; // Outer physical width of Far sensor.
  double m_innerWidthF = 0.0; // Inner width
  double m_lengthF = 0.0;     // Length (radially).
  double m_thicknessF = 0.0;
  double m_radiusF = 0.0;
  double m_activeHalfLengthFar = 0.0; // Active half length
  
  double m_sensorRadius = 0.0; // Physics center (of double wafer)

  double m_sensorOffset = 0.0; // Offset of det element center.
                         // Differs from  m_sensorRadius for truncated middle detector

  const GeoMaterial * m_materialSensor = nullptr;
  const GeoMaterial * m_materialGlass = nullptr;

  GeoPhysVol * m_inactive = nullptr;
  
  const InDetDD::SiDetectorDesign * m_design = nullptr;

  mutable std::atomic_bool m_noElementWarning;
};

#endif // SCT_GEOMODEL_SCT_FWDSENSOR_H



