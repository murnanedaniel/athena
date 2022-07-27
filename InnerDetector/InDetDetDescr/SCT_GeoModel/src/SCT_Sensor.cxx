/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "SCT_GeoModel/SCT_Sensor.h"

#include "SCT_GeoModel/SCT_GeometryManager.h"
#include "SCT_GeoModel/SCT_MaterialManager.h"

#include "SCT_GeoModel/SCT_BarrelModuleParameters.h"

#include "GeoModelKernel/GeoBox.h"
#include "GeoModelKernel/GeoLogVol.h"
#include "GeoModelKernel/GeoFullPhysVol.h"
#include "GeoModelKernel/GeoMaterial.h"

#include "SCT_ReadoutGeometry/SCT_DetectorManager.h"
#include "SCT_ReadoutGeometry/SCT_BarrelModuleSideDesign.h"
#include "InDetReadoutGeometry/SiDetectorElement.h"
#include "ReadoutGeometryBase/InDetDD_Defs.h"
#include "ReadoutGeometryBase/SiCommonItems.h"

#include "GaudiKernel/SystemOfUnits.h"

using namespace InDetDD;

SCT_Sensor::SCT_Sensor(const std::string & name,
                       InDetDD::SCT_DetectorManager* detectorManager,
                       SCT_GeometryManager* geometryManager,
                       SCT_MaterialManager* materials)
  : SCT_UniqueComponentFactory(name, detectorManager, geometryManager, materials),
    m_noElementWarning{true}
{
  getParameters();
  m_logVolume = SCT_Sensor::preBuild();
}


void
SCT_Sensor::getParameters()
{
  
  const SCT_BarrelModuleParameters * parameters = m_geometryManager->barrelModuleParameters();
  m_material  = m_materials->getMaterial(parameters->sensorMaterial());
  m_thickness = parameters->sensorThickness();
  m_length    = 0;
  if (parameters->sensorNumWafers() == 2) {
    m_length = parameters->sensorLength() + parameters->sensorDistCenterToCenter();
  } else {
    m_length = parameters->sensorLength();
  }
  m_width     = parameters->sensorWidth();

}

const GeoLogVol * 
SCT_Sensor::preBuild()
{

  // Build the sensor. Just a simple box.
  const GeoBox * sensorShape = new GeoBox(0.5*m_thickness, 0.5*m_width, 0.5*m_length);
  GeoLogVol * sensorLog = new GeoLogVol(getName(), sensorShape, m_material);

  // Make the moduleside design for this sensor
  makeDesign();

  return sensorLog;
}


void
SCT_Sensor::makeDesign()
{
  //SiDetectorDesign::Axis etaAxis   = SiDetectorDesign::zAxis;
  //SiDetectorDesign::Axis phiAxis   = SiDetectorDesign::yAxis;
  //SiDetectorDesign::Axis depthAxis = SiDetectorDesign::xAxis;

  const SCT_BarrelModuleParameters * parameters = m_geometryManager->barrelModuleParameters();
  
  double stripPitch      = parameters->sensorStripPitch();
  double stripLength     = parameters->sensorStripLength();
  int crystals           = parameters->sensorNumWafers();
  double distCenterToCenter =  parameters->sensorDistCenterToCenter();

  //int cells              = parameters->sensorNumReadoutStrips();
  //int diodes             = parameters->sensorNumStrips();
  //int shift              = parameters->sensorStripShift(); 

  // For consistency with older version we keep shift = 0.
  int cells              = parameters->sensorNumReadoutStrips();
  int diodes             = cells;
  int shift              = 0; 

  double totalDeadLength = distCenterToCenter - stripLength;
    
  double xEtaStripPatternCenter = 0;
  double xPhiStripPatternCenter = 0;

  // Readout direction is in same direction as local phi direction
  bool swapStripReadout = false;

  // The readout side is at the +ve depth direction
  int readoutSide = +1;

  // m_design will be owned and deleted by SCT_DetectorManager
  std::unique_ptr<SCT_BarrelModuleSideDesign> design = std::make_unique<SCT_BarrelModuleSideDesign>(m_thickness,
                                            crystals,
                                            diodes,
                                            cells,
                                            shift,
                                            swapStripReadout,
                                            InDetDD::holes,
                                            stripPitch,
                                            stripLength,
                                            xEtaStripPatternCenter,
                                            xPhiStripPatternCenter,
                                            totalDeadLength,
                                            readoutSide);
  m_design = m_detectorManager->addDesign(std::move(design));

  //
  // Flags to signal if axis can be swapped.
  // For rectangular detector these are all true.
  // This is the default and no action is required. 
  // Can force axes not to be swapped by setting to false.
  // 
  // bool phiSyGaudi::Units::mmetric = true;
  // bool etaSyGaudi::Units::mmetric = true; 
  // bool depthSyGaudi::Units::mmetric = true;
  // m_design->setSyGaudi::Units::mmetry(phiSyGaudi::Units::mmetric, etaSyGaudi::Units::mmetric, depthSyGaudi::Units::mmetric,
  //
}



GeoVPhysVol * 
SCT_Sensor::build(SCT_Identifier id)
{
  GeoFullPhysVol * sensor = new GeoFullPhysVol(m_logVolume); 
  
  // Make detector element and add to collection
  // Only do so if we have a valid id helper.

  //id.print(); // for debugging only

  const SiCommonItems* commonItems =  m_geometryManager->commonItems();

  if (commonItems->getIdHelper()) {

    SiDetectorElement * detElement;

    // detElement will be owned by SCT_DetectorManager
    // and will be deleted in destructor of SiDetectorElementCollection in SCT_DetectorManager
    detElement =  new SiDetectorElement(id.getWaferId(), 
                                        m_design, 
                                        sensor,  
                                        commonItems);
    
    // Add the detector element.
    m_detectorManager->addDetectorElement(detElement);

  } else {
    if (m_noElementWarning) {
      std::cout << "WARNING!!!!: No SCT id helper and so no elements being produced." << std::endl;
      m_noElementWarning = false;
    }
  }
  return sensor;
}
