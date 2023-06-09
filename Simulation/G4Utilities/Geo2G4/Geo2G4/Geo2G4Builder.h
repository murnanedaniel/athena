/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef Geo2G4Builder_H
#define Geo2G4Builder_H

// main builder to create/position all volumes described in a GeoModel Tree

// GeoVPhysVol
#include "Geo2G4/VolumeBuilder.h"
#include "GeoModelKernel/GeoVPhysVol.h"
#include "Geo2G4/LogicalVolume.h"
//#include "Geo2G4/GenericVolumeBuilder.h"

// Typedef
#include "GeoModelUtilities/GeoBorderSurfaceContainer.h"

#include "AthenaKernel/MsgStreamMember.h"

// STL includes
#include <string>
#include <vector>

class GeoMaterial;
class StoreGateSvc;

class Geo2G4Builder {

public:
  // Constructor:
  Geo2G4Builder(std::string detectorName);
  // Destructor:
  ~Geo2G4Builder() {;}

  // Build method - geometry
  LogicalVolume*        BuildTree();

  // Build method - optical surfaces
  void BuildOpticalSurfaces(const GeoBorderSurfaceContainer* surface_container,
                            const OpticalVolumesMap* optical_volumes);

  // Access volume builder:
  VolumeBuilder*        GetVolumeBuilder(std::string);

  HepGeom::Transform3D& GetDetectorTransform() {return motherTransform;}
  /// Log a message using the Athena controlled logging system
  MsgStream& msg( MSG::Level lvl ) const { return m_msg << lvl; }
  /// Check whether the logging system is active at the provided verbosity level
  bool msgLvl( MSG::Level lvl ) const { return m_msg.get().level() <= lvl; }

private:

  // GeoVDetectorManager* theDetectorElement;
  std::string m_detectorName;
  HepGeom::Transform3D motherTransform;
  std::vector<PVConstLink> m_treeTops;
  VolumeBuilder *theBuilder;

  // std::Air in the case when top boolean envelope has to be built
  GeoMaterial* m_matAir;
  StoreGateSvc* m_pDetStore;

  /// Private message stream member
  mutable Athena::MsgStreamMember m_msg;
};

#endif
