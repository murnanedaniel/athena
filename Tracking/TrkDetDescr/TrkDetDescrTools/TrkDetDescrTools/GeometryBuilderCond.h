/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
// GeometryBuilderCond.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef TRKDETDESCRTOOLS_GEOMETRYBUILDERCOND_H
#define TRKDETDESCRTOOLS_GEOMETRYBUILDERCOND_H

// Amg
#include "GeoPrimitives/GeoPrimitives.h"
// Trk
#include "TrkDetDescrInterfaces/IGeometryBuilderCond.h"
#include "TrkDetDescrUtils/GeometrySignature.h"
#include "TrkGeometry/Material.h"
#include "TrkGeometry/TrackingVolumeManipulator.h"
// Gaudi & Athena
#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/ToolHandle.h"

#ifdef TRKDETDESCR_MEMUSAGE
#include "TrkDetDescrUtils/MemoryLogger.h"
#endif

namespace Trk {

class TrackingGeometry;
class TrackingVolume;
class ITrackingVolumeBuilder;
class ITrackingVolumeHelper;
class ITrackingVolumeArrayCreator;

/** @class GeometryBuilderCond

  The Trk::TrackingGeometry Builder for ATLAS Geometry

  It retrieves Trk::TrackingGeometry builders for the subdetectors and joins them together
  to a single Trk::TrackingGeometry.

  @author Andreas.Salzburger@cern.ch
  @author Christos Anastopoulos MT fixes
 */

class GeometryBuilderCond
  : public AthAlgTool
  , public TrackingVolumeManipulator
  , virtual public IGeometryBuilderCond
{

public:
  /** Constructor */
  GeometryBuilderCond(const std::string&, const std::string&, const IInterface*);

  /** Destructor */
  virtual ~GeometryBuilderCond();

  /** AlgTool initialize method */
  virtual StatusCode initialize() override;

  /**
   * TrackingGeometry Interface method - optionally a pointer to Bounds
   * Interface marked as not thread safe
   */
  virtual std::unique_ptr<Trk::TrackingGeometry> trackingGeometry(
    const EventContext& ctx,
    Trk::TrackingVolume* tVol,
    SG::WriteCondHandle<TrackingGeometry>& whandle) const override;

  /** The unique signature */
  virtual GeometrySignature geometrySignature() const override { return Trk::Global; }

private:
  /** TrackingGeometry for ATLAS setup */
  std::unique_ptr<Trk::TrackingGeometry> atlasTrackingGeometry
  (const EventContext& ctx, SG::WriteCondHandle<TrackingGeometry>& whandle) const;

#ifdef TRKDETDESCR_MEMUSAGE
  MemoryLogger m_memoryLogger; //!< in case the memory is logged
#endif

  bool m_createWorld;    //!< Boolean Switch to create World manually
  int m_navigationLevel; //!< NavigationLevel

  std::vector<double> m_worldDimension;          //!< The dimensions of the manually created world
  std::vector<double> m_worldMaterialProperties; //!< The material properties of the created world
  Material m_worldMaterial;                      //!< the world material

  // -------------------------- Tools for geometry building ------------------------------------------------------ //

  ToolHandle<ITrackingVolumeArrayCreator> m_trackingVolumeArrayCreator; //!< Helper Tool to create TrackingVolume Arrays

  ToolHandle<ITrackingVolumeHelper> m_trackingVolumeHelper; //!< Helper Tool to create TrackingVolumes

  ToolHandle<IGeometryBuilderCond> m_inDetGeometryBuilderCond; //!< GeometryBuilderCond for the InnerDetector

  bool m_caloGeometry;                                        //!< switch on TrackingGeometry for the Calorimeters
  ToolHandle<IGeometryBuilderCond> m_caloGeometryBuilderCond; //!< GeometryBuilderCond for the Calorimeters

  bool m_hgtdGeometry;                                        //! switch on TrackingGeometry for HGTD
  ToolHandle<IGeometryBuilderCond> m_hgtdGeometryBuilderCond; //!< GeometryBuilder for the HGTD

  bool m_muonGeometry;                                        //!< GeometryBuilderCond for the Muon System
  ToolHandle<IGeometryBuilderCond> m_muonGeometryBuilderCond; //!< GeometryBuilderCond for the Muon System

  bool m_compactify;        //!< optimize event memory usage: register all surfaces with TG
  bool m_synchronizeLayers; //!< synchronize contained layer dimensions to volumes
};

} // end of namespace

#endif // TRKDETDESCRTOOLS_GEOMETRYBUILDERCOND_H

