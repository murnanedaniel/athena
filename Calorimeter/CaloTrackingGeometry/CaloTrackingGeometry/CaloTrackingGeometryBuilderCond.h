/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
// CaloTrackingGeometryBuilderCond.hm (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef CALORIMETER_CALOTRACKINGGEOMETRYBUILDERCOND_H
#define CALORIMETER_CALOTRACKINGGEOMETRYBUILDERCOND_H

#include "CaloIdentifier/CaloCell_ID.h"
// Gaudi
#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"
// Trk
#include "TrkDetDescrInterfaces/IGeometryBuilderCond.h"
#include "TrkDetDescrUtils/GeometrySignature.h"
#include "TrkDetDescrUtils/LayerIndexSampleMap.h"
// EnvelopeDefinitionService
#include "SubDetectorEnvelopes/IEnvelopeDefSvc.h"
//CaloDDM
#include "CaloDetDescr/CaloDetDescrManager.h"
// STL
#include <vector>

#include <CxxUtils/checker_macros.h>
namespace Trk {
class Material;
class Layer;
class MagneticFieldProperties;
class TrackingGeometry;
class ICaloTrackingVolumeBuilder;
class ITrackingVolumeCreator;
class ITrackingVolumeHelper;
class ITrackingVolumeArrayCreator;
class IMagneticFieldTool;
}

class IEnvelopeDefSvc;

namespace Calo {

/** @class CaloTrackingGeometryBuilderCond

  Retrieves the volume builders from Tile and LAr
  and combines the given volumes to a calorimeter tracking
  geometry.

  It also wraps an inner volume if provided though
  the mehtod interface.

  @author Andreas.Salzburger@cern.ch
  */

// Not thread safe due to mutable
// It is handled via the TrackingGeometryCond
// asssuming that the calo tracking geometry
// Calo usually does not change inside a job
class ATLAS_NOT_THREAD_SAFE CaloTrackingGeometryBuilderCond
  : public AthAlgTool
  , virtual public Trk::IGeometryBuilderCond
{

public:
  /** Constructor */
  CaloTrackingGeometryBuilderCond(const std::string&, const std::string&, const IInterface*);

  /** Destructor */
  virtual ~CaloTrackingGeometryBuilderCond();

  /** AlgTool initailize method.*/
  virtual StatusCode initialize() override;

  /** AlgTool finalize method */
  virtual StatusCode finalize() override;

  /** TrackingGeometry Interface method */
  virtual
  std::unique_ptr<Trk::TrackingGeometry> trackingGeometry(
    const EventContext& ctx,
    const Trk::TrackingVolume* innerVol,
    SG::WriteCondHandle<Trk::TrackingGeometry>& whandle) const override;

  /** The unique signature */
  virtual Trk::GeometrySignature geometrySignature() const override { return Trk::Calo; }

private:
  ToolHandle<Trk::ITrackingVolumeArrayCreator>
    m_trackingVolumeArrayCreator; //!< Helper Tool to create TrackingVolume Arrays

  ToolHandle<Trk::ITrackingVolumeHelper> m_trackingVolumeHelper; //!< Helper Tool to create TrackingVolumes

  ToolHandle<Trk::ITrackingVolumeCreator> m_trackingVolumeCreator; //!< Second helper for volume creation

  ToolHandle<Trk::ICaloTrackingVolumeBuilder> m_lArVolumeBuilder; //!< Volume Builder for the Liquid Argon Calorimeter

  ToolHandle<Trk::ICaloTrackingVolumeBuilder> m_tileVolumeBuilder; //!< Volume Builder for the Tile Calorimeter

  SG::ReadCondHandleKey<CaloDetDescrManager> m_caloMgrKey{ this, "CaloDetDescrManager", "CaloDetDescrManager" };

  Trk::Material* m_caloMaterial; //!< Material properties

  double m_caloEnvelope; //!< Envelope cover for Gap Layers
  // enclosing endcap/cylindervolume
  ServiceHandle<IEnvelopeDefSvc> m_enclosingEnvelopeSvc;

  mutable double m_caloDefaultRadius;      //!< the radius if not built from GeoModel
  mutable double m_caloDefaultHalflengthZ; //!< the halflength in z if not built from GeoModel

  bool m_indexStaticLayers; //!< forces robust indexing for layers

  bool m_recordLayerIndexCaloSampleMap;      //!< for deposition methods
  std::string m_layerIndexCaloSampleMapName; //!< name to record it

  bool m_buildMBTS; //!< MBTS like detectors
  // int				                            m_mbstSurfaceShape;		          //!< MBTS
  // like detectors
  mutable std::vector<double> m_mbtsRadii; //!< MBTS like detectors
  std::vector<double> m_mbtsRadiusGap;     //!< MBTS like detectors
  std::vector<int> m_mbtsPhiSegments;      //!< MBTS like detectors
  std::vector<double> m_mbtsPhiGap;        //!< MBTS like detectors
  std::vector<double> m_mbtsPositionZ;     //!< MBTS like detectors
  std::vector<double> m_mbtsStaggeringZ;   //!< MBTS like detectors

  std::string m_entryVolume; //!< name of the Calo entrance
  std::string m_exitVolume;  //!< name of the Calo container

  mutable RZPairVector m_bpCutouts;

  // MBTS layers
  mutable std::vector<const Trk::Layer*>* m_mbtsNegLayers;
  mutable std::vector<const Trk::Layer*>* m_mbtsPosLayers;
  // ToolHandle<ICaloSurfaceHelper>                m_caloSurfaceHelper;

  /** method to establish a link between the LayerIndex and the CaloCell_ID in an associative container */
  void registerInLayerIndexCaloSampleMap(Trk::LayerIndexSampleMap& licsMAp,
                                         std::vector<CaloCell_ID::CaloSample> ccid,
                                         const Trk::TrackingVolume& vol,
                                         int side = 1) const;

  /** method to build enclosed beam pipe volumes */
  std::pair<const Trk::TrackingVolume*, const Trk::TrackingVolume*> createBeamPipeVolumes(float,
                                                                                          float,
                                                                                          const std::string&,
                                                                                          float&) const;

  /** cleanup of material */
  mutable std::map<const Trk::Material*, bool> m_materialGarbage;
  void throwIntoGarbage(const Trk::Material* mat) const;
};

} // end of namespace

inline void
Calo::CaloTrackingGeometryBuilderCond::throwIntoGarbage(const Trk::Material* mat) const
{
  if (mat)
    m_materialGarbage[mat] = true;
}

#endif // CALORIMETER_CALOTRACKINGGEOMETRYBUILDERCOND_H

