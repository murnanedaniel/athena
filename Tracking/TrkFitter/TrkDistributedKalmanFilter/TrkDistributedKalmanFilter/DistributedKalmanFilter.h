/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

//////////////////////////////////////////////////////////////////
// DistributedKalmanFilter.h
//   Header file for class DistributedKalmanFilter
///////////////////////////////////////////////////////////////////
// (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////
// Author: Dmitry Emeliyanov, RAL
// D.Emeliyanov@rl.ac.uk
///////////////////////////////////////////////////////////////////
#ifndef __TRK_DISTRIBUTED_KALMAN_FILTER__
#define __TRK_DISTRIBUTED_KALMAN_FILTER__

#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/EventContext.h"
#include "GaudiKernel/ToolHandle.h"
#include "TrkEventPrimitives/ParticleHypothesis.h"
#include "TrkFitterInterfaces/ITrackFitter.h"
#include "TrkFitterUtils/FitterTypes.h"

#include "MagFieldConditions/AtlasFieldCacheCondObj.h"
#include "StoreGate/ReadCondHandleKey.h"
#include "TrkToolInterfaces/IRIO_OnTrackCreator.h"
#include "TrkExInterfaces/IExtrapolator.h"

#include <memory>
#include <vector>
class IAlgTool;
class AtlasDetectorID;

namespace Trk {

class Track;            //!> ATLAS standard track class
class TrkBaseNode;      //!> Base class for track filtering nodes
class TrkTrackState;    //!> Simple track state: 5 parameters + covariance
class TrkPlanarSurface; //!> Simple class for planar surface description
class TrackStateOnSurface;
class MeasuredPerigee;

class DistributedKalmanFilter
  : virtual public ITrackFitter
  , public AthAlgTool
{
public:
  // standard AlgTool methods
  DistributedKalmanFilter(const std::string&,
                          const std::string&,
                          const IInterface*);
  virtual ~DistributedKalmanFilter();

  // standard Athena methods
  virtual StatusCode initialize() override;
  virtual StatusCode finalize() override;

  // Filter settings:
  //			RunOutlierRemoval - use logic to remove bad hits
  //			ParticleHypothesis   - allow for multiple scattering and
  //energy loss?


  virtual std::unique_ptr<Track> fit(
    const EventContext& ctx,
    const Track&,
    const RunOutlierRemoval runOutlier = false,
    const ParticleHypothesis matEffects = Trk::nonInteracting) const override;

  virtual std::unique_ptr<Track> fit(const EventContext&,
                                     const Track&,
                                     const PrepRawDataSet&,
                                     const RunOutlierRemoval,
                                     const ParticleHypothesis) const override
  {
    return nullptr;
  }

  // fit a set of PrepRawData objects
  virtual std::unique_ptr<Track> fit(
    const EventContext& ctx,
    const PrepRawDataSet&,
    const TrackParameters& estimatedParametersNearOrigine,
    const RunOutlierRemoval runOutlier = false,
    const ParticleHypothesis matEffects = Trk::nonInteracting) const override;

  virtual std::unique_ptr<Track> fit(
    const EventContext& ctx,
    const MeasurementSet&,
    const TrackParameters& estim,
    const RunOutlierRemoval runo = false,
    const ParticleHypothesis mat = Trk::nonInteracting) const override;

  virtual std::unique_ptr<Track> fit(
    const EventContext& ctx,
    const Track&,
    const MeasurementSet&,
    const RunOutlierRemoval runOutlier = false,
    const ParticleHypothesis matEffects = Trk::nonInteracting) const override;


  virtual std::unique_ptr<Track> fit(const EventContext&,
                                     const Track&,
                                     const Track&,
                                     const RunOutlierRemoval,
                                     const ParticleHypothesis) const override
  {
    return nullptr;
  }

  // fit a set of RIO_OnTrack objects
  std::unique_ptr<Track> fit(
    const EventContext& ctx,
    const RIO_OnTrackSet&,
    const TrackParameters& estimatedParametersNearOrigine,
    const RunOutlierRemoval runOutlier = false,
    const ParticleHypothesis matEffects = Trk::nonInteracting) const;

private:
  using PVPNodes = std::vector<std::unique_ptr<TrkBaseNode> >;
  using PVPSurfaces = std::vector<std::unique_ptr<TrkPlanarSurface> >;
  using PVPTrackStates = std::vector<std::unique_ptr<TrkTrackState> >;

  ///////////////////////////////////////////////////////////////////
  // Private functions:
  ///////////////////////////////////////////////////////////////////

  std::unique_ptr<Trk::TrkTrackState>
  extrapolate(TrkTrackState*,
              TrkPlanarSurface*,
              TrkPlanarSurface*,
              MagField::AtlasFieldCache& fieldCache) const;

  bool runForwardKalmanFilter(PVPNodes& pvpNodes,
                              PVPTrackStates& pvpTrackStates,
                              TrkTrackState*,
                              MagField::AtlasFieldCache& fieldCache) const;
  static void runSmoother(PVPTrackStates& pvpTrackStates) ;
  static int findOutliers(PVPNodes& pvpNodes, double) ;
  static void calculateLRsolution(PVPNodes& pvpNodes) ;
  TrackStateOnSurface* createTrackStateOnSurface(TrkBaseNode*) const;
  void report();
  void report(char fileName[]);
  void getMagneticField(double[3],
                        double*,
                        MagField::AtlasFieldCache& fieldCache) const;


  double integrate(double Rk[5],
                   TrkPlanarSurface* pSB,
                   TrkPlanarSurface* pSE,
                   double* Rf,
                   MagField::AtlasFieldCache& fieldCache) const;

  ///////////////////////////////////////////////////////////////////
  // Private data:
  ///////////////////////////////////////////////////////////////////

  const AtlasDetectorID* m_idHelper;

  ToolHandle<IRIO_OnTrackCreator> m_ROTcreator{
    this,
    "ROTcreator",
    "Trk::RIO_OnTrackCreator/RIO_OnTrackCreator",
    "ROT Creator Tool"
  };
  ToolHandle<IExtrapolator> m_extrapolator{ this,
                                            "ExtrapolatorTool",
                                            "Trk::Extrapolator/Extrapolator",
                                            "Extrapolator Tool" };

  SG::ReadCondHandleKey<AtlasFieldCacheCondObj> m_fieldCacheCondObjInputKey{
    this,
    "AtlasFieldCacheCondObj",
    "fieldCondObj",
    "Name of the Magnetic Field conditions object key"
  };

  // ME temporary fix
  std::vector<double> m_option_sortingRefPoint;
};

inline void
Trk::DistributedKalmanFilter::getMagneticField(
  double gP[3],
  double* pB,
  MagField::AtlasFieldCache& fieldCache) const
{
  pB[0] = 0.0;
  pB[1] = 0.0;
  pB[2] = 0.0;
  double field[3];
  fieldCache.getField(gP, field); // field is returned in kT
  for (int i = 0; i < 3; i++)
    pB[i] = field[i] / Gaudi::Units::kilogauss; // convert to kG
}

} // end of namespace

#endif
