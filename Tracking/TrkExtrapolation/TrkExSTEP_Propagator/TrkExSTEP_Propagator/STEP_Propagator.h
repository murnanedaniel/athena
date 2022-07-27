/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
 */

/////////////////////////////////////////////////////////////////////////////////
//  Header file for class  STEP_Propagator
/////////////////////////////////////////////////////////////////////////////////
// (c) ATLAS Detector software
/////////////////////////////////////////////////////////////////////////////////
// Class for track status propagation through magnetic field
/////////////////////////////////////////////////////////////////////////////////
// Based on RungeKuttaPropagator Version 1.0 21/05/2004 I.Gavrilenko
// Version 0.1 16/2/2005 Esben Lund
// Extended by S.Todorova for propagation to multiple surfaces
/////////////////////////////////////////////////////////////////////////////////

#ifndef STEP_Propagator_H
#define STEP_Propagator_H

#include "AthenaBaseComps/AthAlgTool.h"
#include "AthenaKernel/IAtRndmGenSvc.h"
#include "GaudiKernel/ToolHandle.h"
#include "TrkEventPrimitives/PropDirection.h"
#include "TrkEventPrimitives/SurfaceTypes.h"
#include "TrkExInterfaces/IPropagator.h"
//
#include "TrkExUtils/MaterialInteraction.h"
#include "TrkExUtils/TrackSurfaceIntersection.h"
#include "TrkGeometry/BinnedMaterial.h" //Identified material typedef
#include "TrkMaterialOnTrack/EnergyLoss.h"
#include "TrkParameters/TrackParameters.h" //TrackParameters typedef
//#include "TrkExUtils/ExtrapolationCache.h"
// Amg
#include "EventPrimitives/EventPrimitives.h"
#include "GeoPrimitives/GeoPrimitives.h"
// MagField cache
#include "MagFieldConditions/AtlasFieldCacheCondObj.h"
#include "MagFieldElements/AtlasFieldCache.h"
//
#include "CxxUtils/restrict.h"
#include <list>
#include <vector>

namespace Trk {
class Surface;
class TrackingVolume;
class ITimedMatEffUpdator;
class ScatteringAngles;
class AlignableTrackingVolume;
class ExtrapolationCache;
/**
   @class STEP_Propagator

   STEP (Simultaneous Track and Error Propagation )
   is an algorithm for track parameters propagation through
   magnetic field.

   The algorithm can produce the Jacobian that transports the covariance matrix
   from one set of track parameters at the initial surface to another set of
   track parameters at the destination surface. This is useful for Chi2
   fitting.

   One can choose to perform the transport of the parameters only and omit the transport
   of the associated covariances (propagateParameters).

   The implementation performs the propagation in global coordinates and uses
   Jacobian matrices (see RungeKuttaUtils) for the transformations between the
   global frame and local surface-attached coordinate systems.

   The STEP_Propagator includes material effects in the
   equation of motion and applies corrections to the covariance matrices
   continuously during the parameter transport.  It is designed for the
   propagation  of a particle going through a dense block of material
   (e.g a muon transversing the calorimeter).

   1.The first step of the algorithm is track parameters transformation from
   local presentation for given start surface to global Runge Kutta coordinates.

   2.The second step is propagation through magnetic field with or without jacobian.

   3.Third step is transformation from global Runge Kutta presentation to local
   presentation of given output surface.


   AtaPlane    AtaStraightLine      AtaDisc       AtaCylinder      Perigee
   |               |               |               |              |
   |               |               |               |              |
   V               V               V               V              V
   -----------------------------------------------------------------
   |              Local->Global transformation
   V
   Global position (Runge Kutta presentation)
   |
   |
   Propagation to next surface with or without jacobian
   |
   V              Global->Local transformation
   ----------------------------------------------------------------
   |               |               |               |              |
   |               |               |               |              |
   V               V               V               V              V
   PlaneSurface StraightLineSurface DiscSurface CylinderSurface PerigeeSurface

   For propagation using Runge Kutta method we use global coordinate, direction,
   inverse momentum and Jacobian of transformation. All this parameters we save
   in array P[42].
   /dL0    /dL1    /dPhi   /dThe   /dCM
   X  ->P[0]  dX /   P[ 7]   P[14]   P[21]   P[28]   P[35]
   Y  ->P[1]  dY /   P[ 8]   P[15]   P[22]   P[29]   P[36]
   Z  ->P[2]  dZ /   P[ 9]   P[16]   P[23]   P[30]   P[37]
   Ax ->P[3]  dAx/   P[10]   P[17]   P[24]   P[31]   P[38]
   Ay ->P[4]  dAy/   P[11]   P[18]   P[25]   P[32]   P[39]
   Az ->P[5]  dAz/   P[12]   P[19]   P[26]   P[33]   P[40]
   CM ->P[6]  dCM/   P[13]   P[20]   P[27]   P[34]   P[41]

   where
   in case local presentation

   L0  - first  local coordinate  (surface dependent)
   L1  - second local coordinate  (surface dependent)
   Phi - Azimuthal angle
   The - Polar     angle
   CM  - charge/momentum

   in case global presentation

   X   - global x-coordinate        = surface dependent
   Y   - global y-coordinate        = surface dependent
   Z   - global z-coordinate        = sutface dependent
   Ax  - direction cosine to x-axis = Sin(The)*Cos(Phi)
   Ay  - direction cosine to y-axis = Sin(The)*Sin(Phi)
   Az  - direction cosine to z-axis = Cos(The)
   CM  - charge/momentum            = local CM

   The Runge-Kutta method:

   The equations of motion are solved using an embedded pair of Runge-Kutta formulas.
   This method, Runge-Kutta-Fehlberg, calculates a number of points along the step and
   adds them up in different ways to get two different solutions, of different order, for the
   integration. The higher order solution is used for the propagation and the lower order solution for
   error control. The difference between these solutions is used to estimate the quality of the
   integration (propagation), and to calculate the step size for the next step. If the quality is below
   a given tolerance then the step is rejected and repeated with a shorter step length. This propagator
   uses the TP43 (Tsitouras-Papakostas 4th and 3rd order) Runge-Kutta pair.

   The step size algoritm by L.P.Endresen and J.Myrheim was choosen for its low step rejection and
   effective step size calculation. The low step rejection is achieved by letting the step size
   oscillate around the optimal value instead of repeating steps every time they fall below the
   tolerance level.

   Units are mm, MeV and kiloGauss.

   @author esben.lund@fys.uio.no
**/

class STEP_Propagator final
  : public AthAlgTool
  , virtual public IPropagator
{
  /////////////////////////////////////////////////////////////////////////////////
  // Public methods:
  /////////////////////////////////////////////////////////////////////////////////
public:

  using IPropagator::propagate;

  STEP_Propagator(const std::string&, const std::string&, const IInterface*);

  virtual ~STEP_Propagator();

  /** AlgTool initailize method.*/
  virtual StatusCode initialize() override final;

  /** AlgTool finalize method */
  virtual StatusCode finalize() override final;

  /** Main propagation method NeutralParameters.
   * It is not implemented for STEP propagator.*/
  virtual std::unique_ptr<Trk::NeutralParameters> propagate(
    const Trk::NeutralParameters&,
    const Trk::Surface&,
    Trk::PropDirection,
    const Trk::BoundaryCheck&,
    bool rC = false) const override final;

  /** Propagate parameters and covariance without returning the Jacobian */
  virtual std::unique_ptr<Trk::TrackParameters> propagate(
    const EventContext& ctx,
    const Trk::TrackParameters& trackParameters,
    const Trk::Surface& targetSurface,
    Trk::PropDirection propagationDirection,
    const Trk::BoundaryCheck& boundaryCheck,
    const MagneticFieldProperties& magneticFieldProperties,
    ParticleHypothesis particle,
    bool returnCurv = false,
    const Trk::TrackingVolume* tVol = nullptr) const override final;

  /** Propagate parameters and covariance with search of closest surface */
  virtual std::unique_ptr<Trk::TrackParameters> propagate(
    const EventContext& ctx,
    const Trk::TrackParameters& trackParameters,
    std::vector<Trk::DestSurf>& targetSurfaces,
    Trk::PropDirection propagationDirection,
    const MagneticFieldProperties& magneticFieldProperties,
    ParticleHypothesis particle,
    std::vector<unsigned int>& solutions,
    double& path,
    bool usePathLimit = false,
    bool returnCurv = false,
    const Trk::TrackingVolume* tVol = nullptr) const override final;

  /** Propagate parameters and covariance with search of closest surface
   * time included*/
  virtual std::unique_ptr<Trk::TrackParameters> propagateT(
    const EventContext& ctx,
    const Trk::TrackParameters& trackParameters,
    std::vector<Trk::DestSurf>& targetSurfaces,
    Trk::PropDirection propagationDirection,
    const MagneticFieldProperties& magneticFieldProperties,
    ParticleHypothesis particle,
    std::vector<unsigned int>& solutions,
    Trk::PathLimit& path,
    Trk::TimeLimit& time,
    bool returnCurv,
    const Trk::TrackingVolume* tVol,
    std::vector<Trk::HitInfo>*& hitVector) const override final;

  /** Propagate parameters and covariance with search of closest surface and
   * material collection */
  virtual std::unique_ptr<Trk::TrackParameters> propagateM(
    const EventContext& ctx,
    const Trk::TrackParameters& trackParameters,
    std::vector<Trk::DestSurf>& targetSurfaces,
    Trk::PropDirection propagationDirection,
    const MagneticFieldProperties& magneticFieldProperties,
    ParticleHypothesis particle,
    std::vector<unsigned int>& solutions,
    std::vector<const Trk::TrackStateOnSurface*>*& matstates,
    std::vector<std::pair<std::unique_ptr<Trk::TrackParameters>, int>>*
      intersections,
    double& path,
    bool usePathLimit = false,
    bool returnCurv = false,
    const Trk::TrackingVolume* tVol = nullptr,
    Trk::ExtrapolationCache* = nullptr) const override final;

  /** Propagate parameters and covariance, and return the Jacobian. WARNING:
   * Multiple Scattering is not included in the Jacobian! */
  virtual std::unique_ptr<Trk::TrackParameters> propagate(
    const EventContext& ctx,
    const Trk::TrackParameters& trackParameters,
    const Trk::Surface& targetSurface,
    Trk::PropDirection propagationDirection,
    const Trk::BoundaryCheck& boundaryCheck,
    const MagneticFieldProperties& magneticFieldProperties,
    Trk::TransportJacobian*& jacobian,
    double& pathLimit,
    ParticleHypothesis particle,
    bool returnCurv = false,
    const Trk::TrackingVolume* tVol = nullptr) const override final;

  /** Propagate parameters only */
  virtual std::unique_ptr<Trk::TrackParameters> propagateParameters(
    const EventContext& ctx,
    const Trk::TrackParameters& trackParameters,
    const Trk::Surface& targetSurface,
    Trk::PropDirection propagationDirection,
    const Trk::BoundaryCheck& boundaryCheck,
    const MagneticFieldProperties& magneticFieldProperties,
    ParticleHypothesis particle,
    bool returnCurv = false,
    const Trk::TrackingVolume* tVol = nullptr) const override final;

  /** Propagate parameters and return Jacobian. WARNING: Multiple Scattering is
   * not included in the Jacobian! */
  virtual std::unique_ptr<Trk::TrackParameters> propagateParameters(
    const EventContext& ctx,
    const Trk::TrackParameters& trackParameters,
    const Trk::Surface& targetSurface,
    Trk::PropDirection propagationDirection,
    const Trk::BoundaryCheck& boundaryCheck,
    const MagneticFieldProperties& magneticFieldProperties,
    Trk::TransportJacobian*& jacobian,
    ParticleHypothesis particle,
    bool returnCurv = false,
    const Trk::TrackingVolume* tVol = nullptr) const override final;

  /** Propagate parameters and return path (Similar to propagateParameters */
  virtual const IntersectionSolution* intersect(
    const EventContext& ctx,
    const Trk::TrackParameters& trackParameters,
    const Trk::Surface& targetSurface,
    const Trk::MagneticFieldProperties& magneticFieldProperties,
    ParticleHypothesis particle,
    const Trk::TrackingVolume* tVol = nullptr) const override final;

  /** Intersection and propagation:
   */
  virtual const TrackSurfaceIntersection* intersectSurface(
    const EventContext& ctx,
    const Surface& surface,
    const TrackSurfaceIntersection* trackIntersection,
    const double qOverP,
    const MagneticFieldProperties& mft,
    ParticleHypothesis particle) const override final;

  /** Return a list of positions along the track */
  virtual void globalPositions(
    const EventContext& ctx,
    std::list<Amg::Vector3D>& positionsList,
    const TrackParameters& trackParameters,
    const MagneticFieldProperties& magneticFieldProperties,
    const CylinderBounds& cylinderBounds,
    double maxStepSize,
    ParticleHypothesis particle,
    const Trk::TrackingVolume* tVol = 0) const override final;

  /////////////////////////////////////////////////////////////////////////////////
  // Private methods:
  /////////////////////////////////////////////////////////////////////////////////

  struct Cache
  {
    bool m_detailedElossFlag{ true };
    bool m_solenoid{ false };
    bool m_matPropOK{ true }; //!< Switch for turning off material effects temporarily
    bool m_brem{ false };
    int m_propagateWithPathLimit{};
    size_t m_currentLayerBin{};
    double m_matupd_lastmom{};
    double m_matupd_lastpath{};
    double m_matdump_lastpath{};
    double m_delRad{ 0 };    // deRad/dl;
    double m_delIoni{ 0 };   // deIoni/dl;
    double m_sigmaIoni{ 0 }; // dsigma(ioni)/dl;
    double m_kazL{ 0 };      // kazL constant;
    double m_sigmaRad{ 0 };  // dsigma(rad)/dl;
    // cache for input variance
    double m_inputThetaVariance{};
    double m_stragglingVariance{};

    double m_pathLimit{};
    double m_timeOfFlight{};
    double m_timeStep{};
    double m_particleMass{ 0 }; //!< cache
    double m_charge{};
    double m_combinedThickness{};
    // secondary interactions
    double m_timeIn{};
    double m_bremMom{ 0. };
    double m_bremEmitThreshold{ 0. };
    double m_bremSampleThreshold{ 0. };
    double m_P[45];

    const Trk::BinnedMaterial* m_binMat{ nullptr };
    std::vector<const Trk::TrackStateOnSurface*>* m_matstates{
      nullptr
    }; //!< cache of TrackStateOnSurfaces
    std::vector<std::pair<std::unique_ptr<Trk::TrackParameters>, int>>* m_identifiedParameters{
      nullptr
    };                                                 //!< cache of intersections
    std::vector<Trk::HitInfo>* m_hitVector{ nullptr }; //!< cache of intersections/hit info

    ParticleHypothesis m_particle{};
    const TrackingVolume* m_trackingVolume{ nullptr };
    const Material* m_material{ nullptr };
    Trk::ExtrapolationCache* m_extrapolationCache{
      nullptr
    }; //!< cache for collecting the total X0 ans Elos
    // cache for combined covariance matrix contribution
    AmgSymMatrix(5) m_combinedCovariance;
    // cache for differential covariance matrix contribution ( partial material dump )
    AmgSymMatrix(5) m_covariance;
    Trk::EnergyLoss m_combinedEloss;
    std::vector<std::pair<int, std::pair<double, double>>> m_currentDist;
    MagField::AtlasFieldCache m_fieldCache;

    Cache() { m_currentDist.reserve(100); }
  };

private:
  /////////////////////////////////////////////////////////////////////////////////
  // Main functions for propagation
  /////////////////////////////////////////////////////////////////////////////////
  std::unique_ptr<Trk::TrackParameters> propagateRungeKutta(
    Cache& cache,
    bool errorPropagation,
    const Trk::TrackParameters& trackParameters,
    const Trk::Surface& targetSurface,
    Trk::PropDirection propagationDirection,
    const MagneticFieldProperties& magneticFieldProperties,
    ParticleHypothesis particle,
    const Trk::BoundaryCheck& boundaryCheck,
    double* Jacobian,
    bool returnCurv = false) const;

  /////////////////////////////////////////////////////////////////////////////////
  // Main function for propagation
  // with search of closest surface (ST)
  /////////////////////////////////////////////////////////////////////////////////
  std::unique_ptr<Trk::TrackParameters> propagateRungeKutta(
    Cache& cache,
    bool errorPropagation,
    const Trk::TrackParameters& trackParameters,
    std::vector<DestSurf>& targetSurfaces,
    Trk::PropDirection propagationDirection,
    const MagneticFieldProperties& magneticFieldProperties,
    ParticleHypothesis particle,
    std::vector<unsigned int>& solutions,
    double& path,
    bool returnCurv = false) const;

  /////////////////////////////////////////////////////////////////////////////////
  // Method of the propagation
  /////////////////////////////////////////////////////////////////////////////////
  bool propagateWithJacobian(Cache& cache,
                             bool errorPropagation,
                             Trk::SurfaceType surfaceType,
                             double* targetSurface,
                             double* P,
                             double& path) const;

  ////////////////////////////////////////////////////////////////////////////////
  // Method for propagation with search of closest surface (ST)
  ////////////////////////////////////////////////////////////////////////////////
  bool propagateWithJacobian(Cache& cache,
                             bool errorPropagation,
                             std::vector<DestSurf>& sfs,
                             double* P,
                             Trk::PropDirection propDir,
                             std::vector<unsigned int>& solutions,
                             double& path,
                             double sumPath) const;

  /////////////////////////////////////////////////////////////////////////////////
  // Trajectory model
  /////////////////////////////////////////////////////////////////////////////////
  bool rungeKuttaStep(Cache& cache,
                      bool errorPropagation,
                      double& h,
                      double* P,
                      double* dDir,
                      double* BG1,
                      bool& firstStep,
                      double& distanceStepped) const;

  /////////////////////////////////////////////////////////////////////////////////
  // Get the magnetic field and gradients
  // Input: Globalposition
  // Output: BG, which contains Bx, By, Bz, dBx/dx, dBx/dy, dBx/dz, dBy/dx,
  // dBy/dy, dBy/dz, dBz/dx, dBz/dy, dBz/dz
  /////////////////////////////////////////////////////////////////////////////////
  void getMagneticField(Cache& cache,
                        const Amg::Vector3D& position,
                        bool getGradients,
                        double* ATH_RESTRICT BG) const;

  /////////////////////////////////////////////////////////////////////////////////
  // dg/dlambda for non-electrons (g=dEdX and lambda=q/p).
  /////////////////////////////////////////////////////////////////////////////////
  double dgdlambda(Cache& cache, double l) const;

  /////////////////////////////////////////////////////////////////////////////////
  // Multiple scattering and straggling contributionto the covariance matrix
  // local covariance - to be phased out
  /////////////////////////////////////////////////////////////////////////////////
  void covarianceContribution(Cache& cache,
                              const Trk::TrackParameters* trackParameters,
                              double path,
                              const Trk::TrackParameters* targetParms,
                              AmgSymMatrix(5) * measurementCovariance) const;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Multiple scattering and straggling contributionto the covariance matrix in
  // curvilinear representation
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  void covarianceContribution(Cache& cache,
                              const Trk::TrackParameters* trackParameters,
                              double path,
                              double finalMomentum,
                              AmgSymMatrix(5) * measurementCovariance) const;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // dump material effects
  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  void dumpMaterialEffects(Cache& cache,
                           const Trk::CurvilinearParameters* trackParameters,
                           double path) const;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // update material effects   // to follow change of material
  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  void updateMaterialEffects(Cache& cache,
                             double p,
                             double sinTh,
                             double path) const;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Calculate energy loss in MeV/mm. The radiative effects are scaled by
  // m_radiationScale (1=mean, 0.5=mean(log10), 0.1=mpv)
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  double dEds(Cache& cache, double p) const;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Momentum smearing (simulation mode)
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void smear(Cache& cache,
             double& phi,
             double& theta,
             const Trk::TrackParameters* parms,
             double radDist) const;
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Bremstrahlung (simulation mode)
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void sampleBrem(Cache& cache, double mom) const;

  double m_tolerance; //!< Error tolerance. Low tolerance gives high accuracy
  bool m_materialEffects;    //!< Switch material effects on or off
  bool m_includeBgradients;  //!< Include B-gradients in the error propagation
  bool m_includeGgradient;   //!< Include g-gradient in the error propagation
  double m_momentumCutOff;   //!< Stop propagation below this momentum
  bool m_multipleScattering; //!< Switch multiple scattering on or off
  bool m_energyLoss;
  bool m_detailedEloss;
  bool m_straggling;
  bool m_MPV;
  double m_stragglingScale;
  double m_scatteringScale;
  double m_maxPath;
  double m_maxSteps;
  double m_layXmax;

  // simulation mode
  bool m_simulation;
  ToolHandle<ITimedMatEffUpdator>
    m_simMatUpdator; //!< secondary interactions (brem photon emission)
  /** Random Generator service */
  ServiceHandle<IAtRndmGenSvc> m_rndGenSvc;
  /** Random engine */
  CLHEP::HepRandomEngine* m_randomEngine;
  std::string m_randomEngineName;

  // Read handle for conditions object to get the field cache
  SG::ReadCondHandleKey<AtlasFieldCacheCondObj> m_fieldCacheCondObjInputKey{
    this,
    "AtlasFieldCacheCondObj",
    "fieldCondObj",
    "Name of the Magnetic Field conditions object key"
  };
  void getFieldCacheObject(Cache& cache, const EventContext& ctx) const;
};

}//end of namespace Trk

#endif // STEP_Propagator_H
