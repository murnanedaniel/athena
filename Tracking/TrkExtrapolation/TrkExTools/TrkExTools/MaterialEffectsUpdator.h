/*
   Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
 */

///////////////////////////////////////////////////////////////////
// MaterialEffectsUpdator.h, c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef TRKEXTOOLS_MATERIALEFFECTSUPDATOR_H
#define TRKEXTOOLS_MATERIALEFFECTSUPDATOR_H

#include "TrkExInterfaces/IMaterialEffectsUpdator.h"
// Trk
#include "TrkEventPrimitives/PropDirection.h"
#include "TrkExUtils/MaterialUpdateMode.h"
#include "TrkParameters/TrackParameters.h"
// Gaudi
#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/ToolHandle.h"

#include <string>
#include <vector>

#include <boost/thread/tss.hpp>

#define TRKEXTOOLS_MAXUPDATES 100
#ifndef COVARIANCEUPDATEWITHCHECK
#define COVARIANCEUPDATEWITHCHECK(cov, sign, value)                            \
  cov += (sign > 0 ? value : (value > cov ? 0 : sign * value))
#endif

namespace Trk {

class Layer;
// class TrackParameters;
class MaterialProperties;
class IEnergyLossUpdator;
class IMultipleScatteringUpdator;
class IMaterialMapper;

/** @class MaterialEffectsUpdator

  Point-like (also called surface-based) update of
  TrackParameters and the associated covariance

  Pre/Post/update notation.

  Given a sensitive detector element, it can be
  (e.g silicon sensors) that most of the material
  corresponds to the sensor itself or is located
  behind it.
  
  In this case the material effects would be post-update
  with respect to the measurement update on the given surface,
  e.g in a Kalman filter procedure.
  The same material would  have to be taken into account before, pre-update,
  the measurement update in the smoothing/ backward
  filtering kalman filtering process.

  Pre/Post update make sense in these specific contexts.
  For a uniform detector or non-sensitive material
  a full update is applied.


  @author Andreas.Salzburger@cern.ch
  @author AthenaMT modifications Christos Anastopoulos
  */

class MaterialEffectsUpdator final
  : public AthAlgTool
  , virtual public IMaterialEffectsUpdator
{

public:
  /** AlgTool like constructor */
  MaterialEffectsUpdator(const std::string&,
                         const std::string&,
                         const IInterface*);

  /**Virtual destructor*/
  virtual ~MaterialEffectsUpdator();

  /** AlgTool initailize method.*/
  virtual StatusCode initialize() override;
  /*
   * The concrete cache class for this specialization of the
   * IMaterialEffectsUpdator
   */
  typedef IMaterialEffectsUpdator::ICache ICache;
  /** Updator interface (full update for a layer)
    ---> ALWAYS  pointer to new TrackParameters is returned
    */
  virtual std::unique_ptr<TrackParameters> update(
    ICache& cache,
    const TrackParameters* parm,
    const Layer& sf,
    PropDirection dir = alongMomentum,
    ParticleHypothesis particle = pion,
    MaterialUpdateMode matupmode = addNoise) const override
  {
    return updateImpl(cache, parm, sf, dir, particle, matupmode);
  }
  /** Updator interface (full update for a layer) according to user
    input through MaterialEffectsOnTrack
    ---> ALWAYS pointer to new TrackParameters is returned
    */
  virtual std::unique_ptr<TrackParameters> update(
    ICache& cache,
    const TrackParameters* parm,
    const MaterialEffectsOnTrack& meff,
    Trk::ParticleHypothesis particle = pion,
    MaterialUpdateMode matupmode = addNoise) const override
  {
    return updateImpl(cache, parm, meff, particle, matupmode);
  }

  /** Updator interface (pre-update for a layer):
    ---> ALWAYS pointer to new TrackParameters is returned
    */
  virtual std::unique_ptr<TrackParameters> preUpdate(
    ICache& cache,
    const TrackParameters* parm,
    const Layer& sf,
    PropDirection dir = alongMomentum,
    ParticleHypothesis particle = pion,
    MaterialUpdateMode matupmode = addNoise) const override
  {
    return preUpdateImpl(cache, parm, sf, dir, particle, matupmode);
  }

  /** Updator interface (post-update for a layer):
    ---> ALWAYS pointer to new TrackParameters is returned
    if no postUpdate is to be done : return nullptr
    */
  virtual std::unique_ptr<TrackParameters> postUpdate(
    ICache& cache,
    const TrackParameters& parm,
    const Layer& sf,
    PropDirection dir = alongMomentum,
    ParticleHypothesis particle = pion,
    MaterialUpdateMode matupmode = addNoise) const override final
  {
    return postUpdateImpl(cache, parm, sf, dir, particle, matupmode);
  }
  /** Dedicated Updator interface:-> create new track parameters*/
  virtual std::unique_ptr<TrackParameters> update(
    ICache& cache,
    const TrackParameters& parm,
    const MaterialProperties& mprop,
    double pathcorrection,
    PropDirection dir = alongMomentum,
    ParticleHypothesis particle = pion,
    MaterialUpdateMode matupmode = addNoise) const override final
  {
    return updateImpl(
      cache, parm, mprop, pathcorrection, dir, particle, matupmode);
  }

  /** Validation Action - calls the writing and resetting of the TTree variables
   */
  virtual void validationAction(ICache& cache) const override final
  {
    validationActionImpl(cache);
  }

  /** Only has an effect if m_landauMode == true.
    Resets mutable variables used for non-local calculation of energy loss if
    parm == 0. Otherwise, modifies parm with the final update of the covariance
    matrix*/
  virtual void modelAction(ICache& cache, const TrackParameters* parm = nullptr)
    const override final
  {
    modelActionImpl(cache, parm);
  }

public:
  /*
   * Public methods using the TLS cache.
   */
  virtual std::unique_ptr<TrackParameters> update(
    const TrackParameters* parm,
    const Layer& sf,
    PropDirection dir = alongMomentum,
    ParticleHypothesis particle = pion,
    MaterialUpdateMode matupmode = addNoise) const override final
  {

    ICache& cache = getTLSCache();
    return updateImpl(cache, parm, sf, dir, particle, matupmode);
  }

  virtual std::unique_ptr<TrackParameters> update(
    const TrackParameters* parm,
    const MaterialEffectsOnTrack& meff,
    Trk::ParticleHypothesis particle = pion,
    MaterialUpdateMode matupmode = addNoise) const override final
  {
    ICache& cache = getTLSCache();
    return updateImpl(cache, parm, meff, particle, matupmode);
  }

  virtual std::unique_ptr<TrackParameters> preUpdate(
    const TrackParameters* parm,
    const Layer& sf,
    PropDirection dir = alongMomentum,
    ParticleHypothesis particle = pion,
    MaterialUpdateMode matupmode = addNoise) const override final
  {
    ICache& cache = getTLSCache();
    return preUpdateImpl(cache, parm, sf, dir, particle, matupmode);
  }

  virtual std::unique_ptr<TrackParameters> postUpdate(
    const TrackParameters& parm,
    const Layer& sf,
    PropDirection dir = alongMomentum,
    ParticleHypothesis particle = pion,
    MaterialUpdateMode matupmode = addNoise) const override final
  {
    ICache& cache = getTLSCache();
    return postUpdateImpl(cache, parm, sf, dir, particle, matupmode);
  }

  virtual std::unique_ptr<TrackParameters> update(
    const TrackParameters& parm,
    const MaterialProperties& mprop,
    double pathcorrection,
    PropDirection dir = alongMomentum,
    ParticleHypothesis particle = pion,
    MaterialUpdateMode matupmode = addNoise) const override final
  {
    ICache& cache = getTLSCache();
    return updateImpl(
      cache, parm, mprop, pathcorrection, dir, particle, matupmode);
  }

  virtual void validationAction() const override final
  {
    ICache& cache = getTLSCache();
    validationActionImpl(cache);
  }

  virtual void modelAction(
    const TrackParameters* parm = nullptr) const override final
  {
    ICache& cache = getTLSCache();
    modelActionImpl(cache, parm);
  }

private:
  /* The acutal implementation methods using the tool's
   * concrete  Cache*/
  std::unique_ptr<TrackParameters> updateImpl(
    ICache& cache,
    const TrackParameters* parm,
    const Layer& sf,
    PropDirection dir = alongMomentum,
    ParticleHypothesis particle = pion,
    MaterialUpdateMode matupmode = addNoise) const;

  std::unique_ptr<TrackParameters> updateImpl(
    ICache& cache,
    const TrackParameters* parm,
    const MaterialEffectsOnTrack& meff,
    Trk::ParticleHypothesis particle = pion,
    MaterialUpdateMode matupmode = addNoise) const;

  std::unique_ptr<TrackParameters> preUpdateImpl(
    ICache& cache,
    const TrackParameters* parm,
    const Layer& sf,
    PropDirection dir = alongMomentum,
    ParticleHypothesis particle = pion,
    MaterialUpdateMode matupmode = addNoise) const;

  std::unique_ptr<TrackParameters> postUpdateImpl(
    ICache& cache,
    const TrackParameters& parm,
    const Layer& sf,
    PropDirection dir = alongMomentum,
    ParticleHypothesis particle = pion,
    MaterialUpdateMode matupmode = addNoise) const;

  std::unique_ptr<TrackParameters> updateImpl(
    ICache& cache,
    const TrackParameters* parm,
    const MaterialProperties& mprop,
    double pathcorrection,
    PropDirection dir = alongMomentum,
    ParticleHypothesis particle = pion,
    MaterialUpdateMode matupmode = addNoise) const;

  std::unique_ptr<TrackParameters> updateImpl(
    ICache& cache,
    const TrackParameters& parm,
    const MaterialProperties& mprop,
    double pathcorrection,
    PropDirection dir = alongMomentum,
    ParticleHypothesis particle = pion,
    MaterialUpdateMode matupmode = addNoise) const;

  static void validationActionImpl(ICache& cache);

  static void modelActionImpl(ICache& cache,
                              const TrackParameters* parm = nullptr);

  /** A simple check method for the 'removeNoise' update model */
  bool checkCovariance(AmgSymMatrix(5) & updated) const;
  TrackParameters* finalLandauCovarianceUpdate(
    const TrackParameters* parm) const;

  /* Private Class members*/
  bool m_doCompoundLayerCheck; //!< turn on/off the necessary checks when we may
                               //!< have compound layers
  bool m_doEloss;              //!< steer energy loss On/Off from outside
  bool m_doMs;                 //!< steer multiple scattering On/Off from outside

  bool m_forceMomentum;        //!< Force the momentum to be a specific value
  bool m_xKalmanStraggling;    //!< the momentum Error as calculated in xKalman
  bool m_useMostProbableEloss; //!< use the most probable energy loss

  bool m_msgOutputValidationDirection; //!< validation direction used for screen
                                       //!< output
  bool m_msgOutputCorrections;         //!< screen output of actual corrections

  // ------------ validation variables
  // -------------------------------------------------
  bool m_validationMode;             //!< Switch for validation mode
  bool m_validationIgnoreUnmeasured; //!< Ignore unmeasured TrackParameters
                                     //!< (Navigation!)
  bool m_landauMode;                 //!< If in Landau mode, error propagation is done as for
                                     //!< landaus
  int m_validationDirection;         //!< validation direction
  //  ------------------------------
  double m_momentumCut;    //!< Minimal momentum cut for update
  double m_momentumMax;    //!< Maximal momentum cut for update
  double m_forcedMomentum; //!< Forced momentum value

  ToolHandle<IEnergyLossUpdator> m_eLossUpdator;      //!< AlgoTool for EnergyLoss updates
  ToolHandle<IMultipleScatteringUpdator> m_msUpdator; //!< AlgoTool for MultipleScatterin effects
  ToolHandle<IMaterialMapper> m_materialMapper; //!< the material mapper for recording the layer material

  /*
   * TLS part
   */
  mutable boost::thread_specific_ptr<ICache> m_cache_tls;
  ICache& getTLSCache() const
  {
    ICache* cache = m_cache_tls.get();
    if (!cache) {
      cache = new ICache();
      m_cache_tls.reset(cache);
    }
    return *cache;
  }
};
} // end of namespace

#endif // TRKEXTOOLS_MATERIALEFFECTSUPDATOR_H

