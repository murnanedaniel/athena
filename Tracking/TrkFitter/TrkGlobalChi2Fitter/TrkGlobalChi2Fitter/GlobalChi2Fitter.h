/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#ifndef GLOBALCHI2FITTER_H
#define GLOBALCHI2FITTER_H
//#define GXFDEBUGCODE
#include "TrkDetDescrInterfaces/IMaterialEffectsOnTrackProvider.h"
#include "AthenaBaseComps/AthAlgTool.h"
#include "AthenaBaseComps/AthCheckedComponent.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/EventContext.h"


#include "TrkToolInterfaces/ITrkMaterialProviderTool.h"
#include "TrkToolInterfaces/IResidualPullCalculator.h"
#include "TrkToolInterfaces/IRIO_OnTrackCreator.h"
#include "TrkToolInterfaces/IUpdator.h"
#include "TrkToolInterfaces/IBoundaryCheckTool.h"

#include "TrkExInterfaces/IExtrapolator.h"
#include "TrkExInterfaces/IPropagator.h"
#include "TrkExInterfaces/INavigator.h"
#include "TrkExInterfaces/IMultipleScatteringUpdator.h"
#include "TrkExInterfaces/IEnergyLossUpdator.h"
#include "TrkExInterfaces/IMaterialEffectsUpdator.h"

#include "TrkFitterInterfaces/IGlobalTrackFitter.h"

#include "TrkGlobalChi2Fitter/GXFTrajectory.h"
#include "TrkMaterialOnTrack/MaterialEffectsOnTrack.h"
#include "TrkFitterUtils/FitterStatusCode.h"
#include "TrkEventPrimitives/PropDirection.h"
#include "MagFieldConditions/AtlasFieldCacheCondObj.h"
#include "MagFieldElements/AtlasFieldCache.h"

#include "StoreGate/ReadCondHandleKey.h"
#include "TrkGeometry/TrackingGeometry.h"

#include <memory>
#include <mutex>

/**
 * These headers, as well as other headers in the TrkGlobalChi2Fitter package
 * use modern C++11 memory ownership semantics expressed through smart
 * pointers. This means that ownership is encoded explicitly in the signatures
 * of methods. The general idea is that raw pointers (T *) express non-owning
 * pointers; they imply that one is simply viewing an object, and there is no
 * need to explicitly free them. It is also impossible to transfer ownership of
 * them, as that would be stealing. On the contrary, smart pointers such as
 * std::unique_ptr<T> imply ownership. This gives rise to four general rules:
 *
 * 1. Raw pointers passed as arguments - void f(T *): In this case, the
 *    interface makes a guarantee to the call site: the ownership of the memory
 *    is not changed, and whoever was responsible for the memory remains
 *    responsible for it after the method is called. The callee does not free
 *    or delete the memory, since it does not own it. The method may store an
 *    observing pointer to the argument beyond the lifetime of the callee,
 *    however.
 * 2. Smart pointers are passed as arguments - void f(std::unique_ptr<T>): In
 *    this case, the method assumes full ownership of the pointer. The caller
 *    does no longer own the object and should not attempt to free it, that is
 *    the responsibility of the callee now. It is possible that the callee
 *    passes ownership on to somewhere else, but the beauty if this approach is
 *    that this is abstracted away from the call site, and you do not need to
 *    worry about it as a programmer. In order to call a method like this, the
 *    call site will need to prepare an r-value unique_ptr, which is how the
 *    compiler helps protect us against memory errors in ways that it would not
 *    be able to when using smart pointers. There are four common (and deeply
 *    intertwined, if you think about it) ways of creating an r-value unique
 *    pointer. Firstly, one may pass the return value of an expression
 *    directly. For example, if we have a function std::unique_ptr<T> g(void),
 *    we can call f(g()). Since the r-value generated by g is never turned into
 *    an l-value, we can pass it directly to g. Secondly, we can create a new
 *    object on the heap using std::make_unique, passed directly to the method.
 *    Thirdly, we can turn a raw pointer into a smart pointer using the
 *    std::unique_ptr<T>(T *) constructor, which takes a raw pointer as
 *    argument. Care must be taken that the caller of this constructor has
 *    ownership of the raw pointer, and this should be used primarily as a way
 *    to interface with legacy code where raw pointers have ownership. Finally,
 *    one may std::move an existing l-value std::unique_ptr to turn it into an
 *    r-value to pass to a method. This is probably the most common use case.
 *    Take care that moving a unique pointer nulls it, so the call site should
 *    not attempt to access the unique_ptr after moving it into a method.
 * 3. Raw pointers returned from methods - T * f(void): In this case, there is
 *    no transfer of ownership from the callee to the caller. The caller should
 *    not attempt to delete, free, acquire or pass ownership to the pointer.
 *    Treat it as a normal observant pointer and don't worry about memory
 *    management. Be aware that the raw pointer may outlive the memory it
 *    points to.
 * 4. Smart pointers returned from methods - std::unique_ptr<T> f(void): In
 *    this case, ownership of the memory is transferred from the callee to the
 *    caller, and the caller is responsible for handling the memory management.
 *    It is highly recommended but not enforced to use the smart pointer RAII
 *    mechanisms for this. It is also allowed for the caller to transfer
 *    ownership of the returned pointer, since the interface guarantees that
 *    the callee no longer cares about the object.
 *
 * If you are unfamiliar with smart pointers, this might seem like a whole
 * bunch of unnecessary fanciness, and you may wonder what the purpose of this
 * is. In brief, the idea is to make the interface more predictable for the
 * developer. By encoding the ownership semantics into the function signatures,
 * we make it clear from the signature alone what the ownership mechanics are.
 * Gone are the days of wondering whether a pointer will still exist after
 * passing it to a method, or whether you should be deleting a pointer that was
 * returned from a method. This is all clear from the use of raw and smart
 * pointers now. Of course, it might take some time getting used to reading the
 * signatures at first. In addition, this approach helps the compiler protect
 * us against memory errors, and it encourages us to use smart pointer RAII
 * mechanics to delete and free memory instead of doing it manually, which
 * improves memory safety.
 *
 * While the unique pointer approach provides far, far better memory safety
 * than the traditional raw pointer approach, it is still lacking in many ways
 * which might be striking if you are coming from a safe language such as Rust.
 * Most specifically, the smart pointer approach described above makes no
 * guarantees about the lifetimes of objects. If a method returns a raw viewing
 * pointer to some memory, there is no guarantee that the actual memory it
 * points to will outlive the observing pointer. Unfortunately, there is no
 * robust and performant way around this, and the developer will have to take
 * care to avoid lifetime issues. Caveat emptor.
 *
 * We've made great strides in enforcing a model for memory semantics in the
 * GlobalChi2Fitter code: the number of explicit frees and deletes as well as
 * the use of the new keyword have dropped by some 90%. However, the model is
 * still not enforced perfectly everywhere. In cases where the model is
 * violated, the comments in the header files should inform you of the
 * violations.
 *
 * For more information on smart pointer semantics, please head over to the
 * following websites:
 * * https://timoch.com/blog/2013/04/std-unique_ptr-semantic/
 * * https://en.cppreference.com/w/cpp/memory/unique_ptr
 */

class AtlasDetectorID;


namespace Trk {
  class Track;
  class TransportJacobian;
  class TrackFitInputPreparator;
  class IMagneticFieldTool;
  class MeasuredPerigee;
  class PrepRawDataComparisonFunction;
  class MeasurementBaseComparisonFunction;
  class MaterialEffectsOnTrack;
  class Layer;
  class CylinderLayer;
  class DiscLayer;
  class MagneticFieldProperties;
  class TrackingVolume;
  class Volume;

  class GlobalChi2Fitter: public extends<AthCheckedComponent<AthAlgTool>, IGlobalTrackFitter> {
    struct PropagationResult {
      std::unique_ptr<const TrackParameters> m_parameters;
      std::unique_ptr<TransportJacobian> m_jacobian;
      std::optional<std::vector<std::unique_ptr<const TrackParameters>>> m_preholes;
    };

    /*
     * This struct serves as a very simple container for the five different
     * hole and dead module states that we want to count, which are the dead
     * Pixels, dead SCTs, Pixel holes, SCT holes, and SCT double holes.
     */
    struct TrackHoleCount {
      unsigned int m_pixel_hole = 0;
      unsigned int m_sct_hole = 0;
      unsigned int m_sct_double_hole = 0;
      unsigned int m_pixel_dead = 0;
      unsigned int m_sct_dead = 0;
    };

    struct Cache {
      /*
       * Currently the information about what type of fit is being passed by the
       * presence of a TrackingVolume.
       */
      template <class T>
      static
      void objVectorDeleter(const std::vector<const T *> *ptr) {
        if (ptr) {
          for (const T *elm : *ptr) { delete elm; }
          delete ptr;
        }
      }

      const TrackingGeometry *m_trackingGeometry = nullptr;
      const TrackingVolume *m_caloEntrance = nullptr;
      const TrackingVolume *m_msEntrance = nullptr;

      bool m_calomat, m_extmat;
      bool m_idmat = true;
      bool m_sirecal;
      bool m_getmaterialfromtrack;
      bool m_reintoutl;
      bool m_matfilled = false;
      bool m_acceleration;
      bool m_fiteloss;
      bool m_asymeloss;

      std::vector<double> m_phiweight;
      std::vector<int> m_firstmeasurement;
      std::vector<int> m_lastmeasurement;

      std::vector < const Trk::Layer * >m_negdiscs;
      std::vector < const Trk::Layer * >m_posdiscs;
      std::vector < const Trk::Layer * >m_barrelcylinders;

      bool m_fastmat = true;

      int m_lastiter;
      int m_miniter;

      #ifdef GXFDEBUGCODE
      int m_iterations = 0;
      #endif

      Amg::MatrixX m_derivmat;
      Amg::SymMatrixX m_fullcovmat;

      std::vector< std::unique_ptr< const std::vector < const TrackStateOnSurface *>,
                                    void (*)(const std::vector<const TrackStateOnSurface *> *) > >
        m_matTempStore;

      MagField::AtlasFieldCache m_field_cache;

      FitterStatusCode m_fittercode;

      Cache(const GlobalChi2Fitter *fitter):
        m_calomat(fitter->m_calomat),
        m_extmat(fitter->m_extmat),
        m_sirecal(fitter->m_sirecal),
        m_getmaterialfromtrack(fitter->m_getmaterialfromtrack),
        m_reintoutl(fitter->m_reintoutl),
        m_acceleration(fitter->m_acceleration),
        m_fiteloss(fitter->m_fiteloss),
        m_asymeloss(fitter->m_asymeloss),
        m_miniter(fitter->m_miniter)
      {}

      Cache & operator=(const Cache &) = delete;
    };
  private:

    enum FitterStatusType {
      S_FITS,
      S_SUCCESSFUL_FITS,
      S_MAT_INV_FAIL,
      S_NOT_ENOUGH_MEAS,
      S_PROPAGATION_FAIL,
      S_INVALID_ANGLES,
      S_NOT_CONVERGENT,
      S_HIGH_CHI2,
      S_LOW_MOMENTUM,
      __S_MAX_VALUE
    };

  public:
    GlobalChi2Fitter(
      const std::string &,
      const std::string &,
      const IInterface *
    );

    virtual ~ GlobalChi2Fitter();

    virtual StatusCode initialize() override;
    virtual StatusCode finalize() override;
    /*
     * Bring in default impl with
     * EventContext for now
     */
    using ITrackFitter::fit;

    virtual std::unique_ptr<Track> fit(
      const EventContext& ctx,
      const PrepRawDataSet&,
      const TrackParameters&,
      const RunOutlierRemoval runOutlier = false,
      const ParticleHypothesis matEffects = nonInteracting
      ) const override final;

    virtual std::unique_ptr<Track> fit(
      const EventContext& ctx,
      const Track &,
      const RunOutlierRemoval runOutlier = false,
      const ParticleHypothesis matEffects = nonInteracting
    ) const override final;

    virtual std::unique_ptr<Track> fit(
      const EventContext& ctx,
      const MeasurementSet &,
      const TrackParameters &,
      const RunOutlierRemoval runOutlier = false,
      const ParticleHypothesis matEffects = nonInteracting
    ) const override final;

    virtual std::unique_ptr<Track> fit(
      const EventContext& ctx,
      const Track &,
      const PrepRawDataSet &,
      const RunOutlierRemoval runOutlier = false,
      const ParticleHypothesis matEffects = nonInteracting
    ) const override final;

    virtual std::unique_ptr<Track> fit(
      const EventContext& ctx,
      const Track &,
      const Track &,
      const RunOutlierRemoval runOutlier = false,
      const ParticleHypothesis matEffects = nonInteracting
    ) const override final;

    virtual std::unique_ptr<Track> fit(
      const EventContext& ctx,
      const Track &,
      const MeasurementSet &,
      const RunOutlierRemoval runOutlier = false,
      const ParticleHypothesis matEffects = nonInteracting
    ) const override final;

    virtual Track* alignmentFit(
      AlignmentCache&,
      const Track&,
      const RunOutlierRemoval  runOutlier=false,
      const ParticleHypothesis matEffects=Trk::nonInteracting
    ) const override;

  private:
    static void calculateJac(
      Eigen::Matrix<double, 5, 5> &,
      Eigen::Matrix<double, 5, 5> &,
      int, int
    ) ;

    Track * fitIm(
      const EventContext& ctx,
      Cache & cache,
      const Track & inputTrack,
      const RunOutlierRemoval runOutlier,
      const ParticleHypothesis matEffects
    ) const;

    Track *myfit(
      const EventContext& ctx,
      Cache &,
      GXFTrajectory &,
      const TrackParameters &,
      const RunOutlierRemoval runOutlier = false,
      const ParticleHypothesis matEffects = nonInteracting
    ) const;

    Track *myfit_helper(
      Cache &,
      GXFTrajectory &,
      const TrackParameters &,
      const RunOutlierRemoval runOutlier = false,
      const ParticleHypothesis matEffects = nonInteracting
    ) const;

    Track *mainCombinationStrategy(
      const EventContext& ctx,
      Cache &,
      const Track &,
      const Track &,
      GXFTrajectory &,
      std::vector<MaterialEffectsOnTrack> &
    ) const;

    Track *backupCombinationStrategy(
      const EventContext& ctx,
      Cache &,
      const Track &,
      const Track &,
      GXFTrajectory &,
      std::vector<MaterialEffectsOnTrack> &
    ) const;

    void makeProtoState(
      Cache &,
      GXFTrajectory &,
      const TrackStateOnSurface *,
      int index = -1
     ) const;

    void makeProtoStateFromMeasurement(
      Cache &,
      GXFTrajectory &,
      const MeasurementBase *,
      const TrackParameters * trackpar = nullptr,
      bool isoutlier = false,
      int index = -1
    ) const;

    bool processTrkVolume(
      Cache &,
      const Trk::TrackingVolume * tvol
    ) const;

    /**
     * @brief Find the intersection of a set of track parameters onto a disc
     * surface.
     *
     * Calculates the intersection from a point and momentum in space onto a
     * disc surface which represents a disc-shaped layer in the detector. The
     * position of the intersection can be used to find materials in that layer
     * at that position.
     *
     * @param[in] cache The standard GX2F cache.
     * @param[in] surface The surface to intersect with.
     * @param[in] param1 The main track parameters to calculate the
     * intersection from.
     * @param[in] param2 A secondary set of parameters used for electrons. The
     * purpose of this is not known to us at this time.
     * @param[in] mat A particle hypothesis describing the behaviour of the
     * particle.
     *
     * @returns Nothing if the intersection failed (i.e. there was no
     * intersection), otherwise both an intersection positition as well as the
     * angle of inflection.
     *
     * @note This method can probably be replaced entirely by the straight line
     * intersection method of the appropriate Surface subclass.
     */
    static std::optional<std::pair<Amg::Vector3D, double>> addMaterialFindIntersectionDisc(
      Cache & cache,
      const DiscSurface & surface,
      const TrackParameters & param1,
      const TrackParameters & param2,
      const ParticleHypothesis mat
    ) ;

    /**
     * @brief Find the intersection of a set of track parameters onto a
     * cylindrical surface.
     *
     * See addMaterialFindIntersectionDisc for more information.
     *
     * @note This method can probably be replaced entirely by the straight line
     * intersection method of the appropriate Surface subclass.
     */
    static std::optional<std::pair<Amg::Vector3D, double>> addMaterialFindIntersectionCyl(
      Cache & cache,
      const CylinderSurface & surface,
      const TrackParameters & param1,
      const TrackParameters & param2,
      const ParticleHypothesis mat
    ) ;

    /**
     * @brief Given layer information, probe those layers for scatterers and
     * add them to a track.
     *
     * This is the meat of the pudding, if you will. Given the information that
     * we have about layers, go through them all and find any possible material
     * hits that we need to add to the track.
     *
     * @param[in,out] cache General cache object.
     * @param[in,out] track The track object as it exists now in IR.
     * @param[in] offset The first state after any existing materials.
     * @param[in] layers The list of layers.
     * @param[in] ref1 The first set of reference parameters.
     * @param[in] ref2 The second set of reference parameters.
     * @param[in] mat The particle hypothesis describing the track behaviour.
     *
     * @note Attentive readers may wonder why we pass this function a vector
     * of layers, but not a vector of upstream layers. The reason for this is
     * that the vector of upstream layers is also a member of the cache object.
     */
    void addMaterialUpdateTrajectory(
      Cache & cache,
      GXFTrajectory & track,
      int offset,
      std::vector<std::pair<const Layer *, const Layer *>> & layers,
      const TrackParameters * ref1,
      const TrackParameters * ref2,
      ParticleHypothesis mat
    ) const;

    /**
     * @brief Collect all possible layers that a given track could have passed
     * through.
     *
     * If we are to use layer information to determine possible scatterer hits,
     * we must first gather those layers. That's what this method does. It
     * looks for disc and barrel cylinder layers that the given track might
     * have crossed and collects them into output vectors. One contains layers
     * between states on the track, and the upstream layers lie before the
     * first state of the track.
     *
     * @param[in,out] cache General cache object.
     * @param[out] layers Output vector for layers.
     * @param[out] uplayers Output vector for upstream layers, which lie before
     * the first hit in the track.
     * @param[in] states A list of track states on the track.
     * @param[in] first The first track state.
     * @param[in] last The last track state.
     * @param[in] refpar Reference parameters from which to extrapolate.
     * @param[in] hasmat Are there any existing materials on this track?
     */
    static void addMaterialGetLayers(
      Cache & cache,
      std::vector<std::pair<const Layer *, const Layer *>> & layers,
      std::vector<std::pair<const Layer *, const Layer *>> & uplayers,
      const std::vector<std::unique_ptr<GXFTrackState>> & states,
      GXFTrackState & first,
      GXFTrackState & last,
      const TrackParameters * refpar,
      bool hasmat
    ) ;

    /**
     * @brief A faster strategy for adding scatter material to tracks, works
     * only for inner detector tracks.
     *
     * For every track, we need to add its scatterers. That is to say, we need
     * to determine which bits of non-active material the particle in question
     * may have passed through and add them to the track. This is generally an
     * expensive operation, but we can cut some corners if the track only
     * consists of inner detector hits. Specifically, we can exploit the layer
     * structure of the detector to find possible material hits more quickly
     * and efficiently than using the standard material adding algorithm, which
     * is addMaterial.
     *
     * @param[in,out] cache General cache object, as used everywhere.
     * @param[in,out] trajectory The current state of the track, respresented
     * in the fitter's internal track representation. States may be added to
     * this.
     * @param[in] parameters Starting parameters for the material addition
     * step.
     * @param[in] part Standard representation of particle type, used to
     * determine the behaviour of the particle as it traverses materials.
     */
    void addIDMaterialFast(
      const EventContext& ctx,
      Cache & cache,
      GXFTrajectory & track,
      const TrackParameters * parameters,
      ParticleHypothesis part
    ) const;

    void addMaterial(
      const EventContext& ctx,
      Cache &,
      GXFTrajectory &,
      const TrackParameters *,
      ParticleHypothesis
    ) const;

    /**
     * @brief Helper method which performs an extrapolation with additional
     * logic for hole search.
     *
     * This method is a wrapper around extrapolateStepwise from the
     * extrapolator interface, with the added functionality that it will null
     * any returned track parameters which are on the start and end surface.
     *
     * @param[in] ctx An event context for extrapolation.
     * @param[in] src The track parameters to start extrapolating from.
     * @param[in] dst The track state to extrapolate to.
     * @param[in] propdir The propagation direction.
     * @return A vector of track states, just like normal extrapolation.
     */
    std::vector<std::unique_ptr<const TrackParameters>> holesearchExtrapolation(
      const EventContext & ctx,
      const TrackParameters & src,
      const GXFTrackState & dst,
      PropDirection propdir
    ) const;

    std::unique_ptr<const TrackParameters> makePerigee(
      Cache &,
      const TrackParameters &,
      const ParticleHypothesis
    ) const;

    static void makeTrackFillDerivativeMatrix(
      Cache &,
      GXFTrajectory &
    ) ;

    std::unique_ptr<const TrackParameters> makeTrackFindPerigeeParameters(
      const EventContext &,
      Cache &,
      GXFTrajectory &,
      const ParticleHypothesis
    ) const;

    std::unique_ptr<GXFTrackState> makeTrackFindPerigee(
      const EventContext &,
      Cache &,
      GXFTrajectory &,
      const ParticleHypothesis
    ) const;

    std::unique_ptr<Track> makeTrack(
      const EventContext& ctx,
      Cache &,
      GXFTrajectory &,
      const ParticleHypothesis
    ) const;

    std::unique_ptr<TrackStateOnSurface> makeTSOS(
      GXFTrackState &
    ) const;

    void fillResiduals(
      Cache &,
      GXFTrajectory &,
      int,
      Amg::SymMatrixX &,
      Amg::VectorX &,
      Amg::SymMatrixX &,
      bool &
    ) const;

    void fillDerivatives(
      GXFTrajectory & traj,
      bool onlybrem = false
    ) const;

    FitterStatusCode runIteration(
      const EventContext& ctx,
      Cache &,
      GXFTrajectory &,
      int,
      Amg::SymMatrixX &,
      Amg::VectorX &,
      Amg::SymMatrixX &,
      bool &
    ) const;

    FitterStatusCode updateFitParameters(
      GXFTrajectory &,
      Amg::VectorX &,
      const Amg::SymMatrixX &
    ) const;

    /**
     * @warning This method has some unclear memory ownership mechanics that
     * might not correspond fully with the model described at the beginning of
     * the file. Be aware!
     */
    GXFTrajectory *runTrackCleanerSilicon(
      const EventContext& ctx,
      Cache &,
      GXFTrajectory &,
      Amg::SymMatrixX &,
      Amg::SymMatrixX &,
      Amg::VectorX &,
      bool
    ) const;

    void runTrackCleanerTRT(
      Cache &,
      GXFTrajectory &,
      Amg::SymMatrixX &,
      Amg::VectorX &,
      Amg::SymMatrixX &,
      bool, bool, int
    ) const;

    /**
     * @brief Helper method that encapsulates calls to the propagator tool in
     * the calculateTrackParameters() method.
     *
     * This method encapsulates some of the logic relating to passing or not
     * passing a Jacobian matrix to make the calculateTrackParameters() a lot
     * more readable.
     *
     * For information about parameters see the
     * IPropagator::propagateParameters() method, which this method almost
     * directly wraps.
     */
    PropagationResult calculateTrackParametersPropagateHelper(
      const EventContext &,
      const TrackParameters &,
      const GXFTrackState &,
      PropDirection,
      const MagneticFieldProperties&,
      bool,
      bool
    ) const;

    /**
     * @brief Propagate onto a track state, collecting new track parameters, and
     * optionally the Jacobian and possible holes.
     *
     * This is a helper function for the calculateTrackParameters() method. Its
     * purpose is to propagate from a set of track parameters onto a surface,
     * finding the new track parameters at that surface and optionally the
     * Jacobian and a list of possible holes.
     *
     * This method uses another helper function, aptly called
     * calculateTrackParametersPropagateHelper(), which wraps the actual
     * propagator calls. What calculateTrackParametersPropagate() is call this
     * method once and check the result. If the result is invalid, that is to
     * say we didn't manage to extract a correct set of track parameters, we
     * try again but with the propagation direction flipped.
     *
     * If the calcderiv argument is set, this method will attempt to calculate
     * the Jacobian as well as the new set of track parameters. This involves
     * non-trivial logic which is abstracted away in the underlying helper
     * function.
     *
     * @param[in] ctx An event context.
     * @param[in] prev The origin track parameters to start the propagation.
     * @param[in] ts The destination track state (in GX2F internal form).
     * @param[in] propdir The propagation direction.
     * @param[in] bf The magnetic field properties.
     * @param[in] calcderiv If set, calculate the derivative.
     * @param[in] holesearch If set, search for holes.
     *
     * @return An instance of PropagationResult, which is a struct with three
     * members. Firstly, it contains a unique pointer to a set of track
     * parameters, which are the track parameters at the destination track
     * state following propagation. If these parameters are a nullpointer, that
     * indicates a failure state. Secondly, if requested, the Jacobian is stored
     * in this struct. This may be a nullptr if it was not requested or if the
     * calculation of the Jacobian failed. Thirdly, it contains a vector of
     * possible holes found between the start and end of the propagation. Since
     * these hole states are not always necessary, they are wrapped in a
     * std::optional type.
     */
    PropagationResult calculateTrackParametersPropagate(
      const EventContext &,
      const TrackParameters &,
      const GXFTrackState &,
      PropDirection,
      const MagneticFieldProperties&,
      bool,
      bool
    ) const;

    /**
     * @brief Extracts a collection of track states which are important for
     * hole search.
     *
     * This method helps extract the measurement (and outlier) states from a
     * track. These are the states between which we want to do a hole search,
     * so the result of calling this method can be used as a source of truth
     * for conducting a hole search on the track.
     *
     * This method only returns states between the first and last measurements
     * on the track, which is the region in which we are interested in doing a
     * hole search.
     *
     * As an example, if we denote scatterers as S, and measurements as M, this
     * method would reduce the following track with numbered states:
     *
     * 1 2 3 4 5 6 7 8 9
     * M S S M M S M S S
     *
     * Into a list of references [1, 4, 5, 7].
     *
     * This method ensures that each pair of consecutive states in the return
     * value list is a target for a hole search extrapolation.
     *
     * @param[in] trajectory The trajectory from which to extract states.
     * @return A vector of state references as described above.
     */
    std::vector<std::reference_wrapper<GXFTrackState>> holeSearchStates(
      GXFTrajectory & trajectory
    ) const;

    /**
     * @brief Conduct a hole search between a list of states, possibly reusing
     * existing information.
     *
     * Given a collection of state references, this method will conduct a hole
     * search between consecutive pairs of states, possibly reusing existing
     * information stored in the state data types. The method will check
     * whether the state contains any previous hole search data and use it. If
     * there is no data, it will run additional extrapolations to gather that
     * data. It will then use a helper method to count holes and dead modules
     * and return a total count.
     *
     * In some cases, this method may error. Should this occur, it will return
     * a non-extant value.
     *
     * @param[in] ctx An event context used for extrapolation.
     * @param[in] states A list of states to operate on, using consecutive
     * states as extrapolation regions.
     * @return A list of hole counts if the process succeeded, or a non-extant
     * value in case of an error.
     */
    std::optional<GlobalChi2Fitter::TrackHoleCount> holeSearchProcess(
      const EventContext & ctx,
      const std::vector<std::reference_wrapper<GXFTrackState>> & states
    ) const;

    /**
     * @brief Helper method for the hole search that does the actual counting
     * of holes and dead modules.
     *
     * This is a helper function that does a lot of esoteric and weird things
     * that you most likely won't need to know about. The gist of it is that
     * you pass it a vector of track parameters and a counting object, and it
     * will update those counters according to its analysis of the track
     * parameters.
     *
     * Unfortunately, due to the design of this method, it requires quite a lot
     * of persistent state between invocations for the same track. That's bad
     * design of course, but it is how it is for now. This means that there
     * are quite a few state parameters.
     *
     * @param[in] hc A list of candidate hole track parameters to analyse.
     * @param[in,out] id_set A set of identifiers found to be holes or dead.
     * @param[in,out] sct_set A set of identifiers of SCT holes.
     * @param[in,out] rv The hole count container to update.
     * @param[in] count_holes Holes are counted only if this is enabled.
     * @param[in] count_dead Dead modules are counted only if this is enabled.
     */
    void holeSearchHelper(
      const std::vector<std::unique_ptr<const TrackParameters>> & hc,
      std::set<Identifier> & id_set,
      std::set<Identifier> & sct_set,
      TrackHoleCount & rv,
      bool count_holes,
      bool count_dead
    ) const;

    FitterStatusCode calculateTrackParameters(
      const EventContext& ctx,
      GXFTrajectory&,
      bool) const;

    std::variant<std::unique_ptr<const TrackParameters>, FitterStatusCode> updateEnergyLoss(
      const Surface &,
      const GXFMaterialEffects &,
      const TrackParameters &,
      double,
      int
    ) const;

    static void calculateDerivatives(GXFTrajectory &) ;

    void calculateTrackErrors(GXFTrajectory &, Amg::SymMatrixX &, bool) const;

    std::unique_ptr<TransportJacobian> numericalDerivatives(
      const EventContext& ctx,
      const TrackParameters *,
      const Surface *,
      PropDirection,
      const MagneticFieldProperties&
    ) const;

    virtual int iterationsOfLastFit() const;

    virtual void setMinIterations(int);

    static bool correctAngles(double &, double &) ;

    bool isMuonTrack(const Track &) const;

    void incrementFitStatus(enum FitterStatusType) const;

    /**
     * @brief Initialize a field cache inside a fit cache object.
     *
     * Following the shift from old-style magnetic field services to the new
     * cached implementation for thread safety, we need some additional logic
     * to create a magnetic field cache object and insert it into our fitting
     * cache object for access.
     *
     * @param[in] cache The GX2F cache objects in which to load the magnetic
     * field cache.
     */
    void initFieldCache(
      const EventContext& ctx,
      Cache & cache
    ) const;

    ToolHandle<IRIO_OnTrackCreator> m_ROTcreator {this, "RotCreatorTool", "", ""};
    ToolHandle<IRIO_OnTrackCreator> m_broadROTcreator {this, "BroadRotCreatorTool", "", ""};
    ToolHandle<IUpdator> m_updator {this, "MeasurementUpdateTool", "", ""};
    ToolHandle<IExtrapolator> m_extrapolator {this, "ExtrapolationTool", "Trk::Extrapolator/CosmicsExtrapolator", ""};
    ToolHandle<IMultipleScatteringUpdator> m_scattool {this, "MultipleScatteringTool", "Trk::MultipleScatteringUpdator/AtlasMultipleScatteringUpdator", ""};
    ToolHandle<IEnergyLossUpdator> m_elosstool {this, "EnergyLossTool", "Trk::EnergyLossUpdator/AtlasEnergyLossUpdator", ""};
    ToolHandle<IMaterialEffectsUpdator> m_matupdator {this, "MaterialUpdateTool", "", ""};
    ToolHandle<IPropagator> m_propagator {this, "PropagatorTool", "Trk::StraightLinePropagator/CosmicsPropagator", ""};
    ToolHandle<INavigator> m_navigator {this, "NavigatorTool", "Trk::Navigator/CosmicsNavigator", ""};
    ToolHandle<IResidualPullCalculator> m_residualPullCalculator {this, "ResidualPullCalculatorTool", "Trk::ResidualPullCalculator/ResidualPullCalculator", ""};
    ToolHandle<Trk::ITrkMaterialProviderTool> m_caloMaterialProvider {this, "CaloMaterialProvider", "Trk::TrkMaterialProviderTool/TrkMaterialProviderTool", ""};
    ToolHandle<IMaterialEffectsOnTrackProvider> m_calotool {this, "MuidTool", "Rec::MuidMaterialEffectsOnTrackProvider/MuidMaterialEffectsOnTrackProvider", ""};
    ToolHandle<IMaterialEffectsOnTrackProvider> m_calotoolparam {this, "MuidToolParam", "", ""};
    ToolHandle<IBoundaryCheckTool> m_boundaryCheckTool {this, "BoundaryCheckTool", "", "Boundary checking tool for detector sensitivities" };

    void throwFailedToGetTrackingGeomtry() const;
    const TrackingGeometry* trackingGeometry(Cache& cache,
                                             const EventContext& ctx) const
    {
      if (!cache.m_trackingGeometry)
        cache.m_trackingGeometry = retrieveTrackingGeometry(ctx);
      return cache.m_trackingGeometry;
    }
    const TrackingGeometry* retrieveTrackingGeometry(
      const EventContext& ctx) const
    {
      SG::ReadCondHandle<TrackingGeometry> handle(m_trackingGeometryReadKey,
                                                  ctx);
      if (!handle.isValid()) {
        throwFailedToGetTrackingGeomtry();
      }
      return handle.cptr();
    }

    SG::ReadCondHandleKey<TrackingGeometry> m_trackingGeometryReadKey{
      this,
      "TrackingGeometryReadKey",
      "AtlasTrackingGeometry",
      "Key of the TrackingGeometry conditions data."
    };

    SG::ReadCondHandleKey<AtlasFieldCacheCondObj> m_field_cache_key{
      this,
      "AtlasFieldCacheCondObj",
      "fieldCondObj",
      "Trk::GlobalChi2Fitter field conditions object key"
    };

    const AtlasDetectorID *m_DetID = nullptr;

    Gaudi::Property<bool> m_signedradius {this, "SignedDriftRadius", true};
    Gaudi::Property<bool> m_calomat {this, "MuidMat", false};
    Gaudi::Property<bool> m_extmat {this, "ExtrapolatorMaterial", true};
    Gaudi::Property<bool> m_fillderivmatrix {this, "FillDerivativeMatrix", false};
    Gaudi::Property<bool> m_printderivs {this, "PrintDerivatives", false};
    Gaudi::Property<bool> m_straightlineprop {this, "StraightLine", true};
    Gaudi::Property<bool> m_extensioncuts {this, "TRTExtensionCuts", true};
    Gaudi::Property<bool> m_sirecal {this, "RecalibrateSilicon", false};
    Gaudi::Property<bool> m_trtrecal {this, "RecalibrateTRT", false};
    Gaudi::Property<bool> m_kinkfinding {this, "KinkFinding", false};
    Gaudi::Property<bool> m_decomposesegments {this, "DecomposeSegments", true};
    Gaudi::Property<bool> m_getmaterialfromtrack {this, "GetMaterialFromTrack", true};
    Gaudi::Property<bool> m_domeastrackpar {this, "MeasuredTrackParameters", true};
    Gaudi::Property<bool> m_storemat {this, "StoreMaterialOnTrack", true};
    Gaudi::Property<bool> m_redoderivs {this, "RecalculateDerivatives", false};
    Gaudi::Property<bool> m_reintoutl {this, "ReintegrateOutliers", false};
    Gaudi::Property<bool> m_acceleration {this, "Acceleration", false};
    Gaudi::Property<bool> m_numderiv {this, "NumericalDerivs", false};
    Gaudi::Property<bool> m_fiteloss {this, "FitEnergyLoss", false};
    Gaudi::Property<bool> m_asymeloss {this, "AsymmetricEnergyLoss", true};
    Gaudi::Property<bool> m_useCaloTG {this, "UseCaloTG", false};
    Gaudi::Property<bool> m_rejectLargeNScat {this, "RejectLargeNScat", false};
    Gaudi::Property<bool> m_createSummary {this, "CreateTrackSummary", true};
    Gaudi::Property<bool> m_holeSearch {this, "DoHoleSearch", false};

    Gaudi::Property<double> m_outlcut {this, "OutlierCut", 5.0};
    Gaudi::Property<double> m_p {this, "Momentum", 0.0};
    Gaudi::Property<double> m_chi2cut {this, "TrackChi2PerNDFCut", 1.e15};
    Gaudi::Property<double> m_scalefactor {this, "TRTTubeHitCut", 2.5};
    Gaudi::Property<double> m_minphfcut {this, "MinPHFCut", 0.};

    Gaudi::Property<int> m_maxoutliers {this, "MaxOutliers", 10};
    Gaudi::Property<int> m_maxit {this, "MaxIterations", 30};
    Gaudi::Property<int> m_miniter {this, "MinimumIterations", 1};
    Gaudi::Property<int> m_fixbrem {this, "FixBrem", -1};

    ParticleMasses m_particleMasses;

    /*
     * This little volume defines the inner detector. Its exact size is set at
     * the time that the fitter object is created. We just make one and keep
     * it around to save us a few allocations.
     */
    Trk::Volume m_idVolume;

    /*
     * The following members are mutable. They keep track of the number of
     * fits that have returned with a certain status. Since this must be
     * shared across threads, we protect the array with a mutex, and we mark
     * these members as thread_safe for the ATLAS G++ plugin.
     */
    mutable std::mutex m_fit_status_lock ATLAS_THREAD_SAFE;
    mutable std::array<unsigned int, __S_MAX_VALUE> m_fit_status ATLAS_THREAD_SAFE = {};
  };
}
#endif
