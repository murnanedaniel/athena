/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

/////////////////////////////////////////////////////////////////
// StraightLineSurface.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef TRKSURFACES_STRAIGHTLINESURFACE_H
#define TRKSURFACES_STRAIGHTLINESURFACE_H
// Amg
#include "EventPrimitives/EventPrimitives.h"
#include "GeoPrimitives/GeoPrimitives.h"

// Trk
#include "CxxUtils/CachedValue.h"
#include "TrkDetDescrUtils/SharedObject.h"
#include "TrkParametersBase/ParametersT.h"
#include "TrkSurfaces/CylinderBounds.h"
#include "TrkSurfaces/NoBounds.h"
#include "TrkSurfaces/Surface.h"

class Identifier;
class MsgStream;
template<class SURFACE, class BOUNDS_CNV>
class BoundSurfaceCnv_p1;
template<class SURFACE, class BOUNDS_CNV>
class BoundSurfaceCnv_p2;

namespace Trk {

class TrkDetElementBase;
class LocalParameters;
template<int DIM, class T, class S>
class ParametersT;

/**
 @class StraightLineSurface

 Class for a StraightLineSurface in the ATLAS detector
 to describe dirft tube and straw like detectors.
 It inherits from Surface.

 @author Andreas.Salzburger@cern.ch
 @author Christos Anastopoulos (Thread safety and interface cleanup)
 @author Shaun Roe (interface cleanup)
 
 */

class StraightLineSurface : public Surface
{

public:
  static constexpr SurfaceType staticType = SurfaceType::Line;
  /**Default Constructor - needed for persistency*/
  StraightLineSurface();

  /**Copy constructor*/
  StraightLineSurface(const StraightLineSurface& slsf);

  /**Assignment operator*/
  StraightLineSurface& operator=(const StraightLineSurface& slsf);

  /**Move constructor*/
  StraightLineSurface(StraightLineSurface&& slsf) noexcept = default;

  /**Move Assignment operator*/
  StraightLineSurface& operator=(StraightLineSurface&& slsf) noexcept = default;

  /**Destructor*/
  virtual ~StraightLineSurface() = default;

  /**Constructor from Amg Transform ref (boundless surface)*/
  StraightLineSurface(const Amg::Transform3D& htrans);

  /**Constructor from Amg::Transform and bounds*/
  StraightLineSurface(const Amg::Transform3D& htrans, double radius, double halez);

  /**Constructor from TrkDetElementBase and Element identifier*/
  StraightLineSurface(
    const TrkDetElementBase& detelement,
    const Identifier& id);

  /** Copy constructor with shift*/
  StraightLineSurface(
    const StraightLineSurface& slsf,
    const Amg::Transform3D& transf);

  /** Equality operator*/
  virtual bool operator==(const Surface& sf) const override;

  // Needed to prevent ambiguities with c++20.
  bool operator==(const StraightLineSurface& cf) const;

  /** Implicit constructor*/
  virtual StraightLineSurface* clone() const override final;
  
  /** NVI uniqueClone **/
  std::unique_ptr<StraightLineSurface> uniqueClone() const;

  /** Use the Surface as a ParametersBase constructor, from local parameters -
   * charged */
  virtual Surface::ChargedTrackParametersUniquePtr createUniqueTrackParameters(
    double l1,
    double l2,
    double phi,
    double theta,
    double qop,
    std::optional<AmgSymMatrix(5)> cov = std::nullopt) const override final;

  /** Use the Surface as a ParametersBase constructor, from global parameters -
   * charged*/
  virtual Surface::ChargedTrackParametersUniquePtr createUniqueTrackParameters(
    const Amg::Vector3D& position,
    const Amg::Vector3D& momentum,
    double charge,
    std::optional<AmgSymMatrix(5)> cov = std::nullopt) const override final;

  /** Use the Surface as a ParametersBase constructor, from local parameters -
   * neutral */
  virtual NeutralTrackParametersUniquePtr createUniqueNeutralParameters(
    double l1,
    double l2,
    double phi,
    double theta,
    double qop,
    std::optional<AmgSymMatrix(5)> cov = std::nullopt) const override final;

  /** Use the Surface as a ParametersBase constructor, from global parameters -
   * neutral */
  virtual NeutralTrackParametersUniquePtr createUniqueNeutralParameters(
    const Amg::Vector3D& position,
    const Amg::Vector3D& momentum,
    double charge,
    std::optional<AmgSymMatrix(5)> cov = std::nullopt) const override final;

  /** Use the Surface as a ParametersBase constructor, from local parameters */
  template<int DIM, class T>
  std::unique_ptr<ParametersT<DIM, T, StraightLineSurface>>
  createUniqueParameters(
    double l1,
    double l2,
    double phi,
    double theta,
    double qop,
    std::optional<AmgSymMatrix(DIM)> cov = std::nullopt) const;

  /** Use the Surface as a ParametersBase constructor, from global parameters */
  template<int DIM, class T>
  std::unique_ptr<ParametersT<DIM, T, StraightLineSurface>>
  createUniqueParameters(
    const Amg::Vector3D& position,
    const Amg::Vector3D& momentum,
    double charge,
    std::optional<AmgSymMatrix(DIM)> cov = std::nullopt) const;

  /** Use the Surface as a ParametersBase constructor, from local parameters */
  template<int DIM, class T>
  ParametersT<DIM, T, StraightLineSurface> createParameters(
    double l1,
    double l2,
    double phi,
    double theta,
    double qop,
    std::optional<AmgSymMatrix(DIM)> cov = std::nullopt) const;

  /** Use the Surface as a ParametersBase constructor, from global parameters */
  template<int DIM, class T>
  ParametersT<DIM, T, StraightLineSurface> createParameters(
    const Amg::Vector3D& position,
    const Amg::Vector3D& momentum,
    double charge,
    std::optional<AmgSymMatrix(DIM)> cov = std::nullopt) const;

  /** Return the measurement frame - this is needed for alignment, in particular
     for StraightLine and Perigee Surface
      - the default implementation is the the RotationMatrix3D of the transform
   */
  virtual Amg::RotationMatrix3D measurementFrame(
    const Amg::Vector3D& glopos,
    const Amg::Vector3D& glomom) const override final;

  /** Return the surface type */
  virtual SurfaceType type() const override final;

  /** Specified for StraightLineSurface: LocalToGlobal method without dynamic
   * memory allocation */
  virtual void localToGlobal(
    const Amg::Vector2D& locp,
    const Amg::Vector3D& mom,
    Amg::Vector3D& glob) const override final;

  /** Specified for StraightLineSurface: GlobalToLocal method without dynamic
    memory allocation This method is the true global->local transformation.<br>
    makes use of globalToLocal and indicates the sign of the locR by the given
    momentum

    The calculation of the sign of the radius (or \f$ d_0 \f$) can be done as
    follows:<br> May \f$ \vec d = \vec m - \vec c \f$ denote the difference
    between the center of the line and the global position of the
    measurement/predicted state, then \f$ \vec d \f$ lies within the so called
    measurement plane. The measurement plane is determined by the two orthogonal
    vectors \f$ \vec{measY}= \vec{locZ} \f$ and \f$ \vec{measX} = \vec{measY}
    \times \frac{\vec{p}}{|\vec{p}|} \f$.<br>

    The sign of the radius (\f$ d_{0} \f$ ) is then defined by the projection of
    \f$ \vec{d} \f$ onto \f$ \vec{measX} \f$:<br> \f$ sign = -sign(\vec{d} \cdot
    \vec{measX}) \f$

    \image html SignOfDriftCircleD0.gif
  */
  virtual bool globalToLocal(
    const Amg::Vector3D& glob,
    const Amg::Vector3D& mom,
    Amg::Vector2D& loc) const override final;

  /** Special method for StraightLineSurface - providing a different z estimate
   */
  Amg::Vector3D localToGlobal(
    const Trk::LocalParameters& locpars,
    const Amg::Vector3D& glomom,
    double locZ) const;

  /** Special method for StraightLineSurface - provides the Line direction from
   * cache: speedup */
  const Amg::Vector3D& lineDirection() const;

  /** fast straight line intersection schema - standard: provides closest
     intersection and (signed) path length forceDir is to provide the closest
     forward solution

      b>mathematical motivation:</b>
      Given two lines in parameteric form:<br>
      - @f$ \vec l_{a}(\lambda) = \vec m_a + \lambda \cdot \vec e_{a} @f$ <br>
      - @f$ \vec l_{b}(\mu) = \vec m_b + \mu \cdot \vec e_{b} @f$ <br>
      the vector between any two points on the two lines is given by:
      - @f$ \vec s(\lambda, \mu) = \vec l_{b} - l_{a} = \vec m_{ab} + \mu \cdot
     \vec e_{b} - \lambda \cdot \vec e_{a}
     @f$, <br> when @f$ \vec m_{ab} = \vec m_{b} - \vec m_{a} @f$.<br>
      @f$ \vec s(\lambda_0, \mu_0) @f$  denotes the vector between the two
     closest points <br>
      @f$ \vec l_{a,0} = l_{a}(\lambda_0) @f$ and @f$ \vec l_{b,0} =
     l_{b}(\mu_0) @f$ <br> and is perpenticular to both, @f$ \vec e_{a} @f$ and
     @f$ \vec e_{b} @f$.

      This results in a system of two linear equations:<br>
      - (i) @f$ 0 = \vec s(\lambda_0, \mu_0) \cdot \vec e_a = \vec m_ab \cdot
     \vec e_a + \mu_0 \vec e_a \cdot \vec e_b - \lambda_0 @f$ <br>
      - (ii) @f$ 0 = \vec s(\lambda_0, \mu_0) \cdot \vec e_b = \vec m_ab \cdot
     \vec e_b + \mu_0  - \lambda_0 \vec e_b \cdot \vec e_a @f$ <br>

      Solving (i), (ii) for @f$ \lambda_0 @f$ and @f$ \mu_0 @f$ yields:
      - @f$ \lambda_0 = \frac{(\vec m_ab \cdot \vec e_a)-(\vec m_ab \cdot \vec
     e_b)(\vec e_a \cdot \vec e_b)}{1-(\vec e_a \cdot \vec e_b)^2} @f$ <br>
      - @f$ \mu_0 = - \frac{(\vec m_ab \cdot \vec e_b)-(\vec m_ab \cdot \vec
     e_a)(\vec e_a \cdot \vec e_b)}{1-(\vec e_a \cdot \vec e_b)^2} @f$ <br>
   */
  virtual Intersection straightLineIntersection(
    const Amg::Vector3D& pos,
    const Amg::Vector3D& dir,
    bool forceDir,
    Trk::BoundaryCheck bchk) const override final;

  /** fast straight line distance evaluation to Surface */
  virtual DistanceSolution straightLineDistanceEstimate(
    const Amg::Vector3D& pos,
    const Amg::Vector3D& dir) const override final;

  /** fast straight line distance evaluation to Surface - with bound option*/
  virtual DistanceSolution straightLineDistanceEstimate(
    const Amg::Vector3D& pos,
    const Amg::Vector3D& dir,
    bool Bound) const override final;

  /** the pathCorrection for derived classes with thickness */
  virtual double pathCorrection(const Amg::Vector3D&, const Amg::Vector3D&)
    const override final;

  /** This method checks if the provided GlobalPosition is inside the assigned
    straw radius, but no check is done whether the GlobalPosition is inside
    bounds or not. It overwrites isOnSurface from Base Class as it saves the
    time of sign determination.  */
  virtual bool isOnSurface(
    const Amg::Vector3D& glopo,
    const BoundaryCheck& bchk = true,
    double tol1 = 0.,
    double tol2 = 0.) const override final;

  /**This method returns the bounds of the Surface by reference */
  virtual const SurfaceBounds& bounds() const override final;

  /**This surface calls the iside method of the bouns */
  virtual bool insideBounds(
    const Amg::Vector2D& locpos,
    double tol1 = 0.,
    double tol2 = 0.) const override final;
  virtual bool insideBoundsCheck(
    const Amg::Vector2D& locpos,
    const BoundaryCheck& bchk) const override final;

  /** Return properly formatted class name for screen output */
  virtual std::string name() const override final;

protected: //!< data members
  template<class SURFACE, class BOUNDS_CNV>
  friend class ::BoundSurfaceCnv_p1;
  template<class SURFACE, class BOUNDS_CNV>
  friend class ::BoundSurfaceCnv_p2;
  //!< cache of the line direction (speeds up)
  CxxUtils::CachedValue<Amg::Vector3D> m_lineDirection;
  //!< bounds (shared)
  SharedObject<const CylinderBounds> m_bounds;
  //!< NoBounds as return object when no bounds are declared
  static const NoBounds s_boundless;
};
} // end of namespace

#include "TrkSurfaces/StraightLineSurface.icc"
#endif // TRKSURFACES_STRAIGHTLINESURFACE_H
