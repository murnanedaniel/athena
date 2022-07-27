/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
// CylinderSurface.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef TRKSURFACES_CYLINDERSURFACE_H
#define TRKSURFACES_CYLINDERSURFACE_H

// Trk
#include "TrkDetDescrUtils/SharedObject.h"
#include "TrkEventPrimitives/LocalParameters.h"
#include "TrkEventPrimitives/ParamDefs.h"
#include "TrkParametersBase/ParametersT.h"
#include "TrkSurfaces/CylinderBounds.h"
#include "TrkSurfaces/Surface.h"
// Amg
#include "EventPrimitives/EventPrimitives.h"
#include "GeoPrimitives/GeoPrimitives.h"

class MsgStream;
template<class SURFACE, class BOUNDS_CNV>
class BoundSurfaceCnv_p1;
template<class SURFACE, class BOUNDS_CNV>
class BoundSurfaceCnv_p2;

namespace Trk {

template<int DIM, class T, class S>
class ParametersT;

/**
 @class CylinderSurface
 Class for a CylinderSurface in the ATLAS detector.
 It inherits from Surface.

 The cylinder surface has a special role in the TrackingGeometry,
 since it builds the surfaces of all TrackingVolumes at container level,
 hence, to optimize speed in global to local transformations,
 constructors w/o transform is possible,
 assumint the identity transform to.

 @todo update for new Possibility of CylinderBounds

 @image html CylinderSurface.gif

 @author Andreas.Salzburger@cern.ch
 @author Christos Anastopoulos (Thread safety and interface cleanup)
 @author Shaun Roe (interface cleanup)
 */

class CylinderSurface : public Surface
{

public:
  /** The surface type static constexpr */
  static constexpr SurfaceType staticType = SurfaceType::Cylinder;

  /**Default Constructor*/
  CylinderSurface();

  /**Copy constructor */
  CylinderSurface(const CylinderSurface& csf);

  /**Assignment operator*/
  CylinderSurface& operator=(const CylinderSurface& csf);

  /**Move constructor */
  CylinderSurface(CylinderSurface&& csf) noexcept = default;

  /**Move Assignment operator*/
  CylinderSurface& operator=(CylinderSurface&& csf) noexcept = default;

  /**Destructor*/
  virtual ~CylinderSurface() = default;

  /**Constructor from EigenTransform, radius and halflength*/
  CylinderSurface(const Amg::Transform3D& htrans,
                  double radius,
                  double hlength);

  /**Constructor from EigenTransform, radius halfphi, and halflength*/
  CylinderSurface(const Amg::Transform3D& htrans,
                  double radius,
                  double hphi,
                  double hlength);

  /**Constructor from EigenTransform and CylinderBounds
    - ownership of the bounds is passed */
  CylinderSurface(const Amg::Transform3D& htrans, CylinderBounds* cbounds);

  /**Constructor from EigenTransform from unique_ptr.
     - bounds is not set */
  CylinderSurface(const Amg::Transform3D& htrans);

  /** Constructor from radius and halflength - speed optimized for concentric
   * volumes */
  CylinderSurface(double radius, double hlength);

  /**Constructor from radius halfphi, and halflength - speed optimized for
   * concentric volumes */
  CylinderSurface(double radius, double hphi, double hlength);

  /**Constructor from EigenTransform and CylinderBounds
      - ownership of the bounds is passed
      - speed optimized fron concentric volumes */
  CylinderSurface(CylinderBounds* cbounds);

   /**Copy constructor with shift */
  CylinderSurface(const CylinderSurface& csf, const Amg::Transform3D& transf);

  /**Equality operator*/
  virtual bool operator==(const Surface& sf) const override;

  // Needed to prevent ambiguities with c++20.
  bool operator==(const CylinderSurface& cf) const;

  /**Implicit Constructor*/
  virtual CylinderSurface* clone() const override;
  
  std::unique_ptr<CylinderSurface> uniqueClone() const;
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
  std::unique_ptr<ParametersT<DIM, T, CylinderSurface>> createUniqueParameters(
    double l1,
    double l2,
    double phi,
    double theta,
    double qop,
    std::optional<AmgSymMatrix(DIM)> cov = std::nullopt) const;

  /** Use the Surface as a ParametersBase constructor, from global parameters */
  template<int DIM, class T>
  std::unique_ptr<ParametersT<DIM, T, CylinderSurface>> createUniqueParameters(
    const Amg::Vector3D& position,
    const Amg::Vector3D& momentum,
    double charge,
    std::optional<AmgSymMatrix(DIM)> cov = std::nullopt) const;

  /** Use the Surface as a ParametersBase constructor, from local parameters */
  template<int DIM, class T>
  ParametersT<DIM, T, CylinderSurface> createParameters(
    double l1,
    double l2,
    double phi,
    double theta,
    double qop,
    std::optional<AmgSymMatrix(DIM)> cov = std::nullopt) const;

  /** Use the Surface as a ParametersBase constructor, from global parameters */
  template<int DIM, class T>
  ParametersT<DIM, T, CylinderSurface> createParameters(
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

  /** Returns a global reference point:
    For the Cylinder this is @f$ (R*cos(\phi), R*sin(\phi),0)*transform() @f$
    Where  @f$ \phi @f$ denotes the averagePhi() of the cylinderBounds.
   */
  virtual const Amg::Vector3D& globalReferencePoint() const override final;

  // using from the base class
  using Trk::Surface::normal;
  /**Return method for surface normal information
     at a given local point, overwrites the normal() from base class.*/
  virtual Amg::Vector3D normal(const Amg::Vector2D& locpo) const override final;

  /**Return method for the rotational symmetry axis - the z-Axis of the
   * HepTransform */
  virtual const Amg::Vector3D& rotSymmetryAxis() const;

  /**This method returns the CylinderBounds by reference
   (NoBounds is not possible for cylinder)*/
  virtual const CylinderBounds& bounds() const override final;

  bool hasBounds() const;

  /**This method calls the inside method of CylinderBounds*/
  virtual bool insideBounds(const Amg::Vector2D& locpos,
                            double tol1 = 0.,
                            double tol2 = 0.) const override;

  virtual bool insideBoundsCheck(const Amg::Vector2D& locpos,
                                 const BoundaryCheck& bchk) const override final;

  /** Specialized for CylinderSurface : LocalParameters to Vector2D */
  virtual Amg::Vector2D localParametersToPosition(
    const LocalParameters& locpars) const override final;

  /** Specialized for CylinderSurface : LocalToGlobal method without dynamic
   * memory allocation */
  virtual void localToGlobal(const Amg::Vector2D& locp,
                             const Amg::Vector3D& mom,
                             Amg::Vector3D& glob) const override;

  /** Specialized for CylinderSurface : GlobalToLocal method without dynamic
   * memory allocation - boolean checks if on surface */
  virtual bool globalToLocal(const Amg::Vector3D& glob,
                             const Amg::Vector3D& mom,
                             Amg::Vector2D& loc) const override;

  /** This method returns true if the GlobalPosition is on the Surface for both,
    within or without check of whether the local position is inside boundaries
    or not */
  virtual bool isOnSurface(const Amg::Vector3D& glopo,
                           const BoundaryCheck& bchk = true,
                           double tol1 = 0.,
                           double tol2 = 0.) const override;

  /** fast straight line intersection schema - provides closest intersection and
     (signed) path length

      <b>mathematical motivation:</b>

      The calculation will be done in the 3-dim frame of the cylinder,
      i.e. the symmetry axis of the cylinder is the z-axis, x- and y-axis are
     perpenticular to the the z-axis. In this frame the cylinder is centered
     around the origin. Therefore the two points describing the line have to be
     first recalculated into the new frame. Suppose, this is done, the
     intersection is straight forward:<br> may @f$p_{1}=(p_{1x}, p_{1y},
     p_{1z}), p_{2}=(p_{2x}, p_{2y}, p_{2z}) @f$the two points describing the
     3D-line, then the line in the \f$x-y@f$plane can be written as
      @f$y=kx+d\f$, where @f$k =\frac{p_{2y}-p_{1y}}{p_{2x}-p_{1x}}@f$such as
     @f$d=\frac{p_{2x}p_{1y}-p_{1x}p_{2y}}{p_{2x}-p_{1x}},\f$<br> and intersects
     with the corresponding circle
     @f$x^{2}+y^{2} = R^{2}. @f$<br> The solutions can then be found by a simple
     quadratic equation and reinsertion into the line equation.
  */
  virtual Intersection straightLineIntersection(
    const Amg::Vector3D& pos,
    const Amg::Vector3D& dir,
    bool forceDir = false,
    Trk::BoundaryCheck bchk = false) const override final;

  /** fast distance to Surface */
  virtual DistanceSolution straightLineDistanceEstimate(
    const Amg::Vector3D& pos,
    const Amg::Vector3D& dir) const override;

  /** fast distance to Surface - with bounds directive*/
  virtual DistanceSolution straightLineDistanceEstimate(
    const Amg::Vector3D& pos,
    const Amg::Vector3D& dir,
    bool bound) const override;

  /** the pathCorrection for derived classes with thickness */
  virtual double pathCorrection(const Amg::Vector3D& pos,
                                const Amg::Vector3D& mom) const override;

  /** Return properly formatted class name for screen output */
  virtual std::string name() const override;

protected: //!< data members
  template<class SURFACE, class BOUNDS_CNV>
  friend class ::BoundSurfaceCnv_p1;
  template<class SURFACE, class BOUNDS_CNV>
  friend class ::BoundSurfaceCnv_p2;

  //!< bounds (shared)
  SharedObject<const CylinderBounds> m_bounds;
  //!< The global reference point (== a point on the  surface)
  CxxUtils::CachedUniquePtr<Amg::Vector3D> m_referencePoint;
  //!< The rotational symmetry axis
  CxxUtils::CachedUniquePtr<Amg::Vector3D> m_rotSymmetryAxis;
};

} // end of namespace
#include "TrkSurfaces/CylinderSurface.icc"

#endif // TRKSURFACES_CYLINDERSURFACE_H
