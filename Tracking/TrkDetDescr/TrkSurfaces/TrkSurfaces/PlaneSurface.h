/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
// PlaneSurface.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef TRKSURFACES_PLANESURFACE_H
#define TRKSURFACES_PLANESURFACE_H

// Trk
#include "TrkDetDescrUtils/SharedObject.h"
#include "TrkParametersBase/ParametersT.h"
#include "TrkSurfaces/NoBounds.h"
#include "TrkSurfaces/Surface.h"
#include "TrkSurfaces/SurfaceBounds.h"
// Amg
#include "EventPrimitives/EventPrimitives.h"
#include "GeoPrimitives/GeoPrimitives.h"

class MsgStream;
class Identifier;
template<class SURFACE, class BOUNDS_CNV>
class BoundSurfaceCnv_p1;
template<class SURFACE, class BOUNDS_CNV>
class BoundSurfaceCnv_p2;

namespace Trk {

class LocalDirection;
class LocalParameters;
class TrkDetElementBase;
class RectangleBounds;
class TriangleBounds;
class AnnulusBounds;
class TrapezoidBounds;
class RotatedTrapezoidBounds;
class DiamondBounds;
class EllipseBounds;
class CurvilinearUVT;
template<int DIM, class T, class S>
class ParametersT;

/**
 @class PlaneSurface
 Class for a planaer rectangular or trapezoidal surface in the ATLAS detector.
 It inherits from Surface.

 The Trk::PlaneSurface extends the Surface class with the possibility to convert
 in addition to local to global positions, also local to global direction (vice
 versa). The definition with of a local direciton with respect to a plane can be
 found in the dedicated Trk::LocalDirection class of the TrkEventPrimitives
 package.

 @image html PlaneSurface.gif

 @author Andreas.Salzburger@cern.ch
 @author Christos Anastopoulos (Thread safety and interface cleanup)
 @author Shaun Roe (interface cleanup)
 */

class PlaneSurface : public Surface
{
public:
  /** The surface type static constexpr */
  static constexpr SurfaceType staticType = SurfaceType::Plane;

  /** Default Constructor - needed for persistency*/
  PlaneSurface();

  /** Copy Constructor*/
  PlaneSurface(const PlaneSurface& psf) = default;

  /**Assignment operator*/
  PlaneSurface& operator=(const PlaneSurface& psf) = default;

  /** Move Constructor*/
  PlaneSurface(PlaneSurface&& psf) noexcept = default;

  /**Move assignment operator*/
  PlaneSurface& operator=(PlaneSurface&& psf) noexcept = default;

  /**Destructor*/
  virtual ~PlaneSurface() = default;

  /** Copy Constructor with shift*/
  PlaneSurface(const PlaneSurface& psf, const Amg::Transform3D& transf);

  /** Dedicated Constructor with CurvilinearUVT class */
  PlaneSurface(const Amg::Vector3D& position, const CurvilinearUVT& curvUVT);

  /** Constructor from TrkDetElementBase*/
  PlaneSurface(const TrkDetElementBase& detelement,
               const Amg::Transform3D& transf);
               
  /** Constructor from TrkDetElementBase*/
  PlaneSurface(const TrkDetElementBase& detelement);

  /** Constructor from TrkDetElementBase and Identifier in case one element
   * holds more surfaces*/
  PlaneSurface(const TrkDetElementBase& detelement,
               const Identifier& id,
               const Amg::Transform3D & transf);
               
  /** Constructor from TrkDetElementBase and Identifier in case one element
   * holds more surfaces*/
  PlaneSurface(const TrkDetElementBase& detelement,
               const Identifier& id);

  /** Constructor for planar Surface without Bounds , reference */
  PlaneSurface(const Amg::Transform3D& htrans);

  
  /** Constructor for Rectangular Planes*/
  PlaneSurface(const Amg::Transform3D & htrans, double halephi, double haleta);

  /** Constructor for Trapezoidal Planes*/
  PlaneSurface(const Amg::Transform3D & htrans,
               double minhalephi,
               double maxhalephi,
               double haleta);

  /** Constructor for Planes with provided RectangleBounds - ownership of bounds
   * is passed*/
  PlaneSurface(const Amg::Transform3D & htrans, RectangleBounds* rbounds);

  /** Constructor for Planes with provided TriangleBounds - ownership of bounds
   * is passed*/
  PlaneSurface(const Amg::Transform3D & htrans, TriangleBounds* rbounds);

  /** Constructor for Planes with provided AnnulusBounds - ownership of bounds
   * is passed*/
  PlaneSurface(const Amg::Transform3D& htrans, AnnulusBounds* rbounds);

  /** Constructor for Planes with provided TrapezoidBounds - ownership of bounds
   * is passed*/
  PlaneSurface(const Amg::Transform3D& htrans, TrapezoidBounds* rbounds);

  /** Constructor for Planes with provided RotatedTrapezoidBounds - ownership of
   * bounds is passed*/
  PlaneSurface(const Amg::Transform3D& htrans, RotatedTrapezoidBounds* rbounds);

  /** Constructor for Planes with provided DiamondBounds - ownership of bounds
   * is passed*/
  PlaneSurface(const Amg::Transform3D& htrans, DiamondBounds* rbounds);

  /** Constructor for Planes with provided EllipseBounds - ownership of bounds
   * is passed*/
  PlaneSurface(const Amg::Transform3D& htrans, EllipseBounds* rbounds);

  /** Constructor for Planes with shared object*/
  PlaneSurface(const Amg::Transform3D& htrans,
               Trk::SharedObject<const Trk::SurfaceBounds>& sbounds);

  /**Equality operator*/
  virtual bool operator==(const Surface& sf) const override;

  // Needed to prevent ambiguities with c++20.
  bool operator==(const PlaneSurface& cf) const;

  /**Virtual constructor*/
  virtual PlaneSurface* clone() const override;
  
   /** NVI uniqueClone method */
  std::unique_ptr<PlaneSurface>uniqueClone() const;

  /** Return the surface type */
  virtual SurfaceType type() const override final;

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
    double oop,
    std::optional<AmgSymMatrix(5)> cov = std::nullopt) const override final;

  /** Use the Surface as a ParametersBase constructor, from global parameters
   * - neutral */
  virtual NeutralTrackParametersUniquePtr createUniqueNeutralParameters(
    const Amg::Vector3D& position,
    const Amg::Vector3D& momentum,
    double charge = 0.,
    std::optional<AmgSymMatrix(5)> cov = std::nullopt) const override final;

  /** Use the Surface as a ParametersBase constructor, from local parameters */
  template<int DIM, class T>
  std::unique_ptr<ParametersT<DIM, T, PlaneSurface>> createUniqueParameters(
    double l1,
    double l2,
    double phi,
    double theta,
    double qop,
    std::optional<AmgSymMatrix(DIM)> cov = std::nullopt) const;

  /** Use the Surface as a ParametersBase constructor, from global parameters */
  template<int DIM, class T>
  std::unique_ptr<ParametersT<DIM, T, PlaneSurface>> createUniqueParameters(
    const Amg::Vector3D& position,
    const Amg::Vector3D& momentum,
    double charge,
    std::optional<AmgSymMatrix(DIM)> cov = std::nullopt) const;

  /** Use the Surface as a ParametersBase constructor, from local parameters */
  template<int DIM, class T>
  ParametersT<DIM, T, PlaneSurface> createParameters(
    double l1,
    double l2,
    double phi,
    double theta,
    double qop,
    std::optional<AmgSymMatrix(DIM)> cov = std::nullopt) const;

  /** Use the Surface as a ParametersBase constructor, from global parameters */
  template<int DIM, class T>
  ParametersT<DIM, T, PlaneSurface> createParameters(
    const Amg::Vector3D& position,
    const Amg::Vector3D& momentum,
    double charge,
    std::optional<AmgSymMatrix(DIM)> cov = std::nullopt) const;

  /**This method returns the bounds by reference, static NoBounds in case of no
   * boundaries*/
  virtual const SurfaceBounds& bounds() const override final;

  /**This method calls the inside() method of the Bounds*/
  virtual bool insideBounds(const Amg::Vector2D& locpos,
                            double tol1 = 0.,
                            double tol2 = 0.) const override;

  virtual bool insideBoundsCheck(
    const Amg::Vector2D& locpos,
    const BoundaryCheck& bchk) const override final;

  /** This method returns true if the GlobalPosition is on the Surface for both,
    within or without check of whether the local position is inside boundaries
    or not */
  virtual bool isOnSurface(const Amg::Vector3D& glopo,
                           const BoundaryCheck& bchk = true,
                           double tol1 = 0.,
                           double tol2 = 0.) const override final;

  /** Specified for PlaneSurface: LocalToGlobal method without dynamic memory
   * allocation */
  virtual void localToGlobal(const Amg::Vector2D& locp,
                             const Amg::Vector3D& mom,
                             Amg::Vector3D& glob) const override final;

  /** Specified for PlaneSurface: GlobalToLocal method without dynamic memory
   * allocation - boolean checks if on surface
   */
  virtual bool globalToLocal(const Amg::Vector3D& glob,
                             const Amg::Vector3D& mom,
                             Amg::Vector2D& loc) const override final;

  /** This method transforms a local direction wrt the plane to a global
   * direction */
  void localToGlobalDirection(const Trk::LocalDirection& locdir,
                              Amg::Vector3D& globdir) const;

  /**This method transforms the global direction to a local direction wrt the
   * plane */
  void globalToLocalDirection(const Amg::Vector3D& glodir,
                              Trk::LocalDirection& locdir) const;

  /** fast straight line intersection schema - standard: provides closest
     intersection and (signed) path length forceDir is to provide the closest
     forward solution

      <b>mathematical motivation:</b>

      the equation of the plane is given by: <br>
      @f$ \vec n \cdot \vec x = \vec n \cdot \vec p,@f$ <br>
      where @f$ \vec n = (n_{x}, n_{y}, n_{z})@f$ denotes the normal vector of
     the plane,
      @f$ \vec p = (p_{x}, p_{y}, p_{z})@f$ one specific point on the plane and
     @f$ \vec x = (x,y,z) @f$ all possible points on the plane.<br> Given a line
     with:<br>
      @f$ \vec l(u) = \vec l_{1} + u \cdot \vec v @f$, <br>
      the solution for @f$ u @f$ can be written:
      @f$ u = \frac{\vec n (\vec p - \vec l_{1})}{\vec n \vec v}@f$ <br>
      If the denominator is 0 then the line lies:
      - either in the plane
      - perpenticular to the normal of the plane

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

  /** Return properly formatted class name for screen output */
  virtual std::string name() const override;

protected: //!< data members
  template<class SURFACE, class BOUNDS_CNV>
  friend class ::BoundSurfaceCnv_p1;
  template<class SURFACE, class BOUNDS_CNV>
  friend class ::BoundSurfaceCnv_p2;

  SharedObject<const SurfaceBounds> m_bounds; //!< bounds (shared)
  //!< NoBounds as return object when no bounds are declared
  static const NoBounds s_boundless;
};

} // end of namespace

#include "TrkSurfaces/PlaneSurface.icc"
#endif // TRKSURFACES_PLANESURFACE_H
