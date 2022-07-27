/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
// MaterialEffectsBase.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef TRKMATERIALONTRACK_MATERIALEFFECTSBASE_H
#define TRKMATERIALONTRACK_MATERIALEFFECTSBASE_H

#include "TrkSurfaces/SurfaceHolders.h"
//
#include <bitset>
#include <iosfwd>
#include <memory>

class MsgStream;
class TrackCollectionCnv;
class MaterialEffectsBaseCnv_p1;
class MaterialEffectsBaseCnv_p2;
class MaterialEffectsOnTrackCnv_p2;
class ScatteringAngleOnTrackCnv_p1;

namespace Trk {

/** @brief base class to integrate material effects on Trk::Track in a
    flexible way.
    It holds the pointer to an associated surface (free or from DetStore)
    and the layer thickness which was used in estimated the material
    effects. In case they are external (Calo measurement) they can be 0.
    @author Wolfgang.Liebig <http://consult.cern.ch/xwho/people/>
*/
class MaterialEffectsBase : public SurfacePtrHolderDetEl
{
public:
  enum MaterialEffectsType
  {
    //! contains material effects due to multiple scattering
    ScatteringEffects = 0,
    //! contains energy loss corrections
    // trouble: avoid name clash Trk::Energyloss vs MET::EnergyLoss
    EnergyLossEffects = 1,
    //! contains only thickness, needs M.E.Updator to calculate effects
    MaterialThickness = 2,
    //! contains q/p covariance noise term
    BremPoint = 3,
    //! contains energy loss correction based on Calo measurement
    UsesMeasurement = 4,
    //! contains values obtained by fitting the scatterer or e-loss
    FittedMaterialEffects = 5,
    //! new category
    Unknown = 6,
    NumberOfMaterialEffectsTypes = 7
    // WARNING need to edit MaterialEffectsOnTrack.cxx if these enums change
  };

  //! default constructor for POOL
  MaterialEffectsBase();
  /** @brief base class constructor with information common to all types of
      material effects.
      @param[in] thicknessInX0 is the actually traversed materia t in terms of
     radiation length x0, including all projective and bending corrections.
      @param[in] assocSurf reference to the Surface the material effects are
     expressed at
      @param[in] typePattern bitset to describe and identify the type of
     material effects at base-class level
  */
  MaterialEffectsBase(
    double thicknessInX0,
    const Surface& assocSurf,
    const std::bitset<MaterialEffectsBase::NumberOfMaterialEffectsTypes>&
      typePattern);
  //! destructor. Being virtual forces derived classes to call also this base
  //! destructor
  virtual ~MaterialEffectsBase() = default;

  //! Virtual constructor
  virtual MaterialEffectsBase* clone() const = 0;

  //! NVI uniqueClone
  std::unique_ptr<MaterialEffectsBase> uniqueClone() const
  {
    return std::unique_ptr<MaterialEffectsBase>(clone());
  }

  //! returns the actually traversed material @f$ t/X_0 @f$. Leave 0.0 for
  //! external ME.
  double thicknessInX0() const;

  //! returns the surface to which these m.eff. are associated.
  const Surface& associatedSurface() const;

  /** @brief returns the flags (bits) which types of ME are present
   *
   * Use this method to find out if the ME is of a certain type:
   * i.e. if ( mefot->type(MaterialEffectsBase::EnergyLoss) { //etc }
   *
   * @return true if the MaterialEffectsBase is of this type
   */
  bool type(const MaterialEffectsType& type) const;

  //! returns a string with the type of the object
  std::string dumpType() const;

  //! Interface method for output, can be overloaded by child classes
  virtual MsgStream& dump(MsgStream& sl) const;
  //! Interface method for output, can be overloaded by child classes* */
  virtual std::ostream& dump(std::ostream& sl) const;

protected:
  //! allows POOL converter to recreate transient links to DetStore
  virtual void setValues(const Surface* assocSurface);
  //! copy constructor
  MaterialEffectsBase(const MaterialEffectsBase& rhs) = default;
  //! Assignment operator
  MaterialEffectsBase& operator=(const MaterialEffectsBase& rhs) = default;
  //! move constructor
  MaterialEffectsBase(MaterialEffectsBase&& rhs) = default;
  //! move  Assignment operator
  MaterialEffectsBase& operator=(MaterialEffectsBase&& rhs) = default;

private:
  friend class ::MaterialEffectsBaseCnv_p1;
  friend class ::MaterialEffectsBaseCnv_p2;
  friend class ::MaterialEffectsOnTrackCnv_p2;
  friend class ::TrackCollectionCnv;
  friend class ::ScatteringAngleOnTrackCnv_p1;

  //! @f$ t/X_0    @f$  - the traversed thickness in RadiationLengths
  double m_tInX0{};
  long m_typeFlags{};
};

//! Overload of << operator for MsgStream for debug output
MsgStream&
operator<<(MsgStream& sl, const MaterialEffectsBase& meb);

//! Overload of << operator for std::ostream for debug output
std::ostream&
operator<<(std::ostream& sl, const MaterialEffectsBase& meb);

} // end ns

#include "TrkMaterialOnTrack/MaterialEffectsBase.icc"
#endif // TRKMATERIALONTRACK_MATERIALEFFECTSBASE_H
