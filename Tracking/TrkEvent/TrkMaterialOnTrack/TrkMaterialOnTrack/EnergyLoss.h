/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
// EnergyLoss.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef TRKMATERIALONTRACK_ENERGYLOSS_H
#define TRKMATERIALONTRACK_ENERGYLOSS_H

#include <cassert>
#include <cmath>
#include <iosfwd>

class MsgStream;
class TrackCollectionCnv;

namespace Trk {

/** @brief This class describes energy loss material effects in the ATLAS
    tracking EDM.

    Energy loss through ionisation and/or radiation leads to a change
    (reduction) of the momentum. It uncertainty can be asymmetric in this
    class. The quantity is energy since the calculation from energy to
    momentum can be done better inside the MEFupdators (which know the
    particle hypothesis) than the DetDescr tools.

   @author Common tracking software group

*/
class EnergyLoss
{
  friend class ::TrackCollectionCnv;

public:
  //! default constructor for POOL
  EnergyLoss() = default;
  EnergyLoss(const EnergyLoss&) = default;
  EnergyLoss(EnergyLoss&&) noexcept = default;
  EnergyLoss& operator=(const EnergyLoss&) = default;
  EnergyLoss& operator=(EnergyLoss&&) noexcept = default;
  virtual ~EnergyLoss() = default;
  //! Constructor with @f$\Delta E@f$, @f$\sigma(\Delta E)@f$ and asym. errors
  EnergyLoss(double deltaE,
             double sigmaDeltaE,
             double sigmaMinusDeltaE = 0.0,
             double sigmaPlusDeltaE = 0.0);

  //! Constructor with @f$\Delta E@f$, @f$\sigma(\Delta E)@f$ and component info
  EnergyLoss(double deltaE,
             double sigmaDeltaE,
             double mean_ioni,
             double sigma_ioni,
             double mean_rad,
             double sigma_rad);

  //! Constructor with @f$\Delta E@f$, @f$\sigma(\Delta E)@f$ and component info
  EnergyLoss(double deltaE,
             double sigmaDeltaE,
             double sigmaMinusDeltaE,
             double sigmaPlusDeltaE,
             double mean_ioni,
             double sigma_ioni,
             double mean_rad,
             double sigma_rad,
             double length);

  //! Virtual constructor
  virtual EnergyLoss* clone() const;

  //! returns the @f$ \Delta E @f$
  double deltaE() const;

  //! returns the symmatric error @f$ \sigma(\Delta E) @f$
  double sigmaDeltaE() const;

  //! returns the negative side @f$ \sigma(\Delta E) @f$
  double sigmaMinusDeltaE() const;

  //! returns the positive side @f$ \sigma(\Delta E) @f$
  double sigmaPlusDeltaE() const;

  // access to eloss components
  double meanIoni() const;
  double sigmaIoni() const;
  double meanRad() const;
  double sigmaRad() const;
  double length() const;

  // update from mean values
  void update(double ioni,
              double sigi,
              double rad,
              double sigr,
              bool mpv = false);

  // update
  void update(const EnergyLoss&, bool mpv = false);

  // set
  void set(double eLoss,
           double sigde,
           double ioni,
           double sigi,
           double rad,
           double sigr);

  //! Interface method for output, can be overloaded by child classes
  virtual MsgStream& dump(MsgStream& sl) const;
  //! Interface method for output, can be overloaded by child classes
  virtual std::ostream& dump(std::ostream& sl) const;

private:
  //! @f$ \Delta E @f$        - the estimated or measured energy loss
  double m_deltaE = 0;
  //!< @f$ \sigma(\Delta E) @f$ - error on the energy loss
  double m_sigmaDeltaE = 0;
  //!< @f$ \sigma(\Delta E) @f$ - negative error on the energy loss
  double m_sigmaMinusDeltaE = 0;
  //!< @f$ \sigma(\Delta E) @f$ - positive error on the energy loss
  double m_sigmaPlusDeltaE = 0;
  // additional information about components (cache only, not persistified)
  double m_mean_ioni = 0; // mean value for ionization
  double m_sig_ioni = 0;  // sigma for ionization
  double m_mean_rad = 0;  // mean value for radiation
  double m_sig_rad = 0;   // sigma for radiation
  double m_length = 0;    // 3D length of material
};
//! Overload of << operator for MsgStream for debug output
MsgStream&
operator<<(MsgStream& sl, const EnergyLoss& eloss);

//! Overload of << operator for std::ostream for debug outputstd::ostream&
std::ostream&
operator<<(std::ostream& sl, const EnergyLoss& eloss);

} // end ns

#include "TrkMaterialOnTrack/EnergyLoss.icc"
#endif // TRKMATERIALONTRACK_ENERGYLOSS_H
