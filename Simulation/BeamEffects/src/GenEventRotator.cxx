/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

// class header include
#include "GenEventRotator.h"

// Athena includes
#include "CLHEP/Random/RandGaussZiggurat.h"

// CLHEP includes
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

// HepMC includes
#include "AtlasHepMC/GenEvent.h"
#include "AtlasHepMC/GenParticle.h"

namespace Simulation
{
  /** Constructor **/
  GenEventRotator::GenEventRotator( const std::string& t,
                                    const std::string& n,
                                    const IInterface* p )
    : base_class(t,n,p)
  {
  }

  /** Athena algtool's Hooks */
  StatusCode GenEventRotator::initialize()
  {
    ATH_MSG_VERBOSE("Initializing ...");
    return StatusCode::SUCCESS;
  }

  /** Athena algtool's Hooks */
  StatusCode GenEventRotator::finalize()
  {
    ATH_MSG_VERBOSE("Finalizing ...");
    return StatusCode::SUCCESS;
  }

  StatusCode GenEventRotator::initializeAthenaEvent()
  {
    return StatusCode::SUCCESS;
  }

  StatusCode GenEventRotator::initializeGenEvent(CLHEP::HepLorentzRotation& transform, const EventContext&) const
  {
    // Reset the transformation
    transform = CLHEP::HepLorentzRotation(); //TODO drop this

    // Setting up three rotation matrices, via 3x3 rep. rather than even more messy EulerAngles
    const double sinX = sin(m_xangle*CLHEP::deg);
    const double cosX = cos(m_xangle*CLHEP::deg);
    const CLHEP::HepRotation rotx(CLHEP::HepRep3x3(1.0,0.0,0.0, 0.0,cosX,-sinX, 0.0,sinX,cosX));
    const double sinY = sin(m_yangle*CLHEP::deg);
    const double cosY = cos(m_yangle*CLHEP::deg);
    const CLHEP::HepRotation roty(CLHEP::HepRep3x3(cosY,0.0,sinY, 0.0,1.0,0.0, -sinY,0.0,cosY));
    const double sinZ = sin(m_zangle*CLHEP::deg);
    const double cosZ = cos(m_zangle*CLHEP::deg);
    const CLHEP::HepRotation rotz(CLHEP::HepRep3x3(cosZ,-sinZ,0.0, sinZ,cosZ,0.0, 0.0,0.0,1.0));
    // Combine the three rotations - order of trf execution is unimportant as they are orthogonal.
    const CLHEP::HepRotation rot = rotx*roty*rotz;
    transform.set(rot);
    return StatusCode::SUCCESS;
  }

  /** modifies (displaces) the given GenEvent */
  StatusCode GenEventRotator::manipulate(HepMC::GenEvent& ge, const EventContext& ctx) const
  {
    // Obtain the transformation
    CLHEP::HepLorentzRotation transform = CLHEP::HepLorentzRotation();
    ATH_CHECK( initializeGenEvent(transform, ctx) );

    for(auto particleIter:  ge) {
      rotateParticle(particleIter, transform);
    }
    return StatusCode::SUCCESS;
  }

  void GenEventRotator::rotateParticle(HepMC::GenParticlePtr p,
                                       const CLHEP::HepLorentzRotation& transform) const
  {
    // Apply the same transformation for EVERY HepMC::GenParticle
    const HepMC::FourVector mom = p->momentum();
    CLHEP::HepLorentzVector hv(mom.px(), mom.py(), mom.pz(), mom.e()); //TODO check units
    ATH_MSG_VERBOSE("initial momentum " << hv );
    hv.transform(transform);
    ATH_MSG_VERBOSE("transformed momentum " << hv);
    p->set_momentum(HepMC::FourVector(hv.px(),hv.py(),hv.pz(),hv.e())); //TODO check units
  }

}
