/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

// class header include
#include "GenEventBeamEffectBooster.h"

// Athena includes
#include "AthenaKernel/RNGWrapper.h"
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
  GenEventBeamEffectBooster::GenEventBeamEffectBooster( const std::string& t,
                                                        const std::string& n,
                                                        const IInterface* p )
    : base_class(t,n,p)
  {
  }


  /** Athena algtool's Hooks */
  StatusCode GenEventBeamEffectBooster::initialize()
  {
    ATH_MSG_VERBOSE("Initializing ...");
    // prepare the RandonNumber generation
    ATH_CHECK(m_rndGenSvc.retrieve());
    m_randomEngine = m_rndGenSvc->getEngine(this, m_randomEngineName); //FIXME
    if (!m_randomEngine) {
      ATH_MSG_ERROR("Could not get random number engine from RandomNumberService. Abort.");
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }


  /** Athena algtool's Hooks */
  StatusCode GenEventBeamEffectBooster::finalize()
  {
    ATH_MSG_VERBOSE("Finalizing ...");
    return StatusCode::SUCCESS;
  }

  StatusCode GenEventBeamEffectBooster::initializeAthenaEvent()
  {
    // Refresh properties retrieved from the BeamCondSvc here if
    // necessary.  It may be more thread-safe to retrieve them directly
    // from the BeamCondSvc during initializeGenEvent though?

    // m_sigma_px_b1 = 20.E-6; // angular divergence in x of beam 1 [rad]
    // m_sigma_px_b2 = 21.E-6; // angular divergence in x of beam 2 [rad]
    // m_sigma_py_b1 = 22.E-6; // angular divergence in y of beam 1 [rad]
    // m_sigma_py_b2 = 23.E-6; // angular divergence in y of beam 2 [rad]
    // // Crossing angles, note the convention here is taken from online https://atlas-dcs.cern.ch/index.php?page=LHC_INS::INS_BPM
    // m_xing_x_b1 = 19.E-6;   // half-crossing in x for beam 1 at IP1, i.e. angle between beam 1 and z-axis
    // m_xing_x_b2 = 1.E-6;    // half-crossing in x for beam 2 at IP1, i.e. angle between beam 2 and z-axis
    // m_xing_y_b1 =  150.E-6;  // half-crossing in y for beam 1 at IP1, i.e. angle between beam 1 and z-axis
    // m_xing_y_b2 = -152.E-6; // half-crossing in y for beam 2 at IP1, i.e. angle between beam 2 and z-axis
    // m_dE = 1.1E-4; // Delta_E/E theoretical value
    // m_pbeam1 = 3500.0E3; /// @todo Get from a service
    // m_pbeam2 = 3500.0E3; /// @todo Get from a service
    // m_beam1ParticleMass =  CLHEP::proton_mass_c2;
    // m_beam2ParticleMass =  CLHEP::proton_mass_c2;

    return StatusCode::SUCCESS;
  }

  StatusCode GenEventBeamEffectBooster::initializeGenEvent(CLHEP::HepLorentzRotation& transform, const EventContext& ctx) const
  {
    // Reset the transformation
    transform = CLHEP::HepLorentzRotation(); //TODO drop this

    // Prepare the random engine
    m_randomEngine->setSeed( name(), ctx );
    CLHEP::HepRandomEngine* randomEngine(m_randomEngine->getEngine(ctx));

    // Create beam 1 and 2 momenta, including divergence and crossing-angle effects
    const double px1 = m_pbeam1 * CLHEP::RandGaussZiggurat::shoot(randomEngine, m_xing_x_b2, m_sigma_px_b2);
    const double py1 = m_pbeam1 * CLHEP::RandGaussZiggurat::shoot(randomEngine, m_xing_y_b2, m_sigma_py_b2);
    const double pz1 = CLHEP::RandGaussZiggurat::shoot(randomEngine, sqrt(m_pbeam1*m_pbeam1 - px1*px1 - py1*py1), m_dE);
    const double e1 = sqrt(px1*px1 + py1*py1 + pz1*pz1 + m_beam1ParticleMass*m_beam1ParticleMass);
    CLHEP::HepLorentzVector pp1(px1, py1, pz1, e1);

    const double px2 = -m_pbeam2 * CLHEP::RandGaussZiggurat::shoot(randomEngine, m_xing_x_b1, m_sigma_px_b1); // crossing angle & divergence
    const double py2 = -m_pbeam2 * CLHEP::RandGaussZiggurat::shoot(randomEngine, m_xing_y_b1, m_sigma_py_b1);
    const double pz2 = CLHEP::RandGaussZiggurat::shoot(randomEngine, -sqrt(m_pbeam2*m_pbeam2 - px2*px2 - py2*py2), m_dE); // longitudinal component in + direction energy smeared
    const double e2 = sqrt(px2*px2 + py2*py2 + pz2*pz2 + m_beam2ParticleMass*m_beam2ParticleMass);
    CLHEP::HepLorentzVector pp2(px2, py2, pz2, e2);

    // Now set-up rotation-boost matrix
    const CLHEP::HepLorentzVector psum = pp1 + pp2;
    const CLHEP::HepLorentzVector& dir = pp1;
    // Boost psum back on the direction of dir, adapted from bstback & fromCMframe PYTHIA8, credit to T.Sjostrand
    const double betaX = -psum.x() / psum.t();
    const double betaY = -psum.y() / psum.t();
    const double betaZ = -psum.z() / psum.t();
    const double beta2 = betaX*betaX + betaY*betaY + betaZ*betaZ;
    const double gamma = 1. / sqrt(1. - beta2);
    const double prod1 = betaX * dir.x() + betaY * dir.y() + betaZ * dir.z();
    const double prod2 = gamma * (gamma * prod1 / (1. + gamma) + dir.t());
    // Get the angle of rotation
    const CLHEP::HepLorentzVector back(dir.x() + prod2 * betaX,
                                       dir.y() + prod2 * betaY,
                                       dir.z() + prod2 * betaZ,
                                       gamma * (dir.t() + prod1));
    const double thback = back.theta();
    const double phback = back.phi();
    // Setting up two rotation matrices, via 3x3 rep. rather than even more messy EulerAngles
    const double sph = sin(phback), cph = cos(phback), sth = sin(thback), cth = cos(thback);
    const CLHEP::HepRotation rot1(CLHEP::HepRep3x3(cph, sph, 0.0, -sph, cph, 0.0, 0.0, 0.0, 1.0));
    const CLHEP::HepRotation rot2(CLHEP::HepRep3x3(cth*cph, -sph, sth*cph, cth*sph, cph, sth*sph, -sth, 0.0, cth));
    const CLHEP::HepBoost bst(psum.px()/psum.e(), psum.py()/psum.e(), psum.pz()/psum.e());
    // Combine the two rotations and the boost; matrix multiplication order matters
    const CLHEP::HepRotation rot = rot2*rot1;
    if(m_applyBoost && m_applyDivergence) {
      transform.set(rot,bst);
    }
    else {
      if(m_applyBoost) {
        transform.set(bst);
      }
      if(m_applyDivergence) {
        transform.set(rot);
      }
    }
    return StatusCode::SUCCESS;
  }

  /** modifies (displaces) the given GenEvent */
  StatusCode GenEventBeamEffectBooster::manipulate(HepMC::GenEvent& ge, const EventContext& ctx) const
  {
    // Obtain the transformation
    CLHEP::HepLorentzRotation transform = CLHEP::HepLorentzRotation();
    ATH_CHECK( initializeGenEvent(transform, ctx) );

    for( auto particleIter: ge) {
      this->boostParticle(particleIter, transform);
    }
    return StatusCode::SUCCESS;
  }

  void GenEventBeamEffectBooster::boostParticle(HepMC::GenParticlePtr p,
                                                const CLHEP::HepLorentzRotation& transform) const
  {
    // Apply the same transformation for EVERY HepMC::GenParticle
    const HepMC::FourVector &mom = p->momentum();
    CLHEP::HepLorentzVector hv(mom.px(), mom.py(), mom.pz(), mom.e()); //TODO check units
    ATH_MSG_VERBOSE("BEAMBOOST initial momentum " << hv );
    hv.transform(transform);
    ATH_MSG_VERBOSE("BEAMBOOST transformed momentum " << hv);
    p->set_momentum(HepMC::FourVector(hv.px(),hv.py(),hv.pz(),hv.e())); //TODO check units
  }

}
