/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

// class header include
#include "VertexBeamCondPositioner.h"

// For the speed of light
#include "GaudiKernel/PhysicalConstants.h"

// CLHEP includes
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Transform3D.h"

// RandomNumber generator
#include "AthenaKernel/RNGWrapper.h"
#include "CLHEP/Random/RandGaussZiggurat.h"
#include <cmath>

namespace Simulation
{

  /** Constructor */
  VertexBeamCondPositioner::VertexBeamCondPositioner( const std::string& t,
                                                      const std::string& n,
                                                      const IInterface* p )
    : base_class(t,n,p)
  {
  }

  /** Athena algtool's Hooks */
  StatusCode VertexBeamCondPositioner::initialize()
  {
    ATH_MSG_VERBOSE("Initializing ...");

    // retrieve the random number service
    ATH_CHECK(m_rndGenSvc.retrieve());
    m_randomEngine = m_rndGenSvc->getEngine(this, m_randomEngineName);
    if (!m_randomEngine) {
      ATH_MSG_ERROR("Could not get random number engine from RandomNumberService. Abort.");
      return StatusCode::FAILURE;
    }
    ATH_CHECK(m_beamSpotKey.initialize());
    // everything set up properly
    return StatusCode::SUCCESS;
  }

  /** Athena algtool's Hooks */
  StatusCode VertexBeamCondPositioner::finalize()
  {
    ATH_MSG_VERBOSE("Finalizing ...");
    return StatusCode::SUCCESS;
  }

  /** computes the vertex displacement */
  CLHEP::HepLorentzVector *VertexBeamCondPositioner::generate(const EventContext& ctx) const
  {
    // Prepare the random engine
    m_randomEngine->setSeed( name(), ctx );
    CLHEP::HepRandomEngine* randomEngine(m_randomEngine->getEngine(ctx));
    SG::ReadCondHandle<InDet::BeamSpotData> beamSpotHandle { m_beamSpotKey, ctx };
    // See jira issue ATLASSIM-497 for an explanation of why calling
    // shoot outside the CLHEP::HepLorentzVector constructor is
    // necessary/preferable.
    float vertexX = CLHEP::RandGaussZiggurat::shoot(randomEngine)*beamSpotHandle->beamSigma(0);
    float vertexY = CLHEP::RandGaussZiggurat::shoot(randomEngine)*beamSpotHandle->beamSigma(1);
    float vertexZ = CLHEP::RandGaussZiggurat::shoot(randomEngine)*beamSpotHandle->beamSigma(2);
    // calculate the vertexSmearing
    CLHEP::HepLorentzVector *vertexSmearing =
      new CLHEP::HepLorentzVector( vertexX, vertexY, vertexZ, 0. );

    // (1) code from: Simulation/G4Atlas/G4AtlasUtilities/VertexPositioner.cxx
    const double tx = tan( beamSpotHandle->beamTilt(1) );
    const double ty = tan( beamSpotHandle->beamTilt(0) );

    const double sqrt_abc = sqrt(1. + tx*tx + ty*ty);
    const double sqrt_fgh = sqrt(1. + ty*ty);

    const double a = ty/sqrt_abc;
    const double b = tx/sqrt_abc;
    const double c = 1./sqrt_abc;

    HepGeom::Point3D<double> from1(0,0,1);
    HepGeom::Point3D<double> to1(a,b,c);

    const double f = 1./sqrt_fgh;
    const double g = 0.;
    const double h = -(ty)/sqrt_fgh;

    HepGeom::Point3D<double> from2(1,0,0);
    HepGeom::Point3D<double> to2(f,g,h);

    // first rotation, then translation
    HepGeom::Transform3D transform(
        HepGeom::Rotate3D(from1, from2, to1, to2).getRotation(),
        CLHEP::Hep3Vector( beamSpotHandle->beamPos().x(),
                           beamSpotHandle->beamPos().y(),
                           beamSpotHandle->beamPos().z() )
        );

    // FIXME: don't use endl in MsgStream printouts
    ATH_MSG_VERBOSE("BeamSpotSvc reported beam position as " << beamSpotHandle->beamPos() << std::endl
                    << "\tWidth is (" << beamSpotHandle->beamSigma(0)
                    << ", " << beamSpotHandle->beamSigma(1) << ", "
                    << beamSpotHandle->beamSigma(2) << ")" << std::endl
                    << "\tTilts are " << beamSpotHandle->beamTilt(0) << " and " << beamSpotHandle->beamTilt(1) << std::endl
                    << "\tVertex Position before transform: " << *vertexSmearing);

    // update with the tilt
    *vertexSmearing = transform * HepGeom::Point3D<double>(*vertexSmearing);

    // See if we were asked to do time smearing as well
    if (m_timeSmearing){
      /* This is ballpark code courtesy of Brian Amadio.  He provided some functions based on beam parameters.
         He provided a little trick for pulling out the beam bunch width as well.  Hard coding the crossing angle
         parameter for the time being, as the beam spot service doesn't really provide that yet.  */
      double bunch_length_z = (std::sqrt(2)*beamSpotHandle->beamSigma(2))/0.9; // 0.9 is the crossing angle reduction factor
      //    double tLimit = 2.*(bunch_length_z+bunch_length_z)/Gaudi::Units::c_light;
      //    TF1 func = TF1("func","[0]*exp((-([3]-299792458*x)^2*[2]^2-([3]+299792458*x)^2*[1]^2)/(2*[1]^2*[2]^2))",-1*tLimit,tLimit);
      //    func.SetParameter(0,Gaudi::Units::c_light/(M_PI*bunch_length_z*bunch_length_z));
      //    func.SetParameter(1,bunch_length_z);
      //    func.SetParameter(2,bunch_length_z);
      //    func.SetParameter(3,vertexSmearing->z());
      //    double time_offset = func.GetRandom();

      // Time should be set in units of distance, which is a little funny
      double time_offset = CLHEP::RandGaussZiggurat::shoot(
          randomEngine, vertexSmearing->z()/Gaudi::Units::c_light,
          bunch_length_z/Gaudi::Units::c_light );

      vertexSmearing->setT( vertexSmearing->t() + time_offset*Gaudi::Units::c_light );
    }

    // and return it
    return vertexSmearing;
  }

} // namespace Simulation
