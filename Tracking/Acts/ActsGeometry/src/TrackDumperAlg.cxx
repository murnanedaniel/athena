/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "ActsGeometry/TrackDumperAlg.h"

// ATHENA
#include "GaudiKernel/EventContext.h"
#include "GaudiKernel/ISvcLocator.h"

// ACTS
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsGeometry/ActsGeometryContext.h"

// STL
#include <fstream>
#include <string>

TrackDumperAlg::TrackDumperAlg(const std::string &name,
                                           ISvcLocator *pSvcLocator)
    : AthReentrantAlgorithm(name, pSvcLocator)
{}

StatusCode TrackDumperAlg::initialize() {

  ATH_MSG_DEBUG(name() << "::" << __FUNCTION__);
  ATH_CHECK(m_trackName.initialize());
  ATH_CHECK(m_vertexName.initialize());

  return StatusCode::SUCCESS;
}

StatusCode TrackDumperAlg::execute(const EventContext &ctx) const {

  using namespace Acts::UnitLiterals;
  const auto& gctx = Acts::GeometryContext();

  ATH_MSG_VERBOSE(name() << "::" << __FUNCTION__);

  SG::ReadHandle<xAOD::TrackParticleContainer> tracks(m_trackName, ctx);
  SG::ReadHandle<xAOD::VertexContainer> vertices(m_vertexName, ctx);

  if (!tracks.isValid()){
    ATH_MSG_ERROR("Track particle collection named " << m_trackName.key() << " not found");
    return StatusCode::SUCCESS;
  } else {
    ATH_MSG_DEBUG ("Track particle collection '" << m_trackName.key() << "' retrieved from EventStore.");
  }

  if (!vertices.isValid()){
    ATH_MSG_ERROR("Truth vertex collection named " << m_vertexName.key() << " not found");
    return StatusCode::SUCCESS;
  } else {
    ATH_MSG_DEBUG ("Truth vertex collection '" << m_vertexName.key() << "' retrieved from EventStore.");
  }

  // Output vector to save track parameters and associated xyz covariance block
  std::vector<std::pair<Acts::BoundVector, Acts::ActsSymMatrix<3>>> outTrackInfoVector;

  // Loop over all track particles, convert to Acts bound parameters
  // and convert covariance from bound to free
  for(const auto* trackParticle : *tracks) {
    const auto& params = trackParticle->definingParameters();
    const auto& cov = trackParticle->definingParametersCovMatrix();

    // Convert to Acts bound parameters
    std::shared_ptr<Acts::PerigeeSurface> perigeeSurface =
      Acts::Surface::makeShared<Acts::PerigeeSurface>((trackParticle->perigeeParameters()).associatedSurface().transform());

    Acts::BoundVector actsParams;
    actsParams << params(0), params(1), params(2), params(3), params(4)*1./(1_MeV), 0.;

    // Apply some loose cuts, e.g. on d0, to reduce output size
    if (std::abs(params(0)) > m_d0Cut * 1.0_mm){
      continue;
    }

    Acts::BoundSymMatrix covMat;
    covMat << cov(0,0) , cov(0,1) , cov(0,2) , cov(0,3) , cov(0,4) *1./(1_MeV), 0
    , cov(1,0) , cov(1,1) , cov(1,2) , cov(1,3) , cov(1,4) *1./(1_MeV) , 0
    , cov(2,0) , cov(2,1) , cov(2,2) , cov(2,3) , cov(2,4) *1./(1_MeV) , 0
    , cov(3,0) , cov(3,1) , cov(3,2) , cov(3,3) , cov(3,4) *1./(1_MeV) , 0
    , cov(4,0) *1./(1_MeV) , cov(4,1) *1./(1_MeV) , cov(4,2) *1./(1_MeV) , cov(4,3) *1./(1_MeV) , cov(4,4) *1./(1_MeV*1_MeV), 0
    , 0. , 0. , 0. , 0., 0., 1.;

    Acts::BoundTrackParameters boundParams(perigeeSurface, actsParams, covMat);

    // Convert 6x6 covariance to 8x8 covariance
    auto boundToFree = perigeeSurface->boundToFreeJacobian(gctx, actsParams);
    Acts::FreeSymMatrix freeCov = boundToFree * covMat * boundToFree.transpose();

    outTrackInfoVector.push_back(std::make_pair(actsParams, freeCov.topLeftCorner<3, 3>()));

  }

  std::cout << "vertices: " << (*vertices).size() << std::endl;



  return StatusCode::SUCCESS;
}
