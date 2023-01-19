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

  return StatusCode::SUCCESS;
}

StatusCode TrackDumperAlg::execute(const EventContext &ctx) const {

  ATH_MSG_VERBOSE(name() << "::" << __FUNCTION__);


  std::cout << "HALLO" << std::endl;

  SG::ReadHandle<xAOD::TrackParticleContainer> tracks(m_trackName, ctx);

  if (!tracks.isValid()){
    ATH_MSG_ERROR("Track particle collection named " << m_trackName.key() << " not found");
    return StatusCode::SUCCESS;
  } else {
    ATH_MSG_DEBUG ("Track particle collection '" << m_trackName.key() << "' retrieved from EventStore.");
  }


  for(const auto* trackParticle : *tracks) {
    std::cout << trackParticle->definingParameters().transpose() << std::endl;
    std::cout << trackParticle->definingParametersCovMatrix() << std::endl;
  }

  return StatusCode::SUCCESS;
}
