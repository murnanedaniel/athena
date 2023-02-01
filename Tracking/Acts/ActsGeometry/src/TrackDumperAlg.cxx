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
  ATH_CHECK(m_truthEventName.initialize());
  ATH_CHECK(m_truthPileupEventName.initialize());

  return StatusCode::SUCCESS;
}

StatusCode TrackDumperAlg::execute(const EventContext &ctx) const {

  using namespace Acts::UnitLiterals;
  const auto& gctx = Acts::GeometryContext();

  ATH_MSG_VERBOSE(name() << "::" << __FUNCTION__);

  SG::ReadHandle<xAOD::TrackParticleContainer> tracks(m_trackName, ctx);
  SG::ReadHandle<xAOD::TruthVertexContainer> vertices(m_vertexName, ctx);
  SG::ReadHandle<xAOD::TruthEventContainer> truthEvents(m_truthEventName, ctx);
  SG::ReadHandle<xAOD::TruthPileupEventContainer> truthPileupEvents(m_truthPileupEventName, ctx);

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

  if (!truthEvents.isValid()){
    ATH_MSG_ERROR("Truth event collection named " << m_truthEventName.key() << " not found");
    return StatusCode::SUCCESS;
  } else {
    ATH_MSG_DEBUG ("Truth event collection '" << m_truthEventName.key() << "' retrieved from EventStore.");
  }

  if (!truthPileupEvents.isValid()){
    ATH_MSG_ERROR("Truth PU event collection named " << m_truthPileupEventName.key() << " not found");
    return StatusCode::SUCCESS;
  } else {
    ATH_MSG_DEBUG ("Truth PU event collection '" << m_truthPileupEventName.key() << "' retrieved from EventStore.");
  }

  // Output vector to save track parameters and associated xyz covariance block
  std::vector<std::vector<double>> outTrackInfoVector;

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
    const auto& xyzCov = freeCov.topLeftCorner<3, 3>();

    Acts::FreeVector freeParams = boundToFree * actsParams;
    std::cout << actsParams.transpose() << std::endl;
    std::cout << freeParams.transpose() << std::endl << std::endl;
    std::vector<double> paramsAndCov {freeParams(0), freeParams(1), freeParams(2), freeParams(3), freeParams(4), freeParams(5), freeParams(6), freeParams(7)};

    for (int i = 0; i < xyzCov.rows(); ++i) {
      for (int j = 0; j < xyzCov.cols(); ++j) {
        paramsAndCov.push_back(xyzCov(i, j));
      }
    }
    outTrackInfoVector.push_back(paramsAndCov);
  }

  std::vector<std::vector<double>> truthVtxPositions;

  for (const xAOD::TruthEvent* evt : *truthEvents) {
      // Print hard-scattering info
      const xAOD::TruthVertex* vtx = evt->signalProcessVertex();
      if (vtx){
        std::vector<double> pos{vtx->x(), vtx->y(), vtx->z()};
        truthVtxPositions.push_back(pos);
      }
  }

  for (const auto* evt : *truthPileupEvents) {

    if(evt->nTruthVertices()==0){
      continue;
    }
    auto* vtx = evt->truthVertex(0);
    if(!vtx){
           continue;
        }

        std::vector<double> pos{vtx->x(), vtx->y(), vtx->z()};
        truthVtxPositions.push_back(pos);

        // Make sure the first vtx in that list is in fact a PV with proton(s) as incoming particle(s)
        //
        // bool hasIncomingProton = false;
        // for (unsigned int iPIn = 0; iPIn < vtx->nIncomingParticles(); ++iPIn) {
        //   auto* incomingParticle = vtx->incomingParticle(iPIn);
        //   if (incomingParticle){
        //     int pdgId = incomingParticle->pdgId();
        //     if(pdgId == 2212 || pdgId == 21 || (std::abs(pdgId) > 0 && std::abs(pdgId) < 7)){
        //       hasIncomingProton = true;
        //     }
        //   }
        // }
        // if(!hasIncomingProton){
        //   std::cout << "NO PROTON! nIncomingParticles: " << vtx->nIncomingParticles() << std::endl;
        //   for (unsigned int iPIn = 0; iPIn < vtx->nIncomingParticles(); ++iPIn) {
        //     auto* incomingParticle = vtx->incomingParticle(iPIn);
        //     if (incomingParticle){
        //       int pdgId = incomingParticle->pdgId();
        //       std::cout << "\t" << iPIn << ". pdgId: " << pdgId << std::endl;
        //     }
        //   }
        // }
  }

  std::cout <<"Ntracks: " << outTrackInfoVector.size() << std::endl;
  std::cout << "Ntruthvertices: " << truthVtxPositions.size() << std::endl;



  return StatusCode::SUCCESS;
}
