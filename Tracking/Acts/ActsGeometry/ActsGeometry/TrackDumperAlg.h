/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef ACTSGEOMETRY_TRACKDUMPERALG_H
#define ACTSGEOMETRY_TRACKDUMPERALG_H

// ATHENA
#include "AthenaBaseComps/AthReentrantAlgorithm.h"
#include "CxxUtils/checker_macros.h"
#include "GaudiKernel/ServiceHandle.h"
#include "Gaudi/Property.h"  /*no forward decl: typedef*/
#include "GaudiKernel/ISvcLocator.h"
#include "xAODTracking/TrackParticleContainer.h"

// ACTS
#include "Acts/EventData/TrackParameters.hpp"
// #include "Acts/Geometry/GeometryIdentifier.hpp"

// STL
#include <memory>
#include <vector>
#include <fstream>
#include <mutex>


namespace Acts {
  class TrackingGeometry;
  namespace detail {
    struct Step;
  }
}


class IActsMaterialTrackWriterSvc;

class EventContext;
class IAthRNGSvc;
class IActsExtrapolationTool;
class IActsPropStepRootWriterSvc;

class TrackDumperAlg : public AthReentrantAlgorithm {
public:
  TrackDumperAlg (const std::string& name, ISvcLocator* pSvcLocator);
  StatusCode initialize() override;
  StatusCode execute(const EventContext& ctx) const override;

private:

  SG::ReadHandleKey<xAOD::TrackParticleContainer>  m_trackName{this, "TrackParticles", "InDetTrackParticles", "Collection name for track particles"};
};

#endif // ActsGeometry_ActsExtrapolation_h
