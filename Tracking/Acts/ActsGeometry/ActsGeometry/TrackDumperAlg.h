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
#include "xAODTruth/TruthVertexContainer.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthParticle.h"
#include "xAODTracking/VertexContainer.h"

#include "xAODTruth/TruthEvent.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/TruthEventAuxContainer.h"

#include "xAODTruth/TruthPileupEvent.h"
#include "xAODTruth/TruthPileupEventContainer.h"
#include "xAODTruth/TruthPileupEventAuxContainer.h"



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
  SG::ReadHandleKey<xAOD::TruthVertexContainer>  m_vertexName{this, "TruthVertices", "TruthVertices", "Collection name for truth vertices"};
  SG::ReadHandleKey<xAOD::TruthEventContainer>  m_truthEventName{this, "TruthEvents", "TruthEvents", "Collection name for TruthEvents"};
  SG::ReadHandleKey<xAOD::TruthPileupEventContainer>  m_truthPileupEventName{this, "TruthPileupEvents", "TruthPileupEvents", "Collection name for TruthPileupEvents"};





  DoubleProperty m_d0Cut{this, "d0cut", 5, "Abs d0 cut in mm"};



};

#endif // ActsGeometry_ActsExtrapolation_h
