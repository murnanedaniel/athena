/*
  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
*/

// InDet
#include "InDetTrackingGeometry/RobustTrackingGeometryBuilderCond.h"
// EnvelopeDefinitionService
#include "SubDetectorEnvelopes/IEnvelopeDefSvc.h"
// Trk interfaces
#include "TrkDetDescrInterfaces/ILayerArrayCreator.h"
#include "TrkDetDescrInterfaces/ITrackingVolumeCreator.h"
// Athena
#include "AthenaKernel/IOVInfiniteRange.h"
#include "CxxUtils/checker_macros.h"
// Gaudi
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SystemOfUnits.h"

// constructor
InDet::RobustTrackingGeometryBuilderCond::RobustTrackingGeometryBuilderCond(
    const std::string& t, const std::string& n, const IInterface* p)
    : base_class(t, n, p) {}

// Athena standard methods
// initialize
StatusCode InDet::RobustTrackingGeometryBuilderCond::initialize() {
  // Retrieve the beampipe builders
  // --------------------------------------------------------
  ATH_CHECK(m_beamPipeBuilder.retrieve());
  ATH_MSG_DEBUG("Retrieved tool " << m_beamPipeBuilder);

  // Retrieve the layer builders
  // -----------------------------------------------------------
  ATH_CHECK(m_layerBuilders.retrieve());
  ATH_MSG_DEBUG("Retrieved tool " << m_layerBuilders);

  return InDet::RobustTrackingGeometryBuilderImpl::initialize();
}

std::unique_ptr<Trk::TrackingGeometry>
InDet::RobustTrackingGeometryBuilderCond::trackingGeometry(
    const EventContext& ctx, Trk::TrackingVolume*,
    SG::WriteCondHandle<Trk::TrackingGeometry>& whandle) const {

  return InDet::RobustTrackingGeometryBuilderImpl::trackingGeometryImpl<
      InDet::RobustTrackingGeometryBuilderImpl::Cond>(
      m_layerBuilders, m_beamPipeBuilder, &ctx, &whandle);
}
