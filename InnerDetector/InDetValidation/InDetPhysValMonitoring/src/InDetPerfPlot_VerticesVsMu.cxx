/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

/**
 * @file InDetPerfPlot_VerticesVsMu.cxx
 * @author Simone Pagan Griso
 **/

#include "InDetPerfPlot_VerticesVsMu.h"
#include "xAODTracking/TrackParticle.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTracking/Vertex.h"
#include "xAODTracking/VertexContainer.h"
#include "EventPrimitives/EventPrimitives.h"
#include "EventPrimitives/EventPrimitivesHelpers.h"

InDetPerfPlot_VerticesVsMu::InDetPerfPlot_VerticesVsMu(InDetPlotBase* pParent,
                                                             const std::string& sDir) : InDetPlotBase(pParent, sDir),
  m_vx_n_vs_mu{} {
  // nop
}

void
InDetPerfPlot_VerticesVsMu::initializePlots() {

  book(m_vx_n_vs_mu, "vx_n_vs_mu");

}

void
InDetPerfPlot_VerticesVsMu::fill(const xAOD::VertexContainer& vertices, unsigned int nPU, float weight) {
  // fill number of vertices
  int nGoodVertices(0);

  for (const auto& vtx : vertices.stdcont()) {
    if (vtx->vertexType() == xAOD::VxType::NoVtx) {
      continue; // skip dummy vertex
    }
    nGoodVertices++;
  }
  fillHisto(m_vx_n_vs_mu, nGoodVertices, nPU,weight);
}
