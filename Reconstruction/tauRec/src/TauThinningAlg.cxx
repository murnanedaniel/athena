/*
  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration.
*/

#include "TauThinningAlg.h"
#include "StoreGate/ReadHandle.h"
#include "StoreGate/ThinningHandle.h"
#include "tauRecTools/HelperFunctions.h"

/**
 * @brief Gaudi initialize method.
 */
StatusCode TauThinningAlg::initialize()
{
  ATH_CHECK( m_taus.initialize(m_streamName) );
  ATH_CHECK( m_tauTracks.initialize(m_streamName) );
  ATH_CHECK( m_neutralPFOs.initialize(m_streamName) );
  ATH_CHECK( m_pi0clusters.initialize(m_streamName) );
  ATH_CHECK( m_pi0CellLinks.initialize(m_streamName) );
  ATH_CHECK( m_finalPi0s.initialize(m_streamName) );
  ATH_CHECK( m_shotPFOs.initialize(m_streamName) );
  ATH_CHECK( m_shotclusters.initialize(m_streamName) );
  ATH_CHECK( m_shotCellLinks.initialize(m_streamName) );
  ATH_CHECK( m_hadronicPFOs.initialize(m_streamName) );
  ATH_CHECK( m_secondaryVertices.initialize(m_streamName) );
  ATH_CHECK( m_cells.initialize(m_streamName) );
  ATH_CHECK( m_tauCellLinks.initialize(m_streamName) );

  return StatusCode::SUCCESS;
}

/**
 * @brief Execute the algorithm.
 * @param ctx Current event context.
 */
StatusCode TauThinningAlg::execute (const EventContext& ctx) const
{
  SG::ThinningHandle<xAOD::TauJetContainer> taus (m_taus, ctx);
  taus.thinAll();

  SG::ThinningHandle<xAOD::TauTrackContainer> tauTracks (m_tauTracks, ctx);
  tauTracks.thinAll();

  SG::ThinningHandle<xAOD::PFOContainer> neutralPFOs (m_neutralPFOs, ctx);
  neutralPFOs.thinAll();

  SG::ThinningHandle<xAOD::CaloClusterContainer> pi0clusters (m_pi0clusters, ctx);
  pi0clusters.thinAll();

  SG::ThinningHandle<CaloClusterCellLinkContainer> pi0CellLinks (m_pi0CellLinks, ctx);
  pi0CellLinks.thinAll();

  SG::ThinningHandle<xAOD::PFOContainer> shotPFOs (m_shotPFOs, ctx);
  shotPFOs.thinAll();

  SG::ThinningHandle<xAOD::PFOContainer> hadronicPFOs (m_hadronicPFOs, ctx);
  hadronicPFOs.thinAll();

  SG::ThinningHandle<xAOD::VertexContainer> secondaryVertices (m_secondaryVertices, ctx);
  secondaryVertices.thinAll();

  SG::ThinningHandle<CaloCellContainer> cells (m_cells, ctx);
  cells.thinAll();

  SG::ThinningHandle<CaloClusterCellLinkContainer> tauCellLinks (m_tauCellLinks, ctx);
  tauCellLinks.thinAll();

  // The three following containers didn't exist in r21.
  // Make them optional, so we won't crash processing a r21 ESD.

  std::optional<SG::ThinningHandle<xAOD::ParticleContainer> > finalPi0sOpt;
  if (evtStore()->contains<xAOD::ParticleContainer> (m_finalPi0s.key())) {
    finalPi0sOpt.emplace (m_finalPi0s, ctx);
    finalPi0sOpt->thinAll();
  }

  std::optional<SG::ThinningHandle<xAOD::CaloClusterContainer> > shotclustersOpt;
  if (evtStore()->contains<xAOD::CaloClusterContainer> (m_shotclusters.key())) {
    shotclustersOpt.emplace (m_shotclusters, ctx);
    shotclustersOpt->thinAll();
  }

  std::optional<SG::ThinningHandle<CaloClusterCellLinkContainer> > shotCellLinksOpt;
  if (evtStore()->contains<CaloClusterCellLinkContainer> (m_shotCellLinks.key())) {
    shotCellLinksOpt.emplace (m_shotCellLinks, ctx);
    shotCellLinksOpt->thinAll();
  }

  static const SG::AuxElement::ConstAccessor<char> acc_passThinning("passThinning");

  for (const xAOD::TauJet* tau : *taus) {

    if (!acc_passThinning(*tau)) continue;

    // keep tau
    taus.keep(tau->index());
    
    // keep tau tracks
    for (const xAOD::TauTrack* track : tau->allTracks()) {
      tauTracks.keep(track->index());
    }

    // keep tau cluster cell links and cells within 0.2 of the tau axis
    TLorentzVector tauAxis = tauRecTools::getTauAxis(*tau, m_doVertexCorrection);
    const xAOD::Vertex* tauVertex = tau->vertex();

    auto clusterList = tau->clusters();
    for (const xAOD::IParticle* particle : clusterList) {
      const xAOD::CaloCluster* cluster = static_cast<const xAOD::CaloCluster*>(particle);
      TLorentzVector clusterP4 = cluster->p4();

      // correct the four momentum to point to the tau vertex
      if (tauVertex) {
        xAOD::CaloVertexedTopoCluster vertexedCluster(*cluster, tauVertex->position()); 
        clusterP4 = vertexedCluster.p4();
      }

      if (clusterP4.DeltaR(tauAxis) > 0.2) continue;

      const CaloClusterCellLink* cellLinks = cluster->getCellLinks();
      if (!cellLinks) {
	ATH_MSG_WARNING( "Skipping cluster without cell links." );
	continue;
      }
      
      // cluster cell links
      CaloClusterCellLinkContainer::const_iterator cellLinks_it = std::find(tauCellLinks->begin(), tauCellLinks->end(), cellLinks);
      if(cellLinks_it != tauCellLinks->end()) {
	size_t link_index = std::distance(tauCellLinks->begin(), cellLinks_it);
	tauCellLinks.keep(link_index);
      }
      else {
	ATH_MSG_WARNING( "Could not find cluster cell link in " << m_tauCellLinks.key() << ", skipping cluster." );
	continue;
     }

      // cells
      CaloClusterCellLink::const_iterator it = cellLinks->begin();
      CaloClusterCellLink::const_iterator end = cellLinks->end();
      for (; it != end; ++it) {
	if (it.index() >= cells->size()) {
	  ATH_MSG_WARNING( "Cell index " << it.index() << " is larger than the number of cells in " << m_cells.key() << " (" << cells->size() << ")" );
	  continue;
	}
	cells.keep (it.index());
      }
    }

    // keep neutral PFOs, pi0 clusters, cell links and cells
    for(size_t i=0; i<tau->nNeutralPFOs(); i++) {
      // neutral PFOs
      neutralPFOs.keep(tau->neutralPFO(i)->index());

      // pi0 clusters
      const xAOD::CaloCluster* cluster = tau->neutralPFO(i)->cluster(0);
      pi0clusters.keep(cluster->index());
  
      // pi0 cell links
      const CaloClusterCellLink* cellLinks = cluster->getCellLinks();
      CaloClusterCellLinkContainer::const_iterator cellLinks_it = std::find(pi0CellLinks->begin(), pi0CellLinks->end(), cellLinks);
      if(cellLinks_it != pi0CellLinks->end()) {
	size_t link_index = std::distance(pi0CellLinks->begin(), cellLinks_it);
	pi0CellLinks.keep(link_index);
      }
      else {
	ATH_MSG_WARNING( "Could not find cluster cell link in " << m_pi0CellLinks.key() << ", won't be saved in xAOD." );
      }

      // pi0 cells
      CaloClusterCellLink::const_iterator it = cellLinks->begin();
      CaloClusterCellLink::const_iterator end = cellLinks->end();
      for (; it != end; ++it) {
	if (it.index() >= cells->size()) {
	  ATH_MSG_WARNING( "Cell index " << it.index() << " is larger than the number of cells in " << m_cells.key() << " (" << cells->size() << ")" );
	  continue;
	}
	cells.keep (it.index());
      }
    }

    // keep final pi0s
    if (finalPi0sOpt) {
      for(size_t i=0; i<tau->nPi0s(); i++) {
        finalPi0sOpt->keep(tau->pi0(i)->index());
      }
    }

    // keep shot PFOs, clusters, cell links and cells
    for(size_t i=0; i<tau->nShotPFOs(); i++) {
      // shot PFOs
      shotPFOs.keep(tau->shotPFO(i)->index());

      // shot clusters
      const xAOD::CaloCluster* cluster = tau->shotPFO(i)->cluster(0);
      if (!cluster) continue;
      if (shotclustersOpt) {
        shotclustersOpt->keep(cluster->index());
      }

      // shot cell links
      const CaloClusterCellLink* cellLinks = cluster->getCellLinks();
      if (shotCellLinksOpt) {
        CaloClusterCellLinkContainer::const_iterator cellLinks_it = std::find((*shotCellLinksOpt)->begin(), (*shotCellLinksOpt)->end(), cellLinks);
        if(cellLinks_it != (*shotCellLinksOpt)->end()) {
          size_t link_index = std::distance((*shotCellLinksOpt)->begin(), cellLinks_it);
          shotCellLinksOpt->keep(link_index);
        }
        else {
          ATH_MSG_WARNING( "Could not find cluster cell link in " << m_shotCellLinks.key() << ", won't be saved in xAOD." );
        }
      }

      // shot cells
      CaloClusterCellLink::const_iterator it = cellLinks->begin();
      CaloClusterCellLink::const_iterator end = cellLinks->end();
      for (; it != end; ++it) {
	if (it.index() >= cells->size()) {
	  ATH_MSG_WARNING( "Cell index " << it.index() << " is larger than the number of cells in " << m_cells.key() << " (" << cells->size() << ")" );
	  continue;
	}
	cells.keep (it.index());
      }
    }

    // keep hadronic PFOs
    for(size_t i=0; i<tau->nHadronicPFOs(); i++) {
      hadronicPFOs.keep(tau->hadronicPFO(i)->index());
    }
    
    // keep secondary vertex when present
    if (tau->secondaryVertex() != nullptr) {
      secondaryVertices.keep(tau->secondaryVertex()->index());
    }        
  }

  return StatusCode::SUCCESS;
}
