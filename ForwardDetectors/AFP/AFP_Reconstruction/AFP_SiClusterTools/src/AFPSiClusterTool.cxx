/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

/**
 * @file   AFPSiClusterTool.cxx
 * @author Grzegorz Gach <gach@agh.edu.pl>
 * @date   2017-04-22
 * 
 * @brief  Source file for AFPSiClusterTool used for clustering hits.
 * 
 */


#include "AFP_SiClusterTools/AFPSiClusterTool.h"


AFPSiClusterTool::AFPSiClusterTool (const std::string& type,
                                    const std::string& name,
                                    const IInterface* parent)
  : base_class(type, name, parent)
{
}

std::vector< std::vector<AFPSiClusterLayerBasicObj> > AFPSiClusterTool::initLayers (const std::vector<int>& layersInStations) const
{
  // set number of stations
  std::vector< std::vector<AFPSiClusterLayerBasicObj> > my_layers;
  my_layers.resize(layersInStations.size());

  // initialise number of layers in each station
  int stationID = 0;
  std::vector<std::vector<AFPSiClusterLayerBasicObj> >::iterator layersIter = my_layers.begin();
  for (const int layersN : layersInStations) {
    int layerID = 0;
    std::vector<AFPSiClusterLayerBasicObj>& layersVec = (*layersIter++);
    layersVec.resize(layersN);
    // set station and layer ID numbers
    for (AFPSiClusterLayerBasicObj& theLayer : layersVec) {
      theLayer.setStationID(stationID);
      theLayer.setLayerID(layerID++);
    }
    stationID++;
  } // close for over layersN

  return my_layers;
}

StatusCode AFPSiClusterTool::initialize()
{
  ATH_MSG_DEBUG("Initializing " << name() << "...");

  CHECK( m_AFPSiHitContainerKey.initialize( SG::AllowEmpty ) );

  // retrieve tools
  CHECK( m_clusterAlgToolHandle.retrieve() );
  CHECK( m_rowColToLocalCSToolHandle.retrieve() );
  CHECK( m_siLocAlignDBTool.retrieve() );
  CHECK( m_siGlobAlignDBTool.retrieve() );
  
  // monitoring
  if (!(m_monTool.name().empty())) {
    CHECK( m_monTool.retrieve() );
    ATH_MSG_DEBUG("m_monTool name: " << m_monTool);
  }

  return StatusCode::SUCCESS;
}

StatusCode AFPSiClusterTool::finalize() 
{
  return StatusCode::SUCCESS;
}

StatusCode AFPSiClusterTool::clusterHits(std::unique_ptr<xAOD::AFPSiHitsClusterContainer>& clustr, const EventContext &ctx) const
{
  std::vector< std::vector<AFPSiClusterLayerBasicObj> > my_layers=initLayers(m_numberOfLayersInStations);

  if(clearAllLayers(my_layers).isFailure())
  {
    ATH_MSG_WARNING("clearAllLayers failed");
    return StatusCode::FAILURE;
  }
  
  if(fillLayersWithHits(my_layers, ctx).isFailure())
  {
    ATH_MSG_WARNING("fillLayersWithHits failed");
    return StatusCode::FAILURE;
  }
  
  if(clusterEachLayer(my_layers).isFailure())
  {
    ATH_MSG_WARNING("clusterEachLayer failed");
    return StatusCode::FAILURE;
  }
  
  if(saveToXAOD(clustr, my_layers, ctx).isFailure())
  {
    ATH_MSG_WARNING("saveToXAOD failed");
    return StatusCode::FAILURE;
  }
  
  auto clusterSize = Monitored::Scalar("ClusterSize", clustr->size());
  Monitored::Group( m_monTool, clusterSize);

  return StatusCode::SUCCESS;
}

StatusCode AFPSiClusterTool::clearAllLayers(std::vector< std::vector<AFPSiClusterLayerBasicObj> > &my_layers) const
{
  for (std::vector<AFPSiClusterLayerBasicObj>& station : my_layers)
    for (AFPSiClusterLayerBasicObj& layer : station)
      layer.clearHitsAndClusters();
  
  return StatusCode::SUCCESS;
}

StatusCode AFPSiClusterTool::fillLayersWithHits(std::vector< std::vector<AFPSiClusterLayerBasicObj> > &my_layers, const EventContext &ctx) const
{
  // retrieve hits
  if (m_AFPSiHitContainerKey.empty()) {
    ATH_MSG_DEBUG("AFPSiClusterTool, no input siHitContainer");
    // this is allowed, there might be no AFP data in the input
    return StatusCode::SUCCESS;
  }

  SG::ReadHandle<xAOD::AFPSiHitContainer> siHitContainer( m_AFPSiHitContainerKey, ctx );
  if(!siHitContainer.isValid()) {
     ATH_MSG_WARNING("AFPSiClusterTool failed to retrieve siHitContainer, "<<m_AFPSiHitContainerKey<<", exiting gracefully");
     // unexpected absence of Si hits ?
     return StatusCode::SUCCESS;
  }
  else
  {
    ATH_MSG_DEBUG("AFPSiClusterTool::fillLayersWithHits(), successfully got siHitContainer, there are "<<siHitContainer->size()<<" hits");
    
    auto hitsSize = Monitored::Scalar("HitsSize", siHitContainer->size());
    Monitored::Group( m_monTool, hitsSize);

    try {
      // fill layers with hits
      for(const auto theHit : *siHitContainer)
        my_layers.at(theHit->stationID()).at(theHit->pixelLayerID()).hits().push_back(theHit);
    }
    catch (const std::out_of_range& outOfRange) {
      ATH_MSG_WARNING("Hit with station or pixel ID outside expected range. Aborting pixels clustering.");
      clearAllLayers(my_layers).ignore();
      return StatusCode::FAILURE;
    }
  }

  return StatusCode::SUCCESS;
}

StatusCode AFPSiClusterTool::clusterEachLayer(std::vector< std::vector<AFPSiClusterLayerBasicObj> > &my_layers) const
{
  for (std::vector<AFPSiClusterLayerBasicObj>& station : my_layers)
    for (AFPSiClusterLayerBasicObj& layer : station)
      m_clusterAlgToolHandle->doClustering(layer.hits(), layer.clusters()).ignore();
  
  return StatusCode::SUCCESS;
}

StatusCode AFPSiClusterTool::saveToXAOD(std::unique_ptr<xAOD::AFPSiHitsClusterContainer>& clusterContainer, std::vector< std::vector<AFPSiClusterLayerBasicObj> > &my_layers, const EventContext &ctx) const
{ 

  if (m_AFPSiHitContainerKey.empty()) {
   // this is allowed, there might be no AFP data in the input
    return StatusCode::SUCCESS;
  }

  SG::ReadHandle<xAOD::AFPSiHitContainer> siHitContainer( m_AFPSiHitContainerKey, ctx );
  if(!siHitContainer.isValid())
  {
      // this is allowed, there might be no AFP data in the input
      return StatusCode::SUCCESS;
  }
  
  // fill xAOD container
  
  nlohmann::json dataLA, dataGA;
  bool dataLA_init{false}, dataGA_init{false};
  
  for (std::vector<AFPSiClusterLayerBasicObj>& station : my_layers)
    for (AFPSiClusterLayerBasicObj& layer : station) {
      if(layer.clusters().empty()) continue;
      
      const int stationID = layer.stationID();
      const int layerID = layer.layerID();
      
      if(!dataLA_init)
      {
      	// read from DB only if necessary
      	dataLA=m_siLocAlignDBTool->alignmentData(ctx);
      	dataLA_init=true;
      }
      if(!dataGA_init)
      {
      	// read from DB only if necessary
      	dataGA=m_siGlobAlignDBTool->alignmentData(ctx);
      	dataGA_init=true;
      }
      
      const AFP::SiLocAlignData LA=m_siLocAlignDBTool->alignment(dataLA, stationID, layerID);
      const AFP::SiGlobAlignData GA=m_siGlobAlignDBTool->alignment(dataGA, stationID);
      
      for (const AFPSiClusterBasicObj& theCluster : layer.clusters()) {
        // create xAOD object and set cluster coordinates and errors
        xAOD::AFPSiHitsCluster* xAODCluster = m_rowColToLocalCSToolHandle->newXAODLocal(stationID, layerID, LA, GA, theCluster, clusterContainer);
        
        ATH_MSG_DEBUG("have xAODCluster: pixelLayerID "<<xAODCluster->pixelLayerID()<<", xLocal "<<xAODCluster->xLocal()<<", yLocal "<<xAODCluster->yLocal()<<", zLocal "<<xAODCluster->zLocal()<<", nHits "<<xAODCluster->nHits());
        
        // add links to hits
        for (const xAOD::AFPSiHit* theHit : theCluster.hits()) {
          ElementLink< xAOD::AFPSiHitContainer >* hitLink = new ElementLink< xAOD::AFPSiHitContainer >; // will be taken over by the xAODCluster and deleted
          hitLink->toContainedElement(*siHitContainer, theHit);
          xAODCluster->addHitLink(*hitLink);
        }
      } // end for over layers
    } // end for over stations

  return StatusCode::SUCCESS;
}

