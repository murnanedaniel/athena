/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
// LayerProvider.cxx, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

// STL
#include <sstream>
// Trk include
#include "TrkDetDescrTools/LayerProvider.h"
#include "TrkDetDescrInterfaces/ILayerBuilder.h"
#include "TrkGeometry/Layer.h"
#include "TrkGeometry/CylinderLayer.h"
#include "TrkGeometry/DiscLayer.h"

// constructor
Trk::LayerProvider::LayerProvider(const std::string& t, const std::string& n, const IInterface* p)
: AthAlgTool(t,n,p)
{
    declareInterface<Trk::ILayerProvider>(this);
    
    // Name specification from outside
    declareProperty("LayerBuilder", m_layerBuilder);
    
}

// destructor
Trk::LayerProvider::~LayerProvider()
= default;

// initialize
StatusCode Trk::LayerProvider::initialize()
{
    if (m_layerBuilder.retrieve().isFailure()){
        ATH_MSG_WARNING("Could not retrieve layer builder");
    }
    return StatusCode::SUCCESS;
}

/** LayerBuilder interface method - returning the layers at negative side */
const std::vector<Trk::Layer*> Trk::LayerProvider::negativeLayers() const
{
    // if any cache is there it has to be returned
    if (!m_layerCache.empty()) return m_layerCache;
    // this will fill the cache with positive layers
    return discLayers(-1);
}    

/** LayerBuilder interface method - returning the central layers */
const std::vector<Trk::Layer*> Trk::LayerProvider::centralLayers() const
{
    // central layers
    std::vector<Trk::Layer* > cLayers;
    // retrieving the cylinder layers from the layer builder
    std::unique_ptr<const std::vector<Trk::CylinderLayer*> >
      cylinderLayers = m_layerBuilder->cylindricalLayers();
    // loop over it and push into the return vector;
    if (cylinderLayers){
        for (const auto & cL : (*cylinderLayers))
            cLayers.push_back(cL);
    }
    // and return
    return cLayers;
} 

/** LayerBuilder interface method - returning the layers at negative side */
const std::vector<Trk::Layer*> Trk::LayerProvider::positiveLayers() const
{
    // if any cache is there it has to be returned
    if (!m_layerCache.empty()) return m_layerCache;
    // this will fill the cache with negative layers
    return discLayers(1);
}

/** LayerBuilder interface method - returning the layers at negative side */
const std::vector<Trk::Layer*> Trk::LayerProvider::discLayers(int posneg) const
{
    // get the disc layers
    std::vector <Trk::Layer*>   dLayers;
    // retrieving the cylinder layers from the layer builder
    std::unique_ptr<const std::vector< Trk::DiscLayer*> > discLayers =
      m_layerBuilder->discLayers();
    // loop and fill either cache or dLayers
    if (discLayers){
        // loop over and push into the return/cache vector 
        for (const auto & dL : (*discLayers) ){
            // get the center posituion 
            double zpos = dL->surfaceRepresentation().center().z();
            if (posneg > 0){
                // configured to provide positive and cache negative
                if (zpos > 0.) dLayers.push_back(dL);
                else m_layerCache.push_back(dL);
            } else {
                // configured to provide negative and cache positive
                if (zpos < 0.) dLayers.push_back(dL);
                else m_layerCache.push_back(dL);
            }
        }
    }
    // and return
    return dLayers;
}


// finalize
StatusCode Trk::LayerProvider::finalize()
{
    return StatusCode::SUCCESS;
}

const std::string& Trk::LayerProvider::identification() const
{
    return m_layerBuilder->identification();
}
