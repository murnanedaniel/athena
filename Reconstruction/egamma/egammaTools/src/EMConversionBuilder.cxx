/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

/********************************************************************

NAME:     EMConversionBuilder.cxx
PACKAGE:  offline/Reconstruction/egamma/egammaRec

AUTHORS:  D. Zerwas, B. Lenzi
CREATED:  Jul, 2005
CHANGES:  Mar, 2014 (BL) xAOD migration

PURPOSE:  subAlgorithm which creates an EMConversion object. 

********************************************************************/

// INCLUDE HEADER FILES:

#include "EMConversionBuilder.h"
#include "ConvVxSorter.h"
#include "egammaInterfaces/IEMExtrapolationTools.h"
#include "xAODTracking/VertexContainer.h"
#include "egammaRecEvent/egammaRecContainer.h"
#include "egammaRecEvent/egammaRec.h"

//  END OF HEADER FILES INCLUDE

/////////////////////////////////////////////////////////////////

using CLHEP::GeV;

EMConversionBuilder::EMConversionBuilder(const std::string& type,
					 const std::string& name,
					 const IInterface* parent)
  : egammaBaseTool(type, name, parent),
    m_extrapolationTool("EMExtrapolationTools")
{
  
  // declare interface
  declareInterface<IEMConversionBuilder>(this);
  
  // Name of the input conversion container
  declareProperty("ConversionContainerName",                 
		  m_conversionContainerName="PhotonConversionVertices",
		  "Name of the input conversion container");
	
  // Name of the input egammaRec container
  declareProperty("egammaRecContainerName",                 
		  m_egammaRecContainerName="egammaRecs",
		  "Name of the input egammaRec container");

  // Name of the extrapolation tool
  declareProperty("ExtrapolationTool",
		  m_extrapolationTool,
		  "Handle of the extrapolation tool");
      
  declareProperty("RejectAllTRTConversions", m_rejectAllTRT = false,
		  "Ignore all conversion vertices containing exclusively TRT-only tracks");
    
  declareProperty("MinTRTHits", m_minTRTHits = 0,
		  "minimum number of TRT hits for TRT-only tracks (both single and double track conversion vertices)");
  
  declareProperty("MinSiSingleTrackPt", m_minSiSingleTrackPt = 0*GeV,
		  "minimum pT for Si tracks to be considered single-track conversion vertices");
  
  declareProperty("MinTRTOnlySingleTrackPt", m_minTRTonlySingleTrackPt = 2*GeV,
		  "minimum pT for TRT-only tracks to be considered single-track conversion vertices");
  
  declareProperty("MinTRTOnlyTrackPt", m_minTRTonlyTrackPt = 0*GeV,
		  "minimum pT for each track in TRT-only double-track conversion vertices");
        
  declareProperty("MinSumPt", m_minSumPt = 0*GeV,
		  "minimum sum pT for double track conversion vertices");
}

// =================================================================
// DESTRUCTOR:
EMConversionBuilder::~EMConversionBuilder()
{  
}

// ==================================================================
// INITIALIZE METHOD:  
StatusCode EMConversionBuilder::initialize()
{

  ATH_MSG_DEBUG("Initializing EMConversionBuilder");

  // the extrapolation tool
  if(m_extrapolationTool.retrieve().isFailure()){
    ATH_MSG_ERROR("Cannot retrieve extrapolationTool " << m_extrapolationTool);
    return StatusCode::FAILURE;
  } else {
    ATH_MSG_DEBUG("Retrieved extrapolationTool " << m_extrapolationTool);
  }

  return StatusCode::SUCCESS;
}

// =============================================================
StatusCode EMConversionBuilder::contExecute() 
{

  // retrieve Conversion Container
  const xAOD::VertexContainer* conversions = 0;
  if (evtStore()->contains<xAOD::VertexContainer>(m_conversionContainerName)) {
    StatusCode sc = evtStore()->retrieve(conversions,m_conversionContainerName);
    if(sc.isFailure()) {
      ATH_MSG_WARNING("Could not retrieve Conversion container! EMConversionBuilder will stop.");
      return StatusCode::SUCCESS;
    }
  } else {
    ATH_MSG_DEBUG("Could not find conversion container! Acting as if there are 0 conversions."); 
  }
  
  // retrieve egammaRec container
  const EgammaRecContainer* egammaRecs = 0;
  if (evtStore()->contains<EgammaRecContainer>(m_egammaRecContainerName)) {
    StatusCode sc = evtStore()->retrieve(egammaRecs,m_egammaRecContainerName);
    if(sc.isFailure()) {
      ATH_MSG_WARNING("Could not retrieve egammaRec container! EMConversionBuilder will stop.");
      return StatusCode::SUCCESS;
    }
  } else {
    ATH_MSG_DEBUG("Could not find conversion container! Acting as if there are 0 egammaRecs."); 
  }
  
  
  // Extrapolate each vertex (keeping etatCalo, phiAtCalo)
  // and try to match to each cluster
  
  // xAOD::Vertex does not have method isAvailable for the moment
  // create a lambda function
  auto isAvailable = [](const xAOD::Vertex& vertex, std::string name) { 
    SG::AuxElement::Accessor<float> acc(name, "");
    return acc.isAvailable(vertex); } ; 

  float etaAtCalo, phiAtCalo;
  for (unsigned int iVtx = 0; iVtx < conversions->size(); ++iVtx)
  {
    const xAOD::Vertex *vertex = conversions->at(iVtx);
    
    // Check if vertex was already decorated with etaAtCalo, phiAtCalo  
    if (isAvailable(*vertex, "etaAtCalo") &&
        isAvailable(*vertex, "phiAtCalo") )
    {
      etaAtCalo = vertex->auxdata<float>("etaAtCalo");
      phiAtCalo = vertex->auxdata<float>("phiAtCalo");
    }
    // check extrapolation, skip vertex in case of failure
    else if (!m_extrapolationTool->getEtaPhiAtCalo(vertex, &etaAtCalo, &phiAtCalo))
      continue;
    
    for (auto& egRec : *egammaRecs)
    {
      const xAOD::CaloCluster *cluster = egRec->caloCluster();
      if (!m_extrapolationTool->matchesAtCalo(cluster, vertex, etaAtCalo, phiAtCalo))
        continue;
      const ElementLink< xAOD::VertexContainer > vertexLink( *conversions, iVtx );
      
      // If this is the best (or the first) vertex, push front and keep deltaEta, deltaPhi
      if (!egRec->getNumberOfVertices() || ConvVxSorter()(*vertex, *egRec->vertex()))
      {
        egRec->pushFrontVertex( vertexLink );
        egRec->setDeltaEtaVtx( cluster->etaBE(2) - etaAtCalo );
        egRec->setDeltaPhiVtx( m_phiHelper.diff(cluster->phiBE(2), phiAtCalo) );
      }
      else // Not the best vertex, push back
        egRec->pushBackVertex( vertexLink );       
    }
    
  }
  
  return StatusCode::SUCCESS;
}

// =============================================================
StatusCode EMConversionBuilder::hltExecute(egammaRec* egRec, const xAOD::VertexContainer* conversions)
{

  if (!egRec || !conversions)
  {
    ATH_MSG_WARNING("trackExecute: NULL pointer to egammaRec or VertexContainer");
    return StatusCode::SUCCESS;
  }
  
  
  // Extrapolate each vertex (keeping etatCalo, phiAtCalo)
  // and try to match to each cluster
  
  // xAOD::Vertex does not have method isAvailable for the moment
  // create a lambda function
  auto isAvailable = [](const xAOD::Vertex& vertex, std::string name) { 
    SG::AuxElement::Accessor<float> acc(name, "");
    return acc.isAvailable(vertex); } ; 

  float etaAtCalo, phiAtCalo;
  for (unsigned int iVtx = 0; iVtx < conversions->size(); ++iVtx)
  {
    const xAOD::Vertex *vertex = conversions->at(iVtx);
    
    // Check if vertex was already decorated with etaAtCalo, phiAtCalo  
    if (isAvailable(*vertex, "etaAtCalo") &&
        isAvailable(*vertex, "phiAtCalo") )
    {
      etaAtCalo = vertex->auxdata<float>("etaAtCalo");
      phiAtCalo = vertex->auxdata<float>("phiAtCalo");
    }
    // check extrapolation, skip vertex in case of failure
    else if (!m_extrapolationTool->getEtaPhiAtCalo(vertex, &etaAtCalo, &phiAtCalo))
      continue;
    
    const xAOD::CaloCluster *cluster = egRec->caloCluster();
    if (!m_extrapolationTool->matchesAtCalo(cluster, vertex, etaAtCalo, phiAtCalo))
      continue;
    const ElementLink< xAOD::VertexContainer > vertexLink( *conversions, iVtx );
    
    // If this is the best (or the first) vertex, push front and keep deltaEta, deltaPhi
    if (!egRec->getNumberOfVertices() || ConvVxSorter()(*vertex, *egRec->vertex()))
      {
        egRec->pushFrontVertex( vertexLink );
        egRec->setDeltaEtaVtx( cluster->etaBE(2) - etaAtCalo );
        egRec->setDeltaPhiVtx( m_phiHelper.diff(cluster->phiBE(2), phiAtCalo) );
      }
    else // Not the best vertex, push back
      egRec->pushBackVertex( vertexLink );       
  }
  
  return StatusCode::SUCCESS;
}


// ==================================================================
// FINALIZE METHOD:  
StatusCode EMConversionBuilder::finalize()
{

  return StatusCode::SUCCESS;

}
