/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#include "PhysValCluster.h"
#include "xAODCaloEvent/CaloClusterContainer.h"

PhysValCluster::PhysValCluster (const std::string& type, const std::string& name, const IInterface* parent ) : ManagedMonitorToolBase( type, name, parent ) {
  declareProperty( "ClusterContainerName", m_clusterContainerName="");
}

PhysValCluster::~PhysValCluster() {}

StatusCode PhysValCluster::initialize(){
  ATH_CHECK(ManagedMonitorToolBase::initialize());
  return StatusCode::SUCCESS;
}


StatusCode PhysValCluster::bookHistograms(){

  std::string theName = "TopoClusters/TopoClusters_JetETMiss/JetETMiss_"+m_clusterContainerName;

  m_clusterValidationPlots.reset(new ClusterValidationPlots(0,theName, theName));
  m_clusterValidationPlots->setDetailLevel(100);
  m_clusterValidationPlots->initialize();
  std::vector<HistData> hists = m_clusterValidationPlots->retrieveBookedHistograms();

  for (auto hist : hists) {
    ATH_CHECK(regHist(hist.first,hist.second,all));
  }
  
  return StatusCode::SUCCESS;     
   
}

StatusCode PhysValCluster::fillHistograms(){

  const xAOD::CaloClusterContainer* theClusterContainer = NULL;

  ATH_CHECK( evtStore()->retrieve(theClusterContainer,m_clusterContainerName) );

  if (!theClusterContainer){
    if (msgLvl(MSG::WARNING))  msg(MSG::WARNING) << " Have NULL pointer to xAOD::PFOContainer " << endreq;
    return StatusCode::SUCCESS;
  }
  
  xAOD::CaloClusterContainer::const_iterator firstCluster = theClusterContainer->begin();
  xAOD::CaloClusterContainer::const_iterator lastCluster = theClusterContainer->end();

  for (; firstCluster != lastCluster; ++firstCluster) {
    const xAOD::CaloCluster* theCluster = *firstCluster;
    m_clusterValidationPlots->fill(*theCluster);
  }

  return StatusCode::SUCCESS;

}

StatusCode PhysValCluster::procHistograms(){
   return StatusCode::SUCCESS;
}
