/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#include "tauMonitoring/tauMonitorAlgorithm.h"



tauMonitorAlgorithm::tauMonitorAlgorithm( const std::string& name, ISvcLocator* pSvcLocator )
:AthMonitorAlgorithm(name,pSvcLocator)
,m_doRandom(true)
{
  declareProperty("TauRecContainer", m_TauContainerKey="TauJets");
  declareProperty("etaMin", m_etaMin=-1.);
  declareProperty("etaMax", m_etaMax=2.6);
  declareProperty("kinGroupName", m_kinGroupName="tauMonKinGroupBA");

}


tauMonitorAlgorithm::~tauMonitorAlgorithm() {}


StatusCode tauMonitorAlgorithm::initialize() {
    using namespace Monitored;

    ATH_CHECK( m_TauContainerKey.initialize() );

    return AthMonitorAlgorithm::initialize();
}


StatusCode tauMonitorAlgorithm::fillHistograms( const EventContext& ctx ) const {
    using namespace Monitored;


    auto eta = Monitored::Scalar<float>("eta");
    auto phi = Monitored::Scalar<float>("phi");
    auto pt = Monitored::Scalar<float>("pt");
    auto charge = Monitored::Scalar<int>("charge");
    auto nTracks = Monitored::Scalar<int>("nTracks");
    auto ntaus = Monitored::Scalar<int>("ntaus");
    auto nClusters = Monitored::Scalar<int>("nClusters");
    auto lb = Monitored::Scalar<int>("lb",0);



    SG::ReadHandle<xAOD::TauJetContainer> taus(m_TauContainerKey, ctx);
    if (! taus.isValid() ) {
      ATH_MSG_ERROR("evtStore() does not contain tau Collection with name "<< m_TauContainerKey);
      return StatusCode::FAILURE;
    }

    ntaus   = taus->size();

    for (const auto& tau : *taus) {
      // do stuff with taus
      eta = tau->eta();
      phi = tau->phi();
      pt = tau->pt();
      charge = tau->charge();
      nTracks = tau->nTracks();
      nClusters = tau->detail<int>(xAOD::TauJetParameters::numTopoClusters) ;
      lb = GetEventInfo(ctx)->lumiBlock();
      ANA_MSG_INFO( "groupName = " << m_kinGroupName );


      if (m_etaMin < eta && eta < m_etaMax){
      	fill(m_kinGroupName, eta, phi, pt, charge, nTracks, ntaus, nClusters,lb);
      }



    }
    

    return StatusCode::SUCCESS;
}
