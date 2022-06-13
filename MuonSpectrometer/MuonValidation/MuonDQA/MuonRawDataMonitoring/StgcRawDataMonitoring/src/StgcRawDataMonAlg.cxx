/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Package : sTGCRawDataMonAlg
// Author: Sebastian Fuenzalida Garrido
// Local supervisor: Edson Carquin Lopez
// Technical supervisor: Gerardo Vasquez
//
// DESCRIPTION:
// Subject: sTGC --> sTGC raw data monitoring
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "MuonReadoutGeometry/MuonDetectorManager.h"
#include "MuonReadoutGeometry/MuonStation.h"
#include "MuonReadoutGeometry/sTgcReadoutElement.h"
#include "MuonDQAUtils/MuonChamberNameConverter.h"
#include "MuonDQAUtils/MuonChambersRange.h"
#include "MuonCalibIdentifier/MuonFixedId.h"

#include "StgcRawDataMonitoring/StgcRawDataMonAlg.h"

#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTracking/TrackParticle.h"
#include "xAODTracking/TrackingPrimitives.h"
#include "MuonPrepRawData/MuonPrepDataContainer.h"
#include "MuonRIO_OnTrack/sTgcClusterOnTrack.h"
#include "AthenaMonitoring/AthenaMonManager.h"
#include "MuonPrepRawData/sTgcPrepData.h"
#include "MuonSegment/MuonSegment.h"


/////////////////////////////////////////////////////////////////////////////
// *********************************************************************
// Public Methods
// ********************************************************************* 
/////////////////////////////////////////////////////////////////////////////

StgcRawDataMonAlg::StgcRawDataMonAlg( const std::string& name, ISvcLocator* pSvcLocator ) : AthMonitorAlgorithm(name,pSvcLocator)	      
{
  //Declare the property 
}


StatusCode StgcRawDataMonAlg::initialize()
{   
  // init message stream
  ATH_MSG_DEBUG("initialize sTGCRawDataMonAlg" );
  ATH_MSG_DEBUG("******************" );
  ATH_MSG_DEBUG("doSTGCESD: " << m_doSTGCESD );
  ATH_MSG_DEBUG("******************" );
  
  ATH_CHECK(AthMonitorAlgorithm::initialize());
  ATH_CHECK(m_DetectorManagerKey.initialize());
  ATH_CHECK(m_idHelperSvc.retrieve());

  ATH_MSG_INFO(" Found the MuonIdHelperSvc " );
  ATH_CHECK(m_muonKey.initialize());
  ATH_CHECK(m_sTGCContainerKey.initialize() );
  
  ATH_MSG_DEBUG(" end of initialize " );
  ATH_MSG_INFO("sTGCRawDataMonAlg initialization DONE " );

  return StatusCode::SUCCESS;
} 

StatusCode StgcRawDataMonAlg::fillHistograms(const EventContext& ctx) const
{  
  const xAOD::TrackParticleContainer *meTPContainer = nullptr;
  ATH_CHECK(evtStore() -> retrieve(meTPContainer,"ExtrapolatedMuonTrackParticles" ));

  std::vector<const Muon::sTgcPrepData*> prep_data_vector;
  SG::ReadHandle<Muon::sTgcPrepDataContainer> sTGC_container(m_sTGCContainerKey,ctx);

  if (m_doSTGCESD) 
    {
      Histograms::sTGCSummaryHistogramStruct summaryPlots[2][2][4];
            
      for(const Muon::sTgcPrepDataCollection* coll : *sTGC_container)
	{
	  for (const Muon::sTgcPrepData* prd : *coll)
	    {
	      prep_data_vector.push_back(prd);
	      if (m_do_sTgc_overview) fillsTGCOverviewHistograms(prep_data_vector, prd);
	      fillsTGCSummaryHistograms(prd, summaryPlots);
	    }
	}
    }
    
  prep_data_vector.clear();
  
  return StatusCode::SUCCESS;
}

void StgcRawDataMonAlg::fillsTGCOverviewHistograms(const std::vector<const Muon::sTgcPrepData*> &prd, const Muon::sTgcPrepData *stg_object) const 
{    
  auto charge_all = Monitored::Collection("charge_all", prd, [] (const Muon::sTgcPrepData *aux) 
					  {
					    auto charge_var = aux -> charge(); 
					    return charge_var;
					  });
  auto numberofstrips_percluster = Monitored::Collection("numberofstrips_percluster", prd, [] (const Muon::sTgcPrepData *aux) {const std::vector<Identifier> &stripIds = aux -> rdoList(); return stripIds.size();});
  
  fill("sTGCMonitor", charge_all, numberofstrips_percluster);
  
  std::vector<short int> strip_times_target = stg_object -> stripTimes();
  std::vector<int> strip_charges_target = stg_object -> stripCharges();
  std::vector<short unsigned int> strip_number_target = stg_object -> stripNumbers();
  std::vector<int> strip_statEta_strip_target = stg_object -> stripCharges();

  auto time_all = Monitored::Collection("time_all", prd, [] (const Muon::sTgcPrepData *aux) {return aux -> time();});      
  auto strip_times = Monitored::Collection("strip_times", strip_times_target);
  auto strip_charges = Monitored::Collection("strip_charges", strip_charges_target);
  auto strip_number = Monitored::Collection("strip_number", strip_number_target);

  fill("sTGCMonitor", time_all, strip_times, strip_charges, strip_number);

  auto x_mon = Monitored::Collection("x_mon", prd, [] (const Muon::sTgcPrepData *aux) {Amg::Vector3D pos = aux -> globalPosition(); return pos.x();});
  auto y_mon = Monitored::Collection("y_mon", prd, [] (const Muon::sTgcPrepData *aux) {Amg::Vector3D pos = aux -> globalPosition(); return pos.y();});
  auto z_mon = Monitored::Collection("z_mon", prd, [] (const Muon::sTgcPrepData *aux) {Amg::Vector3D pos = aux -> globalPosition(); return pos.z();});
  auto R_mon = Monitored::Collection("R_mon", prd, [] (const Muon::sTgcPrepData *aux) {Amg::Vector3D pos = aux -> globalPosition(); return std::hypot(pos.x(), pos.y());});

  fill("sTGCMonitor", x_mon, y_mon, z_mon, R_mon);
}

void StgcRawDataMonAlg::fillsTGCSummaryHistograms(const Muon::sTgcPrepData *stg_object, Histograms::sTGCSummaryHistogramStruct (&vects)[2][2][4]) const
{
  Identifier Id    = stg_object -> identify();
  int stationPhi   = m_idHelperSvc -> stgcIdHelper().stationPhi(Id);
  int stationEta   = m_idHelperSvc -> stgcIdHelper().stationEta(Id);
  int multiplet    = m_idHelperSvc -> stgcIdHelper().multilayer(Id);
  int gas_gap      = m_idHelperSvc -> stgcIdHelper().gasGap(Id);  
  int channel_type = m_idHelperSvc -> stgcIdHelper().channelType(Id);
  int charge       = stg_object -> charge();
 
  int stationName  = m_idHelperSvc -> stgcIdHelper().stationName(Id);

  int isector      = (stationName == 58) ? 1 : 0;
  int iside        = (stationEta > 0) ? 1 : 0;

  int stationPhiComplete = get_sectorPhi_from_stationPhi_stName(stationPhi, GeometricSectors::sTGC_Sector[isector]);
 
  auto &Vectors = vects[iside][multiplet - 1][gas_gap - 1];
  Vectors.strip_charges_vec = stg_object -> stripCharges();  
  Vectors.stationEta_perPhi_vec.push_back(stationEta);
  Vectors.strip_numbers_perPhi_vec = stg_object -> stripNumbers();
  Vectors.charge_perPhi_vec.push_back(charge);
  Vectors.charge_vec.push_back(charge);
  Vectors.stationPhi_vec.push_back(stationPhiComplete);
  Vectors.stationEta_vec.push_back(stationEta);

  if (channel_type == 0)
    {
      auto charge_perLayer_pad_ = Monitored::Collection("charge_pad_" + GeometricSectors::sTGC_Side[iside] + "_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), Vectors.charge_vec);
      auto stationPhi_ = Monitored::Collection("sector_" + GeometricSectors::sTGC_Side[iside] + "_phi_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), Vectors.stationPhi_vec);
      auto stationEta_ = Monitored::Collection("sector_" + GeometricSectors::sTGC_Side[iside] + "_eta_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), Vectors.stationEta_vec);
      fill("sTGC_sideGroup" + GeometricSectors::sTGC_Side[iside], charge_perLayer_pad_, stationPhi_, stationEta_);    
    }
  
  else if (channel_type == 1)
    {      
      auto charge_perLayer_strip_ = Monitored::Collection("charge_strip_" + GeometricSectors::sTGC_Side[iside] + "_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), Vectors.charge_vec);
      auto stationPhi_ = Monitored::Collection("sector_" + GeometricSectors::sTGC_Side[iside] + "_phi_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), Vectors.stationPhi_vec);
      auto stationEta_ = Monitored::Collection("sector_" + GeometricSectors::sTGC_Side[iside] + "_eta_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), Vectors.stationEta_vec); 
      auto stationEta_perPhi_ = Monitored::Collection("sector_" + GeometricSectors::sTGC_Side[iside] + "_eta_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap) + "_stationPhi_" + std::to_string(stationPhiComplete), Vectors.stationEta_perPhi_vec);
      auto stripNumber_perLayer_perPhi_strip_ = Monitored::Collection("stripNumber_strip_" + GeometricSectors::sTGC_Side[iside] + "_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap) + "_stationPhi_" + std::to_string(stationPhiComplete), Vectors.strip_numbers_perPhi_vec);
      auto charge_perLayer_perPhi_strip_ = Monitored::Collection("charge_strip_" + GeometricSectors::sTGC_Side[iside] + "_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap) + "_stationPhi_" + std::to_string(stationPhiComplete), Vectors.charge_perPhi_vec);
      fill("sTGC_sideGroup" + GeometricSectors::sTGC_Side[iside], charge_perLayer_strip_, stationPhi_, stationEta_, stationEta_perPhi_, stripNumber_perLayer_perPhi_strip_, charge_perLayer_perPhi_strip_);    
    }
  
  else if (channel_type == 2)
    {
      auto charge_perLayer_wire_ = Monitored::Collection("charge_wire_" + GeometricSectors::sTGC_Side[iside] + "_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), Vectors.charge_vec);
      auto stationPhi_ = Monitored::Collection("sector_" + GeometricSectors::sTGC_Side[iside] + "_phi_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), Vectors.stationPhi_vec);
      auto stationEta_ = Monitored::Collection("sector_" + GeometricSectors::sTGC_Side[iside] + "_eta_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), Vectors.stationEta_vec);
      fill("sTGC_sideGroup" + GeometricSectors::sTGC_Side[iside], charge_perLayer_wire_, stationPhi_, stationEta_);    
    }

  Vectors.strip_charges_vec.clear();
  Vectors.charge_vec.clear();
  Vectors.stationPhi_vec.clear();
  Vectors.stationEta_vec.clear();
  Vectors.stationEta_perPhi_vec.clear();
  Vectors.strip_numbers_perPhi_vec.clear();
  Vectors.charge_perPhi_vec.clear();
}

StatusCode StgcRawDataMonAlg::fillsTGCHistograms( const Muon::sTgcPrepData* ) const
{
  return StatusCode::SUCCESS;
}

