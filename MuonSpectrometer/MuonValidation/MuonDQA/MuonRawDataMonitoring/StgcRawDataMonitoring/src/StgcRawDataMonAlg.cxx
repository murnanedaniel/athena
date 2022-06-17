/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////////////////////////////////////////
// Package : sTgcRawDataMonAlg
// Author: Sebastian Fuenzalida Garrido
// Local supervisor: Edson Carquin Lopez
// Technical supervisor: Gerardo Vasquez
//
// DESCRIPTION:
// Subject: sTgc --> sTgc raw data monitoring
///////////////////////////////////////////////////////////////////////////////////////////////////////

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

sTgcRawDataMonAlg::sTgcRawDataMonAlg( const std::string& name, ISvcLocator* pSvcLocator ) : AthMonitorAlgorithm(name,pSvcLocator)	      
{
  //Declare the property 
}


StatusCode sTgcRawDataMonAlg::initialize()
{   
  ATH_CHECK(AthMonitorAlgorithm::initialize());
  ATH_CHECK(m_DetectorManagerKey.initialize());
  ATH_CHECK(m_idHelperSvc.retrieve());
  
  ATH_CHECK(m_muonKey.initialize());
  ATH_CHECK(m_meTrkKey.initialize());
  ATH_CHECK(m_sTgcContainerKey.initialize());
  
  return StatusCode::SUCCESS;
} 

StatusCode sTgcRawDataMonAlg::fillHistograms(const EventContext& ctx) const
{  
  //const xAOD::TrackParticleContainer *meTPContainer = nullptr;
  //ATH_CHECK(evtStore() -> retrieve(meTPContainer,"ExtrapolatedMuonTrackParticles" ));

  std::vector<const Muon::sTgcPrepData*> prep_data_vector;
  
  SG::ReadHandle<xAOD::TrackParticleContainer> meTPContainer{m_meTrkKey, ctx};
  SG::ReadHandle<Muon::sTgcPrepDataContainer> sTgc_container(m_sTgcContainerKey, ctx);

  if (m_dosTgcESD) 
    {
      if (m_dosTgcOverview)
	{
	  Histograms::sTgcSummaryHistogramStruct summaryPlots[2][2][4];
            
	  for(const Muon::sTgcPrepDataCollection* coll : *sTgc_container)
	    {
	      for (const Muon::sTgcPrepData* prd : *coll)
		{
		  prep_data_vector.push_back(prd);
		  fillsTgcOverviewHistograms(prep_data_vector, prd);
		  fillsTgcSummaryHistograms(prd, summaryPlots);
		}
	    }
	}
    }
     
  return StatusCode::SUCCESS;
}

void sTgcRawDataMonAlg::fillsTgcOverviewHistograms(const std::vector<const Muon::sTgcPrepData*> &prd, const Muon::sTgcPrepData *sTgc_object) const 
{    
  auto charge_all = Monitored::Collection("charge_all", prd, [] (const Muon::sTgcPrepData *aux) {return aux -> charge();});
  auto numberofstrips_percluster = Monitored::Collection("numberofstrips_percluster", prd, [] (const Muon::sTgcPrepData *aux) {const std::vector<Identifier> &stripIds = aux -> rdoList(); return stripIds.size();});
  
  fill("sTgcMonitor", charge_all, numberofstrips_percluster);
  
  std::vector<short int> strip_times_target = sTgc_object-> stripTimes();
  std::vector<int> strip_charges_target = sTgc_object-> stripCharges();
  std::vector<short unsigned int> strip_number_target = sTgc_object-> stripNumbers();
  std::vector<int> strip_statEta_strip_target = sTgc_object-> stripCharges();

  auto time_all = Monitored::Collection("time_all", prd, [] (const Muon::sTgcPrepData *aux) {return aux -> time();});      
  auto strip_times = Monitored::Collection("strip_times", strip_times_target);
  auto strip_charges = Monitored::Collection("strip_charges", strip_charges_target);
  auto strip_number = Monitored::Collection("strip_number", strip_number_target);

  fill("sTgcMonitor", time_all, strip_times, strip_charges, strip_number);

  auto x_mon = Monitored::Collection("x_mon", prd, [] (const Muon::sTgcPrepData *aux) {Amg::Vector3D pos = aux -> globalPosition(); return pos.x();});
  auto y_mon = Monitored::Collection("y_mon", prd, [] (const Muon::sTgcPrepData *aux) {Amg::Vector3D pos = aux -> globalPosition(); return pos.y();});
  auto z_mon = Monitored::Collection("z_mon", prd, [] (const Muon::sTgcPrepData *aux) {Amg::Vector3D pos = aux -> globalPosition(); return pos.z();});
  auto R_mon = Monitored::Collection("R_mon", prd, [] (const Muon::sTgcPrepData *aux) {Amg::Vector3D pos = aux -> globalPosition(); return std::hypot(pos.x(), pos.y());});

  fill("sTgcMonitor", x_mon, y_mon, z_mon, R_mon);
}

void sTgcRawDataMonAlg::fillsTgcSummaryHistograms(const Muon::sTgcPrepData *sTgc_object, Histograms::sTgcSummaryHistogramStruct (&vects)[2][2][4]) const
{
  Identifier Id    = sTgc_object-> identify();
  int stationPhi   = m_idHelperSvc -> stgcIdHelper().stationPhi(Id);
  int stationEta   = m_idHelperSvc -> stgcIdHelper().stationEta(Id);
  int multiplet    = m_idHelperSvc -> stgcIdHelper().multilayer(Id);
  int gas_gap      = m_idHelperSvc -> stgcIdHelper().gasGap(Id);  
  int channel_type = m_idHelperSvc -> stgcIdHelper().channelType(Id);
  int charge       = sTgc_object-> charge();
  int iside        = (stationEta > 0) ? 1 : 0;
  std::string stationName = m_idHelperSvc -> stgcIdHelper().stationNameString(m_idHelperSvc -> stgcIdHelper().stationName(Id));
  int stationPhiComplete = get_sectorPhi_from_stationPhi_stName(stationPhi, stationName);
 
  auto &Vectors = vects[iside][multiplet - 1][gas_gap - 1];
  Vectors.strip_charges_vec = sTgc_object-> stripCharges();  
  Vectors.stationEta_perPhi_vec.push_back(stationEta);
  Vectors.strip_numbers_perPhi_vec = sTgc_object-> stripNumbers();
  Vectors.charge_perPhi_vec.push_back(charge);
  Vectors.charge_vec.push_back(charge);
  Vectors.stationPhi_vec.push_back(stationPhiComplete);
  Vectors.stationEta_vec.push_back(stationEta);

  if (channel_type == 0)
    {
      auto charge_perLayer_pad_ = Monitored::Collection("charge_pad_" + GeometricSectors::sTgc_Side[iside] + "_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), Vectors.charge_vec);
      auto stationPhi_ = Monitored::Collection("sector_" + GeometricSectors::sTgc_Side[iside] + "_phi_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), Vectors.stationPhi_vec);
      auto stationEta_ = Monitored::Collection("sector_" + GeometricSectors::sTgc_Side[iside] + "_eta_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), Vectors.stationEta_vec);
      fill("sTgc_sideGroup" + GeometricSectors::sTgc_Side[iside], charge_perLayer_pad_, stationPhi_, stationEta_);    
    }
  
  else if (channel_type == 1)
    {      
      auto charge_perLayer_strip_ = Monitored::Collection("charge_strip_" + GeometricSectors::sTgc_Side[iside] + "_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), Vectors.charge_vec);
      auto stationPhi_ = Monitored::Collection("sector_" + GeometricSectors::sTgc_Side[iside] + "_phi_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), Vectors.stationPhi_vec);
      auto stationEta_ = Monitored::Collection("sector_" + GeometricSectors::sTgc_Side[iside] + "_eta_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), Vectors.stationEta_vec); 
      auto stationEta_perPhi_ = Monitored::Collection("sector_" + GeometricSectors::sTgc_Side[iside] + "_eta_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap) + "_stationPhi_" + std::to_string(stationPhiComplete), Vectors.stationEta_perPhi_vec);
      auto stripNumber_perLayer_perPhi_strip_ = Monitored::Collection("stripNumber_strip_" + GeometricSectors::sTgc_Side[iside] + "_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap) + "_stationPhi_" + std::to_string(stationPhiComplete), Vectors.strip_numbers_perPhi_vec);
      auto charge_perLayer_perPhi_strip_ = Monitored::Collection("charge_strip_" + GeometricSectors::sTgc_Side[iside] + "_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap) + "_stationPhi_" + std::to_string(stationPhiComplete), Vectors.charge_perPhi_vec);
      fill("sTgc_sideGroup" + GeometricSectors::sTgc_Side[iside], charge_perLayer_strip_, stationPhi_, stationEta_, stationEta_perPhi_, stripNumber_perLayer_perPhi_strip_, charge_perLayer_perPhi_strip_);    
    }
  
  else if (channel_type == 2)
    {
      auto charge_perLayer_wire_ = Monitored::Collection("charge_wire_" + GeometricSectors::sTgc_Side[iside] + "_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), Vectors.charge_vec);
      auto stationPhi_ = Monitored::Collection("sector_" + GeometricSectors::sTgc_Side[iside] + "_phi_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), Vectors.stationPhi_vec);
      auto stationEta_ = Monitored::Collection("sector_" + GeometricSectors::sTgc_Side[iside] + "_eta_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), Vectors.stationEta_vec);
      fill("sTgc_sideGroup" + GeometricSectors::sTgc_Side[iside], charge_perLayer_wire_, stationPhi_, stationEta_);    
    }

  Vectors.strip_charges_vec.clear();
  Vectors.charge_vec.clear();
  Vectors.stationPhi_vec.clear();
  Vectors.stationEta_vec.clear();
  Vectors.stationEta_perPhi_vec.clear();
  Vectors.strip_numbers_perPhi_vec.clear();
  Vectors.charge_perPhi_vec.clear();
}

