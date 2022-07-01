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

#include "StgcRawDataMonitoring/StgcRawDataMonAlg.h"
//Stgc cxx includes

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
  ATH_CHECK(m_idHelperSvc.retrieve());
  ATH_CHECK(m_sTgcContainerKey.initialize());
  ATH_CHECK(m_muonKey.initialize());
  ATH_CHECK(m_meTrkKey.initialize());
  return StatusCode::SUCCESS;
} 

StatusCode sTgcRawDataMonAlg::fillHistograms(const EventContext& ctx) const
{  
  int lumiblock = -1;
  lumiblock = GetEventInfo(ctx)->lumiBlock();
  SG::ReadHandle<Muon::sTgcPrepDataContainer> sTgc_container(m_sTgcContainerKey, ctx);
  ATH_CHECK(sTgc_container.isValid());
  
  if (m_dosTgcESD && m_dosTgcOverview)  
    {
      for(const Muon::sTgcPrepDataCollection* coll : *sTgc_container)
	{	  
	  for (const Muon::sTgcPrepData* prd : *coll)
	    {
	      fillsTgcOverviewHistograms(prd, *coll);
	      fillsTgcSummaryHistograms(prd);
	    }
	}
    }
  SG::ReadHandle<xAOD::TrackParticleContainer> meTPContainer{m_meTrkKey,ctx};
    if (!meTPContainer.isValid()) {
	ATH_MSG_FATAL("Nope. Could not retrieve "<<m_meTrkKey.fullKey());
	return StatusCode::FAILURE;
    }
    clusterFromTrack(meTPContainer.cptr(),lumiblock);
   
  return StatusCode::SUCCESS;
}

void sTgcRawDataMonAlg::fillsTgcOverviewHistograms(const Muon::sTgcPrepData *sTgc_object, const Muon::MuonPrepDataCollection<Muon::sTgcPrepData> &prd) const 
{   
  auto charge_all = Monitored::Collection("charge_all", prd, [] (const Muon::sTgcPrepData *aux) 
					  {
					    return aux -> charge();
					  });
  
  auto numberofstrips_percluster = Monitored::Collection("numberofstrips_percluster", prd, [] (const Muon::sTgcPrepData *aux) 
							 {
							   const std::vector<Identifier> &stripIds = aux -> rdoList(); 
							   return stripIds.size();
							 });
  
  fill("sTgcMonitor", charge_all, numberofstrips_percluster);
  
  std::vector<short int> strip_times_target = sTgc_object-> stripTimes();
  std::vector<int> strip_charges_target = sTgc_object-> stripCharges();
  std::vector<short unsigned int> strip_number_target = sTgc_object-> stripNumbers();

  auto time_all = Monitored::Collection("time_all", prd, [] (const Muon::sTgcPrepData *aux) 
					{
					  return aux -> time();
					});
      
  auto strip_times = Monitored::Collection("strip_times", strip_times_target);
  auto strip_charges = Monitored::Collection("strip_charges", strip_charges_target);
  auto strip_number = Monitored::Collection("strip_number", strip_number_target);

  fill("sTgcMonitor", time_all, strip_times, strip_charges, strip_number);

  auto x_mon = Monitored::Collection("x_mon", prd, [] (const Muon::sTgcPrepData *aux) 
				     {
				       Amg::Vector3D pos = aux -> globalPosition(); 
				       return pos.x();
				     });
  
  auto y_mon = Monitored::Collection("y_mon", prd, [] (const Muon::sTgcPrepData *aux) 
				     {
				       Amg::Vector3D pos = aux -> globalPosition(); 
				       return pos.y();
				     });
  
  auto z_mon = Monitored::Collection("z_mon", prd, [] (const Muon::sTgcPrepData *aux) 
				     {
				       Amg::Vector3D pos = aux -> globalPosition(); 
				       return pos.z();
				     });
  
  auto R_mon = Monitored::Collection("R_mon", prd, [] (const Muon::sTgcPrepData *aux) 
				     {
				       Amg::Vector3D pos = aux -> globalPosition(); 
				       return std::hypot(pos.x(), pos.y());
				     });

  fill("sTgcMonitor", x_mon, y_mon, z_mon, R_mon);
}

void sTgcRawDataMonAlg::fillsTgcSummaryHistograms(const Muon::sTgcPrepData *sTgc_object) const
{
  Identifier Id    = sTgc_object   -> identify();
  int stationPhi   = m_idHelperSvc -> stgcIdHelper().stationPhi(Id);
  int stationEta   = m_idHelperSvc -> stgcIdHelper().stationEta(Id);
  int multiplet    = m_idHelperSvc -> stgcIdHelper().multilayer(Id);
  int gas_gap      = m_idHelperSvc -> stgcIdHelper().gasGap(Id);  
  int channel_type = m_idHelperSvc -> stgcIdHelper().channelType(Id);
  int charge       = sTgc_object   -> charge();
  int iside        = (stationEta > 0) ? 1 : 0;
  std::string stationName = m_idHelperSvc -> stgcIdHelper().stationNameString(m_idHelperSvc -> stgcIdHelper().stationName(Id));
  int stationPhiComplete = get_sectorPhi_from_stationPhi_stName(stationPhi, stationName);
 
  std::vector<int> strip_charges_vec = sTgc_object -> stripCharges();  
  std::vector<short unsigned int> strip_numbers_perPhi_vec = sTgc_object -> stripNumbers();

  std::vector<int> charge_vec;
  std::vector<int> stationPhi_vec;
  std::vector<int> stationEta_vec;
  std::vector<int> stationEta_perPhi_vec;

  charge_vec.push_back(charge);
  stationPhi_vec.push_back(stationPhiComplete);
  stationEta_vec.push_back(stationEta);
  stationEta_perPhi_vec.push_back(stationEta);
  
  if (channel_type == 0)
    {
      auto charge_perLayer_pad_ = Monitored::Collection("charge_pad_" + GeometricSectors::sTgc_Side[iside] + "_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), charge_vec);
      auto stationPhi_ = Monitored::Collection("sector_" + GeometricSectors::sTgc_Side[iside] + "_phi_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), stationPhi_vec);
      auto stationEta_ = Monitored::Collection("sector_" + GeometricSectors::sTgc_Side[iside] + "_eta_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), stationEta_vec);
      fill("sTgc_sideGroup" + GeometricSectors::sTgc_Side[iside], charge_perLayer_pad_, stationPhi_, stationEta_);    
    }
  
  else if (channel_type == 1)
    {      
      auto charge_perLayer_strip_ = Monitored::Collection("charge_strip_" + GeometricSectors::sTgc_Side[iside] + "_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), charge_vec);
      auto stationPhi_ = Monitored::Collection("sector_" + GeometricSectors::sTgc_Side[iside] + "_phi_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), stationPhi_vec);
      auto stationEta_ = Monitored::Collection("sector_" + GeometricSectors::sTgc_Side[iside] + "_eta_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), stationEta_vec); 
      auto stationEta_perPhi_ = Monitored::Collection("sector_" + GeometricSectors::sTgc_Side[iside] + "_eta_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap) + "_stationPhi_" + std::to_string(stationPhiComplete), stationEta_perPhi_vec);
      auto stripNumber_perLayer_perPhi_strip_ = Monitored::Collection("stripNumber_strip_" + GeometricSectors::sTgc_Side[iside] + "_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap) + "_stationPhi_" + std::to_string(stationPhiComplete), strip_numbers_perPhi_vec);
      auto charge_perLayer_perPhi_strip_ = Monitored::Collection("charge_strip_" + GeometricSectors::sTgc_Side[iside] + "_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap) + "_stationPhi_" + std::to_string(stationPhiComplete), charge_vec);
      fill("sTgc_sideGroup" + GeometricSectors::sTgc_Side[iside], charge_perLayer_strip_, stationPhi_, stationEta_, stationEta_perPhi_, stripNumber_perLayer_perPhi_strip_, charge_perLayer_perPhi_strip_);    
    }
  
  else if (channel_type == 2)
    {
      auto charge_perLayer_wire_ = Monitored::Collection("charge_wire_" + GeometricSectors::sTgc_Side[iside] + "_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), charge_vec);
      auto stationPhi_ = Monitored::Collection("sector_" + GeometricSectors::sTgc_Side[iside] + "_phi_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), stationPhi_vec);
      auto stationEta_ = Monitored::Collection("sector_" + GeometricSectors::sTgc_Side[iside] + "_eta_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), stationEta_vec);
      fill("sTgc_sideGroup" + GeometricSectors::sTgc_Side[iside], charge_perLayer_wire_, stationPhi_, stationEta_);    
    }
}
void sTgcRawDataMonAlg::clusterFromTrack(const xAOD::TrackParticleContainer*  muonContainer, int lb) const
{
        std::vector<int> sector_CSide_eta1_ontrack;
        std::vector<int> stationPhi_CSide_eta1_ontrack;
	for(const xAOD::TrackParticle* meTP : *muonContainer) {
		
		if(!meTP) continue;

		auto eta_trk = Monitored::Scalar<float>("eta_trk", meTP->eta());
                auto phi_trk = Monitored::Scalar<float>("phi_trk", meTP->phi());
                fill("sTgcMonitor",eta_trk,phi_trk);
                const Trk::Track* meTrack = meTP->track();
                if(!meTrack) continue;
		// get the vector of measurements on track
		const DataVector<const Trk::MeasurementBase>* meas = meTrack->measurementsOnTrack();
                bool isMM=false;
                for(const Trk::MeasurementBase* it : *meas) {
			const Trk::RIO_OnTrack* rot = dynamic_cast<const Trk::RIO_OnTrack*>(it);
                        if(!rot) continue;
			Identifier rot_id = rot->identify();
                        if(!m_idHelperSvc->isMM(rot_id)) continue;
			isMM=true;
                        const Muon::MMClusterOnTrack* cluster = dynamic_cast<const Muon::MMClusterOnTrack*>(rot);
			if(!cluster) continue;

			std::string stName = m_idHelperSvc->mmIdHelper().stationNameString(m_idHelperSvc->mmIdHelper().stationName(rot_id));
			int stEta          = m_idHelperSvc->mmIdHelper().stationEta(rot_id);
			int stPhi          = m_idHelperSvc->mmIdHelper().stationPhi(rot_id);
			int multi          = m_idHelperSvc->mmIdHelper().multilayer(rot_id);
			int gap            = m_idHelperSvc->mmIdHelper().gasGap(rot_id);
			int ch             = m_idHelperSvc->mmIdHelper().channel(rot_id);
                        // MMS and MML phi sectors
			//				int phisec = (stNumber%2==0) ? 1 : 0;
			int sectorPhi = get_sectorPhi_from_stationPhi_stName(stPhi,stName); // 1->16
			int PCB = get_PCB_from_channel(ch);
			int iside = (stEta > 0) ? 1 : 0;
                        const int pcb_counter = 5;
                        if(stEta==-1) {
			 stationPhi_CSide_eta1_ontrack.push_back(sectorPhi);
                    sector_CSide_eta1_ontrack.push_back(get_bin_for_occ_CSide_pcb_eta1_hist(stEta, multi, gap, PCB));
                         auto stationPhi_cSide_eta1_ontrack = Monitored::Collection("stationPhi_CSide_eta1_ontrack",stationPhi_CSide_eta1_ontrack);
                         auto sector_cSide_eta1_ontrack = Monitored::Collection("sector_CSide_eta1_ontrack",sector_CSide_eta1_ontrack); 
                      
                       fill("sTgcMonitor",stationPhi_cSide_eta1_ontrack,sector_cSide_eta1_ontrack);
			}


                        
                 

                }        
                     
                          
	        
        }
        
}
