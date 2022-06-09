/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Package : sTGCRawDataMonAlg
// Author: Sebastian Fuenzalida 
// Local supervisor: Edson Carquin 
// Technical supervisor: Gerardo Vasquez
//
// DESCRIPTION:
// Subject: sTGC-->Offline Muon Data Quality
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

/*
namespace GeometricSectors
{
  static const std::array<std::string,2> sTGC_Side = {"CSide", "ASide"};
  static const std::array<std::string,2> sTGC_Sector = {"S", "L"};
  static const std::array<std::string,2> EtaSector = {"1","2"};
}
*/

  
//static constexpr double const toDeg = 180/M_PI;

//1e=1.6X10-4 fC                                                             
//static constexpr double conversion_charge=1.6E-04;

  
//struct sTGCOverviewHistogramStruct {
//std::vector<int> charge_all;
//std::vector<int> numberofstrips_percluster;
/*
    std::vector<int> statEta_strip
    std::vector<int> time_all;
    std::vector<int> strip_number;
    std::vector<short int> strip_times;
    std::vector<int> strip_charges;
    std::vector<float> R_mon;
    std::vector<float> z_mon;
    std::vector<float> x_mon;
    std::vector<float> y_mon;
    std::vector<int> stationPhi_ASide_eta1_ontrack;
    std::vector<int> stationPhi_ASide_eta2_ontrack;
    std::vector<int> stationPhi_CSide_eta1_ontrack;
    std::vector<int> stationPhi_CSide_eta2_ontrack;
    std::vector<int> sector_ASide_eta2_ontrack;
    std::vector<int> sector_CSide_eta1_ontrack;
    std::vector<int> sector_ASide_eta1_ontrack;
    std::vector<int> sector_CSide_eta2_ontrack;
    std::vector<int> sector_lb_ASide_eta1_ontrack;
    std::vector<int> sector_lb_ASide_eta2_ontrack;
    std::vector<int> sector_lb_CSide_eta1_ontrack;
    std::vector<int> sector_lb_CSide_eta2_ontrack;
    std::vector<int> stationPhi_ASide_eta1;
    std::vector<int> stationPhi_ASide_eta2;
    std::vector<int> stationPhi_CSide_eta1;
    std::vector<int> stationPhi_CSide_eta2;
    std::vector<int> sector_CSide_eta2;
    std::vector<int> sector_ASide_eta2;
    std::vector<int> sector_CSide_eta1;
    std::vector<int> sector_ASide_eta1;
    std::vector<int> sector_lb_ASide_eta1;
    std::vector<int> sector_lb_ASide_eta2;
std::vector<int> sector_lb_CSide_eta1;
std::vector<int> sector_lb_CSide_eta2;
    */
//};

//struct sTGCSummaryHistogramStruct {
    
/*
    std::vector<int> strip_number;
    std::vector<int> sector_strip;
    std::vector<int> charge;
    std::vector<short int> strip_times;
    std::vector<int> strip_charges;
    std::vector<float> x_ontrack;
    std::vector<float> y_ontrack;
*/
//std::vector<float> residuals;
//};

//}

/////////////////////////////////////////////////////////////////////////////
// *********************************************************************
// Public Methods
// ********************************************************************* 

StgcRawDataMonAlg::StgcRawDataMonAlg( const std::string& name, ISvcLocator* pSvcLocator ) : AthMonitorAlgorithm(name,pSvcLocator)
											    //m_muonSelectionTool("CP::MuonSelectionTool/MuonSelectionTool")
											    //m_sTGCContainerKey("STGC_Measurements")
{
  //Declare the property 
  //declareProperty("sTGCPrepDataContainerName",m_sTGCContainerKey);
}

/*---------------------------------------------------------*/
StatusCode StgcRawDataMonAlg::initialize()
/*---------------------------------------------------------*/
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
  //const Muon::sTgcPrepData *stg_object;
  
  //int lumiblock = -1;
  //int charge = -1;
  
  //lumiblock = GetEventInfo(ctx)->lumiBlock();

  //ATH_MSG_INFO("sTGCRawDataMonAlg::sTGC RawData Monitoring Histograms being filled" );

  SG::ReadHandle<Muon::sTgcPrepDataContainer> sTGC_container(m_sTGCContainerKey,ctx);

  //int counter_1 = 0;
  //int counter_2 = 0;
  //ATH_MSG_INFO("****** sTGCContainer->size() : " << sTGC_container->size());

  if (m_doSTGCESD) 
    {
      Histograms::sTGCSummaryHistogramStruct summaryPlots[2][2][4];
      //sTGCOverviewHistogramStruct overviewPlots;
      //sTGCSummaryHistogramStruct summaryPlots[2][2][8][2][2][4];
      
      //loop in sTGCPrepDataContainer
      
      
      
      /*
      for(const Muon::sTgcPrepDataCollection* coll : *sTGC_container)
	{
	  for (const Muon::sTgcPrepData* prd : *coll)
	    {
	      //ATH_CHECK(fillsTGCOverviewVects(prd, overviewPlots));
	      //ATH_MSG_INFO("fillsTGCOverviewVects was successful");
	      //ATH_CHECK(fillsTGCSummaryVects(prd, summaryPlots));
	      //ATH_MSG_INFO("fillsTGCSummaryVects was successful");
	      //ATH_CHECK(fillsTGCHistograms(prd));
	      //ATH_MSG_INFO("fillsTGCHistograms was successful");
	      prep_data_vector.push_back(prd);
	      //++counter_2;
	      //ATH_MSG_INFO("COUNTER2: " << counter_2);
	    }

	  //++counter_1;
	  //ATH_MSG_INFO("COUNTER1: " << counter_1);
	}
      */

      //const std::vector<int> strip_charges_vector; 
      
      for(const Muon::sTgcPrepDataCollection* coll : *sTGC_container)
	{
	  for (const Muon::sTgcPrepData* prd : *coll)
	    {
	      prep_data_vector.push_back(prd);
	      if (m_do_sTgc_overview) fillsTGCOverviewHistograms(prep_data_vector, prd);
	      //clusterFromTrack(meTPContainer, lumiblock, prd);
	      fillsTGCSummaryHistograms(prd, summaryPlots);
	      //strip_charges_vector = prd -> stripCharges();
	      //if (prd -> stripCharges().size() != 0) fillsTGCSummaryHistograms(prd, summaryPlots); // Put a condition to call this function. Size strip charge diff
	      //fillsTGCSummaryHistograms(prd, summaryPlots);
	    }
	}

      //ATH_MSG_INFO("fillsTGCOverviewHistograms was successful");
      
      //ATH_CHECK( fillsTGCSummaryHistograms(summaryPlots) );
      
      //ATH_MSG_INFO("fillsTGCSummaryHistograms was successful");
      
      
      /*
      for(const Muon::sTgcPrepDataCollection* coll : *sTGC_container)
	{
	  for (const Muon::sTgcPrepData* prd : *coll)
	    {
	      clusterFromTrack(meTPContainer, lumiblock, prd);
	    }
	}
      */
      
    }
  
  //ATH_MSG_INFO("fill Histograms was successful");
  
  prep_data_vector.clear();
  
  return StatusCode::SUCCESS;
}


//StatusCode StgcRawDataMonAlg::fillsTGCOverviewVects( const Muon::sTgcPrepData* prd, sTGCOverviewHistogramStruct& vects ) const {
  
  //const std::vector<Identifier>& stripIds = prd->rdoList();
  // number of strips in this cluster (cluster size)
  //unsigned int nStrips = stripIds.size();  
  // Returns the charge (number of electrons) converted in fC
  //int charge=prd->charge()*conversion_charge;

  //ATH_MSG_INFO("Now on fillsTGCOverviewVects function");
  /*
  Identifier Id = prd->identify();
                 
  const std::vector<uint16_t>& stripNumbers=prd->stripNumbers();
 
  std::string stName   = m_idHelperSvc->stgcIdHelper().stationNameString(m_idHelperSvc->stgcIdHelper().stationName(Id));
  int gas_gap          = m_idHelperSvc->stgcIdHelper().gasGap(Id);
  int stationNumber    = m_idHelperSvc->stgcIdHelper().stationName(Id);
  int stationEta       = m_idHelperSvc->stgcIdHelper().stationEta(Id);
  int stationPhi       = m_idHelperSvc->stgcIdHelper().stationPhi(Id);
  int multiplet        = m_idHelperSvc->stgcIdHelper().multilayer(Id);
  int channel          = m_idHelperSvc->stgcIdHelper().channel(Id);

  ATH_MSG_INFO("stName: " << stName << " stationPhi: " << stationPhi);
  */  
  // Returns the time (in ns)
  //int drift_time=prd->time();
  // Returns the microTPC angle (radians converted in degrees)
  //float mu_TPC_angle=prd->angle()*toDeg;
  // Returns the microTPC chisq Prob.
  //float mu_TPC_chi2=prd->chisqProb();
  //const std::vector<short int>& strip_times=prd->stripTimes();
  //const std::vector<int>& strip_charges=prd->stripCharges();

  //Amg::Vector3D pos    = prd->globalPosition();
  
  //float R=std::hypot(pos.x(),pos.y());

  //int PCB=get_PCB_from_channel(channel);

  //MM gaps are back to back, so the direction of the drift (time) is different for the even and odd gaps -> flip for the even gaps

  //if (gas_gap % 2 == 0) { mu_TPC_angle = -mu_TPC_angle; } 
  
  //vects.charge_all.push_back(charge);
  //vects.numberofstrips_percluster.push_back(nStrips);
  /*
  vects.time_all.push_back(drift_time);
  //vects.strip_times.insert(strip_times);
  //vects.strip_charges.insert(strip_charges);
  vects.x_mon.push_back(pos.x());
  vects.y_mon.push_back(pos.y());
  vects.z_mon.push_back(pos.z());
vects.R_mon.push_back(R);
  */
    

  /*
  //MMS and MML phi sectors
  int phisec=0;
  if (stationNumber%2 == 0) phisec=1;
    
  //16 phi sectors, 8 stationPhi times 2 stName, MMS and MML
  int sectorPhi=get_sectorPhi_from_stationPhi_stName(stationPhi,stName);  

  ATH_MSG_INFO("After get_sectorPhi_from_stationPhi_stName function" << vects.charge_all.size());
  
  //Occupancy plots with PCB granularity further divided for each eta sector: -2, -1, 1, 2
  
  //Filling Vectors for stationEta=-2
  if (stationEta==-2){
    vects.sector_CSide_eta2.push_back(get_bin_for_occ_CSide_pcb_eta2_hist(stationEta,multiplet,gas_gap,PCB));
    vects.stationPhi_CSide_eta2.push_back(sectorPhi);
    vects.sector_lb_CSide_eta2.push_back(get_bin_for_occ_lb_CSide_pcb_eta2_hist(stationEta,multiplet,gas_gap,PCB,phisec));
  }
  //Filling Vectors for stationEta=-1
  else if (stationEta==-1){
    vects.sector_CSide_eta1.push_back(get_bin_for_occ_CSide_pcb_eta1_hist(stationEta,multiplet,gas_gap,PCB));
    vects.stationPhi_CSide_eta1.push_back(sectorPhi);
    vects.sector_lb_CSide_eta1.push_back(get_bin_for_occ_lb_CSide_pcb_eta1_hist(stationEta,multiplet,gas_gap,PCB,phisec));
  }
  //Filling Vectors for stationEta=1
  else if (stationEta==1){
    vects.sector_ASide_eta1.push_back(get_bin_for_occ_ASide_pcb_eta1_hist(stationEta,multiplet,gas_gap,PCB));
    vects.stationPhi_ASide_eta1.push_back(sectorPhi);
    vects.sector_lb_ASide_eta1.push_back(get_bin_for_occ_lb_ASide_pcb_eta1_hist(stationEta,multiplet,gas_gap,PCB,phisec));
  }
  //Filling Vectors for stationEta=2
  else {
    vects.sector_ASide_eta2.push_back(get_bin_for_occ_ASide_pcb_eta2_hist(stationEta,multiplet,gas_gap,PCB));
    vects.stationPhi_ASide_eta2.push_back(sectorPhi);
    vects.sector_lb_ASide_eta2.push_back(get_bin_for_occ_lb_ASide_pcb_eta2_hist(stationEta,multiplet,gas_gap,PCB,phisec));
  }

  ATH_MSG_INFO("Checkpoint 3");

  //loop on each strip                                                                                       
  
  int sIdx = 0; // index-counter for the vector of Id's                                
  for (const Identifier& id : stripIds){    
    
    std::string stName_strip   = m_idHelperSvc->stgcIdHelper().stationNameString(m_idHelperSvc->stgcIdHelper().stationName(id));
    int stationEta_strip       = m_idHelperSvc->stgcIdHelper().stationEta(id);
    ATH_MSG_INFO("sIdx: " << sIdx);    
    vects.statEta_strip.push_back(stationEta_strip);
    ATH_MSG_INFO("Checkpoint 4 "<< stripNumbers.size());
    if(stripNumbers.size()!=0) vects.strip_number.push_back(stripNumbers[sIdx]);
    ATH_MSG_INFO("Checkpoint 5 ");

  }
  */
  //ATH_MSG_INFO("Running of fillsTGCOverviewVects was successful");

  //return StatusCode::SUCCESS;
//}


void StgcRawDataMonAlg::fillsTGCOverviewHistograms(const std::vector<const Muon::sTgcPrepData*> &prd, const Muon::sTgcPrepData *stg_object) const 
{    
  //1e=1.6X10-4 fC                                                                                          
  //static constexpr double conversion_charge=1.6E-04;
 
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

  /*
  ATH_MSG_INFO("Strip Times: " << stg_object -> stripTimes());
  ATH_MSG_INFO("Strip Number: " << stg_object -> stripNumbers());
  ATH_MSG_INFO("Strip Charges: " << stg_object -> stripCharges());
  ATH_MSG_INFO("Charge: " << stg_object -> charge());
  */

  auto time_all = Monitored::Collection("time_all", prd, [] (const Muon::sTgcPrepData *aux) {return aux -> time();});      
  auto strip_times = Monitored::Collection("strip_times", strip_times_target);
  auto strip_charges = Monitored::Collection("strip_charges", strip_charges_target);
  auto strip_number = Monitored::Collection("strip_number", strip_number_target);
  //auto statEta_strip = Monitored::Collection("statEta_strip", );

  fill("sTGCMonitor", time_all, strip_times, strip_charges, strip_number);

  auto x_mon = Monitored::Collection("x_mon", prd, [] (const Muon::sTgcPrepData *aux) {Amg::Vector3D pos = aux -> globalPosition(); return pos.x();});
  auto y_mon = Monitored::Collection("y_mon", prd, [] (const Muon::sTgcPrepData *aux) {Amg::Vector3D pos = aux -> globalPosition(); return pos.y();});
  auto z_mon = Monitored::Collection("z_mon", prd, [] (const Muon::sTgcPrepData *aux) {Amg::Vector3D pos = aux -> globalPosition(); return pos.z();});
  auto R_mon = Monitored::Collection("R_mon", prd, [] (const Muon::sTgcPrepData *aux) {Amg::Vector3D pos = aux -> globalPosition(); return std::hypot(pos.x(), pos.y());});

  fill("sTGCMonitor", x_mon, y_mon, z_mon, R_mon);
  
  // We need to see carefully this part!

  /*
  Identifier Id = stg_object -> identify();
  int gas_gap    = m_idHelperSvc -> stgcIdHelper().gasGap(Id);
  int multiplet    = m_idHelperSvc -> stgcIdHelper().multilayer(Id);
  int stationEta    = m_idHelperSvc -> stgcIdHelper().stationEta(Id);
  int stationPhi     = m_idHelperSvc -> stgcIdHelper().stationPhi(Id);
  int channel        = m_idHelperSvc -> stgcIdHelper().channel(Id);
  int PCB = get_PCB_from_channel(channel);
  int iside = (stationEta > 0) ? 1 : 0;
  int channel_type       = m_idHelperSvc -> stgcIdHelper().channelType(Id);
  
  
  ATH_MSG_INFO("In fillsTGCOverviewHistograms");
  ATH_MSG_INFO("Side: " << iside << ", PCB: " << PCB << ", station phi: " << stationPhi << ", station eta: " << stationEta << ", multiplet: " << multiplet << ", gas gap: " << gas_gap << ", strip charge: " << stg_object -> stripCharges() << ", strip number: " << stg_object -> stripNumbers() << ", channel type: " << channel_type << ", charge: " << stg_object -> charge());
  

  const int pcb_counter = 5;
  int PCBeta12 = (std::abs(stationEta) == 2) ? (PCB + pcb_counter) : PCB;

  auto lb_mon = Monitored::Scalar<int>("lb_mon", lb);

  std::vector<int> sector_lb_vec;
  sector_lb_vec.push_back(get_bin_for_occ_lb_pcb_hist(multiplet, gas_gap, PCBeta12));

  for (int statPhi = 0; statPhi < 16; ++statPhi)
    {
      for (int iside = 0; iside < 2; ++iside)
	{
	  auto sector_lb = Monitored::Collection("sector_lb_" + GeometricSectors::sTGC_Side[iside] + "_phi" + std::to_string(statPhi + 1), sector_lb_vec);
	  std::string sTGC_sideGroup = "sTGC_sideGroup" + GeometricSectors::sTGC_Side[iside];
	  fill(sTGC_sideGroup, lb_mon, sector_lb);
	}
    }

  sector_lb_vec.clear();

  */

  //  auto &Vectors = vects[iside][stationPhi - 1][stationEta][multiplet - 1][gas_gap - 1];
  //Vectors.strip_charges_vec = stg_object -> stripCharges();


  /*
  auto time_all = Monitored::Collection("time_all", vects.time_all);
  auto strip_times = Monitored::Collection("strip_times", vects.strip_times);
  auto strip_charges = Monitored::Collection("strip_charges", vects.strip_charges);
  auto strip_number = Monitored::Collection("strip_number", vects.strip_number);
  auto statEta_strip = Monitored::Collection("statEta_strip", vects.statEta_strip);

  fill("sTGCMonitor",time_all,strip_times,strip_charges,strip_number,statEta_strip);

  auto x_mon = Monitored::Collection("x_mon", vects.x_mon);
  auto y_mon = Monitored::Collection("y_mon", vects.y_mon);
  auto z_mon = Monitored::Collection("z_mon", vects.z_mon);
  auto R_mon = Monitored::Collection("R_mon", vects.R_mon);

  fill("sTGCMonitor",x_mon,y_mon,z_mon,R_mon);
    
  auto lb_mon = Monitored::Scalar<int>("lb_mon", lb);
    
  auto sector_lb_CSide_eta2 = Monitored::Collection("sector_lb_CSide_eta2",vects.sector_lb_CSide_eta2);
  auto sector_lb_CSide_eta1 = Monitored::Collection("sector_lb_CSide_eta1",vects.sector_lb_CSide_eta1);
  auto sector_lb_ASide_eta2 = Monitored::Collection("sector_lb_ASide_eta2",vects.sector_lb_ASide_eta2);
  auto sector_lb_ASide_eta1 = Monitored::Collection("sector_lb_ASide_eta1",vects.sector_lb_ASide_eta1);
    
  fill("sTGCMonitor",lb_mon,sector_lb_CSide_eta2,sector_lb_CSide_eta1,sector_lb_ASide_eta1,sector_lb_ASide_eta2);
    
  auto sector_CSide_eta2 = Monitored::Collection("sector_CSide_eta2",vects.sector_CSide_eta2);
  auto sector_CSide_eta1 = Monitored::Collection("sector_CSide_eta1",vects.sector_CSide_eta1);
  auto sector_ASide_eta1 = Monitored::Collection("sector_ASide_eta1",vects.sector_ASide_eta1);
  auto sector_ASide_eta2 = Monitored::Collection("sector_ASide_eta2",vects.sector_ASide_eta2);
  auto stationPhi_CSide_eta1 = Monitored::Collection("stationPhi_CSide_eta1",vects.stationPhi_CSide_eta1);
  auto stationPhi_CSide_eta2 = Monitored::Collection("stationPhi_CSide_eta2",vects.stationPhi_CSide_eta2);
  auto stationPhi_ASide_eta1 = Monitored::Collection("stationPhi_ASide_eta1",vects.stationPhi_ASide_eta1);
  auto stationPhi_ASide_eta2 = Monitored::Collection("stationPhi_ASide_eta2",vects.stationPhi_ASide_eta2);
  
  fill("sTGCMonitor",sector_CSide_eta1,sector_CSide_eta2,sector_ASide_eta1,sector_ASide_eta2,stationPhi_CSide_eta1,stationPhi_CSide_eta2,stationPhi_ASide_eta1,stationPhi_ASide_eta2);
  */
}


//StatusCode StgcRawDataMonAlg::fillsTGCSummaryVects( const Muon::sTgcPrepData* prd, sTGCSummaryHistogramStruct (&vects)[2][2][8][2][2][4]) const{
/*
  Identifier Id = prd->identify();
  const std::vector<Identifier>& stripIds = prd->rdoList();

  std::string stName   = m_idHelperSvc->stgcIdHelper().stationNameString(m_idHelperSvc->stgcIdHelper().stationName(Id));
  int thisStationNumber    = m_idHelperSvc->stgcIdHelper().stationName(Id);
  int thisStationEta       = m_idHelperSvc->stgcIdHelper().stationEta(Id);
  int thisStationPhi       = m_idHelperSvc->stgcIdHelper().stationPhi(Id);
  int thisMultiplet        = m_idHelperSvc->stgcIdHelper().multilayer(Id);
  int thisGasgap          = m_idHelperSvc->stgcIdHelper().gasGap(Id);
  int thisCharge=prd->charge()*conversion_charge;

  //float thisMu_TPC_angle=prd->angle()*toDeg;
    
  //if ( thisGasgap % 2 == 0 ) { thisMu_TPC_angle = - thisMu_TPC_angle; }
    
  //MMS and MML phi sectors
  int phisec=0;
  if (thisStationNumber%2 == 0) phisec=1;
  
  //CSide and ASide
  int iside=0;
  if(thisStationEta>0) iside=1;
  
  //2 eta sectors depending on Eta=+-1 (0) and +-2 (1)
  int sectorEta=get_sectorEta_from_stationEta(thisStationEta);

  auto& Vectors = vects[iside][phisec][thisStationPhi-1][sectorEta][thisMultiplet-1][thisGasgap-1];
  
  //Vectors.mu_TPC_angle.push_back(thisMu_TPC_angle);
  Vectors.charge.push_back(thisCharge);
  
  //loop on strips
  int sIdx = 0;
  const std::vector<uint16_t>& stripNumbers=prd->stripNumbers();
  
  for ( const Identifier& id : stripIds){
    
  int stationEta       = m_idHelperSvc->stgcIdHelper().stationEta(id);
int gas_gap          = m_idHelperSvc->stgcIdHelper().gasGap(Id);
int multiplet        = m_idHelperSvc->stgcIdHelper().multilayer(Id); 

//    Filling Vectors for both sides, considering each strip  
if(stripNumbers.size()!=0)Vectors.strip_number.push_back(stripNumbers[sIdx]);
if(iside==1)    Vectors.sector_strip.push_back(get_bin_for_occ_ASide_hist(stationEta,multiplet,gas_gap));
if(iside==0)    Vectors.sector_strip.push_back(get_bin_for_occ_CSide_hist(stationEta,multiplet,gas_gap));
    
}
*/
//return StatusCode::SUCCESS;
//}




void StgcRawDataMonAlg::fillsTGCSummaryHistograms(const Muon::sTgcPrepData *stg_object, Histograms::sTGCSummaryHistogramStruct (&vects)[2][2][4]) const
{
  
  Identifier Id = stg_object -> identify();
  
  int stationPhi     = m_idHelperSvc -> stgcIdHelper().stationPhi(Id);
  int stationEta    = m_idHelperSvc -> stgcIdHelper().stationEta(Id);
  int multiplet    = m_idHelperSvc -> stgcIdHelper().multilayer(Id);
  int gas_gap    = m_idHelperSvc -> stgcIdHelper().gasGap(Id);
  
  int iside = (stationEta > 0) ? 1 : 0;
  
  //int channel        = m_idHelperSvc -> stgcIdHelper().channel(Id);
  //int PCB = get_PCB_from_channel(channel);
  int channel_type       = m_idHelperSvc -> stgcIdHelper().channelType(Id);
 
  int charge = stg_object -> charge();
  
  /*
  if (channel_type == 1)
    {
      ATH_MSG_INFO("In fillsTGCSummaryHistograms");
      ATH_MSG_INFO("Side: " << iside << ", PCB: " << PCB << ", station phi: " << stationPhi << ", station eta: " << stationEta << ", multiplet: " << multiplet << ", gas gap: " << gas_gap << ", strip charge: " << stg_object -> stripCharges() << ", strip number: " << stg_object -> stripNumbers() << ", channel type: " << channel_type << ", global charge: " << charge);
    }
  */
  
  auto &Vectors = vects[iside][multiplet - 1][gas_gap - 1];
  Vectors.strip_charges_vec = stg_object -> stripCharges();
    
  Vectors.stationEta_perPhi_vec.push_back(stationEta);
  Vectors.strip_numbers_perPhi_vec = stg_object -> stripNumbers();
  Vectors.charge_perPhi_vec.push_back(charge);

  Vectors.charge_vec.push_back(charge);
  Vectors.stationPhi_vec.push_back(stationPhi);
  Vectors.stationEta_vec.push_back(stationEta);

  
  if (channel_type == 0)
    {
      //ATH_MSG_INFO("Pad");
      //ATH_MSG_INFO("Side: " << GeometricSectors::sTGC_Side[iside] << ", station phi: " << stationPhi << ", station eta: " << stationEta << ", multiplet: " << multiplet << ", gas gap: " << gas_gap << ", strip charge: " << stg_object -> stripCharges() << ", strip number: " << stg_object -> stripNumbers() << ", channel type: " << channel_type << ", global charge: " << charge);
      
      auto charge_perLayer_pad_ = Monitored::Collection("charge_pad_" + GeometricSectors::sTGC_Side[iside] + "_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), Vectors.charge_vec);
      auto stationPhi_ = Monitored::Collection("sector_" + GeometricSectors::sTGC_Side[iside] + "_phi_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), Vectors.stationPhi_vec);
      auto stationEta_ = Monitored::Collection("sector_" + GeometricSectors::sTGC_Side[iside] + "_eta_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), Vectors.stationEta_vec);
      fill("sTGC_sideGroup" + GeometricSectors::sTGC_Side[iside], charge_perLayer_pad_, stationPhi_, stationEta_);    
    }

  
  else if (channel_type == 1)
    {      
      //ATH_MSG_INFO("Strip");
      //ATH_MSG_INFO("Side: " << GeometricSectors::sTGC_Side[iside] << ", station phi: " << stationPhi << ", station eta: " << stationEta << ", multiplet: " << multiplet << ", gas gap: " << gas_gap << ", strip charge: " << stg_object -> stripCharges() << ", strip number: " << stg_object -> stripNumbers() << ", channel type: " << channel_type << ", global charge: " << charge);
      
      auto charge_perLayer_strip_ = Monitored::Collection("charge_strip_" + GeometricSectors::sTGC_Side[iside] + "_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), Vectors.charge_vec);
      auto stationPhi_ = Monitored::Collection("sector_" + GeometricSectors::sTGC_Side[iside] + "_phi_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), Vectors.stationPhi_vec);
      auto stationEta_ = Monitored::Collection("sector_" + GeometricSectors::sTGC_Side[iside] + "_eta_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), Vectors.stationEta_vec); 
      auto stationEta_perPhi_ = Monitored::Collection("sector_" + GeometricSectors::sTGC_Side[iside] + "_eta_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap) + "_stationPhi_" + std::to_string(stationPhi), Vectors.stationEta_perPhi_vec);
      auto stripNumber_perLayer_perPhi_strip_ = Monitored::Collection("stripNumber_strip_" + GeometricSectors::sTGC_Side[iside] + "_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap) + "_stationPhi_" + std::to_string(stationPhi), Vectors.strip_numbers_perPhi_vec);
      auto charge_perLayer_perPhi_strip_ = Monitored::Collection("charge_strip_" + GeometricSectors::sTGC_Side[iside] + "_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap) + "_stationPhi_" + std::to_string(stationPhi), Vectors.charge_perPhi_vec);
      fill("sTGC_sideGroup" + GeometricSectors::sTGC_Side[iside], charge_perLayer_strip_, stationPhi_, stationEta_, stationEta_perPhi_, stripNumber_perLayer_perPhi_strip_, charge_perLayer_perPhi_strip_);    
    }
  
  else if (channel_type == 2)
    {
      //ATH_MSG_INFO("Wire");
      //ATH_MSG_INFO("Side: " << GeometricSectors::sTGC_Side[iside]  << ", station phi: " << stationPhi << ", station eta: " << stationEta << ", multiplet: " << multiplet << ", gas gap: " << gas_gap << ", strip charge: " << stg_object -> stripCharges() << ", strip number: " << stg_object -> stripNumbers() << ", channel type: " << channel_type << ", global charge: " << charge);
      
      auto charge_perLayer_wire_ = Monitored::Collection("charge_wire_" + GeometricSectors::sTGC_Side[iside] + "_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), Vectors.charge_vec);
      auto stationPhi_ = Monitored::Collection("sector_" + GeometricSectors::sTGC_Side[iside] + "_phi_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), Vectors.stationPhi_vec);
      auto stationEta_ = Monitored::Collection("sector_" + GeometricSectors::sTGC_Side[iside] + "_eta_multiplet_" + std::to_string(multiplet) + "_gasgap_" + std::to_string(gas_gap), Vectors.stationEta_vec);
      fill("sTGC_sideGroup" + GeometricSectors::sTGC_Side[iside], charge_perLayer_wire_, stationPhi_, stationEta_);    
    }

    
  /*
  for (int side = 0; side < 2; ++side)
    {
      std::string sTGC_sideGroup = "sTGC_sideGroup" + GeometricSectors::sTGC_Side[side];
      
      for (int statPhi = 0; statPhi < 8; ++statPhi)
	{
	  for (int statEta = 0; statEta < 4; ++statEta)
	    {
	      for (int multip = 0; multip < 2; ++multip)
		{
		  for (int gasgap = 0; gasgap < 4; ++gasgap)
		    {
		      auto &Vectors = vects[side][statPhi][statEta][multip][gasgap];
		      auto stationPhi = Monitored::Collection("sector_" + GeometricSectors::sTGC_Side[side] + "_phi_" + std::to_string(statPhi + 1), Vectors.stationPhi_vec);
		      auto stationEta = Monitored::Collection("sector_" + GeometricSectors::sTGC_Side[side] + "_eta_" + std::to_string(statEta + 1), Vectors.stationEta_vec);
		      
		      if (channel_type == 1)
			{
			  auto charge_perLayer_strip = Monitored::Collection("strip_charge_" + GeometricSectors::sTGC_Side[side] + "_phi_" + std::to_string(statPhi + 1) + "_eta_" + std::to_string(statEta + 1) + "_multiplet_" + std::to_string(multip + 1) + "_gasgap_" + std::to_string(gasgap + 1), Vectors.charge_vec);
			  fill(sTGC_sideGroup, stationPhi, stationEta, charge_perLayer_strip);    
			}

		    }
		}
	    }
	}	
    }
  
*/

  /*
  auto stationEta_ASide_eta3    = Monitored::Collection("sector_ASide_eta_3", Vectors.stationEta_ASide_eta3_vec);
  auto stationEta_CSide_eta3    = Monitored::Collection("sector_CSide_eta_3", Vectors.stationEta_CSide_eta3_vec);
  
  fill("sTGCMonitor", stationEta_ASide_eta3, stationEta_CSide_eta3);
  */
  Vectors.strip_charges_vec.clear();
  Vectors.charge_vec.clear();
  Vectors.stationPhi_vec.clear();
  Vectors.stationEta_vec.clear();
  
  Vectors.stationEta_perPhi_vec.clear();
  Vectors.strip_numbers_perPhi_vec.clear();
  Vectors.charge_perPhi_vec.clear();

  
  //return StatusCode::SUCCESS;
  
  //Vectors.stationEta_ASide_eta3_vec.clear();
  //Vectors.stationEta_CSide_eta3_vec.clear();

  /*
  for (int i = 0; i < Vectors.strip_charges_vec.size(); ++i)
    {
      ATH_MSG_INFO("Entries strip charge vector:" << Vectors.strip_charges_vec.at(i));
    }
  */


  

  //auto summary_plots = Monitored::Collection("sector_lb_" + GeometricSectors::sTGC_Side[iside] "_phi" + std::to_string(statPhi + 1), sector_lb_vec);
  
  //ATH_MSG_INFO("Lumiblock:" << lb_mon);

  
  //ATH_MSG_INFO("Lumiblock:" << lb_mon);



/*  
    
  for (int iside=0;iside<2;iside++){
    std::string sTGC_sideGroup = "sTGC_sideGroup"+sTGC_Side[iside];
    
    for (int isector=0;isector<2;isector++){
      for( int statPhi=0; statPhi<8; statPhi++) {
      for( int statEta=0; statEta<2; statEta++) {
        for( int multiplet=0; multiplet<2; multiplet++) {
	    for( int gas_gap=0; gas_gap<4; gas_gap++) {
	          
	          auto& Vectors = vects[iside][isector][statPhi][statEta][multiplet][gas_gap];
                          
		        auto sector_strip=Monitored::Collection("sector_strip_"+sTGC_Side[iside]+"_"+sTGC_Sector[isector]+"_phi"+std::to_string(statPhi+1),Vectors.sector_strip);
			      auto strip_number = Monitored::Collection("strip_number_"+sTGC_Side[iside]+"_"+sTGC_Sector[isector]+"_phi"+std::to_string(statPhi+1), Vectors.strip_number);
			            
			            fill(sTGC_sideGroup,strip_number,sector_strip);

				          auto charge_perLayer = Monitored::Collection("charge_"+sTGC_Side[iside]+"_sector_"+sTGC_Sector[isector]+"_phi"+std::to_string(statPhi+1)+"_stationEta"+EtaSector[statEta]+"_multiplet"+std::to_string(multiplet+1)+"_gas_gap"+std::to_string(gas_gap+1), Vectors.charge);
					        auto strip_times = Monitored::Collection("strip_times_"+sTGC_Side[iside]+"_sector_"+sTGC_Sector[isector]+"_phi"+std::to_string(statPhi+1)+"_stationEta"+EtaSector[statEta]+"_multiplet"+std::to_string(multiplet+1)+"_gas_gap"+std::to_string(gas_gap+1),Vectors.strip_charges);

						      fill(sTGC_sideGroup,charge_perLayer,strip_times);
						            
						          }
							    }
							    }
      }
    }
  }
*/
  
  
    
}




StatusCode StgcRawDataMonAlg::fillsTGCHistograms( const Muon::sTgcPrepData* ) const{

  return StatusCode::SUCCESS;
}


/*
  for (int iside=0;iside<2;iside++){
  std::string sTGC_sideGroup = "sTGC_sideGroup"+sTGC_Side[iside];
    for( int statPhi=0; statPhi<16; statPhi++) {
        for( int statEta=0; statEta<2; statEta++) {
	      for( int multiplet=0; multiplet<2; multiplet++) {
	      for( int gas_gap=0; gas_gap<4; gas_gap++) {
	        auto& vects=summaryPlots_full[iside][statPhi][statEta][multiplet][gas_gap];;
		auto residuals_gap = Monitored::Collection("residuals_"+sTGC_Side[iside]+"_phi"+std::to_string(statPhi+1)+"_stationEta"+EtaSector[statEta]+"_multiplet"+std::to_string(multiplet+1)+"_gas_gap"+std::to_string(gas_gap+1),vects.residuals);
std::string sTGC_GapGroup = "sTGC_GapGroup"+std::to_string(gas_gap+1);

fill(sTGC_sideGroup,residuals_gap);
}
}}}}




for (const Trk::TrackStateOnSurface* trkState: *meTrack->trackStateOnSurfaces()) {
  
if(!(trkState)) continue;
Identifier surfaceId = (trkState)->surface().associatedDetectorElementIdentifier();
if(!m_idHelperSvc->issTgc(surfaceId)) continue;
  
const Amg::Vector3D& pos    = (trkState)->trackParameters()->position();
int stEta= m_idHelperSvc->stgcIdHelper().stationEta(surfaceId);
int multi = m_idHelperSvc->stgcIdHelper().multilayer(surfaceId);
int gap=  m_idHelperSvc->stgcIdHelper().gasGap(surfaceId);

//CSide and ASide                                                                                  
int iside=0;
if(stEta>0) iside=1;

auto& Vectors = summaryPlots[iside][multi-1][gap-1];  

//Filling x-y position vectors using the trackStateonSurface 
Vectors.x_ontrack.push_back(pos.x());
Vectors.y_ontrack.push_back(pos.y());
    
}
} // if meTrack
} // if muon
} //loop on muonContainer

auto& vects=overviewPlots;

auto stationPhi_CSide_eta1_ontrack = Monitored::Collection("stationPhi_CSide_eta1_ontrack",vects.stationPhi_CSide_eta1_ontrack);
auto stationPhi_CSide_eta2_ontrack = Monitored::Collection("stationPhi_CSide_eta2_ontrack",vects.stationPhi_CSide_eta2_ontrack);
auto stationPhi_ASide_eta1_ontrack = Monitored::Collection("stationPhi_ASide_eta1_ontrack",vects.stationPhi_ASide_eta1_ontrack);
auto stationPhi_ASide_eta2_ontrack = Monitored::Collection("stationPhi_ASide_eta2_ontrack",vects.stationPhi_ASide_eta2_ontrack);
auto sector_ASide_eta1_ontrack = Monitored::Collection("sector_ASide_eta1_ontrack",vects.sector_ASide_eta1_ontrack);
auto sector_ASide_eta2_ontrack = Monitored::Collection("sector_ASide_eta2_ontrack",vects.sector_ASide_eta2_ontrack);
auto sector_CSide_eta2_ontrack = Monitored::Collection("sector_CSide_eta2_ontrack",vects.sector_CSide_eta2_ontrack);                                                                                         
auto sector_CSide_eta1_ontrack = Monitored::Collection("sector_CSide_eta1_ontrack",vects.sector_CSide_eta1_ontrack);   

auto lb_ontrack = Monitored::Scalar<int>("lb_ontrack", lb);

auto sector_lb_CSide_eta2_ontrack = Monitored::Collection("sector_lb_CSide_eta2_ontrack",vects.sector_lb_CSide_eta2_ontrack);
auto sector_lb_CSide_eta1_ontrack = Monitored::Collection("sector_lb_CSide_eta1_ontrack",vects.sector_lb_CSide_eta1_ontrack);
auto sector_lb_ASide_eta2_ontrack = Monitored::Collection("sector_lb_ASide_eta2_ontrack",vects.sector_lb_ASide_eta2_ontrack);
auto sector_lb_ASide_eta1_ontrack = Monitored::Collection("sector_lb_ASide_eta1_ontrack",vects.sector_lb_ASide_eta1_ontrack);

fill("sTGCMonitor",stationPhi_CSide_eta1_ontrack,stationPhi_CSide_eta2_ontrack,stationPhi_ASide_eta1_ontrack,stationPhi_ASide_eta2_ontrack,sector_CSide_eta1_ontrack,sector_CSide_eta2_ontrack,sector_ASide_eta1_ontrack,sector_ASide_eta2_ontrack,sector_lb_CSide_eta2_ontrack,sector_lb_CSide_eta1_ontrack,sector_lb_ASide_eta2_ontrack,sector_lb_ASide_eta1_ontrack,lb_ontrack);
    
for (int iside=0;iside<2;iside++){
std::string sTGC_sideGroup = "sTGC_sideGroup"+sTGC_Side[iside];
for( int multiplet=0; multiplet<2; multiplet++) {
for( int gas_gap=0; gas_gap<4; gas_gap++) {
      
auto& Vectors = summaryPlots[iside][multiplet][gas_gap];

auto x_ontrack = Monitored::Collection("x_"+sTGC_Side[iside]+"_multiplet"+std::to_string(multiplet+1)+"_gas_gap_"+std::to_string(gas_gap+1)+"_ontrack", Vectors.x_ontrack);
auto y_ontrack = Monitored::Collection("y_"+sTGC_Side[iside]+"_multiplet"+std::to_string(multiplet+1)+"_gas_gap_"+std::to_string(gas_gap+1)+"_ontrack", Vectors.y_ontrack);
      
fill(sTGC_sideGroup,x_ontrack,y_ontrack);
}          
}
}
*/           
  








/*


void StgcRawDataMonAlg::clusterFromTrack(const xAOD::TrackParticleContainer *muonContainer, int lb, const Muon::sTgcPrepData *stg_object) const
{
  std::vector<float> x_ontrack_vec;
  std::vector<float> y_ontrack_vec;

  Amg::Vector3D pos = stg_object -> globalPosition();

  //ATH_MSG_INFO("pos.x() = " << pos.x());
  //ATH_MSG_INFO("pos.y() = " << pos.y());

  x_ontrack_vec.push_back(pos.x());
  y_ontrack_vec.push_back(pos.y());

  //ATH_MSG_INFO("pos.x() vec = " << x_ontrack_vec.size());
  //ATH_MSG_INFO("pos.y() vec = " << y_ontrack_vec.size());

  for(int iside = 0; iside < 2; ++iside) 
    {
      for(int multiplet = 0; multiplet < 2; ++multiplet) 
	{
	  for(int gas_gap = 0; gas_gap < 4; ++gas_gap) 
	    {
	      std::string sTGC_sideGroup = "sTGC_sideGroup" + GeometricSectors::sTGC_Side[iside];

	      auto x_ontrack = Monitored::Collection("x_" + GeometricSectors::sTGC_Side[iside] + "_multiplet" + std::to_string(multiplet + 1) + "_gas_gap_" + std::to_string(gas_gap + 1) + "_ontrack", x_ontrack_vec);
	      auto y_ontrack = Monitored::Collection("y_" + GeometricSectors::sTGC_Side[iside] + "_multiplet" + std::to_string(multiplet + 1) + "_gas_gap_" + std::to_string(gas_gap + 1) + "_ontrack", y_ontrack_vec);
	      fill(sTGC_sideGroup, x_ontrack, y_ontrack);
	    }
	}
    }
  
  x_ontrack_vec.clear();
  y_ontrack_vec.clear();
}
*/

/*
void StgcRawDataMonAlg::clusterFromTrack(const xAOD::TrackParticleContainer*  muonContainer, int lb) const{

  //sTGCSummaryHistogramStruct summaryPlots[2][2][4];
  //sTGCSummaryHistogramStruct summaryPlots_full[2][16][2][2][4];
  //sTGCOverviewHistogramStruct overviewPlots;

  int nmu=0;

  for (const xAOD::TrackParticle* meTP  : *muonContainer){

    if (meTP) {
      nmu++;
      auto eta_trk = Monitored::Scalar<float>("eta_trk", meTP->eta());
      auto phi_trk = Monitored::Scalar<float>("phi_trk", meTP->phi());

      // retrieve the original track                                                       
      const Trk::Track* meTrack = meTP->track();
      if (meTrack) {
        // get the vector of measurements on track                                                           

	const DataVector<const Trk::MeasurementBase>* meas = meTrack->measurementsOnTrack();

	
	for(const Trk::MeasurementBase* it: *meas){
	  
	  const Trk::RIO_OnTrack* rot = dynamic_cast<const Trk::RIO_OnTrack*>(it);
	      
	    if (rot) {
	        Identifier rot_id = rot->identify();
		    if (m_idHelperSvc->issTgc(rot_id)) {
		          const Muon::sTgcClusterOnTrack* cluster = dynamic_cast<const Muon::sTgcClusterOnTrack*>(rot);
			              
			        if (cluster) {
				
				std::string stName   = m_idHelperSvc->stgcIdHelper().stationNameString(m_idHelperSvc->stgcIdHelper().stationName(rot_id));
				int stNumber    = m_idHelperSvc->stgcIdHelper().stationName(rot_id);
				int stEta= m_idHelperSvc->stgcIdHelper().stationEta(rot_id);
				int stPhi= m_idHelperSvc->stgcIdHelper().stationPhi(rot_id);
				int multi = m_idHelperSvc->stgcIdHelper().multilayer(rot_id);
				int gap=  m_idHelperSvc->stgcIdHelper().gasGap(rot_id);
				int ch=  m_idHelperSvc->stgcIdHelper().channel(rot_id);
                        
				//MMS and MML phi sectors                                                                    
				int phisec=0;
				if (stNumber%2 == 0) phisec=1;
				
				int sectorPhi=get_sectorPhi_from_stationPhi_stName(stPhi,stName);
				
				int PCB=get_PCB_from_channel(ch);
				
				
				auto& vects=overviewPlots;
				
				//Occupancy plots with PCB granularity further divided for each eta sector: -2, -1, 1, 2  
//Filling Vectors for stationEta=-1 - cluster on track
if (stEta==-1){
vects.stationPhi_CSide_eta1_ontrack.push_back(sectorPhi);
vects.sector_CSide_eta1_ontrack.push_back(get_bin_for_occ_CSide_pcb_eta1_hist(stEta,multi,gap,PCB));
vects.sector_lb_CSide_eta1_ontrack.push_back(get_bin_for_occ_lb_CSide_pcb_eta1_hist(stEta,multi,gap,PCB,phisec));
}
//Filling Vectors for stationEta=-2 - cluster on track 
 else if (stEta==-2){
vects.stationPhi_CSide_eta2_ontrack.push_back(sectorPhi);
vects.sector_CSide_eta2_ontrack.push_back(get_bin_for_occ_CSide_pcb_eta2_hist(stEta,multi,gap,PCB));
vects.sector_lb_CSide_eta2_ontrack.push_back(get_bin_for_occ_lb_CSide_pcb_eta2_hist(stEta,multi,gap,PCB,phisec));
}
//Filling Vectors for stationEta=1 - cluster on track 
 else if (stEta==1){
vects.stationPhi_ASide_eta1_ontrack.push_back(sectorPhi);
vects.sector_ASide_eta1_ontrack.push_back(get_bin_for_occ_ASide_pcb_eta1_hist(stEta,multi,gap,PCB));
vects.sector_lb_ASide_eta1_ontrack.push_back(get_bin_for_occ_lb_ASide_pcb_eta1_hist(stEta,multi,gap,PCB,phisec));
}
//Filling Vectors for stationEta=2 - cluster on track 
 else {
vects.stationPhi_ASide_eta2_ontrack.push_back(sectorPhi);
vects.sector_ASide_eta2_ontrack.push_back(get_bin_for_occ_ASide_pcb_eta2_hist(stEta,multi,gap,PCB));
vects.sector_lb_ASide_eta2_ontrack.push_back(get_bin_for_occ_lb_ASide_pcb_eta2_hist(stEta,multi,gap,PCB,phisec));
    
}


float x =cluster->localParameters()[Trk::loc1] ;

for (const Trk::TrackStateOnSurface* trkState: *meTrack->trackStateOnSurfaces()) {

if(!(trkState)) continue;
Identifier surfaceId = (trkState)->surface().associatedDetectorElementIdentifier();
if(!m_idHelperSvc->issTgc(surfaceId)) continue;

int trk_stEta= m_idHelperSvc->stgcIdHelper().stationEta(surfaceId);
int trk_stPhi= m_idHelperSvc->stgcIdHelper().stationPhi(surfaceId);
int trk_multi = m_idHelperSvc->stgcIdHelper().multilayer(surfaceId);
int trk_gap=  m_idHelperSvc->stgcIdHelper().gasGap(surfaceId);
if(  trk_stPhi==stPhi  and trk_stEta==stEta and trk_multi==multi and trk_gap==gap ){
double x_trk = trkState->trackParameters()->parameters()[Trk::loc1];
double y_trk = trkState->trackParameters()->parameters()[Trk::locY];
        
int stPhi16=0;
//int stPhi16 = (stName == "L") ?  (2*stPhi - 1) : (2*stPhi);
if (stName=="L")   stPhi16=2*stPhi-1;
        
if (stName=="S")  stPhi16=2*stPhi;

int iside=0;
if(stEta>0) iside=1;

float stereo_angle=      0.02618;
//    float stereo_correction=stereo_angle*y_trk;
float stereo_correction=sin(stereo_angle)*y_trk;
if(multi==1 && gap<3) stereo_correction=0;
if(multi==2 && gap>2) stereo_correction=0;
if(multi==1 && gap==3 ) stereo_correction*=-1;
if(multi==2 && gap==1 ) stereo_correction*=-1;
if(multi==1 && gap<3) stereo_angle=0;
if(multi==2 && gap>2) stereo_angle=0;
float res_stereo = (x - x_trk)*cos(stereo_angle) - stereo_correction;
auto residual_mon = Monitored::Scalar<float>("residual", res_stereo);
auto stPhi_mon = Monitored::Scalar<float>("stPhi_mon",stPhi16);
fill("sTGCMonitor",residual_mon, eta_trk, phi_trk, stPhi_mon);
//int abs_stEta= get_sectorEta_from_stationEta(stEta);
//auto& vectors = summaryPlots_full[iside][stPhi16-1][abs_stEta-1][multi-1][gap-1];
//vectors.residuals.push_back(res_stereo);
}
}
} //if cluster
} //isMM
} // if rot
} // loop on meas
}
}
}
}
*/
/*
  for (int iside=0;iside<2;iside++){
    std::string sTGC_sideGroup = "sTGC_sideGroup"+sTGC_Side[iside];
      for( int statPhi=0; statPhi<16; statPhi++) {
          for( int statEta=0; statEta<2; statEta++) {
	        for( int multiplet=0; multiplet<2; multiplet++) {
		for( int gas_gap=0; gas_gap<4; gas_gap++) {
		  auto& vects=summaryPlots_full[iside][statPhi][statEta][multiplet][gas_gap];;
		    auto residuals_gap = Monitored::Collection("residuals_"+sTGC_Side[iside]+"_phi"+std::to_string(statPhi+1)+"_stationEta"+EtaSector[statEta]+"_multiplet"+std::to_string(multiplet+1)+"_gas_gap"+std::to_string(gas_gap+1),vects.residuals);
		      std::string sTGC_GapGroup = "sTGC_GapGroup"+std::to_string(gas_gap+1);

		        fill(sTGC_sideGroup,residuals_gap);
			}
			      }}}}




			      for (const Trk::TrackStateOnSurface* trkState: *meTrack->trackStateOnSurfaces()) {
			          
			        if(!(trkState)) continue;
				  Identifier surfaceId = (trkState)->surface().associatedDetectorElementIdentifier();
				    if(!m_idHelperSvc->issTgc(surfaceId)) continue;
				        
				      const Amg::Vector3D& pos    = (trkState)->trackParameters()->position();
				        int stEta= m_idHelperSvc->stgcIdHelper().stationEta(surfaceId);
					  int multi = m_idHelperSvc->stgcIdHelper().multilayer(surfaceId);
					    int gap=  m_idHelperSvc->stgcIdHelper().gasGap(surfaceId);

					      //CSide and ASide                                                                                  
					        int iside=0;
						  if(stEta>0) iside=1;

						    auto& Vectors = summaryPlots[iside][multi-1][gap-1];  

						      //Filling x-y position vectors using the trackStateonSurface 
						        Vectors.x_ontrack.push_back(pos.x());
							  Vectors.y_ontrack.push_back(pos.y());
							        
							  }
      } // if meTrack
    } // if muon
  } //loop on muonContainer

*/

/*
  auto& vects=overviewPlots;

  auto stationPhi_CSide_eta1_ontrack = Monitored::Collection("stationPhi_CSide_eta1_ontrack",vects.stationPhi_CSide_eta1_ontrack);
  auto stationPhi_CSide_eta2_ontrack = Monitored::Collection("stationPhi_CSide_eta2_ontrack",vects.stationPhi_CSide_eta2_ontrack);
  auto stationPhi_ASide_eta1_ontrack = Monitored::Collection("stationPhi_ASide_eta1_ontrack",vects.stationPhi_ASide_eta1_ontrack);
  auto stationPhi_ASide_eta2_ontrack = Monitored::Collection("stationPhi_ASide_eta2_ontrack",vects.stationPhi_ASide_eta2_ontrack);
  auto sector_ASide_eta1_ontrack = Monitored::Collection("sector_ASide_eta1_ontrack",vects.sector_ASide_eta1_ontrack);
  auto sector_ASide_eta2_ontrack = Monitored::Collection("sector_ASide_eta2_ontrack",vects.sector_ASide_eta2_ontrack);
  auto sector_CSide_eta2_ontrack = Monitored::Collection("sector_CSide_eta2_ontrack",vects.sector_CSide_eta2_ontrack);                                                                                         
  auto sector_CSide_eta1_ontrack = Monitored::Collection("sector_CSide_eta1_ontrack",vects.sector_CSide_eta1_ontrack);   

  auto lb_ontrack = Monitored::Scalar<int>("lb_ontrack", lb);

  auto sector_lb_CSide_eta2_ontrack = Monitored::Collection("sector_lb_CSide_eta2_ontrack",vects.sector_lb_CSide_eta2_ontrack);
  auto sector_lb_CSide_eta1_ontrack = Monitored::Collection("sector_lb_CSide_eta1_ontrack",vects.sector_lb_CSide_eta1_ontrack);
  auto sector_lb_ASide_eta2_ontrack = Monitored::Collection("sector_lb_ASide_eta2_ontrack",vects.sector_lb_ASide_eta2_ontrack);
  auto sector_lb_ASide_eta1_ontrack = Monitored::Collection("sector_lb_ASide_eta1_ontrack",vects.sector_lb_ASide_eta1_ontrack);

  fill("sTGCMonitor",stationPhi_CSide_eta1_ontrack,stationPhi_CSide_eta2_ontrack,stationPhi_ASide_eta1_ontrack,stationPhi_ASide_eta2_ontrack,sector_CSide_eta1_ontrack,sector_CSide_eta2_ontrack,sector_ASide_eta1_ontrack,sector_ASide_eta2_ontrack,sector_lb_CSide_eta2_ontrack,sector_lb_CSide_eta1_ontrack,sector_lb_ASide_eta2_ontrack,sector_lb_ASide_eta1_ontrack,lb_ontrack);
    
  for (int iside=0;iside<2;iside++){
    std::string sTGC_sideGroup = "sTGC_sideGroup"+sTGC_Side[iside];
    for( int multiplet=0; multiplet<2; multiplet++) {
      for( int gas_gap=0; gas_gap<4; gas_gap++) {
      
auto& Vectors = summaryPlots[iside][multiplet][gas_gap];

auto x_ontrack = Monitored::Collection("x_"+sTGC_Side[iside]+"_multiplet"+std::to_string(multiplet+1)+"_gas_gap_"+std::to_string(gas_gap+1)+"_ontrack", Vectors.x_ontrack);
auto y_ontrack = Monitored::Collection("y_"+sTGC_Side[iside]+"_multiplet"+std::to_string(multiplet+1)+"_gas_gap_"+std::to_string(gas_gap+1)+"_ontrack", Vectors.y_ontrack);
      
fill(sTGC_sideGroup,x_ontrack,y_ontrack);
      }          
    }
  }
*/           
 
 
/*
void StgcRawDataMonAlg::clusterFromTrack(const xAOD::TrackParticleContainer*  muonContainer, int lb) const
{
  
  for(const xAOD::TrackParticle* meTP : *muonContainer) {
    
    if(!meTP) continue;

    auto eta_trk = Monitored::Scalar<float>("eta_trk", meTP->eta());
    auto phi_trk = Monitored::Scalar<float>("phi_trk", meTP->phi());

    //retrieve the original track
    const Trk::Track* meTrack = meTP->track();
    if(!meTrack) continue;
    // get the vector of measurements on track
    const DataVector<const Trk::MeasurementBase>* meas = meTrack->measurementsOnTrack();

    for(const Trk::MeasurementBase* it : *meas) {
      const Trk::RIO_OnTrack* rot = dynamic_cast<const Trk::RIO_OnTrack*>(it);
      if(!rot) continue;
      Identifier rot_id = rot->identify();
      if(!m_idHelperSvc->isMM(rot_id)) continue;
      const Muon::sTgcClusterOnTrack* cluster = dynamic_cast<const Muon::sTgcClusterOnTrack*>(rot);
      if(!cluster) continue;

      std::string stName = m_idHelperSvc->stgcIdHelper().stationNameString(m_idHelperSvc->mmIdHelper().stationName(rot_id));
      int stEta          = m_idHelperSvc->stgcIdHelper().stationEta(rot_id);
      int stPhi          = m_idHelperSvc->stgcIdHelper().stationPhi(rot_id);
      int multi          = m_idHelperSvc->stgcIdHelper().multilayer(rot_id);
      int gap            = m_idHelperSvc->stgcIdHelper().gasGap(rot_id);
      int ch             = m_idHelperSvc->stgcIdHelper().channel(rot_id);

      // MMS and MML phi sectors
      //int phisec = (stNumber%2==0) ? 1 : 0;
      int sectorPhi = get_sectorPhi_from_stationPhi_stName(stPhi,stName); // 1->16
      int PCB = get_PCB_from_channel(ch);
      int iside = (stEta > 0) ? 1 : 0;

      //auto& vects = overviewPlots;
      //auto& thisSect = occupancyPlots[sectorPhi-1][iside];
      //auto& vect = sumPlots[iside][sectorPhi-1][std::abs(stEta)-1][multi-1][gap-1];

      const Muon::sTgcPrepData* prd = cluster->prepRawData();
      const std::vector<Identifier>& stripIds = prd->rdoList();
      unsigned int csize = stripIds.size();
      const std::vector<uint16_t>& stripNumbers = prd->stripNumbers();
      float charge = prd->charge();
      std::vector<short int> s_times = prd->stripTimes();

      //cts.numberofstrips_percluster.push_back(csize);
      //cts.charge_all.push_back(charge);

      //float c_time = 0;
      //for(unsigned int sIdx=0; sIdx<stripIds.size(); ++sIdx)
      //{
      // vects.strp_times.push_back(s_times.at(sIdx));
      // vect.strip_number.push_back(stripNumbers[sIdx]);
      //  vect.strp_times.push_back(s_times.at(sIdx));
      //  c_time += s_times.at(sIdx);
      //}
      //c_time /= s_times.size();
      //vects.cl_times.push_back(c_time);
      //vect.cl_times.push_back(c_time);
      //vect.charge.push_back(charge);


      
      const int pcb_counter = 5; // index for the PCBs from 1 to 8 as done globally for the two detector components (abs(eta)=1 and abs(eta)=2)
      int PCBeta12 = (std::abs(stEta) == 2) ? (PCB + pcb_counter) : PCB;
      thisSect.sector_lb_ontrack.push_back(get_bin_for_occ_lb_pcb_hist(multi,gap,PCBeta12));
      
      // Occupancy plots with PCB granularity further divided for each eta sector: -2, -1, 1, 2
      // Filling Vectors for stationEta=-1 - cluster on track
      if(stEta==-1) {
      vects.stationPhi_CSide_eta1_ontrack.push_back(sectorPhi);
      vects.sector_CSide_eta1_ontrack.push_back(get_bin_for_occ_CSide_pcb_eta1_hist(stEta, multi, gap, PCB));
      }
      // Filling Vectors for stationEta=-2 - cluster on track
      else if(stEta==-2) {
      vects.stationPhi_CSide_eta2_ontrack.push_back(sectorPhi);
      vects.sector_CSide_eta2_ontrack.push_back(get_bin_for_occ_CSide_pcb_eta2_hist(stEta,multi,gap,PCB));
      }
      // Filling Vectors for stationEta=1 - cluster on track
      else if(stEta==1) {
      vects.stationPhi_ASide_eta1_ontrack.push_back(sectorPhi);
      vects.sector_ASide_eta1_ontrack.push_back(get_bin_for_occ_ASide_pcb_eta1_hist(stEta,multi,gap,PCB));
    }
      // Filling Vectors for stationEta=2 - cluster on track
      else {
      vects.stationPhi_ASide_eta2_ontrack.push_back(sectorPhi);
      vects.sector_ASide_eta2_ontrack.push_back(get_bin_for_occ_ASide_pcb_eta2_hist(stEta,multi,gap,PCB));
    }
      
      
      float x = cluster->localParameters()[Trk::loc1];
      for(const Trk::TrackStateOnSurface* trkState : *meTrack->trackStateOnSurfaces()) {
	
	if(!(trkState)) continue;
	Identifier surfaceId = (trkState)->surface().associatedDetectorElementIdentifier();
	if(!m_idHelperSvc->issTgc(surfaceId)) continue;
	int trk_stEta = m_idHelperSvc->stgcIdHelper().stationEta(surfaceId);
	int trk_stPhi = m_idHelperSvc->stgcIdHelper().stationPhi(surfaceId);
	int trk_multi = m_idHelperSvc->stgcIdHelper().multilayer(surfaceId);
	int trk_gap   = m_idHelperSvc->stgcIdHelper().gasGap(surfaceId);

	if( (trk_stPhi == stPhi) && (trk_stEta == stEta) && (trk_multi == multi) && (trk_gap == gap)) {
	  double x_trk = trkState->trackParameters()->parameters()[Trk::loc1];
	  int sectorPhi = get_sectorPhi_from_stationPhi_stName(trk_stPhi,stName); // 1->16
	  int side = (stEta > 0) ? 1 : 0;
	  float res_stereo = (x - x_trk);
	  if(m_do_stereoCorrection) {
	    float stereo_angle = ((multi == 1 && gap < 3) || (multi == 2 && gap > 2)) ? 0 : 0.02618;
	    double y_trk = trkState->trackParameters()->parameters()[Trk::locY];
	    float stereo_correction = ( (multi == 1 && gap < 3) || (multi == 2 && gap > 2) ) ? 0 : ( ((multi == 1 && gap == 3) || (multi == 2 && gap ==1 )) ? (-std::sin(stereo_angle)*y_trk) : std::sin(stereo_angle)*y_trk );
	    res_stereo = (x - x_trk)*std::cos(stereo_angle) - stereo_correction;
	  }
	  auto residual_mon = Monitored::Scalar<float>("residual", res_stereo);
	  auto stPhi_mon = Monitored::Scalar<float>("stPhi_mon",sectorPhi);
	  fill("sTGCMonitor", residual_mon, eta_trk, phi_trk, stPhi_mon);
	  //int abs_stEta = get_sectorEta_from_stationEta(stEta); // 0 or 1
	  //auto& vectors = summaryPlots_full[side][sectorPhi-1][abs_stEta][multi-1][gap-1];
	  //vectors.residuals.push_back(res_stereo);
	}
      }//TrackStates
    }
  } // loop on meas
}

*/
