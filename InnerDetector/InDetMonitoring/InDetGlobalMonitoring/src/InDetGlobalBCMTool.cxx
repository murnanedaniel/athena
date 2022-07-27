/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

/** @file InDetGlobalBCMTool.cxx
 * Implementation of inner detector beam conditions monitoring tool
 *
 * @author Heidi Sandaker <Heidi.Sandaker@cern.ch> @n
 * Alex Kastanas <Alex.Kastanas@cern.ch> @n
 * Jennifer Jentzsch <Jennifer.Jentzsch@cern.ch>      
 *
 *
 *********************************************************************************/

//Local
#include "InDetGlobalBCMTool.h"
//Framework
#include "GaudiKernel/IInterface.h"
#include "GaudiKernel/StatusCode.h"
#include "StoreGate/ReadHandle.h"
//Root
#include "TH1F.h"
#include "TH2F.h"

//Primary Vertex
#include "TrkEventPrimitives/ParamDefs.h"
#include "VxVertex/VxCandidate.h"
#include "VxVertex/VxTrackAtVertex.h"
#include "InDetGlobalPrimaryVertexMonTool.h"
//Pixel and SCT stuff
#include "InDetRawData/PixelRDORawData.h"

//Standard c++
#include <string>
#include <list>
#include <vector>


#define BCM_LVL1A 18

//Some variables
static const unsigned int bcid_start 				= 0; //to choose the first bcid displayed
static const unsigned int bcids_displayed 			= 3570; //to choose the number of bcids displayed when bcid as x-axis
static const signed int ecr_start                             = 0;    //to chose the first ECR displayed
static const unsigned int ecrs_displayed                      = 256;  //to chose the number of ecrs displayed when ECR as X-axis
static const unsigned int bc_readout			= 31;   // length of read out per L1A, at the moment 31 BCs re read out
static const unsigned int lb_start 				= 0;
static const unsigned int lb_max				= 5000;

// the declare property can be used to make variables accessible through python
InDetGlobalBCMTool::InDetGlobalBCMTool(
				       const std::string & type, 
				       const std::string & name,
				       const IInterface* parent)
  : InDetGlobalMotherMonTool(type, name, parent),
    m_detector("ID"),
    m_pulse_width_all{nullptr,nullptr},
    m_hits_lvl1a{nullptr,nullptr},
    m_hits_bcid{nullptr, nullptr},
    m_hitdistribution{nullptr, nullptr},
    m_pulse_position_gen{nullptr, nullptr},
    m_pulse_width_gen{nullptr, nullptr},
    m_strange_signals{nullptr, nullptr},
    m_highocc(nullptr),
    m_deltat_vs_hits{nullptr,nullptr},
    m_deltat_aligned_vs_hits{nullptr,nullptr},
    m_deltat_vs_bcid{nullptr,nullptr},
    m_deltat_aligned_vs_bcid{nullptr,nullptr},
    m_deltat_vs_lb{nullptr,nullptr},
    m_deltat_aligned_vs_lb{nullptr,nullptr},
    m_deltat_vs_ecr{nullptr,nullptr},
    m_deltat_aligned_vs_ecr{nullptr,nullptr},
    m_deltat_vs_PrVertex{nullptr,nullptr},
    m_deltat_vs_pixhits{nullptr,nullptr},
    m_deltat_vs_pixhitsEC{nullptr,nullptr},
    m_deltat_vs_pixhitsBR{nullptr,nullptr},
    m_deltat_vs_scthits{nullptr,nullptr},
    m_deltat_vs_scthitsBR{nullptr,nullptr},
    m_sct_vs_pix_col{nullptr,nullptr},
    m_sct_vs_pix_bkg{nullptr,nullptr}
{
  declareProperty("Detector", m_detector); 
  declareProperty("histFolder",m_histFolder= "InDetGlobal/PrimaryVertex");
}

StatusCode InDetGlobalBCMTool::initialize()
{
    ATH_CHECK(m_vxContainerName.initialize());
    ATH_CHECK(m_eventInfoKey.initialize());

    return StatusCode::SUCCESS;
}

//---------------------------------------------------------
StatusCode InDetGlobalBCMTool::bookHistogramsRecurrent()
{
    MonGroup monGr_shift ( this, "InDetGlobal/BCM", run, ATTRIB_UNMANAGED);
    MonGroup monGr_GlobalHistograms ( this, "InDetGlobal/BCM/GlobalHistograms", run, ATTRIB_UNMANAGED);
    MonGroup monGr_SideHistograms ( this, "InDetGlobal/BCM/SideHistograms", run, ATTRIB_UNMANAGED);
    MonGroup monGr_DetectorHistograms ( this, "InDetGlobal/BCM/DetectorHistograms", run, ATTRIB_UNMANAGED);

    bool status = true;

    if (newRunFlag()){
    // Example of plot registration per new run
    status &= registerHist(monGr_shift,m_nExamplePlot = new TH1F("m_nExample","Example plot BCM",5,0,5));

    std::string name,title,station_name,gain,pulse;
    
    /*************************************************************************************************************
     * Register Monitoring Histograms
     *************************************************************************************************************/

    name = "HitsVsLvl1AAll";
    title = "Hits vs LVL1A All";
    status &= registerHist(monGr_GlobalHistograms,m_hits_lvl1a_all = new TH1F(name.c_str(),title.c_str(),32,0,32));
    m_hits_lvl1a_all->GetXaxis()->SetTitle("lvl1a [25 ns]");
    m_hits_lvl1a_all->GetYaxis()->SetTitle("# of hits");

    name = "ChannelVsLvl1a";
    title = "Channel vs LVL1 A";
    status &= registerHist(monGr_GlobalHistograms,m_ChannelVsLvl1a = new TH2F(name.c_str(), title.c_str(), 64, 0, 64, 16, 0, 16));
    m_ChannelVsLvl1a->GetXaxis()->SetTitle("LVL1 A [25 ns]");
    m_ChannelVsLvl1a->GetYaxis()->SetTitle("Channel #");
    
    name = "ChannelVsBCID";
    title = "Channel vs BCID";
    status &= registerHist(monGr_GlobalHistograms,m_ChannelVsBCID = new TH2F(name.c_str(),title.c_str(), bcids_displayed,bcid_start, bcid_start+bcids_displayed, 16, 0,16));
    m_ChannelVsBCID->GetXaxis()->SetTitle("BCID [25 ns]");
    m_ChannelVsBCID->GetYaxis()->SetTitle("Channel #");
    
    name = "ChannelVsECR";
    title = "Channel vs ECR";
    status &= registerHist(monGr_GlobalHistograms,m_ChannelVsECR = new TH2F(name.c_str(), title.c_str(), 256, 0, 256, 16, 0, 16));
    m_ChannelVsECR->GetXaxis()->SetTitle("ECR [5 s]");
    m_ChannelVsECR->GetYaxis()->SetTitle("Channel #");

    name = "NumberOfEvents";
    title = "Number of monitored events";
    status &= registerHist(monGr_GlobalHistograms,m_NumberOfEvents = new TH1F(name.c_str(), title.c_str(), 3,0,3));
    m_NumberOfEvents->GetXaxis()->SetTitle("sourceID [0=BCM, 1=LowHorizontal, 2=LowVertical]");
    m_NumberOfEvents->GetYaxis()->SetTitle("# of events");    

    name = "AbortFraction";
    title = "Abort Fraction";
    status &= registerHist(monGr_GlobalHistograms,m_AbortFraction = new TH1F(name.c_str(), title.c_str(), 101, 0, 101));
    m_AbortFraction->GetXaxis()->SetTitle("Abort Fraction %");
    m_AbortFraction->GetYaxis()->SetTitle("# of Hits");

    name = "AbortFractionROD0";
    title = "Abort Fraction ROD0";
    status &= registerHist(monGr_GlobalHistograms,m_AbortFractionROD0 = new TH1F(name.c_str(), title.c_str(), 101, 0, 101));
    m_AbortFractionROD0->GetXaxis()->SetTitle("Abort Fraction %");
    m_AbortFractionROD0->GetYaxis()->SetTitle("# of Hits");

    name = "AbortFractionROD1";
    title = "Abort Fraction ROD1";
    status &= registerHist(monGr_GlobalHistograms,m_AbortFractionROD1 = new TH1F(name.c_str(), title.c_str(), 101, 0, 101));
    m_AbortFractionROD1->GetXaxis()->SetTitle("Abort Fraction %");
    m_AbortFractionROD1->GetYaxis()->SetTitle("# of Hits");

    name = "AbortFractionVsBCID";
    title = "Abort Fraction Vs BCID";
    status &= registerHist(monGr_GlobalHistograms,m_AbortFractionVsBCID = new TH1F(name.c_str(), title.c_str(), 3563, 0, 3563));
    m_AbortFractionVsBCID->GetXaxis()->SetTitle("BCID [25 ns]");
    m_AbortFractionVsBCID->GetYaxis()->SetTitle("Abort Fraction %");

    name = "AbortFractionROD0VsBCID";
    title = "Abort Fraction ROD0 Vs BCID";
    status &= registerHist(monGr_GlobalHistograms,m_AbortFractionROD0VsBCID = new TH1F(name.c_str(), title.c_str(), 3563, 0, 3563));
    m_AbortFractionROD0VsBCID->GetXaxis()->SetTitle("BCID [25 ns]");
    m_AbortFractionROD0VsBCID->GetYaxis()->SetTitle("Abort Fraction %");

    name = "AbortFractionROD1VsBCID";
    title = "Abort Fraction ROD1 Vs BCID";
    status &= registerHist(monGr_GlobalHistograms,m_AbortFractionROD1VsBCID = new TH1F(name.c_str(), title.c_str(), 3563, 0, 3563));
    m_AbortFractionROD1VsBCID->GetXaxis()->SetTitle("BCID [25 ns]");
    m_AbortFractionROD1VsBCID->GetYaxis()->SetTitle("Abort Fraction %");

    name = "AbortFractionVsECR";
    title = "Abort Fraction Vs ECR";
    status &= registerHist(monGr_GlobalHistograms,m_AbortFractionVsECR = new TH1F(name.c_str(), title.c_str(), 256, 0, 256));
    m_AbortFractionVsECR->GetXaxis()->SetTitle("ECR [5 s]");
    m_AbortFractionVsECR->GetYaxis()->SetTitle("Abort Fraction %");

    name = "AbortFractionVsLB";
    title = "Abort Fraction Vs LB";
    status &= registerHist(monGr_GlobalHistograms,m_AbortFractionVsLB = new TH2F(name.c_str(), title.c_str(), lb_max, lb_start, lb_max, 100, 0, 1));
    m_AbortFractionVsLB->GetXaxis()->SetTitle("LB");
    m_AbortFractionVsLB->GetYaxis()->SetTitle("Abort Fraction %");

    name = "AbortFractionROD0VsECR";
    title = "Abort Fraction ROD0 Vs ECR";
    status &= registerHist(monGr_GlobalHistograms,m_AbortFractionROD0VsECR = new TH1F(name.c_str(), title.c_str(), 256, 0, 256));
    m_AbortFractionROD0VsECR->GetXaxis()->SetTitle("ECR [5 s]");
    m_AbortFractionROD0VsECR->GetYaxis()->SetTitle("Abort Fraction %");

    name = "AbortFractionROD1VsECR";
    title = "Abort Fraction ROD1 Vs ECR";
    status &= registerHist(monGr_GlobalHistograms,m_AbortFractionROD1VsECR = new TH1F(name.c_str(), title.c_str(), 256, 0, 256));
    m_AbortFractionROD1VsECR->GetXaxis()->SetTitle("BCID [5 s]");
    m_AbortFractionROD1VsECR->GetYaxis()->SetTitle("Abort Fraction %");

    name = "timewalkAll";
    title = "timewalk All";
    status &= registerHist(monGr_GlobalHistograms,m_timewalk_all = new TH2F(name.c_str(),title.c_str(),64, 0, 64, 32, 0, 32));
    m_timewalk_all->GetXaxis()->SetTitle("pulse pos in time bins [25/64 ns]");
    m_timewalk_all->GetYaxis()->SetTitle("pulse width in time bins [25/64 ns]");

    name = "NumberHighOccupancyEventsVsLB";
    title = "Number of High Occupancy Events Vs LB";
    status &= registerHist(monGr_GlobalHistograms, m_highocc = new TH1F(name.c_str(),title.c_str(), lb_max, lb_start, lb_max));
    m_highocc->GetXaxis()->SetTitle("LB number");
    m_highocc->GetYaxis()->SetTitle("Number of High Occupancy Events");


    // station_nr=0 for A-side, 1 for C-side
    // gain_value=0 if low , 1 if high
    // pulse_nr=0 for pulse1, 1 for pulse2
    
    for (unsigned int gain_value=0;gain_value<2;gain_value++)
      {
	if (gain_value==0){gain = "Low";} 
	else gain="High";   

	name = "PulsePositionAll";
	title = "Pulse Position All";
	status &= registerHist(monGr_GlobalHistograms,m_pulse_position_all[gain_value] = new TH1F(name.c_str(),title.c_str(),64, 0, 64));
	m_pulse_position_all[gain_value]->GetXaxis()->SetTitle("time bins [25/64 ns]");
	m_pulse_position_all[gain_value]->GetYaxis()->SetTitle("# of hits");

	name = "PulseWidthAll";
	title = "Pulse Width All";
	status &= registerHist(monGr_GlobalHistograms,m_pulse_width_all[gain_value] = new TH1F(name.c_str(),title.c_str(),32, 0, 32));
	m_pulse_width_all[gain_value]->GetXaxis()->SetTitle("time bins [25/64 ns]");
	m_pulse_width_all[gain_value]->GetYaxis()->SetTitle("# of hits");
	
	name = "HitDistribution" + gain + "Gain";
	title = "Hits vs Channel " + gain + " Gain";
	status &= registerHist(monGr_GlobalHistograms,m_hitdistribution[gain_value] = new TH1F(name.c_str(),title.c_str(),8,0,8));
	m_hitdistribution[gain_value]->GetXaxis()->SetTitle("detector");
	m_hitdistribution[gain_value]->GetYaxis()->SetTitle("# of hits");

	name = "HitsVsLvl1A" + gain + "Gain";
	title = "Hits vs LVL1 A  " + gain + " Gain";
	status &= registerHist(monGr_GlobalHistograms,m_hits_lvl1a[gain_value] = new TH1F(name.c_str(),title.c_str(),32,0,32));
	m_hits_lvl1a[gain_value]->GetXaxis()->SetTitle("lvl1a [25 ns]");
	m_hits_lvl1a[gain_value]->GetYaxis()->SetTitle("# of hits");
       
	name = "HitsVsBCID" + gain + "Gain";
	title = "Hits vs BCID  " + gain + " Gain";
	status &= registerHist(monGr_GlobalHistograms,m_hits_bcid[gain_value] = new TH1F(name.c_str(),title.c_str(), bcids_displayed, bcid_start, bcid_start+bcids_displayed));
	m_hits_bcid[gain_value]->GetXaxis()->SetTitle("BCID [25 ns]");
	m_hits_bcid[gain_value]->GetYaxis()->SetTitle("# of hits");\
	
	name = "StrangeSignals" + gain + "Gain";
	title = "Strange Signals  " + gain + " Gain";
	status &= registerHist(monGr_GlobalHistograms,m_strange_signals[gain_value] = new TH1F(name.c_str(),title.c_str(),16,0,16));
	m_strange_signals[gain_value]->GetXaxis()->SetTitle("Channel");
	m_strange_signals[gain_value]->GetYaxis()->SetTitle("# of hits");

	name = "PulsePosition" + gain + "Gain";
	title = "Pulse Position  " + gain + " Gain";
	status &= registerHist(monGr_GlobalHistograms,m_pulse_position_gen[gain_value] = new TH1F(name.c_str(),title.c_str(),64, 0, 64));
	m_pulse_position_gen[gain_value]->GetXaxis()->SetTitle("time bins [25/64 ns]");
	m_pulse_position_gen[gain_value]->GetYaxis()->SetTitle("# of hits");

	name = "PulseWidth" + gain + "Gain";
	title = "Pulse Width  " + gain + " Gain";
	status &= registerHist(monGr_GlobalHistograms,m_pulse_width_gen[gain_value] = new TH1F(name.c_str(),title.c_str(),32, 0, 32));
	m_pulse_width_gen[gain_value]->GetXaxis()->SetTitle("time bins [25/64 ns]");
	m_pulse_width_gen[gain_value]->GetYaxis()->SetTitle("# of hits");

	name = "Deltat" + gain + "Gain";
	title = "#Delta t " + gain + " Gain";
	status &= registerHist(monGr_GlobalHistograms, m_deltat_vs_hits[gain_value] = new TH1F(name.c_str(),title.c_str(),50,-25,25));
	m_deltat_vs_hits[gain_value]->GetXaxis()->SetTitle("#Delta t [ns]");
	m_deltat_vs_hits[gain_value]->GetYaxis()->SetTitle("# of hits");
	
	name = "DeltatAligned" + gain + "Gain";
	title = "#Delta t Aligned " + gain + " Gain";
	status &= registerHist(monGr_GlobalHistograms, m_deltat_aligned_vs_hits[gain_value] = new TH1F(name.c_str(),title.c_str(),50,-25,25));
	m_deltat_aligned_vs_hits[gain_value]->GetXaxis()->SetTitle("#Delta t [ns]");
	m_deltat_aligned_vs_hits[gain_value]->GetYaxis()->SetTitle("# of hits");
	
	name = "DeltatVsBCID" + gain + "Gain";
	title = "#Delta t vs BCID " + gain + " Gain";
	status &= registerHist(monGr_GlobalHistograms, m_deltat_vs_bcid[gain_value] = new TH2F(name.c_str(),title.c_str(), bcids_displayed,bcid_start, bcid_start+bcids_displayed,50,-25,25));
	m_deltat_vs_bcid[gain_value]->GetXaxis()->SetTitle("BCID [25 ns]");
	m_deltat_vs_bcid[gain_value]->GetYaxis()->SetTitle("#Delta t [ns]");

	name = "DeltatAlignedVsBCID" + gain + "Gain";
	title = "#Delta t Aligned vs BCID " + gain + " Gain";
	status &= registerHist(monGr_GlobalHistograms, m_deltat_aligned_vs_bcid[gain_value] = new TH2F(name.c_str(),title.c_str(), bcids_displayed,bcid_start, bcid_start+bcids_displayed,50,-25,25));
	m_deltat_aligned_vs_bcid[gain_value]->GetXaxis()->SetTitle("BCID [25 ns]");
	m_deltat_aligned_vs_bcid[gain_value]->GetYaxis()->SetTitle("#Delta t [ns]");

	name = "DeltatVsLB" + gain + "Gain";
	title = "#Delta t vs LB " + gain + " Gain";
	status &= registerHist(monGr_GlobalHistograms, m_deltat_vs_lb[gain_value] = new TH2F(name.c_str(),title.c_str(), lb_max, lb_start, lb_max, 50,-25,25));
	m_deltat_vs_lb[gain_value]->GetXaxis()->SetTitle("LB number");
	m_deltat_vs_lb[gain_value]->GetYaxis()->SetTitle("#Delta t [ns]");

	name = "DeltatAlignedVsLB" + gain + "Gain";
	title = "#Delta t Aligned vs LB " + gain + " Gain";
	status &= registerHist(monGr_GlobalHistograms, m_deltat_aligned_vs_lb[gain_value] = new TH2F(name.c_str(),title.c_str(), lb_max, lb_start, lb_max, 50,-25,25));
	m_deltat_aligned_vs_lb[gain_value]->GetXaxis()->SetTitle("LB number");
	m_deltat_aligned_vs_lb[gain_value]->GetYaxis()->SetTitle("#Delta t [ns]");

	name = "DeltatVsECR" + gain + "Gain";
	title = "#Delta t vs ECR " + gain + " Gain";
	status &= registerHist(monGr_GlobalHistograms, m_deltat_vs_ecr[gain_value] = new TH2F(name.c_str(),title.c_str(), ecrs_displayed,ecr_start, ecr_start+ecrs_displayed,50,-25,25));
	m_deltat_vs_ecr[gain_value]->GetXaxis()->SetTitle("ECR");
	m_deltat_vs_ecr[gain_value]->GetYaxis()->SetTitle("#Delta t [ns]");

	name = "DeltatAlignedVsECR" + gain + "Gain";
	title = "#Delta t Aligned vs ECR " + gain + " Gain";
	status &= registerHist(monGr_GlobalHistograms, m_deltat_aligned_vs_ecr[gain_value] = new TH2F(name.c_str(),title.c_str(), ecrs_displayed,ecr_start, ecr_start+ecrs_displayed,50,-25,25));
	m_deltat_aligned_vs_ecr[gain_value]->GetXaxis()->SetTitle("ECR");
	m_deltat_aligned_vs_ecr[gain_value]->GetYaxis()->SetTitle("#Delta t [ns]");

	name = "DeltatVsPrVertex" + gain + "Gain";
	title = "#Delta t vs Primary Vertex " + gain + " Gain";
	status &= registerHist(monGr_GlobalHistograms, m_deltat_vs_PrVertex[gain_value] = new TH2F(name.c_str(),title.c_str(),100,-200,200,50,-25,25));
	m_deltat_vs_PrVertex[gain_value]->GetXaxis()->SetTitle("Primary Vertex z position [mm] ");
	m_deltat_vs_PrVertex[gain_value]->GetYaxis()->SetTitle("#Delta t [ns]");

	name = "DeltatVsPixelHits" + gain + "Gain";
	title = "#Delta t vs NumberPixelHits " + gain + " Gain";
	status &= registerHist(monGr_GlobalHistograms, m_deltat_vs_pixhits[gain_value] = new TH2F(name.c_str(),title.c_str(), 100, 0, 10000, 50,-25,25));
	m_deltat_vs_pixhits[gain_value]->GetXaxis()->SetTitle("Number Pixel Hits");
	m_deltat_vs_pixhits[gain_value]->GetYaxis()->SetTitle("#Delta t [ns]");
	  
	name = "DeltatVsPixelHitsBR" + gain + "Gain";
	title = "#Delta t vs NumberPixelHits in Barrel " + gain + " Gain";
	status &= registerHist(monGr_GlobalHistograms, m_deltat_vs_pixhitsBR[gain_value] = new TH2F(name.c_str(),title.c_str(), 100, 0, 10000, 50,-25,25));
	m_deltat_vs_pixhitsBR[gain_value]->GetXaxis()->SetTitle("Number Barrel Pixel Hits");
	m_deltat_vs_pixhitsBR[gain_value]->GetYaxis()->SetTitle("#Delta t [ns]");
	  
	name = "DeltatVsPixelHitsEC" + gain + "Gain";
	title = "#Delta t vs NumberPixelHits in Endcap " + gain + " Gain";
	status &= registerHist(monGr_GlobalHistograms, m_deltat_vs_pixhitsEC[gain_value] = new TH2F(name.c_str(),title.c_str(), 100, 0, 10000, 50,-25,25));
	m_deltat_vs_pixhitsEC[gain_value]->GetXaxis()->SetTitle("Number Endcap Pixel Hits");
	m_deltat_vs_pixhitsEC[gain_value]->GetYaxis()->SetTitle("#Delta t [ns]");
	  
	name = "DeltatVsSctHits" + gain + "Gain";
	title = "#Delta t vs NumberSctHits " + gain + " Gain";
	status &= registerHist(monGr_GlobalHistograms, m_deltat_vs_scthits[gain_value] = new TH2F(name.c_str(),title.c_str(), 100, 0, 20000, 50,-25,25));
	m_deltat_vs_scthits[gain_value]->GetXaxis()->SetTitle("Number Sct Hits");
	m_deltat_vs_scthits[gain_value]->GetYaxis()->SetTitle("#Delta t [ns]");
	  
	name = "DeltatVsSctHitsBR" + gain + "Gain";
	title = "#Delta t vs NumberSctHits in Barrel " + gain + " Gain";
	status &= registerHist(monGr_GlobalHistograms, m_deltat_vs_scthitsBR[gain_value] = new TH2F(name.c_str(),title.c_str(), 100, 0, 20000, 50,-25,25));
	m_deltat_vs_scthitsBR[gain_value]->GetXaxis()->SetTitle("Number Barrel Sct Hits");
	m_deltat_vs_scthitsBR[gain_value]->GetYaxis()->SetTitle("#Delta t [ns]");
	  
	name = "DeltatVsSctHitsEC" + gain + "Gain";
	title = "#Delta t vs NumberSctHits in Endcap " + gain + " Gain";
	status &= registerHist(monGr_GlobalHistograms, m_deltat_vs_scthitsEC[gain_value] = new TH2F(name.c_str(),title.c_str(), 100, 0, 20000, 50,-25,25));
	m_deltat_vs_scthitsEC[gain_value]->GetXaxis()->SetTitle("Number Endcap Sct Hits");
	m_deltat_vs_scthitsEC[gain_value]->GetYaxis()->SetTitle("#Delta t [ns]");

	name = "SctVsPixHitsCol" + gain + "Gain";
	title = "#NumberSctHits vs NumberPixHits w/ Collison delta t " + gain + " Gain";
	status &= registerHist(monGr_GlobalHistograms, m_sct_vs_pix_col[gain_value] = new TH2F(name.c_str(),title.c_str(), 100, 0, 20000, 100, 0, 10000));
	m_sct_vs_pix_col[gain_value]->GetXaxis()->SetTitle("Number Pixel Hits");
	m_sct_vs_pix_col[gain_value]->GetYaxis()->SetTitle("Number Sct Hits");
	  
	name = "SctVsPixHitsBkg" + gain + "Gain";
	title = "#NumberSctHits vs NumberPixHits w/ Background delta t " + gain + " Gain";
	status &= registerHist(monGr_GlobalHistograms, m_sct_vs_pix_bkg[gain_value] = new TH2F(name.c_str(),title.c_str(), 100, 0, 20000, 100, 0, 10000));
	m_sct_vs_pix_bkg[gain_value]->GetXaxis()->SetTitle("Number Pixel Hits");
	m_sct_vs_pix_bkg[gain_value]->GetYaxis()->SetTitle("Number Sct Hits");
	  
	for (unsigned int pulse_nr=0;pulse_nr<2;pulse_nr++) {

	  if (pulse_nr==0){pulse = "1";} 
	  else pulse="2";   
	    
	  name="Pulse"+pulse+"Width"+gain+"Gain";
	  title="Pulse"+ pulse +"Width" + gain  + "Gain";
	  status &= registerHist(monGr_GlobalHistograms, m_pulse_width[gain_value][pulse_nr] = new TH1F(name.c_str(),title.c_str(),32,0,32));
	  m_pulse_width[gain_value][pulse_nr]->GetXaxis()->SetTitle("time bins [25/32 ns]");
	  m_pulse_width[gain_value][pulse_nr]->GetYaxis()->SetTitle("# of hits");

	  name= "Pulse"+pulse+"Position"+ gain+"Gain";
	  title = "Pulse " + pulse + " Position  " + gain + " Gain";
	  status &= registerHist(monGr_GlobalHistograms, m_pulse_position[gain_value][pulse_nr] = new TH1F(name.c_str(),title.c_str(),64,0,64));
	  m_pulse_position[gain_value][pulse_nr]->GetXaxis()->SetTitle("time bins [25/64 ns]");
	  m_pulse_position[gain_value][pulse_nr]->GetYaxis()->SetTitle("# of hits");

	  name = "HitDistributionPulse"+pulse+gain+"Gain";
	  title = "Hitdistribution Pulse " + pulse + " " + gain + " Gain";
	  status &= registerHist(monGr_GlobalHistograms, m_hitdistribution_pulse[gain_value][pulse_nr] = new TH1F(name.c_str(),title.c_str(),8,0,8));
	  m_hitdistribution_pulse[gain_value][pulse_nr]->GetXaxis()->SetTitle("detector");
	  m_hitdistribution_pulse[gain_value][pulse_nr]->GetYaxis()->SetTitle("# of hits");

	}//end of pulse loop
      }//end of gain loop
    
    /**************************
     * Register Side Histograms
     ***************************/  
    
    for(unsigned int station_nr=0;station_nr<2;station_nr++)
      {
	if (station_nr==0) station_name="A-side";
	else station_name="C-side";
	
	for (unsigned int gain_value=0;gain_value<2;gain_value++)
	  {
	    if (gain_value==0){gain="Low";} 
	    else gain="High";   
	        
	    name="HitsVsLvl1A"+station_name+"Gain_gain";
	    title = "" + station_name + " Hits vs LVL1 A  " + gain + " Gain";
	    status &= registerHist(monGr_SideHistograms,m_stat_lvl1a[station_nr][gain_value] = new TH1F(name.c_str(),title.c_str(),32,0,32));
	    m_stat_lvl1a[station_nr][gain_value]->GetXaxis()->SetTitle("lvl1a [25 ns]");
	    m_stat_lvl1a[station_nr][gain_value]->GetYaxis()->SetTitle("# of hits");
	        
	    name="PulsePosition"+station_name+"Gain_"+gain;
	    title = "" + station_name + " Pulse Position  " + gain + " Gain";
	    status &= registerHist(monGr_SideHistograms, m_stat_pulse_position_gen[station_nr][gain_value] = new TH1F(name.c_str(),title.c_str(),64,0,64));
	    m_stat_pulse_position_gen[station_nr][gain_value]->GetXaxis()->SetTitle("time bins [25/64 ns]");
	    m_stat_pulse_position_gen[station_nr][gain_value]->GetYaxis()->SetTitle("# of hits");
	        
	    name="PulseWidth"+station_name+"Gain_"+gain;
	    title = "" + station_name + " Pulse Width  " + gain + " Gain";
	    status &= registerHist(monGr_SideHistograms,m_stat_pulse_width_gen[station_nr][gain_value] = new TH1F(name.c_str(),title.c_str(),32,0,32));
	    m_stat_pulse_width_gen[station_nr][gain_value]->GetXaxis()->SetTitle("time bins [25/32 ns]");
	    m_stat_pulse_width_gen[station_nr][gain_value]->GetYaxis()->SetTitle("# of hits");

	    for (unsigned int pulse_nr=0; pulse_nr<2;pulse_nr++)
	      {
		if (pulse_nr==0){pulse="1";} 
		else pulse="2";   
		
		name="Pulse"+station_name+"Position"+pulse+"Gain_"+gain;
		title=station_name+" Pulse "+pulse+" Position "+gain+" Gain";
		status &= registerHist(monGr_SideHistograms, m_stat_pulse_position[station_nr][gain_value][pulse_nr] = new TH1F(name.c_str(),title.c_str(),64,0,64));
		m_stat_pulse_position[station_nr][gain_value][pulse_nr]->GetXaxis()->SetTitle("time bins [25/64 ns]");
		m_stat_pulse_position[station_nr][gain_value][pulse_nr]->GetYaxis()->SetTitle("# of hits");
		
		name="Pulse"+station_name+"Width"+pulse+"Gain_"+gain;
		title=station_name+" Pulse "+pulse+" Width "+gain+" Gain";
		status &= registerHist(monGr_SideHistograms,m_stat_pulse_width[station_nr][gain_value][pulse_nr] = new TH1F(name.c_str(),title.c_str(),32,0,32));
		m_stat_pulse_width[station_nr][gain_value][pulse_nr]->GetXaxis()->SetTitle("time bins [25/32 ns]");
		m_stat_pulse_width[station_nr][gain_value][pulse_nr]->GetYaxis()->SetTitle("# of hits");
			
	      }//end of pulse_nr loop
	  }//end of gain loop
      }//end of side loop 
    
    /*************************
     * Register Detector Histograms
     *************************/
    std::string detector_name;
    
    for (unsigned int detector=0;detector<8;detector++) {

      if (detector<4) station_name="A-side";
      else station_name="C-side";
      
      switch(detector) {
      case 0:
	detector_name="+X";
	break;
      case 1:
	detector_name="+Y";
	break;
      case 2:     
	detector_name="-X";
	break;
      case 3:
	detector_name="-Y";
	break;
      case 4:
	detector_name="+X";
	break;
      case 5:
	detector_name="+Y";
	break;
      case 6:
	detector_name="-X";
	break;
      case 7:
	detector_name="-Y";
	break;
      }
      
      for (unsigned int gain_value=0;gain_value<2;gain_value++)
	{
	  if (gain_value==0){gain="Low";} 
	  else gain="High";   
	    
	  name="HitsVsLvl1A"+station_name+"Gain_"+detector_name+gain;
	  title= station_name+ " "+ detector_name +" Hits vs LVL1 A "+ gain +" Gain";
	  status &= registerHist(monGr_DetectorHistograms,m_det_lvl1a[detector][gain_value] = new TH1F(name.c_str(),title.c_str(),32,0,32));
	  m_det_lvl1a[detector][gain_value]->GetXaxis()->SetTitle("lvl1a [25 ns]");
	  m_det_lvl1a[detector][gain_value]->GetYaxis()->SetTitle("# of hits");
	    
	  name="PulsePosition"+station_name+"Gain_"+detector_name+gain;
	  title=station_name+" "+detector_name+" Pulse Position "+gain +" Gain";
	  status &= registerHist(monGr_DetectorHistograms, m_det_pulse_position_gen[detector][gain_value] = new TH1F(name.c_str(),title.c_str(),64,0,64));
	  m_det_pulse_position_gen[detector][gain_value]->GetXaxis()->SetTitle("time bins [25/64 ns]");
	  m_det_pulse_position_gen[detector][gain_value]->GetYaxis()->SetTitle("# of hits");
	    
	  name="PulseWidth"+station_name+"Gain_"+detector_name+gain;
	  title=station_name+" "+detector_name+" Pulse Width "+gain+" Gain";
	  status &= registerHist(monGr_DetectorHistograms,m_det_pulse_width_gen[detector][gain_value] = new TH1F(name.c_str(),title.c_str(),32,0,32));
	  m_det_pulse_width_gen[detector][gain_value]->GetXaxis()->SetTitle("time bins [25/32 ns]");
	  m_det_pulse_width_gen[detector][gain_value]->GetYaxis()->SetTitle("# of hits");
	  
	  name="timewalk"+station_name+"_"+detector_name+gain;
	  title=station_name+" "+detector_name+" timewalk "+ gain +" Gain";
	  status &= registerHist(monGr_GlobalHistograms,m_det_timewalk[detector][gain_value] = new TH2F(name.c_str(),title.c_str(),64, 0, 64, 32, 0, 32));
	  m_det_timewalk[detector][gain_value]->GetXaxis()->SetTitle("pulse pos in time bins [25/64 ns]");
	  m_det_timewalk[detector][gain_value]->GetYaxis()->SetTitle("pulse width in time bins [25/64 ns]");    
	  
	  for (unsigned int pulse_nr=0;pulse_nr<2;pulse_nr++) {
	    if (pulse_nr==0){pulse="1";} 
	    else pulse="2";   
	        
	    name="Pulse"+station_name+"Position"+detector_name+"Gain_"+pulse+gain;
	    title=station_name+" "+ detector_name +" Pulse "+pulse+" Position "+ gain +" Gain";
	    status &= registerHist(monGr_DetectorHistograms, m_det_pulse_position[detector][gain_value][pulse_nr] = new TH1F(name.c_str(),title.c_str(),64,0,64));
	    m_det_pulse_position[detector][gain_value][pulse_nr]->GetXaxis()->SetTitle("time bins [25/64 ns]");
	    m_det_pulse_position[detector][gain_value][pulse_nr]->GetYaxis()->SetTitle("# of hits");
	        
	    name="Pulse"+station_name+"Width"+detector_name+"Gain_"+pulse+gain;
	    title=station_name + " " + detector_name + " Pulse " + pulse + "  Width " + gain + " Gain";
	    status &= registerHist(monGr_DetectorHistograms,m_det_pulse_width[detector][gain_value][pulse_nr] = new TH1F(name.c_str(),title.c_str(),32,0,33));
	    m_det_pulse_width[detector][gain_value][pulse_nr]->GetXaxis()->SetTitle("time bins [25/32 ns]");
	    m_det_pulse_width[detector][gain_value][pulse_nr]->GetYaxis()->SetTitle("# of hits");
	  }//end of pulse loop
	}//end of gain loop 
    }//end of detector loop
  } 

  if (status) return StatusCode::SUCCESS;
  else return StatusCode::FAILURE;
}

//---------------------------------------------------------
/***************************************************************************
 * At first a method to fill the global and the pulse 1 and 2 histos
 **************************************************************************/
void FillGlobalHistos(TH2F* h1, TH2F* h2, TH2F* h3, TH1F* h4, TH1F* h5, int pulsewid, int channel, int bcid, int ecr, int lvl1a){
  if (pulsewid!=0){
    h1->Fill(lvl1a, channel);
    h2->Fill(bcid, channel);
    h3->Fill(ecr, channel);
    h4->Fill(lvl1a);
    h5->Fill(bcid);
  }
}

void FillPulseHistos(TH1F* h1, TH1F* h2, TH1F* h3, TH1F* h4, int pulsepos, int pulsewid, int tmp = 0, TH1F* h5 = nullptr){
  if (pulsewid!=0){
    h1->Fill(pulsepos);
    h2->Fill(pulsewid);
    h3->Fill(pulsepos);
    h4->Fill(pulsewid);
    if (h5 != nullptr){
      h5->Fill(tmp);
    }
  }
}
/****************************************************************************
 * Some methods one needs for the calculation of Delta t
 ***************************************************************************/
deltat_data::deltat_data():
    bcid(0),
    ecr(0),
    position(0),
    detector(0),
    lvl1a(0)
{
}

bool deltat_data::operator<(const deltat_data &data) const{
  return (bcid<data.bcid);
}

bool bcid_select(const deltat_data &data){
  return (data.bcid <data.bcid_max);
}

StatusCode InDetGlobalBCMTool::fillHistograms(){

  //Filling of histograms (loop over collections) :
  
  if      (m_detector == "sct") { m_nExampleInt = 2; }   
  else if (m_detector == "trt") { m_nExampleInt = 3; }
  else if (m_detector == "pixel"){m_nExampleInt = 1; }
  else if (m_detector == "ID" ) { m_nExampleInt = 4; }
  else { m_nExampleInt = 0; }
  
  m_nExamplePlot->Fill(m_nExampleInt);

  //retrieve LB number for plots vs LB
  SG::ReadHandle<xAOD::EventInfo> evtInfo(m_eventInfoKey,Gaudi::Hive::currentContext());
  if (!evtInfo.isValid()) {
    m_current_LB = 0;
    if ( msgLvl(MSG::WARNING) ){
      msg(MSG::WARNING) << "Could not retrieve the event information container" << endmsg;
    }
  } else {
      m_current_LB = evtInfo->lumiBlock();
  }

  //lists for calculation of deltat  
  std::vector<std::list<deltat_data> > positions_A(2);
  std::vector<std::list<deltat_data> > positions_C(2);

  //double z to store Primary Vertex z position
  double z = 1001; //that it can not mix with data if primary vertex calculation fails, z_max = 1000
  
  // Primary vertex monitoring
  SG::ReadHandle<VxContainer> vxContainer(m_vxContainerName);
  if (evtStore()->contains<VxContainer>(m_vxContainerName.key())) {
    if (!vxContainer.isValid()) {
      ATH_MSG_DEBUG ("Could not retrieve primary vertex container with key "+m_vxContainerName.key());
      return StatusCode::SUCCESS;
    } else {
	  
	int nPriVtx = 0;
	int nPileupVtx = 0;
	
	if (vxContainer->size() == 1) {
	    if ( msgLvl(MSG::DEBUG ) ) {
		msg(MSG::DEBUG) <<  "vxContainer size = 1 ->contains only dummy vertex" << endmsg; 
	    }
	}
	else{
	    for (VxContainer::const_iterator vxIter = vxContainer->begin(); vxIter != vxContainer->end(); ++vxIter)
	    {
		// Count different types of vertices
		if ((*vxIter)->vertexType() == Trk::PriVtx) nPriVtx++;
		if ((*vxIter)->vertexType() == Trk::PileUp) nPileupVtx++;
		
		// Select primary vertex
		if ((*vxIter)->vertexType() != Trk::PriVtx) continue;
		z = (*vxIter)->recVertex().position().z();
	    }
	}
    }
  } else {
      ATH_MSG_DEBUG ("StoreGate doesn't contain primary vertex container with key "+m_vxContainerName.key());
  }

  
  if ( m_BCM_RDO != nullptr && !m_BCM_RDO->empty() ){
    BCM_RDO_Container::const_iterator BCM_RDO_itr     = m_BCM_RDO->begin();
    BCM_RDO_Container::const_iterator BCM_RDO_itr_end = m_BCM_RDO->end();

    bool ROD0 = true; //bools for m_NumberOfEvents histo
    bool ROD1 = true;

    //int abortcountHGRod0 = 0; //FIXME abortFraction calculator
    //int abortcountLGRod0 = 0;
    //int abortcountHGRod1 = 0;
    //int abortcountLGRod1 = 0;

    int nROD0HitLG[bc_readout] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int nROD0HitHG[bc_readout] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int nROD1HitLG[bc_readout] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int nROD1HitHG[bc_readout] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int nROD0BCID[bc_readout]  = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int nROD1BCID[bc_readout]  = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    unsigned int ecr = -1;
    int ROD = -1;

    int BC_counter = -1;
    int channel_counter = -1;

    for (; BCM_RDO_itr != BCM_RDO_itr_end; ++BCM_RDO_itr) { // loops over 16 channels (HR counter 0 to 15)
      channel_counter++;
      if ( !(*BCM_RDO_itr)->empty()) {
	BCM_RDO_Collection::const_iterator RDO_element        = (*BCM_RDO_itr)->begin();
	BCM_RDO_Collection::const_iterator RDO_element_last   = (*BCM_RDO_itr)->end();
    	BC_counter = -1;
	
	for (; RDO_element != RDO_element_last; ++RDO_element){ // loops over 31 BCs read out per L1A
          ++BC_counter;
	  if (*RDO_element == nullptr)
	  {
	      ATH_MSG_WARNING ("NULL pointer!");
	      continue;
	  }
	    
	  int bcm_lvl1a           = (*RDO_element)->getLVL1A();
	  int bcm_channel         = (*RDO_element)->getChannel();
	  int bcm_bcid            = (*RDO_element)->getBCID();
	  int bcm_pulse1position  = (*RDO_element)->getPulse1Position(); 
	  int bcm_pulse1width     = (*RDO_element)->getPulse1Width();
	  int bcm_pulse2position  = (*RDO_element)->getPulse2Position();
	  int bcm_pulse2width     = (*RDO_element)->getPulse2Width();




	  /*********************************
	  *Filling Abort Fraction arrays
	  **********************************/
	  
	  if (channel_counter == 0 || channel_counter == 2 || channel_counter == 4 || channel_counter == 6  || channel_counter ==  9 || channel_counter ==  11 || channel_counter ==  13 || channel_counter == 15)  ROD = 0;
	  else ROD = 1;
	  if (bcm_pulse1width>0 || bcm_pulse2width>0){
	    if (ROD == 0 && channel_counter < 8 && nROD0HitLG[BC_counter]<3) nROD0HitLG[BC_counter]++;
	    if (ROD == 1 && channel_counter < 8 && nROD1HitLG[BC_counter]<3) nROD1HitLG[BC_counter]++;
	    if (ROD == 0) nROD0BCID[BC_counter] = bcm_bcid;
	    if (ROD == 0 && channel_counter > 7 && nROD0HitHG[BC_counter]<3) nROD0HitHG[BC_counter]++;
	    if (ROD == 1 && channel_counter > 7 && nROD1HitHG[BC_counter]<3) nROD1HitHG[BC_counter]++;
	    if (ROD == 1) nROD1BCID[BC_counter] = bcm_bcid;
	  }








	
		  
	  /*********************************
	   *Filling NumberOfEvent histograms
	   *********************************/
	  if (ROD0 && (bcm_channel == 0 || bcm_channel == 2 || bcm_channel == 4 || bcm_channel == 6 || bcm_channel == 9 || bcm_channel == 11 || bcm_channel == 13 || bcm_channel == 15)){
	    m_NumberOfEvents->Fill(1.5);
	    m_NumberOfEvents->Fill(0.5);
	    ROD0 = false; // FIXME seems to be a bug: first loop is over channel, so alternating RODs, bit bool ROD0 (and ROD1) is being set to false after filling first time! as well filled, if all BC was empty!
	  }
	  if (ROD1 && (bcm_channel == 1 || bcm_channel == 3 || bcm_channel == 5 || bcm_channel == 7 || bcm_channel == 8 || bcm_channel == 10 || bcm_channel == 12 || bcm_channel == 14)){
	    m_NumberOfEvents->Fill(2.5);
	    m_NumberOfEvents->Fill(0.5);
	    ROD1 = false;
	  }
	    
	  if ( bcm_pulse1width !=0 || bcm_pulse2width !=0){
	        
	    // strange signals: You see something where there should be nothing, because the pulse widths are zero. // FIXME condition never fulfillable
	    if (( bcm_pulse1width == 0 && bcm_pulse2width == 0) &&  (bcm_pulse1position !=0 || bcm_pulse2position !=0)){
	      if (bcm_channel < 8) m_strange_signals[0]->Fill(bcm_channel);
	      if (bcm_channel > 7) m_strange_signals[1]->Fill(bcm_channel);
	    }
	        
	    //Event information to get ECR
	    ecr = evtInfo.isValid() ? (evtInfo->extendedLevel1ID() & 0xff000000) >> 24 : 0;
	        
	    /******************************************************************************************
	     *Global Histograms
	     ******************************************************************************************/

	    if ( bcm_channel < 8 ){
	      if (bcm_pulse1width != 0) m_hitdistribution[0]->Fill(bcm_channel);
	      if (bcm_pulse2width != 0) m_hitdistribution[0]->Fill(bcm_channel);
	      FillGlobalHistos(m_ChannelVsLvl1a, m_ChannelVsBCID, m_ChannelVsECR, m_hits_lvl1a[0], m_hits_bcid[0], bcm_pulse1width, bcm_channel, bcm_bcid, ecr, bcm_lvl1a);
	      FillGlobalHistos(m_ChannelVsLvl1a, m_ChannelVsBCID, m_ChannelVsECR, m_hits_lvl1a[0], m_hits_bcid[0], bcm_pulse2width, bcm_channel, bcm_bcid, ecr, bcm_lvl1a);
	      FillPulseHistos(m_pulse_position_gen[0], m_pulse_width_gen[0], m_pulse_position[0][0], m_pulse_width[0][0], bcm_pulse1position, bcm_pulse1width, bcm_channel,  m_hitdistribution_pulse[0][0]);
	      FillPulseHistos(m_pulse_position_gen[0], m_pulse_width_gen[0], m_pulse_position[0][1], m_pulse_width[0][1], bcm_pulse2position, bcm_pulse2width, bcm_channel,  m_hitdistribution_pulse[0][1]);
	    }
	    if (bcm_channel > 7){
	      if (bcm_pulse1width != 0) m_hitdistribution[1]->Fill(bcm_channel-8);
	      if (bcm_pulse2width != 0) m_hitdistribution[1]->Fill(bcm_channel-8);
	      FillGlobalHistos(m_ChannelVsLvl1a, m_ChannelVsBCID, m_ChannelVsECR, m_hits_lvl1a[1], m_hits_bcid[1], bcm_pulse1width, bcm_channel, bcm_bcid, ecr, bcm_lvl1a);
	      FillGlobalHistos(m_ChannelVsLvl1a, m_ChannelVsBCID, m_ChannelVsECR, m_hits_lvl1a[1], m_hits_bcid[1], bcm_pulse2width, bcm_channel, bcm_bcid, ecr, bcm_lvl1a);
	      FillPulseHistos(m_pulse_position_gen[1], m_pulse_width_gen[1], m_pulse_position[1][0], m_pulse_width[1][0], bcm_pulse1position, bcm_pulse1width, bcm_channel-8, m_hitdistribution_pulse[1][0]);
	      FillPulseHistos(m_pulse_position_gen[1], m_pulse_width_gen[1], m_pulse_position[1][1], m_pulse_width[1][1], bcm_pulse2position, bcm_pulse2width, bcm_channel-8, m_hitdistribution_pulse[1][1]);
	    }
	          
	    /************************************************************************************************
	     *Side Histograms
	     *************************************************************************************************/
	          
	    // A Side Low Gain
	    if (bcm_channel < 4){
	      deltat_data s1a;
	      s1a.bcid = bcm_bcid;
	      s1a.position = bcm_pulse1position;
	      s1a.ecr = ecr;
	      s1a.detector = bcm_channel%4;
              s1a.lvl1a = bcm_lvl1a;
	      positions_A[0].push_back(s1a);
	      if (bcm_pulse2width != 0){
		deltat_data s2a;
		s2a.bcid = bcm_bcid;
		s2a.position = bcm_pulse2position;
		s2a.ecr = ecr;
	        s2a.detector = bcm_channel%4;
		s2a.lvl1a = bcm_lvl1a;
		positions_A[0].push_back(s2a);
	      }
	      FillPulseHistos(m_stat_pulse_position_gen[0][0], m_stat_pulse_width_gen[0][0], m_stat_pulse_position[0][0][0], m_stat_pulse_width[0][0][0], bcm_pulse1position, bcm_pulse1width, bcm_lvl1a, m_stat_lvl1a[0][0]);
	      FillPulseHistos(m_stat_pulse_position_gen[0][0], m_stat_pulse_width_gen[0][0], m_stat_pulse_position[0][0][1], m_stat_pulse_width[0][0][1], bcm_pulse2position, bcm_pulse2width, bcm_lvl1a, m_stat_lvl1a[0][0]);
	    }

	    // A Side High Gain
	    if (bcm_channel > 7 && bcm_channel < 12){
	      deltat_data s1a;
	      s1a.bcid = bcm_bcid;
	      s1a.position = bcm_pulse1position;
	      s1a.ecr = ecr;
	      s1a.detector = bcm_channel%4;
              s1a.lvl1a = bcm_lvl1a;
	      positions_A[1].push_back(s1a);
	      if (bcm_pulse2width != 0){
		deltat_data s2a;
		s2a.bcid = bcm_bcid;
		s2a.position = bcm_pulse2position;
		s2a.ecr = ecr;
	        s2a.detector = bcm_channel%4;
		s2a.lvl1a = bcm_lvl1a;
		positions_A[1].push_back(s2a);
	      }
	      FillPulseHistos(m_stat_pulse_position_gen[0][1], m_stat_pulse_width_gen[0][1], m_stat_pulse_position[0][1][0], m_stat_pulse_width[0][1][0], bcm_pulse1position, bcm_pulse1width, bcm_lvl1a, m_stat_lvl1a[0][1]);
	      FillPulseHistos(m_stat_pulse_position_gen[0][1], m_stat_pulse_width_gen[0][1], m_stat_pulse_position[0][1][1], m_stat_pulse_width[0][1][1], bcm_pulse2position, bcm_pulse2width, bcm_lvl1a, m_stat_lvl1a[0][1]);
	    }
	    // C Side Low Gain
	    if (bcm_channel > 3 && bcm_channel < 8){
	      deltat_data s1c;
	      s1c.bcid = bcm_bcid;
	      s1c.position = bcm_pulse1position;
	      s1c.ecr = ecr;
	      s1c.detector = bcm_channel%4;
              s1c.lvl1a = bcm_lvl1a;
	      positions_C[0].push_back(s1c);
	      if (bcm_pulse2width != 0){
		deltat_data s2c;
		s2c.bcid = bcm_bcid;
		s2c.position = bcm_pulse2position;
		s2c.ecr = ecr;
	        s2c.detector = bcm_channel%4;
		s2c.lvl1a = bcm_lvl1a;
		positions_C[0].push_back(s2c);
	      }
	      FillPulseHistos(m_stat_pulse_position_gen[1][0], m_stat_pulse_width_gen[1][0], m_stat_pulse_position[1][0][0], m_stat_pulse_width[1][0][0], bcm_pulse1position, bcm_pulse1width, bcm_lvl1a, m_stat_lvl1a[1][0]);
	      FillPulseHistos(m_stat_pulse_position_gen[1][0], m_stat_pulse_width_gen[1][0], m_stat_pulse_position[1][0][1], m_stat_pulse_width[1][0][1], bcm_pulse2position, bcm_pulse2width, bcm_lvl1a, m_stat_lvl1a[1][0]);
	    }
	    // C Side High Gain
	    if (bcm_channel > 11){
	      deltat_data s1c;
	      s1c.bcid = bcm_bcid;
	      s1c.position = bcm_pulse1position;
	      s1c.ecr = ecr;
	      s1c.detector = bcm_channel%4;
              s1c.lvl1a = bcm_lvl1a;
	      positions_C[1].push_back(s1c);
	      if (bcm_pulse2width != 0){
		deltat_data s2c;
		s2c.bcid = bcm_bcid;
		s2c.position = bcm_pulse2position;
		s2c.ecr = ecr;
	        s2c.detector = bcm_channel%4;
		s2c.lvl1a = bcm_lvl1a;
		positions_C[1].push_back(s2c);
	      }
	      FillPulseHistos(m_stat_pulse_position_gen[1][1], m_stat_pulse_width_gen[1][1], m_stat_pulse_position[1][1][0], m_stat_pulse_width[1][1][0], bcm_pulse1position, bcm_pulse1width, bcm_lvl1a, m_stat_lvl1a[1][1]);
	      FillPulseHistos(m_stat_pulse_position_gen[1][1], m_stat_pulse_width_gen[1][1], m_stat_pulse_position[1][1][1], m_stat_pulse_width[1][1][1], bcm_pulse2position, bcm_pulse2width, bcm_lvl1a, m_stat_lvl1a[1][1]);
	    }
	          
	    /****************************************************************************
	     *And filling the detector histograms
	     ****************************************************************************/
	    switch (bcm_channel){
	    case 0:
	      FillPulseHistos(m_det_pulse_position_gen[0][0], m_det_pulse_width_gen[0][0], m_det_pulse_position[0][0][0], m_det_pulse_width[0][0][0], bcm_pulse1position, bcm_pulse1width, bcm_lvl1a, m_det_lvl1a[0][0]);
	      FillPulseHistos(m_det_pulse_position_gen[0][0], m_det_pulse_width_gen[0][0], m_det_pulse_position[0][0][1], m_det_pulse_width[0][0][1], bcm_pulse2position, bcm_pulse2width, bcm_lvl1a, m_det_lvl1a[0][0]);
	      break;
	    case 1:
	      FillPulseHistos(m_det_pulse_position_gen[1][0], m_det_pulse_width_gen[1][0], m_det_pulse_position[1][0][0], m_det_pulse_width[1][0][0], bcm_pulse1position, bcm_pulse1width, bcm_lvl1a, m_det_lvl1a[1][0]);
	      FillPulseHistos(m_det_pulse_position_gen[1][0], m_det_pulse_width_gen[1][0], m_det_pulse_position[1][0][1], m_det_pulse_width[1][0][1], bcm_pulse2position, bcm_pulse2width, bcm_lvl1a, m_det_lvl1a[1][0]);
	      break;
	    case 2:
	      FillPulseHistos(m_det_pulse_position_gen[2][0], m_det_pulse_width_gen[2][0], m_det_pulse_position[2][0][0], m_det_pulse_width[2][0][0], bcm_pulse1position, bcm_pulse1width, bcm_lvl1a, m_det_lvl1a[2][0]);
	      FillPulseHistos(m_det_pulse_position_gen[2][0], m_det_pulse_width_gen[2][0], m_det_pulse_position[2][0][1], m_det_pulse_width[2][0][1], bcm_pulse2position, bcm_pulse2width, bcm_lvl1a, m_det_lvl1a[2][0]);
	      break;
	    case 3:
	      FillPulseHistos(m_det_pulse_position_gen[3][0], m_det_pulse_width_gen[3][0], m_det_pulse_position[3][0][0], m_det_pulse_width[3][0][0], bcm_pulse1position, bcm_pulse1width, bcm_lvl1a, m_det_lvl1a[3][0]);
	      FillPulseHistos(m_det_pulse_position_gen[3][0], m_det_pulse_width_gen[3][0], m_det_pulse_position[3][0][1], m_det_pulse_width[3][0][1], bcm_pulse2position, bcm_pulse2width, bcm_lvl1a, m_det_lvl1a[3][0]);
	      break;
	    case 4:
	      FillPulseHistos(m_det_pulse_position_gen[4][0], m_det_pulse_width_gen[4][0], m_det_pulse_position[4][0][0], m_det_pulse_width[4][0][0], bcm_pulse1position, bcm_pulse1width, bcm_lvl1a, m_det_lvl1a[4][0]);
	      FillPulseHistos(m_det_pulse_position_gen[4][0], m_det_pulse_width_gen[4][0], m_det_pulse_position[4][0][1], m_det_pulse_width[4][0][1], bcm_pulse2position, bcm_pulse2width, bcm_lvl1a, m_det_lvl1a[4][0]);
	      break;
	    case 5:
	      FillPulseHistos(m_det_pulse_position_gen[5][0], m_det_pulse_width_gen[5][0], m_det_pulse_position[5][0][0], m_det_pulse_width[5][0][0], bcm_pulse1position, bcm_pulse1width, bcm_lvl1a, m_det_lvl1a[5][0]);
	      FillPulseHistos(m_det_pulse_position_gen[5][0], m_det_pulse_width_gen[5][0], m_det_pulse_position[5][0][1], m_det_pulse_width[5][0][1], bcm_pulse2position, bcm_pulse2width, bcm_lvl1a, m_det_lvl1a[5][0]);
	      break;
	    case 6:
	      FillPulseHistos(m_det_pulse_position_gen[6][0], m_det_pulse_width_gen[6][0], m_det_pulse_position[6][0][0], m_det_pulse_width[6][0][0], bcm_pulse1position, bcm_pulse1width, bcm_lvl1a, m_det_lvl1a[6][0]);
	      FillPulseHistos(m_det_pulse_position_gen[6][0], m_det_pulse_width_gen[6][0], m_det_pulse_position[6][0][1], m_det_pulse_width[6][0][1], bcm_pulse2position, bcm_pulse2width, bcm_lvl1a, m_det_lvl1a[6][0]);
	      break;
	    case 7:
	      FillPulseHistos(m_det_pulse_position_gen[7][0], m_det_pulse_width_gen[7][0], m_det_pulse_position[7][0][0], m_det_pulse_width[7][0][0], bcm_pulse1position, bcm_pulse1width, bcm_lvl1a, m_det_lvl1a[7][0]);
	      FillPulseHistos(m_det_pulse_position_gen[7][0], m_det_pulse_width_gen[7][0], m_det_pulse_position[7][0][1], m_det_pulse_width[7][0][1], bcm_pulse2position, bcm_pulse2width, bcm_lvl1a, m_det_lvl1a[7][0]);
	      break;
	    case 8:
	      FillPulseHistos(m_det_pulse_position_gen[0][1], m_det_pulse_width_gen[0][1], m_det_pulse_position[0][1][0], m_det_pulse_width[0][1][0], bcm_pulse1position, bcm_pulse1width, bcm_lvl1a, m_det_lvl1a[0][1]);
	      FillPulseHistos(m_det_pulse_position_gen[0][1], m_det_pulse_width_gen[0][1], m_det_pulse_position[0][1][1], m_det_pulse_width[0][1][1], bcm_pulse2position, bcm_pulse2width, bcm_lvl1a, m_det_lvl1a[0][1]);
	      break;
	    case 9:
	      FillPulseHistos(m_det_pulse_position_gen[1][1], m_det_pulse_width_gen[1][1], m_det_pulse_position[1][1][0], m_det_pulse_width[1][1][0], bcm_pulse1position, bcm_pulse1width, bcm_lvl1a, m_det_lvl1a[1][1]);
	      FillPulseHistos(m_det_pulse_position_gen[1][1], m_det_pulse_width_gen[1][1], m_det_pulse_position[1][1][1], m_det_pulse_width[1][1][1], bcm_pulse2position, bcm_pulse2width, bcm_lvl1a, m_det_lvl1a[1][1]);
	      break;
	    case 10:
	      FillPulseHistos(m_det_pulse_position_gen[2][1], m_det_pulse_width_gen[2][1], m_det_pulse_position[2][1][0], m_det_pulse_width[2][1][0], bcm_pulse1position, bcm_pulse1width, bcm_lvl1a, m_det_lvl1a[2][1]);
	      FillPulseHistos(m_det_pulse_position_gen[2][1], m_det_pulse_width_gen[2][1], m_det_pulse_position[2][1][1], m_det_pulse_width[2][1][1], bcm_pulse2position, bcm_pulse2width, bcm_lvl1a, m_det_lvl1a[2][1]);
	      break;
	    case 11:
	      FillPulseHistos(m_det_pulse_position_gen[3][1], m_det_pulse_width_gen[3][1], m_det_pulse_position[3][1][0], m_det_pulse_width[3][1][0], bcm_pulse1position, bcm_pulse1width, bcm_lvl1a, m_det_lvl1a[3][1]);
	      FillPulseHistos(m_det_pulse_position_gen[3][1], m_det_pulse_width_gen[3][1], m_det_pulse_position[3][1][1], m_det_pulse_width[3][1][1], bcm_pulse2position, bcm_pulse2width, bcm_lvl1a, m_det_lvl1a[3][1]);
	      break;
	    case 12:
	      FillPulseHistos(m_det_pulse_position_gen[4][1], m_det_pulse_width_gen[4][1], m_det_pulse_position[4][1][0], m_det_pulse_width[4][1][0], bcm_pulse1position, bcm_pulse1width, bcm_lvl1a, m_det_lvl1a[4][1]);
	      FillPulseHistos(m_det_pulse_position_gen[4][1], m_det_pulse_width_gen[4][1], m_det_pulse_position[4][1][1], m_det_pulse_width[4][1][1], bcm_pulse2position, bcm_pulse2width, bcm_lvl1a, m_det_lvl1a[4][1]);
	      break;
	    case 13:
	      FillPulseHistos(m_det_pulse_position_gen[5][1], m_det_pulse_width_gen[5][1], m_det_pulse_position[5][1][0], m_det_pulse_width[5][1][0], bcm_pulse1position, bcm_pulse1width, bcm_lvl1a, m_det_lvl1a[5][1]);
	      FillPulseHistos(m_det_pulse_position_gen[5][1], m_det_pulse_width_gen[5][1], m_det_pulse_position[5][1][1], m_det_pulse_width[5][1][1], bcm_pulse2position, bcm_pulse2width, bcm_lvl1a, m_det_lvl1a[5][1]);
	      break;
	    case 14:
	      FillPulseHistos(m_det_pulse_position_gen[6][1], m_det_pulse_width_gen[6][1], m_det_pulse_position[6][1][0], m_det_pulse_width[6][1][0], bcm_pulse1position, bcm_pulse1width, bcm_lvl1a, m_det_lvl1a[6][1]);
	      FillPulseHistos(m_det_pulse_position_gen[6][1], m_det_pulse_width_gen[6][1], m_det_pulse_position[6][1][1], m_det_pulse_width[6][1][1], bcm_pulse2position, bcm_pulse2width, bcm_lvl1a, m_det_lvl1a[6][1]);
	      break;
	    case 15:
	      FillPulseHistos(m_det_pulse_position_gen[7][1], m_det_pulse_width_gen[7][1], m_det_pulse_position[7][1][0], m_det_pulse_width[7][1][0], bcm_pulse1position, bcm_pulse1width, bcm_lvl1a, m_det_lvl1a[7][1]);
	      FillPulseHistos(m_det_pulse_position_gen[7][1], m_det_pulse_width_gen[7][1], m_det_pulse_position[7][1][1], m_det_pulse_width[7][1][1], bcm_pulse2position, bcm_pulse2width, bcm_lvl1a, m_det_lvl1a[7][1]);
	      break;
	    }
	  }
	}// End of for BCM_RDO_element
      }
    } // End of for BCM_RDO_itr 
    

    /********************************************    
     *count pixel hits
     ********************************************/    
    unsigned int n_pix_hits[3] = {0};

    if(m_pixRdoContainer != nullptr){
        PixelRDO_Container::const_iterator colNextpix = m_pixRdoContainer->begin();
        PixelRDO_Container::const_iterator colNextpix_end = m_pixRdoContainer->end();
        for( ; colNextpix != colNextpix_end; ++colNextpix){
            const InDetRawDataCollection<PixelRDORawData>* PIX_Collection(*colNextpix);
	    if(!PIX_Collection) continue;

            Identifier id = PIX_Collection->identify();
            switch ( m_pixelID->barrel_ec( id )  )
                    {
                    case 0:   // Barrel
                        n_pix_hits[1] += PIX_Collection->size();
                        break;
                    case -2:  // ECA
                        n_pix_hits[0] += PIX_Collection->size();
                        break;
                    case 2:   // ECC
                        n_pix_hits[2] += PIX_Collection->size();
                        break;
                    }
        }
    }


    /********************************************    
     *count sct hits
     ********************************************/    
    unsigned int n_sct_hits[3] = {0};

    if(m_sctRdoContainer != nullptr){
        SCT_RDO_Container::const_iterator colNextsct = m_sctRdoContainer->begin();
        SCT_RDO_Container::const_iterator colNextsct_end = m_sctRdoContainer->end();
        for( ; colNextsct != colNextsct_end; ++colNextsct){
            const InDetRawDataCollection<SCT_RDORawData>* SCT_Collection(*colNextsct);
	    if(!SCT_Collection) continue;

            Identifier id = SCT_Collection->identify();
            switch ( m_sctID->barrel_ec( id )  )
                    {
                    case 0:   // Barrel
                        n_sct_hits[1] += SCT_Collection->size();
                        break;
                    case -2:  // ECA
                        n_sct_hits[0] += SCT_Collection->size();
                        break;
                    case 2:   // ECC
                        n_sct_hits[2] += SCT_Collection->size();
                        break;
                    }
        }
    }

    /********************************************    
     *Delta t histograms
     ********************************************/
    for(unsigned int gain =0;gain<2;gain++){
     
      positions_A[gain].sort();
      positions_C[gain].sort();

      /*
	Code to calculate the weight for the event. Could be useful in the future
	but taken out now as it is unused at the moment, causing warnings during
	compilation.
	
      float weight = 1.0;
      float ndt = 0;

      for (std::list<deltat_data>::iterator it_a =positions_A[gain].begin();it_a!=positions_A[gain].end();it_a++)
          for (std::list<deltat_data>::iterator it_c =positions_C[gain].begin();it_c!=positions_C[gain].end();it_c++)
              if ( (*it_a).bcid == (*it_c).bcid && (*it_a).lvl1a == BCM_LVL1A && (*it_c).lvl1a == BCM_LVL1A )
                  ndt++;

      if ( ndt )
	  weight = 1 / ndt;
      */
      float deltat=350.0;//so that it can't mix with normal data if calculation of deltat fails
      while (!positions_A[gain].empty()){
	unsigned int bcid=(positions_A[gain].front()).bcid;
	unsigned int detector_a=(positions_A[gain].front()).detector;


	for (std::list<deltat_data>::iterator it_c =positions_C[gain].begin();it_c!=positions_C[gain].end();++it_c){
	  if(bcid<(*it_c).bcid) continue;
	  if (bcid==(*it_c).bcid) { // i.e. (positions_A[gain].front()).bcid == (positions_C[gain].begin()).bcid
	    int deltatbins=(*it_c).position-(positions_A[gain].front()).position;
	    deltat = deltatbins/64.0*25.0;
	    m_deltat_vs_hits[gain]->Fill(deltat);
	    m_deltat_vs_bcid[gain]->Fill((*it_c).bcid, deltat);
	    m_deltat_vs_lb[gain]->Fill(m_current_LB, deltat);
	    m_deltat_vs_ecr[gain]->Fill((*it_c).ecr, deltat);
	    m_deltat_vs_PrVertex[gain]->Fill(z, deltat);
	    m_deltat_vs_pixhits[gain]->Fill(n_pix_hits[0]+n_pix_hits[1]+n_pix_hits[2], deltat);
	    m_deltat_vs_pixhitsBR[gain]->Fill(n_pix_hits[1], deltat);
	    m_deltat_vs_pixhitsEC[gain]->Fill(n_pix_hits[0]+n_pix_hits[2], deltat);
	    m_deltat_vs_scthits[gain]->Fill(n_sct_hits[0]+n_sct_hits[1]+n_sct_hits[2], deltat);
	    m_deltat_vs_scthitsBR[gain]->Fill(n_sct_hits[1], deltat);
	    m_deltat_vs_scthitsEC[gain]->Fill(n_sct_hits[0]+n_sct_hits[2], deltat);
	    if (abs(deltat) < 6.25) m_sct_vs_pix_col[gain]->Fill(n_sct_hits[0]+n_sct_hits[1]+n_sct_hits[2], n_pix_hits[0]+n_pix_hits[1]+n_pix_hits[2]);
	    else m_sct_vs_pix_bkg[gain]->Fill(n_sct_hits[0]+n_sct_hits[1]+n_sct_hits[2], n_pix_hits[0]+n_pix_hits[1]+n_pix_hits[2]);
	    if (detector_a==(*it_c).detector) {
	      m_deltat_aligned_vs_ecr[gain]->Fill(ecr,deltat);
	      m_deltat_aligned_vs_bcid[gain]->Fill(bcid,deltat);
	      m_deltat_aligned_vs_lb[gain]->Fill(m_current_LB, deltat);
	      m_deltat_aligned_vs_hits[gain]->Fill(deltat);
	    }
	  }//end of bcid if
	}//end of c loop
    // Need to set bcid_max before using bcid_select method
    for (deltat_data& data : positions_C[gain]) {
      data.bcid_max = (positions_A[gain].front()).bcid;
    }
	positions_C[gain].remove_if(bcid_select);
	positions_A[gain].pop_front();
      }//end of a loop
    }//end of gain loop


	  /*********************************
	  *Filling Abort Fraction arrays
	  **********************************/


    for(unsigned int bc_counter = 0; bc_counter < bc_readout; bc_counter++){
      double nROD0abortfraction = 100*((nROD0HitLG[bc_counter]*11)+nROD0HitHG[bc_counter])/36.0;
      double nROD1abortfraction = 100*((nROD1HitLG[bc_counter]*11)+nROD1HitHG[bc_counter])/36.0;
      if (nROD0abortfraction != 0 || nROD1abortfraction !=0) {
	double nabortfraction = (nROD0abortfraction+nROD1abortfraction)/2;
	m_AbortFraction->Fill(nabortfraction);
	m_AbortFractionVsECR->Fill(ecr, nabortfraction);
	m_AbortFractionVsLB->Fill(m_current_LB, nabortfraction);
	if (nROD0BCID[bc_counter]==nROD1BCID[bc_counter]) {
	  m_AbortFractionVsBCID->Fill(nROD0BCID[bc_counter], nabortfraction);
	}
      }
      if (nROD0abortfraction != 0) {
	m_AbortFractionROD0->Fill(nROD0abortfraction);
	m_AbortFractionROD0VsECR->Fill(ecr,nROD0abortfraction);
	m_AbortFractionROD0VsBCID->Fill(nROD0BCID[bc_counter],nROD0abortfraction);
      }
      if (nROD1abortfraction != 0) {
	m_AbortFractionROD1->Fill(nROD1abortfraction);
	m_AbortFractionROD1VsECR->Fill(ecr,nROD1abortfraction);
	m_AbortFractionROD1VsBCID->Fill(nROD1BCID[bc_counter],nROD1abortfraction);
      }
    }
	  
    
  }// end if bcm_rdo not empty 
  return StatusCode::SUCCESS;
}
