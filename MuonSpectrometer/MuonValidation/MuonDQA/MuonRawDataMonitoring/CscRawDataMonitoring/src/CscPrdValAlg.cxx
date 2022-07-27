/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
/*   NAME    : CscPrdValAlg.cxx
 *   PACKAGE : MuonRawDataMonitoring/CscRawDataMonitoring
 *   PURPOSE : CSC PrepRawData (PRD) monitoring
 *   AUTHOR  : Ken Johns, Xiaowen Lei (U.Arizona)
 *
 *   MODIFIED: N. Benekos(Illinois)
 *             V. Kaushik(U.Arizona) : 2009-05-12
 *             J. Veatch(U.Arizona)  : 2012-02-08
 */
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

// Athena include(s)
#include "CscPrdValAlg.h"
#include "MuonPrepRawData/CscStripPrepDataCollection.h"
#include "MuonPrepRawData/CscStripPrepData.h"

// ROOT include(s)
#include "TClass.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace Muon;

namespace CscPrdBins {
  void PrdBinLabels(TH1 *h, int side) {
    h->GetXaxis()->SetTitle("");
    h->GetXaxis()->SetLabelSize(0.03);
    if(side == -1) {
      for(size_t j=6; j<86; j++) {
        if( j%5 != 0 ) {
          float xmid = h->GetBinLowEdge(j) + h->GetBinWidth(j);
          xmid -= 1.0; xmid *= -1.0;
          int seclay = int(ceil(xmid*10)); // get sector/layer
          int sec = seclay/10;
          int lay = (seclay%10)/2 ;
          lay = (4 - lay) + 1;
          h->GetXaxis()->SetBinLabel(j,Form("%c%02d:%d",(sec%2==0?'S':'L'),sec,lay));
        } // end for
      } // end if
    } else if (side == 1) {
      for(size_t j=6; j<86; j++) {
        if( j%5 != 0 ) {
          float xmid = h->GetBinLowEdge(j) + h->GetBinWidth(j);
          int seclay = int(ceil(xmid*10)); // get sector/layer
          int sec = seclay/10;
          int lay = (seclay%10)/2 ;
          h->GetXaxis()->SetBinLabel(j,Form("%c%02d:%d",(sec%2==0?'S':'L'),sec,lay));
        }
      } // end for
    } // end else

  } // BinLabels

} // end namespace CscBins


//
// constructor ----------------------------------------------------------------
//
CscPrdValAlg::CscPrdValAlg(const std::string & type, const std::string & name, const IInterface* parent) : 
  ManagedMonitorToolBase(type, name, parent),
  m_stripFitter(name),
  m_cscprd_oviewEA(nullptr),
  m_cscprd_oviewEC(nullptr)
{

  declareProperty("CSCStripFitter", m_stripFitter);
  declareProperty("CSCPrepRawDataPath", m_cscPRDPath = "Muon/MuonRawDataMonitoring/CSC/PRD");
  declareProperty("NoiseCutADC", m_cscNoiseCut = 50);
  declareProperty("MapYXandRZ", m_mapxyrz = false);

  m_cscGenPath = m_cscPRDPath.substr(0,m_cscPRDPath.find("CSC"));

  initHistograms();

} // end CscPrdValAlg::CscPrdValAlg constructor


//
// destructor ----------------------------------------------------------------
//
CscPrdValAlg::~CscPrdValAlg() {
  if(m_cscprd_oviewEA) {
    delete m_cscprd_oviewEA;
    m_cscprd_oviewEA = nullptr;
  }
  if(m_cscprd_oviewEC) {
    delete m_cscprd_oviewEC;
    m_cscprd_oviewEC = nullptr;
  }
  ATH_MSG_DEBUG( "CscPrdValAlg: in destructor" );
}


//
// initialize ----------------------------------------------------------------
//
StatusCode CscPrdValAlg::initialize() {
  ATH_MSG_INFO( "CscPrdValAlg: in initialize" );
  ATH_CHECK(ManagedMonitorToolBase::initialize());
  ATH_CHECK(m_idHelperSvc.retrieve());
  ATH_CHECK(m_stripFitter.retrieve());
  ATH_MSG_DEBUG( "CscPrdValAlg " << name() << ": retrieved " << m_stripFitter );
  ATH_CHECK(m_cscPrdKey.initialize());
  ATH_CHECK(m_eventInfo.initialize());
  return StatusCode::SUCCESS;
} // end CscPrdValAlg::initialize()


//
// initHistograms ----------------------------------------------------------------
//
void CscPrdValAlg::initHistograms() {

  m_h2csc_prd_hitmap = nullptr;
  m_h2csc_prd_hitmap_signal = nullptr;
  m_h2csc_prd_hitmap_noise = nullptr;

  m_h2csc_prd_hitmap_norm = nullptr;
  m_h2csc_prd_hitmap_norm_signal = nullptr;
  m_h2csc_prd_hitmap_norm_noise = nullptr;

  m_h2csc_prd_hitmap_signal_EA = nullptr;
  m_h1csc_prd_hitmap_signal_EA_count = nullptr;
  m_h1csc_prd_hitmap_signal_EA_occupancy = nullptr;

  m_h2csc_prd_hitmap_norm_signal_EA = nullptr;

  m_h2csc_prd_hitmap_signal_EC = nullptr;
  m_h1csc_prd_hitmap_signal_EC_count = nullptr;
  m_h1csc_prd_hitmap_signal_EC_occupancy = nullptr;

  m_h2csc_prd_hitmap_norm_signal_EC = nullptr;

  m_h2csc_prd_occvslb_EA = nullptr;
  m_h2csc_prd_occvslb_EC = nullptr;

  m_h2csc_prd_eta_vs_phi_hitmap = nullptr;
  m_h2csc_prd_r_vs_z_hitmap = nullptr;
  m_h2csc_prd_y_vs_x_hitmap = nullptr;

  m_h2csc_prd_phicluswidth = nullptr;
  m_h2csc_prd_phicluswidth_signal = nullptr;
  m_h2csc_prd_phicluswidth_noise = nullptr;

  m_h2csc_prd_etacluswidth = nullptr;
  m_h2csc_prd_etacluswidth_signal = nullptr;
  m_h2csc_prd_etacluswidth_noise = nullptr;

  m_h2csc_prd_phicluscount = nullptr;
  m_h2csc_prd_phicluscount_signal = nullptr;
  m_h2csc_prd_phicluscount_noise = nullptr;

  m_h2csc_prd_etacluscount = nullptr;
  m_h2csc_prd_etacluscount_signal = nullptr;
  m_h2csc_prd_etacluscount_noise = nullptr;


  m_h1csc_prd_maxdiffamp = nullptr;

  // Correlation plots
  m_h2csc_prd_eta_vs_phi_cluscount = nullptr;
  m_h2csc_prd_eta_vs_phi_cluscount_signal = nullptr;
  m_h2csc_prd_eta_vs_phi_cluscount_noise = nullptr;

  m_h2csc_prd_eta_vs_phi_cluswidth = nullptr;
  m_h2csc_prd_eta_vs_phi_cluswidth_signal = nullptr;
  m_h2csc_prd_eta_vs_phi_cluswidth_noise = nullptr;

}

//
// bookPrdHistograms ----------------------------------------------------------------
//
void CscPrdValAlg::bookPrdHistograms() {

  m_regHShift.clear();
  m_regHExpert.clear();
  m_regHOviewEA.clear();
  m_regHOviewEC.clear();

  /// ******************** DO NOT MODIFY (begin) *********************************************** ///
  int nxbins  = 242;  // 192 bins for precision, 48 for transverse strips + 2 extra
  float nxmin = -49.; // -1 -> -48 (for transverse)
  float nxmax = 193.; // 1 -> 192 (for precision)

  int nybins  = 175;  // 32 chambers (16 per side x 5 layers per chamber) + 1 set extra
  float nymin = -17.; //
  float nymax = 18.;  //

  int nybinsEA  = 90;  // 16 chambers (x 5 layers per chamber) + 1 set extra
  float nyminEA = 0.;  //
  float nymaxEA = 18.; //

  int nybinsEC  = 85;   // 16 chambers (x 5 layers per chamber) + 1 set extra
  float nyminEC = -17.; //
  float nymaxEC = 0.;   //

  int nxbinsLB  = 2510; // 2500 LB (x 16 chambers x 5 layers per chamber) + 10 extra
  float nxminLB = -10.;   //
  float nxmaxLB = 2500.;  //
  /// ******************** DO NOT MODIFY (end) *********************************************** ///

  // book histograms

  // strip hitmap
  m_h2csc_prd_hitmap = new TH2F("h2csc_prd_hitmap", "Hit Occupancy;channel;[sector] + [0.2 #times layer]",
      nxbins,nxmin,nxmax,nybins,nymin,nymax);
  m_h2csc_prd_hitmap_noise = new TH2F("h2csc_prd_hitmap_noise", "Noise Occupancy;channel;[sector] + [0.2 #times layer]",
      nxbins,nxmin,nxmax,nybins,nymin,nymax);
  m_h2csc_prd_hitmap_signal = new TH2F("h2csc_prd_hitmap_signal", "Signal Occupancy;channel;[sector] + [0.2 #times layer]",
      nxbins,nxmin,nxmax,nybins,nymin,nymax);

  m_h2csc_prd_hitmap_norm = new TH2F("h2csc_prd_hitmap_norm", "Normalized Hit Occupancy;channel;[sector] + [0.2 #times layer]",
      nxbins,nxmin,nxmax,nybins,nymin,nymax);
  m_h2csc_prd_hitmap_norm_noise = new TH2F("h2csc_prd_hitmap_norm_noise", "Normalized Noise Occupancy;channel;[sector] + [0.2 #times layer]",
      nxbins,nxmin,nxmax,nybins,nymin,nymax);
  m_h2csc_prd_hitmap_norm_signal = new TH2F("h2csc_prd_hitmap_norm_signal", "Normalized Signal Occupancy;channel;[sector] + [0.2 #times layer]",
      nxbins,nxmin,nxmax,nybins,nymin,nymax);

  m_h2csc_prd_hitmap_signal_EA = new TH2F("h2csc_prd_hitmap_signal_EA", "EndCap A: Signal Occupancy;channel;[sector] + [0.2 #times layer]",
      nxbins,nxmin,nxmax,nybinsEA,nyminEA,nymaxEA);
  m_h1csc_prd_hitmap_signal_EA_count = new TH1F("h1csc_prd_hitmap_signal_EA_count", "EndCap A: Signal Occupancy;channel;entries/channel;",
      nxbins,nxmin,nxmax);
  m_h1csc_prd_hitmap_signal_EA_occupancy = new TH1F("h1csc_prd_hitmap_signal_EA_occupancy", "EndCap A: Signal Occupancy;[sector] + [0.2 #times layer];entries/layer",
      nybinsEA,nyminEA,nymaxEA);
  CscPrdBins::PrdBinLabels(m_h1csc_prd_hitmap_signal_EA_occupancy,1);

  m_h2csc_prd_hitmap_norm_signal_EA = new TH2F("h2csc_prd_hitmap_norm_signal_EA", "EndCap A: Signal Occupancy;channel;[sector] + [0.2 #times layer]",
      nxbins,nxmin,nxmax,nybinsEA,nyminEA,nymaxEA);

  m_h2csc_prd_hitmap_signal_EC = new TH2F("h2csc_prd_hitmap_signal_EC", "EndCap C: Signal Occupancy;channel;[sector] + [0.2 #times layer]",
      nxbins,nxmin,nxmax,nybinsEC,nyminEC,nymaxEC);
  m_h1csc_prd_hitmap_signal_EC_count = new TH1F("h1csc_prd_hitmap_signal_EC_count", "EndCap C: Signal Occupancy;channel;entries/channel;",
      nxbins,nxmin,nxmax);
  m_h1csc_prd_hitmap_signal_EC_occupancy = new TH1F("h1csc_prd_hitmap_signal_EC_occupancy", "EndCap C: Signal Occupancy;[sector] + [0.2 #times layer];entries/layer",
      nybinsEC,nyminEC,nymaxEC);
  CscPrdBins::PrdBinLabels(m_h1csc_prd_hitmap_signal_EC_occupancy,-1);

  m_h2csc_prd_hitmap_norm_signal_EC = new TH2F("h2csc_prd_hitmap_norm_signal_EC", "EndCap C: Signal Occupancy;channel;[sector] + [0.2 #times layer]",
      nxbins,nxmin,nxmax,nybinsEC,nyminEC,nymaxEC);

  m_h2csc_prd_occvslb_EA = new TH2F("h2csc_prd_occvslb_EA", "EndCap A: Layer Signal Occupancy Per LB;LB;[sector] + [0.2 #times layer]",
      nxbinsLB,nxminLB,nxmaxLB,nybinsEA,nyminEA,nymaxEA);
  m_h2csc_prd_occvslb_EC = new TH2F("h2csc_prd_occvslb_EC", "EndCap C: Layer Signal Occupancy Per LB;LB;[sector] + [0.2 #times layer]",
      nxbinsLB,nxminLB,nxmaxLB,nybinsEC,nyminEC,nymaxEC);

  // PRD eta-cluster width
  m_h2csc_prd_etacluswidth = new TH2F("h2csc_prd_etacluswidth",
      "PRD precision-cluster width;no.of strips;[sector] + [0.2 #times layer]",
      192,0,192,nybins,nymin,nymax);
  m_h2csc_prd_etacluswidth_noise = new TH2F("h2csc_prd_etacluswidth_noise",
      "PRD precision-cluster width;no.of strips;[sector] + [0.2 #times layer]",
      192,0,192,nybins,nymin,nymax);
  m_h2csc_prd_etacluswidth_signal = new TH2F("h2csc_prd_etacluswidth_signal",
      "PRD precision-cluster width;no.of strips;[sector] + [0.2 #times layer]",
      192,0,192,nybins,nymin,nymax);

  // PRD phi-cluster width
  m_h2csc_prd_phicluswidth = new TH2F("h2csc_prd_phicluswidth",
      "PRD precision-cluster width;no.of strips;[sector] + [0.2 #times layer]",
      48,0,48,nybins,nymin,nymax);
  m_h2csc_prd_phicluswidth_noise = new TH2F("h2csc_prd_phicluswidth_noise",
      "PRD precision-cluster width;no.of strips;[sector] + [0.2 #times layer]",
      48,0,48,nybins,nymin,nymax);
  m_h2csc_prd_phicluswidth_signal = new TH2F("h2csc_prd_phicluswidth_signal",
      "PRD precision-cluster width;no.of strips;[sector] + [0.2 #times layer]",
      48,0,48,nybins,nymin,nymax);

  // PRD eta-cluster count
  m_h2csc_prd_etacluscount = new TH2F("h2csc_prd_etacluscount",
      "PRD precision-cluster count;no.of clusters;[sector] + [0.2 #times layer]",
      20,0,20,nybins,nymin,nymax);
  m_h2csc_prd_etacluscount_noise = new TH2F("h2csc_prd_etacluscount_noise",
      "PRD precision-cluster count;no.of clusters;[sector] + [0.2 #times layer]",
      20,0,20,nybins,nymin,nymax);
  m_h2csc_prd_etacluscount_signal = new TH2F("h2csc_prd_etacluscount_signal",
      "PRD precision-cluster count;no.of clusters;[sector] + [0.2 #times layer]",
      20,0,20,nybins,nymin,nymax);

  // PRD phi-cluster count
  m_h2csc_prd_phicluscount = new TH2F("h2csc_prd_phicluscount",
      "PRD transverse-cluster count;no.of clusters;[sector] + [0.2 #times layer]",
      20,0,20,nybins,nymin,nymax);
  m_h2csc_prd_phicluscount_noise = new TH2F("h2csc_prd_phicluscount_noise",
      "PRD transverse-cluster count;no.of clusters;[sector] + [0.2 #times layer]",
      20,0,20,nybins,nymin,nymax);
  m_h2csc_prd_phicluscount_signal = new TH2F("h2csc_prd_phicluscount_signal",
      "PRD transverse-cluster count;no.of clusters;[sector] + [0.2 #times layer]",
      20,0,20,nybins,nymin,nymax);

  // correlation histograms
  m_h2csc_prd_eta_vs_phi_cluscount = new TH2F("h2csc_prd_eta_vs_phi_cluscount",
      "Eta vs. Phi Cluster count correlation;#varphi-cluster count;#eta-cluster count",100,0,100,100,0,100);
  m_h2csc_prd_eta_vs_phi_cluscount_noise = new TH2F("h2csc_prd_eta_vs_phi_cluscount_noise",
      "Eta vs. Phi Noise-Cluster count correlation;#varphi-cluster count;#eta-cluster count",100,0,100,100,0,100);
  m_h2csc_prd_eta_vs_phi_cluscount_signal = new TH2F("h2csc_prd_eta_vs_phi_cluscount_signal",
      "Eta vs. Phi Signal-Cluster count correlation;#varphi-cluster count;#eta-cluster count",100,0,100,100,0,100);

  m_h2csc_prd_eta_vs_phi_cluswidth = new TH2F("h2csc_prd_eta_vs_phi_cluswidth",
      "Eta vs. Phi Cluster width correlation;#varphi-cluster width;#eta-cluster width",100,0,100,100,0,100);

  if(m_mapxyrz) {
    m_h2csc_prd_r_vs_z_hitmap = new TH2F("h2csc_prd_r_vs_z_hitmap",
        "R vs. Z Cluster hitmap;z(mm);R(mm)",200, -10000., 10000., 40, 0., 4000);

    m_h2csc_prd_y_vs_x_hitmap = new TH2F("h2csc_prd_y_vs_x_hitmap",
        "Y vs. X Cluster hitmap;x(mm);y(mm)",100, -5000., 5000., 100, -5000., 5000);
  }

  //m_regHShift.push_back(m_h1csc_prd_maxdiffamp);            // shift
  m_regHExpert.push_back(m_h2csc_prd_hitmap);                 // expert
  m_regHShift.push_back(m_h2csc_prd_hitmap_signal);           // shift, overview (dq-flag)
  m_regHExpert.push_back(m_h2csc_prd_hitmap_noise);           // expert

  m_regHExpert.push_back(m_h2csc_prd_hitmap_norm);                 // expert
  m_regHShift.push_back(m_h2csc_prd_hitmap_norm_signal);           // shift, overview (dq-flag)
  m_regHExpert.push_back(m_h2csc_prd_hitmap_norm_noise);           // expert

  m_regHOviewEA.push_back(m_h2csc_prd_hitmap_signal_EA);
  m_regHOviewEA.push_back(m_h1csc_prd_hitmap_signal_EA_count);
  m_regHOviewEA.push_back(m_h1csc_prd_hitmap_signal_EA_occupancy);

  m_regHOviewEA.push_back(m_h2csc_prd_hitmap_norm_signal_EA);

  m_regHOviewEC.push_back(m_h2csc_prd_hitmap_signal_EC);
  m_regHOviewEC.push_back(m_h1csc_prd_hitmap_signal_EC_count);
  m_regHOviewEC.push_back(m_h1csc_prd_hitmap_signal_EC_occupancy);

  m_regHOviewEC.push_back(m_h2csc_prd_hitmap_norm_signal_EC);

  m_regHOviewEA.push_back(m_h2csc_prd_occvslb_EA);
  m_regHOviewEC.push_back(m_h2csc_prd_occvslb_EC);

  m_regHExpert.push_back(m_h2csc_prd_etacluswidth);           // expert
  m_regHExpert.push_back(m_h2csc_prd_etacluswidth_signal);    // expert
  m_regHExpert.push_back(m_h2csc_prd_etacluswidth_noise);     // expert

  m_regHExpert.push_back(m_h2csc_prd_phicluswidth);           // expert
  m_regHExpert.push_back(m_h2csc_prd_phicluswidth_signal);    // expert
  m_regHExpert.push_back(m_h2csc_prd_phicluswidth_noise);     // expert

  m_regHExpert.push_back(m_h2csc_prd_etacluscount);           // expert
  m_regHShift.push_back(m_h2csc_prd_etacluscount_signal);     // shift
  m_regHExpert.push_back(m_h2csc_prd_etacluscount_noise);     // expert

  m_regHExpert.push_back(m_h2csc_prd_phicluscount);           // expert
  m_regHShift.push_back(m_h2csc_prd_phicluscount_signal);     // shift
  m_regHExpert.push_back(m_h2csc_prd_phicluscount_noise);     // expert

  m_regHExpert.push_back(m_h2csc_prd_eta_vs_phi_cluscount);           // expert
  m_regHExpert.push_back(m_h2csc_prd_eta_vs_phi_cluscount_signal);    // expert
  m_regHExpert.push_back(m_h2csc_prd_eta_vs_phi_cluscount_noise);     // expert

  m_regHExpert.push_back(m_h2csc_prd_eta_vs_phi_cluswidth);       // expert
  if(m_mapxyrz) {
    m_regHShift.push_back(m_h2csc_prd_r_vs_z_hitmap);               // shift
    m_regHShift.push_back(m_h2csc_prd_y_vs_x_hitmap);               // shift
  }

}

//
// bookHistograms ----------------------------------------------------------------
//
StatusCode CscPrdValAlg::bookHistograms() {
  ATH_MSG_DEBUG( "CscPrdValAlg: in bookHistograms" );

  StatusCode sc = StatusCode::SUCCESS;

  bookPrdHistograms();

  //declare a group of histograms
  MonGroup cscprd_shift( this, m_cscPRDPath+"/Shift", run, ATTRIB_MANAGED );
  MonGroup cscprd_expert( this, m_cscPRDPath+"/Expert", run, ATTRIB_MANAGED );

  // register shift histograms with service
  for(size_t j = 0; j < m_regHShift.size(); j++) {
    sc=cscprd_shift.regHist(m_regHShift[j]);
    if(sc.isFailure()) {
      ATH_MSG_ERROR ( "Failed to register shift histogram \"" << m_regHShift[j]->GetName() << "\""  );
      return sc;
    } // end if
  } // end for shift

  // register expert histograms with service
  for(size_t j = 0; j < m_regHExpert.size(); j++) {
    sc=cscprd_expert.regHist(m_regHExpert[j]);
    if(sc.isFailure()) {
      ATH_MSG_ERROR ( "Failed to register expert histogram \"" << m_regHExpert[j]->GetName() << "\""  );
      return sc;
    } // end if
  } // end for expert

  // register overview histograms for EA
  std::vector<TH1 *>::iterator iT;
  m_cscprd_oviewEA = new MonGroup( this, m_cscGenPath+"CSC/Overview/CSCEA/PRD", run, ATTRIB_MANAGED );
  iT = m_regHOviewEA.begin();
  ATH_MSG_DEBUG (  "Found " << m_regHOviewEA.size() << " CSCEA Overview Histograms " );
  for (; iT != m_regHOviewEA.end(); ++iT) {
    ATH_MSG_DEBUG ( "Registering CSCEA Overview Histogram: " << (*iT)->GetName() );
    sc = m_cscprd_oviewEA->regHist(*iT);
    ATH_MSG_DEBUG ( "Successfully Registered CSCEA Overview Histogram: " << (*iT)->GetName() );
    if ( sc.isFailure() ) {
      ATH_MSG_ERROR (  "Cannot register overview histogram for Endcap A: " << (*iT)->GetName() );
      return sc;
    }
  }

  // register overview histograms for EC
  m_cscprd_oviewEC = new MonGroup( this, m_cscGenPath+"CSC/Overview/CSCEC/PRD", run, ATTRIB_MANAGED );
  iT = m_regHOviewEC.begin();
  ATH_MSG_DEBUG (  "Found " << m_regHOviewEC.size() << " CSCEC Overview Histograms " );
  for (; iT != m_regHOviewEC.end(); ++iT) {
    ATH_MSG_DEBUG ( "Registering CSCEC Overview Histogram: " << (*iT)->GetName() );
    sc = m_cscprd_oviewEC->regHist(*iT);
    ATH_MSG_DEBUG ( "Successfully Registered CSCEC Overview Histogram: " << (*iT)->GetName() );
    if ( sc.isFailure() ) {
      ATH_MSG_ERROR (  "Cannot register overview histogram for Endcap C: " << (*iT)->GetName() );
      return sc;
    }
  }


  //} // end if isNewRun

  return sc;
} // End CscPrdValAlg::bookHistograms

//
// fillHistograms ----------------------------------------------------------------
//
StatusCode CscPrdValAlg::fillHistograms()  {

  // Part 1: Get the messaging service, print where you are
  ATH_MSG_DEBUG( "CscPrdValAlg: in fillHistograms" );

  SG::ReadHandle<CscStripPrepDataContainer> CscPRD(m_cscPrdKey);

  //Get lumiblock info
  StatusCode sc = fillLumiBlock();
  if( sc.isFailure() ){
    ATH_MSG_ERROR( "Could Not Get LumiBlock Info" );
    return sc;
  }

  // ==============================================================================
  // Field           Range               Notes
  // ==============================================================================
  // StationName     unsigned integer    maps to "CSS", "CSL", etc.
  // StationEta      [-1,1]              -1 for backward, 1 for forward endcap
  // StationPhi      [1,8]               increases with Phi
  // Technology      [1]                 maps to "CSC"
  // ChamberLayer    [1,2]               increases with |Z|
  // WireLayer       [1,4]               increases with |Z|
  // MeasuresPhi     [0,1]               0 if measures R, 1 if measures Phi
  // Strip           [1,n]               increases with R   for MeasuresPhi=0
  //                                     increases with Phi for MeasuresPhi=1
  // ==============================================================================



  // Begin Event ==================================================
  ATH_MSG_DEBUG ( " BEGIN  EVENT ========================================== "  );
  ATH_MSG_DEBUG(" Size of PRD Container  : " << CscPRD->size());

  for (CscStripPrepDataContainer::const_iterator it = CscPRD->begin(); it != CscPRD->end(); ++it) {
    const CscStripPrepDataCollection *prd = *it;
    ATH_MSG_DEBUG ( " Size of Collection     : " << prd->size()  );
    size_t noStrips = prd->size();  // no. of strips in this cluster = m_stripIds.size()
    size_t nEtaClusWidthCnt[5], nPhiClusWidthCnt[5];    // cluster position in each phi-layer
    int clusCount[33][9], sigclusCount[33][9];
    for(size_t kl = 0; kl < 33; kl++ ) {
      for(size_t km = 0; km < 9; km++ ) {
        if(kl == 0 && km < 5) {
          nEtaClusWidthCnt[km] = 0;
          nPhiClusWidthCnt[km] = 0;
        }
        clusCount[kl][km] = 0;
        sigclusCount[kl][km] = 0;
      } // end loop over km
    } // end loop over kl

    // loop over PRD-clusters
    // Loop over strip id's vector -- this is just one strip even though its a vector of ID's
    ATH_MSG_DEBUG ( " BEGIN Loop over Strips ========================================== "  );
    for (CscStripPrepDataCollection::const_iterator ic = (*it)->begin(); ic != (*it)->end(); ++ic) { // for-loop over PRD collection
      const CscStripPrepData& praw = **ic;

      // Identify the PRD cluster
      Identifier prawId = praw.identify();
      int stationName = m_idHelperSvc->cscIdHelper().stationName(prawId);
      std::string stationString = m_idHelperSvc->cscIdHelper().stationNameString(stationName);
      int chamberType = stationString == "CSS" ? 0 : 1;
      int stationEta  = m_idHelperSvc->cscIdHelper().stationEta(prawId);
      int stationPhi  = m_idHelperSvc->cscIdHelper().stationPhi(prawId);
      int wireLayer   = m_idHelperSvc->cscIdHelper().wireLayer(prawId);
      int measuresPhi = m_idHelperSvc->cscIdHelper().measuresPhi(prawId);
      int stripId     = m_idHelperSvc->cscIdHelper().strip(prawId);

      int sectorNo  = stationEta * (2 * stationPhi - chamberType);

      // compute the indices to store cluster count
      int ns = sectorNo < 0 ? sectorNo*(-1) : sectorNo+16; // [-16 -> -1] shifted to [1 -> 16] and [+1 -> +16] shifted to [+17 -> +32]
      int nl = (measuresPhi ? wireLayer : wireLayer+4);  // [ 1 -> 4] (phi-layers) and [5 -> 8] (eta-layers)

      clusCount[ns][nl]++;

      // indices for ns = [+1 -> +32]; 32 places (index '0' is not counted); allocated 33 places
      // indices for nl = [+1 -> +8]; 8 places (index '0' is not counted); allocated 9 places
      ATH_MSG_DEBUG(" ns = " << ns << "\tm_nl = " << nl << "\tm_sec = " << sectorNo << "\t m_lay= " 
          << wireLayer << "\t strip = " << stripId << "\tmPhi = " << measuresPhi);


      // y-axis fill value
      // sector# +2 layer 1 maps to +2 + 0.2*(1-1) + 0.1 = +2.1
      // sector# +2 layer 2 maps to +2 + 0.2*(2-1) + 0.1 = +2.3
      // sector# +2 layer 3 maps to +2 + 0.2*(3-1) + 0.1 = +2.5
      // sector# +2 layer 4 maps to +2 + 0.2*(4-1) + 0.1 = +2.7
      float secLayer = sectorNo + 0.2 * (wireLayer - 1) + 0.1;
      int xfac = measuresPhi ? -1 : 1;        // [-1 -> -48] / [+1 -> +192]

      // x-axis fill value
      float spid = stripId * xfac;  
      m_h2csc_prd_hitmap->Fill(spid, secLayer);
      measuresPhi ? m_h2csc_prd_phicluswidth->Fill(noStrips,secLayer) : m_h2csc_prd_etacluswidth->Fill(noStrips,secLayer);

      if(m_mapxyrz) {
        float x = praw.globalPosition().x();
        float y = praw.globalPosition().y();
        float z = praw.globalPosition().z();
        float r = sqrt(x*x + y*y);
        ATH_MSG_DEBUG(" prd x = " << x << "\t y = " << y << "\t z = " << z );
        m_h2csc_prd_r_vs_z_hitmap->Fill(z,r);
        m_h2csc_prd_y_vs_x_hitmap->Fill(y,x);
      } // end if(m_mapxyrz)

      // Fit this strip and get Charge (in units of: # of electrons)
      ICscStripFitter::Result res;
      res = m_stripFitter->fit(praw);

      ATH_MSG_DEBUG ( "Strip q +- dq = " << res.charge  << " +- " << res.dcharge << "\t t +- dt = "
          << res.time << " +- " <<  res.dtime << "\t w +- dw = " << res.width << " +- "
          << res.dwidth << "\t status= " << res.status << "\t chisq= " << res.chsq);


      // determine of the cluster is a noise/signal cluster Max_Delta_ADC > NoiseCut
      float kiloele = 1.0e-3; // multiply # of electrons by this number to get kiloElectrons (1 ke = 1 ADC)
      float qstripADC = res.charge * kiloele;

      // By default res.status = -1
      // if strip fit is success res.status = 0
      // If fit fails use the peak sample. In this case res.status = 1

      bool signal = ((qstripADC > m_cscNoiseCut) && (res.status >= 0)) ? true : false;

      // increment the signal-cluster count
      if(signal) {
        sigclusCount[ns][nl]++;
        measuresPhi ? nPhiClusWidthCnt[wireLayer]++ : nEtaClusWidthCnt[wireLayer]++ ;
        m_h2csc_prd_hitmap_signal->Fill(spid,secLayer);
        if(stationEta == -1) {
          m_h2csc_prd_hitmap_signal_EC->Fill(spid,secLayer);
          m_h1csc_prd_hitmap_signal_EC_count->Fill(spid);
          m_h1csc_prd_hitmap_signal_EC_occupancy->Fill(secLayer);
          m_h2csc_prd_occvslb_EC->Fill(m_lumiblock,secLayer);
        } // end if(stationEta == -1)
        else {
          m_h2csc_prd_hitmap_signal_EA->Fill(spid,secLayer);
          m_h1csc_prd_hitmap_signal_EA_count->Fill(spid);
          m_h1csc_prd_hitmap_signal_EA_occupancy->Fill(secLayer);
          m_h2csc_prd_occvslb_EA->Fill(m_lumiblock,secLayer);
        } // end else if(stationEta == -1)
        if(!measuresPhi) {
          m_h2csc_prd_etacluswidth_signal->Fill(noStrips,secLayer);
        } // end if(!measuresPhi)
        else {
          m_h2csc_prd_phicluswidth_signal->Fill(noStrips,secLayer);
        } // end else if(!measuresPhi)
      } // end if(signal)
      else {
        m_h2csc_prd_hitmap_noise->Fill(spid,secLayer);
        if(!measuresPhi) {
          m_h2csc_prd_etacluswidth_noise->Fill(noStrips,secLayer);
        } // end if(!measuresPhi)
        else {
          m_h2csc_prd_phicluswidth_noise->Fill(noStrips,secLayer);
        } // end else if(!measuresPhi)
      } // end else if(signal)


    } // end for-loop over PRD collection
    ATH_MSG_DEBUG ( " End loop over PRD collection======================" );

    for(size_t lcnt = 1; lcnt < 5; lcnt++ ) {
      m_h2csc_prd_eta_vs_phi_cluswidth->Fill(nPhiClusWidthCnt[lcnt],nEtaClusWidthCnt[lcnt]);
    } // end loop over lcnt

    int numeta = 0, numphi = 0;
    int numetasignal = 0, numphisignal = 0;
    for(int kl = 1; kl < 33; kl++ ) {

      for(int km = 1; km < 9; km++ ) {
        int lay = (km > 4 && km < 9) ? km-4 : km;  // 1,2,3,4 (phi-layers)     5-4, 6-4, 7-4, 8-4 (eta-layers)
        bool mphi = (km > 0 && km < 5) ? true : false; // 1,2,3,4 (phi-layers) 5,6,7,8 (eta-layers)
        std::string wlay = mphi ? "Phi-Layer " : "Eta-Layer: ";

        int count = clusCount[kl][km];
        int scount = sigclusCount[kl][km];


        if(count) {
          float secLayer = kl-16 + 0.2 * (lay - 1) + 0.1;
          if(mphi) {
            m_h2csc_prd_phicluscount->Fill(count,secLayer); // all phi-cluster counts
            numphi += count;
            if(scount) {
              numphisignal += scount;
              m_h2csc_prd_phicluscount_signal->Fill(scount,secLayer); // signal phi-cluster count
              m_h2csc_prd_phicluscount_noise->Fill((count-scount),secLayer); // noise phi-cluster count
            } // end if(scount) 
            else {
              m_h2csc_prd_phicluscount_noise->Fill(count,secLayer); // noise phi-cluster count
            } // end else if(scount)
          } // end if(mphi)
          else {
            m_h2csc_prd_etacluscount->Fill(count,secLayer);
            numeta += count;
            if(scount) {
              numetasignal += scount;
              m_h2csc_prd_etacluscount_signal->Fill(scount,secLayer); // signal eta-cluster count
              m_h2csc_prd_etacluscount_noise->Fill((count-scount),secLayer); // noise eta-cluster count
            } // end if(scount)
            else {
              m_h2csc_prd_etacluscount_noise->Fill(count,secLayer); // noise eta-cluster count
            } // end else if(scount)
          } // end else if(mphi)
          ATH_MSG_DEBUG ( wlay << "Counts sec: [" << kl-16 << "]\tlayer: [" << km << "] = " <<
              secLayer << "\t = " << count << "\t" << scount);
        } // end if(count)
      } // end for km
    } // end for kl

    m_h2csc_prd_eta_vs_phi_cluscount->Fill(numphi,numeta);
    m_h2csc_prd_eta_vs_phi_cluscount_signal->Fill(numphisignal,numetasignal);
    m_h2csc_prd_eta_vs_phi_cluscount_noise->Fill(numphi-numphisignal, numeta-numetasignal);

  } // end for-loop over container
  ATH_MSG_DEBUG ( " End EVENT======================" );

  ATH_MSG_DEBUG( "CscPrdValAlg: fillHistograms reports success" );

  return sc;
} // end CscPrdValAlg::fillHistograms()

  //
  // checkHists ----------------------------------------------------------------
  //
  StatusCode CscPrdValAlg::checkHists(bool /* fromFinalize */) {

    ATH_MSG_DEBUG( "CscPrdValAlg: in checkHists" );

    return StatusCode::SUCCESS;
  }

  //
  // fillLumiBlock ----------------------------------------------------------------
  //
  StatusCode CscPrdValAlg::fillLumiBlock(){

    m_lumiblock = -1;

    SG::ReadHandle<xAOD::EventInfo> evt(m_eventInfo);

    m_lumiblock = evt->lumiBlock();

    return StatusCode::SUCCESS;

  }
