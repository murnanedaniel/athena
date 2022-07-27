/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

//*****************************************************************************
//  Filename : TileAANtuple.cxx
//  Author   : Alexander Solodkov
//  Created  : April, 2010
//
//  DESCRIPTION:
//     Algorithm which puts in ntuple all raw data (digits) and results of Opt filter
//
//  HISTORY:
//     29-Apr-2010 First version - simplified version of TileTBAANtuple
//
//  BUGS:
//
//*****************************************************************************

//Tile includes
#include "TileRec/TileAANtuple.h"
#include "TileIdentifier/TileHWID.h"
#include "TileCalibBlobObjs/TileCalibUtils.h"
#include "TileConditions/TileCablingService.h"
#include "TileConditions/ITileBadChanTool.h"
#include "TileConditions/TileInfo.h"
#include "TileDetDescr/TileDetDescrManager.h"
#include "TileConditions/TileCondToolEmscale.h"
#include "TileEvent/TileDigitsContainer.h"
#include "TileEvent/TileBeamElemContainer.h"
#include "TileEvent/TileRawChannelContainer.h"
#include "TileEvent/TileContainer.h"
#include "TileEvent/TileLaserObject.h"
#include "TileEvent/TileMuonReceiverContainer.h"
#include "TileByteStream/TileHid2RESrcID.h"
#include "TileIdentifier/TileTBFrag.h"
#include "TileL2Algs/TileL2Builder.h"

// Calo includes
#include "CaloDetDescr/CaloDetDescrElement.h"
#include "CaloDetDescr/MbtsDetDescrManager.h"
#include "Identifier/IdentifierHash.h"
#include "CaloIdentifier/TileID.h"

//Atlas include
#include "AthenaKernel/errorcheck.h"
#include "xAODEventInfo/EventInfo.h"
#include "StoreGate/ReadHandleKey.h"
#include "ByteStreamCnvSvcBase/ROBDataProviderSvc.h"

#include "eformat/ROBFragment.h"
#include "eformat/FullEventFragment.h"

// Gaudi includes
#include "GaudiKernel/ITHistSvc.h"
#include "GaudiKernel/ThreadLocalContext.h"

#include "TTree.h"
#include "TFile.h"
#include <iomanip>
#include "boost/date_time/local_time/local_time.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"
#include <iostream>
#include <sstream>

#define CLEAR(array)				\
memset(array,0,sizeof(array))

#define CLEAR1(array)				\
memset(array,-1,sizeof(array))

#define CLEAR2(array,size)			\
  memset(array,0,sizeof(array)/size)

#define CLEAR3(array,size)			\
  memset(array,-1,sizeof(array)/size)

// clear for MF arrays with two gains
#define CLEAR4(array,size)                      \
memset(array,0,sizeof(*array)*N_ROS2*N_MODULES*N_CHANS*m_nSamples/size)

// clear for sample arrays with two gains
#define CLEAR5(array,size)                      \
memset(array,-1,sizeof(*array)*N_ROS2*N_MODULES*N_CHANS*m_nSamples/size)

// clear for TMDB sample arrays with one gain
#define CLEAR6(array)                           \
memset(array,0,sizeof(*array)*N_ROS*N_MODULES*N_TMDBCHANS*m_nSamples)

#define NAME1(s1)				\
s1.c_str()

#define NAME2(s1,s2)				\
(s1+s2).c_str()

#define NAME3(s1,s2,s3)				\
(s1+s2+s3).c_str()

#define NAME5(s1,s2,s3,s4,s5)			\
(s1+s2+s3+s4+s5).c_str()

#define sample_ind(r,m,c,i) (((r*N_MODULES + m)*N_CHANS + c)*m_nSamples) + i

#define sample_ind_TMDB(r,m,c,i) (((r*N_MODULES + m)*N_TMDBCHANS + c)*m_nSamples) + i

// Constructor & deconstructor
/** @class TileAANtuple
 *  @brief class to produce TileCal commissioning ntuples
 */


TileAANtuple::TileAANtuple(std::string name, ISvcLocator* pSvcLocator)
: AthAlgorithm(name, pSvcLocator)
, m_evTime(0)
, m_run(0)
, m_evt(0)
, m_lumiBlock(0)
, m_HHMMSS(0)
, m_dateTime()
, m_trigType(0)
, m_dspFlags(0)
, m_l1ID()
, m_l1Type()
, m_evBCID()
, m_evType()
, m_cispar()
, m_las_version(0)
, m_las_BCID(0)
, m_las_Filt(0)
, m_las_ReqAmp(0)
, m_las_MeasAmp(0)
, m_las_Temperature(0)
, m_arrays (std::make_unique<Arrays>())
, m_qdctimeout(0)
, m_tdctimeout(0)
, m_daqtype(0)
, m_nBadDr(0)
, m_nBadHV(0)
, m_nBadDCS(0)
, m_nBadDB(0)
, m_nBadTotal(0)
, m_rchUnit(TileRawChannelUnit::MegaElectronVolts)
, m_dspUnit(TileRawChannelUnit::ADCcounts)
, m_ntuplePtr(0)
, m_DCSntuplePtr(0)
, m_thistSvc("THistSvc", name)
, m_tileID(0)
, m_tileHWID(0)
, m_cabling(0)
, m_tileMgr(0)
, m_tileBadChanTool("TileBadChanTool")
, m_tileToolEmscale("TileCondToolEmscale")
, m_l2Builder()
, m_sumEt_xx()
, m_sumEz_xx()
, m_sumE_xx()
, m_sumEt_yy()
, m_sumEz_yy()
, m_sumE_yy()
, m_sumEt_zz()
, m_sumEz_zz()
, m_sumE_zz()
, m_bad()
{
  declareProperty("TileCondToolEmscale", m_tileToolEmscale);
  declareProperty("TileDigitsContainer", m_digitsContainerKey = "TileDigitsCnt");
  declareProperty("TileDigitsContainerFlt", m_fltDigitsContainerKey = "" /* "TileDigitsFlt" */);
  declareProperty("TileBeamElemContainer", m_beamElemContainerKey = "TileBeamElemCnt");
  declareProperty("TileRawChannelContainer", m_rawChannelContainerKey = "TileRawChannelCnt");
  declareProperty("TileRawChannelContainerFit", m_fitRawChannelContainerKey = "");      //
  declareProperty("TileRawChannelContainerFitCool", m_fitcRawChannelContainerKey = ""); // don't create
  declareProperty("TileRawChannelContainerOpt", m_optRawChannelContainerKey = "");      // by default
  declareProperty("TileRawChannelContainerQIE", m_qieRawChannelContainerKey = "");      // processed QIE data
  declareProperty("TileRawChannelContainerOF1", m_of1RawChannelContainerKey = "");      //
  declareProperty("TileRawChannelContainerDsp", m_dspRawChannelContainerKey = "");      //
  declareProperty("TileRawChannelContainerMF", m_mfRawChannelContainerKey = "");      //
  declareProperty("TileRawChannelContainerWiener", m_wienerRawChannelContainerKey = "");//
  declareProperty("TileMuRcvRawChannelContainer", m_tileMuRcvRawChannelContainerKey = "MuRcvRawChCnt");// TMDB
  declareProperty("TileMuRcvDigitsContainer", m_tileMuRcvDigitsContainerKey = "MuRcvDigitsCnt");// TMDB
  declareProperty("TileMuRcvContainer", m_tileMuRcvContainerKey = "TileMuRcvCnt");// TMDB
  declareProperty("TileLaserObject", m_laserObjectKey = "" /* "TileLaserObj" */);       //
  declareProperty("TileL2Cnt", m_l2CntKey = "TileL2Cnt");
  declareProperty("CalibrateEnergy", m_calibrateEnergy = true);
  declareProperty("UseDspUnits", m_useDspUnits = false);
  declareProperty("OfflineUnits", m_finalUnit = TileRawChannelUnit::MegaElectronVolts);
  declareProperty("CalibMode", m_calibMode = false);
  declareProperty("CompareMode", m_compareMode = false);
  declareProperty("BSInput", m_bsInput = true);
  declareProperty("PMTOrder", m_pmtOrder = false);

  declareProperty("StreamName", m_streamName = "AANT");
  declareProperty("NTupleID", m_ntupleID = "h2000");
  declareProperty("TreeSize", m_treeSize = 16000000000LL);
  
  declareProperty("CheckDCS",m_checkDCS = false);
  declareProperty("DCSBranches",m_DCSBranches = 111111111);

  declareProperty("SkipEvents", m_skipEvents = 0);
  declareProperty("NSamples", m_nSamples=7);
  declareProperty("Reduced", m_reduced=false);
  declareProperty("CompressionSettings", m_compressSettings = -1);

  m_evtNr = -1;
}

TileAANtuple::~TileAANtuple() {
}

/// Alg standard interface function
StatusCode TileAANtuple::initialize() {
  ATH_MSG_INFO( "Initialization started");

  //=== get TileCablingSvc
  ATH_CHECK( m_cablingSvc.retrieve() );

  // find TileCablingService
  m_cabling = TileCablingService::getInstance();

  // retrieve TileDetDescr Manager det store
  ATH_CHECK( detStore()->retrieve(m_tileMgr) );
  
  // retrieve TileID helper from det store
  ATH_CHECK( detStore()->retrieve(m_tileID) );
  ATH_CHECK( detStore()->retrieve(m_tileHWID) );
  
  //=== get TileDCSTool
  if (m_checkDCS) {
    ATH_CHECK( m_tileDCS.retrieve() );
  } else {
    m_tileDCS.disable();
  }
  
  //=== get TileBadChanTool
  ATH_CHECK( m_tileBadChanTool.retrieve() );
  
  //=== get TileCondToolEmscale
  ATH_CHECK( m_tileToolEmscale.retrieve() );
  
  //=== get TileL2Builder
  if (m_compareMode) {
    ATH_CHECK( m_l2Builder.retrieve() );
  }

  ATH_CHECK( m_DQstatusKey.initialize() );

  int sample_size = N_ROS2*N_MODULES*N_CHANS*m_nSamples;
  int sample_TMDB_size = N_ROS*N_MODULES*N_TMDBCHANS*m_nSamples;
  m_arrays->m_sample = (short *) malloc(sample_size*sizeof(short));
  m_arrays->m_sampleFlt = (short *) malloc(sample_size*sizeof(short));
  m_arrays->m_sampleTMDB = (unsigned char *) malloc(sample_TMDB_size*sizeof(unsigned char));

  ATH_CHECK( m_beamElemContainerKey.initialize(m_bsInput) );
  ATH_CHECK( m_digitsContainerKey.initialize(SG::AllowEmpty) );
  ATH_CHECK( m_fltDigitsContainerKey.initialize(SG::AllowEmpty) );
  ATH_CHECK( m_laserObjectKey.initialize(SG::AllowEmpty) );
  ATH_CHECK( m_tileMuRcvContainerKey.initialize(SG::AllowEmpty) );
  ATH_CHECK( m_tileMuRcvDigitsContainerKey.initialize(SG::AllowEmpty) );
  ATH_CHECK( m_tileMuRcvRawChannelContainerKey.initialize(SG::AllowEmpty) );
  ATH_CHECK( m_mfRawChannelContainerKey.initialize(SG::AllowEmpty) );
  ATH_CHECK( m_rawChannelContainerKey.initialize(SG::AllowEmpty) );
  ATH_CHECK( m_fitRawChannelContainerKey.initialize(SG::AllowEmpty) );
  ATH_CHECK( m_fitcRawChannelContainerKey.initialize(SG::AllowEmpty) );
  ATH_CHECK( m_optRawChannelContainerKey.initialize(SG::AllowEmpty) );
  ATH_CHECK( m_qieRawChannelContainerKey.initialize(SG::AllowEmpty) );
  ATH_CHECK( m_dspRawChannelContainerKey.initialize(SG::AllowEmpty) );
  ATH_CHECK( m_of1RawChannelContainerKey.initialize(SG::AllowEmpty) );
  ATH_CHECK( m_wienerRawChannelContainerKey.initialize(SG::AllowEmpty) );
  ATH_CHECK( m_l2CntKey.initialize(m_compareMode) );
  
  ATH_MSG_INFO( "initialization completed" ) ;
  return StatusCode::SUCCESS;
}


StatusCode TileAANtuple::ntuple_initialize(const EventContext& ctx,
                                           const TileDQstatus& DQstatus)
{
  if (m_bsInput) {
    ATH_CHECK( m_robSvc.retrieve() );
    ATH_CHECK( m_decoder.retrieve() );
    const TileHid2RESrcID* hid2re = m_decoder->getHid2re();
    m_ROBID.push_back( hid2re->getRobFromFragID(DIGI_PAR_FRAG) );
    m_ROBID.push_back( hid2re->getRobFromFragID(LASER_OBJ_FRAG) );
  }

  uint32_t calib = DQstatus.calibMode();
  bool calibMode  = (calib == 1);
  if ( calibMode != m_calibMode && calib!=0xFFFFFFFF ) {
    ATH_MSG_INFO( "Calib mode from data is " << calibMode );
    ATH_MSG_INFO( "  Overwriting calib mode " );
    m_calibMode = calibMode;
  }
  
  if (m_finalUnit < TileRawChannelUnit::ADCcounts
      || m_finalUnit > TileRawChannelUnit::OnlineMegaElectronVolts) {
    
    m_finalUnit = -1;
    if ( !m_useDspUnits && m_calibrateEnergy ) {
      m_useDspUnits = true;
      ATH_MSG_INFO( "Final offline units are not set, will use DSP units" );
    }
  }
  
  if ( !m_calibrateEnergy && m_useDspUnits) {
    ATH_MSG_INFO( "calibrateEnergy is disabled, don't want to use DSP units" );
    m_useDspUnits = false;
  }
  
  ATH_MSG_INFO( "calibMode " << m_calibMode );
  ATH_MSG_INFO( "calibrateEnergy " << m_calibrateEnergy );
  ATH_MSG_INFO( "offlineUnits " << m_finalUnit );
  ATH_MSG_INFO( "useDspUnits " << m_useDspUnits );
  
  // set event number to 0 before first event
  m_evtNr = 0;
  
  ATH_CHECK( m_thistSvc.retrieve() );

  if (m_compressSettings >= 0) {
    ATH_CHECK( m_fileMgr.retrieve() );
  }
  
  if(initNTuple(ctx).isFailure()) {
    ATH_MSG_ERROR( " Error during ntuple initialization" );
  }
  
  ATH_MSG_INFO( "ntuple initialization completed" );
  return StatusCode::SUCCESS;
}


StatusCode TileAANtuple::execute() {
  const EventContext& ctx = Gaudi::Hive::currentContext();
  const TileDQstatus* DQstatus = SG::makeHandle (m_DQstatusKey, ctx).get();

  if (m_evtNr < 0) {
    if (ntuple_initialize(ctx, *DQstatus).isFailure()) {
      ATH_MSG_ERROR( "ntuple_initialize failed" );
    }
  }
  
  if (m_evtNr%1000 == 0) {
    ATH_MSG_INFO( m_evtNr << " events processed so far" );
  }
  
  if (ntuple_clear().isFailure()) {
    ATH_MSG_ERROR( "ntuple_clear failed" );
  }
  
  bool empty = true;

  // store BeamElements
  if (!m_beamElemContainerKey.key().empty()) {
    empty &= storeBeamElements(*DQstatus).isFailure();
  }
  
  //store Laser Object
  if (!m_laserObjectKey.empty()) {
    empty &= storeLaser(ctx).isFailure();
  }
  
  // store TileDigits
  empty &= storeDigits(ctx, m_fltDigitsContainerKey,m_arrays->m_sampleFlt,m_arrays->m_gainFlt,false).isFailure();
  empty &= storeDigits(ctx, m_digitsContainerKey,   m_arrays->m_sample,   m_arrays->m_gain,   true ).isFailure();
  
  // store TileRawChannels
  // start from DSP channels - so we can find out what is the DSP units
  empty &= storeRawChannels(ctx,   m_dspRawChannelContainerKey,    m_arrays->m_eDsp,     m_arrays->m_tDsp,     m_arrays->m_chi2Dsp,    m_arrays->m_pedDsp,    true ).isFailure();
  empty &= storeRawChannels(ctx,   m_rawChannelContainerKey,       m_arrays->m_ene,      m_arrays->m_time,     m_arrays->m_chi2,       m_arrays->m_ped,       false).isFailure();
  empty &= storeMFRawChannels(ctx, m_mfRawChannelContainerKey,     m_arrays->m_eMF,      m_arrays->m_tMF,      m_arrays->m_chi2MF,     m_arrays->m_pedMF,     false).isFailure();
  empty &= storeRawChannels(ctx,   m_fitRawChannelContainerKey,    m_arrays->m_eFit,     m_arrays->m_tFit,     m_arrays->m_chi2Fit,    m_arrays->m_pedFit,    false).isFailure();
  empty &= storeRawChannels(ctx,   m_fitcRawChannelContainerKey,   m_arrays->m_eFitc,    m_arrays->m_tFitc,    m_arrays->m_chi2Fitc,   m_arrays->m_pedFitc,   false).isFailure();
  empty &= storeRawChannels(ctx,   m_optRawChannelContainerKey,    m_arrays->m_eOpt,     m_arrays->m_tOpt,     m_arrays->m_chi2Opt,    m_arrays->m_pedOpt,    false).isFailure();
  empty &= storeRawChannels(ctx,   m_qieRawChannelContainerKey,    m_arrays->m_eQIE,     m_arrays->m_tQIE,     m_arrays->m_chi2QIE,    m_arrays->m_pedQIE,    false).isFailure();
  empty &= storeRawChannels(ctx,   m_of1RawChannelContainerKey,    m_arrays->m_eOF1,     m_arrays->m_tOF1,     m_arrays->m_chi2OF1,    m_arrays->m_pedOF1,    false).isFailure();
  empty &= storeRawChannels(ctx,   m_wienerRawChannelContainerKey, m_arrays->m_eWiener,  m_arrays->m_tWiener,  m_arrays->m_chi2Wiener, m_arrays->m_pedWiener, false).isFailure();
  
  // store TMDB data
  //
  empty &= storeTMDBDecision(ctx).isFailure();
  empty &= storeTMDBDigits(ctx).isFailure();
  empty &= storeTMDBRawChannel(ctx).isFailure();

  m_evTime = 0;

  if (m_bsInput) {
    const eformat::FullEventFragment<const uint32_t*>* event = nullptr;
    const eformat::ROBFragment<const uint32_t*>* robFrag = nullptr;
    event = m_robSvc->getEvent();
    std::vector<const ROBDataProviderSvc::ROBF*> robf;
    // keep pointer to whole event and to CIS PAR frag internally
    m_robSvc->getROBData(m_ROBID, robf);
    robFrag = (robf.size() > 0 ) ? robf[0] : nullptr;
    if (event) {
      m_evTime = event->bc_time_seconds();
      if ( robFrag ) {
        // Store ROD header info from collection
        int rod = N_RODS-1;
        m_l1ID[rod]   = robFrag->rod_lvl1_id();
        m_l1Type[rod] = robFrag->rod_lvl1_trigger_type();
        m_evType[rod] = robFrag->rod_detev_type();
        m_evBCID[rod] = robFrag->rod_bc_id();
        if (m_trigType == 0) m_trigType = -m_l1Type[rod]; // make negative to distinguish from TileCal internal trig types
      }
    }
  }

  m_lumiBlock = -1; // placeholder
  
  //Get run and event numbers
  m_run = ctx.eventID().run_number();
  m_evt = ctx.eventID().event_number();
    
  if ( ctx.eventID().lumi_block() ){
    m_lumiBlock = ctx.eventID().lumi_block();
  }
    
  //Get timestamp of the event
  if (ctx.eventID().time_stamp() > 0) {
    m_evTime = ctx.eventID().time_stamp();
  }
  
  if (m_evTime>0) {
    using namespace boost::local_time;
    using namespace boost::posix_time;
    
    /*
     // just an example how to read file with time zones
     tz_database tz_db;
     try {
     tz_db.load_from_file("../data/date_time_zonespec.csv");
     time_zone_ptr gva_tz = tz_db.time_zone_from_region("Europe/Zurich");
     }catch(data_not_accessible dna) {
     std::cerr << "Error with time zone data file: " << dna.what() << std::endl;
     //exit(EXIT_FAILURE);
     }catch(bad_field_count bfc) {
     std::cerr << "Error with time zone data file: " << bfc.what() << std::endl;
     //exit(EXIT_FAILURE);
     }
     */
    //"Europe/Zurich","CET","CET","CEST","CEST","+01:00:00","+01:00:00","-1;0;3","+02:00:00","-1;0;10","+03:00:00"
    static const time_zone_ptr gva_tz(new posix_time_zone((std::string)"CET+01CEST01:00:00,M3.5.0/02:00:00,M10.5.0/03:00:00"));
    local_date_time gva_time(from_time_t(m_evTime),gva_tz);
    
    //std::ostringstream otime;
    //otime << gva_time; // time in the format YYYY-MMM-DD HH:MM:SS TZ
    //strncpy(m_dateTime,otime.str().c_str(),31);
    
    //time_duration hms(gva_time.time_of_day()); - will give time of the day in GMT
    //m_HHMMSS = hms.hours()*10000+hms.minutes()*100+hms.seconds();
    
    struct tm gva_tm(to_tm(gva_time));
    strftime(m_dateTime, 32, "%Y-%b-%d %H:%M:%S %Z", &gva_tm);
    m_HHMMSS = gva_tm.tm_hour*10000+gva_tm.tm_min*100+gva_tm.tm_sec;
    
    // the code below is only valid when running at CERN (in Geneva time zone)
    //struct tm *time = localtime((time_t*)(&m_evTime));
    //m_HHMMSS = time->tm_hour*10000+time->tm_min*100+time->tm_sec;
    //strftime(m_dateTime, 32, "%Y-%m-%d %H:%M:%S %Z", time);
    
  } else {
    m_HHMMSS = -1;
    m_dateTime[0] = '\0'; // empty string
  }
  
  // store DCS data
  if (m_checkDCS) {
    empty &= storeDCS().isFailure();
  }
 
  if (empty) {
    ATH_MSG_WARNING( "Some problems in execute - ntuple was not filled at all" );
  }
  
  // FIRST 4 EVENTS ARE SKIPPED TO RETRIEVE LASER PEDESTALS
  if (m_ntuplePtr && m_evtNr >= m_skipEvents){
    m_ntuplePtr->Fill();
  } // IF
  
  ++m_evtNr;
  
  // Execution completed.
  
  ATH_MSG_DEBUG( "execute() completed successfully" );
  return StatusCode::SUCCESS;
}


//
// Here the LASER object is opened and corresponding variable are stored
//
StatusCode TileAANtuple::storeLaser (const EventContext& ctx) {
  
  ATH_MSG_DEBUG("TileAANtuple::storeLaser()");
  const char* gainnames[2]  = {"LG","HG"};
  
  const TileLaserObject* laserObj = SG::makeHandle(m_laserObjectKey, ctx).get();
  
  m_las_BCID = laserObj->getBCID();
  
  m_las_Filt = laserObj->getFiltNumber();
  m_las_ReqAmp = laserObj->getDiodeCurrOrd();
  m_las_MeasAmp = laserObj->getDiodeCurrMeas();
  m_las_Temperature = laserObj->getPumpDiodeTemp();
  ATH_MSG_VERBOSE( "Laser BCID " << m_las_BCID
                   << " Filt " << m_las_Filt
                   << " ReqAmp " << m_las_ReqAmp
                   << " MeasAmp " << m_las_MeasAmp
                   << " Temp " << m_las_Temperature );
  
  ATH_MSG_DEBUG("LASER"<<(laserObj->isLASERII()?"II":"I")<<" VERSION IS " << laserObj->getVersion());
  
  if(laserObj->isLASERII()){
    m_qdctimeout = laserObj->getQDCTimeout();
    m_tdctimeout = laserObj->getTDCTimeout();
    m_daqtype = laserObj->getDaqType();
    if (msgLvl(MSG::DEBUG)) {
      msg(MSG::DEBUG) << "DAQ Type    " << m_daqtype << endmsg;
      msg(MSG::DEBUG) << "QDC TimeOut " << m_qdctimeout << endmsg;
      msg(MSG::DEBUG) << "TDC TimeOut " << m_tdctimeout << endmsg;
    }

    // RETRIEVE SIGNAL IN ADC COUNTS
    for(int chan=0;chan<28;++chan){
      int ch=chan>>1;
      int gn=chan&1;
      // MONITORING DIODES
      m_arrays->m_chan[chan] = laserObj->getDiodeADC(ch,gn);
      ATH_MSG_DEBUG("LASERII CHANNEL " << ch << " ("<<gainnames[gn]<<") " << m_arrays->m_chan[chan]);
    } // FOR

    for(int chan=28;chan<32;++chan){
      int ch=(chan-28)>>1;
      int gn=chan&1;
      // MONITORING PMTS
      m_arrays->m_chan[chan] = laserObj->getPMADC(ch,gn);
      ATH_MSG_DEBUG("LASERII PMT " << ch << " ("<<gainnames[gn]<<") " << m_arrays->m_chan[chan]);
    } // FOR
    
    // RETRIEVE PEDESTALS IF NOT ALREADY SET
    for(int chan=0;chan<32;++chan){
      int ch=chan>>1;
      int gn=chan&1;
      if(laserObj->isSet(ch, gn, 0) && laserObj->getMean (ch,gn,0)>0) m_arrays->m_chan_Ped[chan]    = laserObj->getMean (ch,gn,0);
      if(laserObj->isSet(ch, gn, 2) && laserObj->getMean (ch,gn,2)>0) m_arrays->m_chan_Led[chan]    = laserObj->getMean (ch,gn,2);
      if(laserObj->isSet(ch, gn, 3) && laserObj->getMean (ch,gn,3)>0) m_arrays->m_chan_Alpha[chan]  = laserObj->getMean (ch,gn,3);
      if(laserObj->isSet(ch, gn, 1) && laserObj->getMean (ch,gn,1)>0) m_arrays->m_chan_Lin[chan]    = laserObj->getMean (ch,gn,1);
      if(laserObj->isSet(ch, gn, 0) && laserObj->getSigma(ch,gn,0)>0) m_arrays->m_chan_SPed[chan]   = laserObj->getSigma(ch,gn,0);
      if(laserObj->isSet(ch, gn, 2) && laserObj->getSigma(ch,gn,2)>0) m_arrays->m_chan_SLed[chan]   = laserObj->getSigma(ch,gn,2);
      if(laserObj->isSet(ch, gn, 3) && laserObj->getSigma(ch,gn,3)>0) m_arrays->m_chan_SAlpha[chan] = laserObj->getSigma(ch,gn,3);
      if(laserObj->isSet(ch, gn, 1) && laserObj->getSigma(ch,gn,1)>0) m_arrays->m_chan_SLin[chan]   = laserObj->getSigma(ch,gn,1);
      
      // DEBUG OUTPUT
      if (msgLvl(MSG::DEBUG)) {
        msg(MSG::DEBUG) << gainnames[gn] << " CHAN " << ch << " SIG= " << m_arrays->m_chan[chan] << endmsg;
        msg(MSG::DEBUG) << gainnames[gn] << " CHAN " << ch << " PED= " << m_arrays->m_chan_Ped[chan]   << "+/-" << m_arrays->m_chan_SPed[chan]   << " ( " << laserObj->isSet(ch, gn, 0) << " ) " << endmsg;
        msg(MSG::DEBUG) << gainnames[gn] << " CHAN " << ch << " PED= " << m_arrays->m_chan_Lin[chan]   << "+/-" << m_arrays->m_chan_SLin[chan]   << " ( " << laserObj->isSet(ch, gn, 1) << " ) " << endmsg;
        msg(MSG::DEBUG) << gainnames[gn] << " CHAN " << ch << " LED= " << m_arrays->m_chan_Led[chan]   << "+/-" << m_arrays->m_chan_SLed[chan]   << " ( " << laserObj->isSet(ch, gn, 2) << " ) " << endmsg;
        msg(MSG::DEBUG) << gainnames[gn] << " CHAN " << ch << " ALP= " << m_arrays->m_chan_Alpha[chan] << "+/-" << m_arrays->m_chan_SAlpha[chan] << " ( " << laserObj->isSet(ch, gn, 3) << " ) " << endmsg;
      } // IF
    } // FOR
  } // IF
  else{
    for (unsigned int gn=0; gn<TileLaserObject::nbGains; ++gn) {
      for (unsigned int i=0; i<TileLaserObject::nbPmts; ++i) {
        m_arrays->m_las_PMT_ADC[gn][i] = laserObj->getPMADC(i,gn);
        m_arrays->m_las_PMT_TDC[gn][i] = laserObj->getTDC(i,gn);
        m_arrays->m_las_PMT_Ped[gn][i] = laserObj->getPMPedestal(i,gn);
        m_arrays->m_las_PMT_Ped_RMS[gn][i] = laserObj->getPMSigmaPedestal(i,gn);
        ATH_MSG_VERBOSE( "LasPMT" << i << " g " << gn
                         << " adc " << m_arrays->m_las_PMT_ADC[gn][i]
                         << " ped " << m_arrays->m_las_PMT_Ped[gn][i]
                         << " rms " << m_arrays->m_las_PMT_Ped_RMS[gn][i]
                         << " tdc " << m_arrays->m_las_PMT_TDC[gn][i] );
      } // FOR
      
      for (unsigned int i=0; i<14; ++i) {
        m_arrays->m_las_D_ADC[gn][i] = laserObj->getDiodeADC(i,gn);
        m_arrays->m_las_D_Ped[gn][i] = laserObj->getDiodePedestal(i,gn);
        m_arrays->m_las_D_Ped_RMS[gn][i] = laserObj->getDiodeSigmaPedestal(i,gn);
        m_arrays->m_las_D_Alpha[gn][i] = laserObj->getAlpha(i,gn);
        m_arrays->m_las_D_Alpha_RMS[gn][i] = laserObj->getSigmaAlpha(i,gn);
        m_arrays->m_las_D_AlphaPed[gn][i] = laserObj->getPedestalAlpha(i,gn);
        m_arrays->m_las_D_AlphaPed_RMS[gn][i] = laserObj->getSigmaPedAlpha(i,gn);
        
        ATH_MSG_VERBOSE( "LasD" << i << " g " << gn
                         << " adc " << m_arrays->m_las_D_ADC[gn][i]
                         << " ped " << m_arrays->m_las_D_Ped[gn][i]
                         << " rms " << m_arrays->m_las_D_Ped_RMS[gn][i]
                         << " alp " << m_arrays->m_las_D_Alpha[gn][i]
                         << " rms " << m_arrays->m_las_D_Alpha_RMS[gn][i]
                         << " ape " << m_arrays->m_las_D_AlphaPed[gn][i]
                         << " rms " << m_arrays->m_las_D_AlphaPed_RMS[gn][i] );
      } // FOR
    } // FOR
  } // ELSE
  
  return StatusCode::SUCCESS;
}

StatusCode TileAANtuple::storeBeamElements(const TileDQstatus& DQstatus) {
  
  const uint32_t* cispar = DQstatus.cispar();
  
  uint32_t oldval = 0;
  int last = 0;
  for(int i = 0; i< N_CISPAR; ++i) {
    m_cispar[i] = cispar[i];
    if (msgLvl(MSG::VERBOSE)) {
      if (oldval != cispar[i]) {
        if (last < i-1) {
          ATH_MSG_VERBOSE( "cispar[" << last << ".." << i-1 << "] = "
                           << oldval  );
        } else if (last == i-1) {
          ATH_MSG_VERBOSE( "cispar[" << last << "] = " << oldval  );
        }
        last = i;
        oldval = cispar[i];
      }
    }
  }
  
  if (msgLvl(MSG::VERBOSE)) {
    if (last < N_CISPAR-1) {
      ATH_MSG_VERBOSE( "cispar[" << last << ".." << N_CISPAR-1 << "] = "
                       << oldval  );
    } else {
      ATH_MSG_VERBOSE( "cispar[" << last << "] = "<< oldval  );
    }
  }
  
  m_trigType = cispar[12];
  
  return StatusCode::SUCCESS;
}


/**
 /// Fill ntuple with data from TRC.
 */
StatusCode
TileAANtuple::storeRawChannels(const EventContext& ctx
                               , const SG::ReadHandleKey<TileRawChannelContainer>& containerKey
                               , float ene[N_ROS2][N_MODULES][N_CHANS]
                               , float time[N_ROS2][N_MODULES][N_CHANS]
                               , float chi2[N_ROS2][N_MODULES][N_CHANS]
                               , float ped[N_ROS2][N_MODULES][N_CHANS]
                               , bool fillAll)
{
  if (containerKey.empty()) {// empty name, nothing to do
    return StatusCode::FAILURE;
  }
  
  // get named container
  const TileRawChannelContainer* rcCnt =
    SG::makeHandle (containerKey, ctx).get();
  ATH_MSG_VERBOSE( "Container ID " << containerKey.key() );
  
  TileRawChannelUnit::UNIT rChUnit = rcCnt->get_unit();
  ATH_MSG_VERBOSE( "RawChannel unit is " << rChUnit );
  
  bool dspCont = ( rChUnit >= TileRawChannelUnit::OnlineADCcounts );
  if (dspCont) { // this is container with DSP results
    m_dspUnit = rChUnit;
    m_dspFlags = rcCnt->get_bsflags() >> 16;
    ATH_MSG_VERBOSE( "DSP flag is 0x" << MSG::hex << m_dspFlags << MSG::dec
                    << " DSP unit is " << m_dspUnit);
    
  } else if ((m_useDspUnits || m_finalUnit >= TileRawChannelUnit::OnlineADCcounts)
             && rChUnit != TileRawChannelUnit::ADCcounts) {
    ATH_MSG_ERROR( "RawChannel units are not ADC counts, can't apply DSP-like calibration" );
    return StatusCode::FAILURE;
  }
  
  if (m_calibrateEnergy) {
    if (m_useDspUnits) { // calibrate a-la online
      m_rchUnit = m_dspUnit;
    } else { // convert to final units
      m_rchUnit = (TileRawChannelUnit::UNIT)m_finalUnit;
    }
  } else {
    m_rchUnit = rChUnit;
  }
  ATH_MSG_VERBOSE( "Final RawChannel unit is " << m_rchUnit );
  
  std::vector<float> sumE(3);
  float E[48];
  int gain[48];
  if (m_compareMode && dspCont) memset(m_bad,0,sizeof(m_bad));
  
  // Get iterator for all TRCColl in TRCCont
  TileRawChannelContainer::const_iterator itColl = (*rcCnt).begin();
  TileRawChannelContainer::const_iterator itCollEnd = (*rcCnt).end();
  
  TileRawChannelCollection::const_iterator it, itEnd;
  
  // Go through all TileRawChannelCollections
  for(; itColl != itCollEnd; ++itColl) {
    int fragId = (*itColl)->identify();
    int drawerIdx = TileCalibUtils::getDrawerIdxFromFragId(fragId);
    int drawer = fragId & 0x3F;
    int ROS = (fragId>>8);
    int rosI = ROS-1;
    int rosL = rosI;
    int rosH = rosI + N_ROS;
    
    ATH_MSG_VERBOSE( "TRC ("<< containerKey.key()
                    <<") Event# "<< m_evtNr
                    << " Frag id 0x" << MSG::hex << fragId << MSG::dec
                    << " ROS " << ROS
                    << " drawer " << drawer );
    
    // go through all TileRawChannels in collection
    it = (*itColl)->begin();
    itEnd = (*itColl)->end();
    
    int cmpCounter = 0;
    if (m_compareMode) {
      memset(E, 0, sizeof(E));
      memset(gain, 0, sizeof(gain));
    }
    
    for(; it != itEnd; ++it) {
      const TileRawChannel* rch = (*it);
      
      HWIdentifier hwid = rch->adc_HWID();
      
      // determine channel
      int channel = m_tileHWID->channel(hwid);
      // convert channel number to PMT number if needed
      if (m_pmtOrder) channel = digiChannel2PMT(ROS,channel);
      
      // determine gain and set ros index accordingly
      int adc = m_tileHWID->adc(hwid);
      if (m_calibMode) {
        if (m_compareMode) {
          ++cmpCounter;
          if(cmpCounter>48) rosI = rosH;
          else              rosI = rosL;
        } else {
          if(adc == 1) rosI = rosH;
          else         rosI = rosL;
        }
      }
      
      /// final calibration
      float energy =  rch->amplitude();
      if (m_rchUnit != rChUnit) {
        if (m_rchUnit < TileRawChannelUnit::OnlineOffset)
          energy = m_tileToolEmscale->channelCalib(drawerIdx, channel, adc, energy, rChUnit, m_rchUnit);
        else
          energy = m_tileToolEmscale->channelCalibOnl(drawerIdx, channel, adc, energy, m_rchUnit);
      }
      
      ene[rosI][drawer][channel] = energy;
      time[rosI][drawer][channel] = rch->time();
      chi2[rosI][drawer][channel] = rch->quality();
      ped[rosI][drawer][channel] = rch->pedestal();
      if (m_arrays->m_gain[rosI][drawer][channel] < 0)
        m_arrays->m_gain[rosI][drawer][channel] = adc;
      
      if (m_compareMode) { // filling array for SumEt calculations
        E[channel] = energy;
        gain[channel] = adc;
        if (dspCont) { // use bad flag from DSP container only
          m_bad[rosL][drawer][channel] = (rch->quality()>15.99);
          //} else {
          //m_bad[rosL][drawer][channel] = m_tileBadChanTool->getAdcStatus(drawerIdx, channel, adc).isBad();
        }
      }
      
      if (msgLvl(MSG::VERBOSE)) {
        int index,pmt;
        rch->cell_ID_index(index,pmt);
        ATH_MSG_VERBOSE( "TRC ch " << channel
                         << " gain " << adc
                         << " type " << std::min(index,0)
                         << " ene=" << energy
                         << " time=" << rch->time()
                         << " chi2=" << rch->quality()
                         << " ped=" << rch->pedestal()  );
      }
    }
    
    if (fillAll) {
      
      m_arrays->m_ROD_GlobalCRC[rosL][drawer] = (*itColl)->getFragGlobalCRC() & 1;
      m_arrays->m_ROD_BCID[rosL][drawer] = (*itColl)->getFragDSPBCID();
      m_arrays->m_ROD_DMUMask[rosL][drawer][0] = (*itColl)->getFragRODChipMask();
      m_arrays->m_ROD_DMUMask[rosL][drawer][1] = (*itColl)->getFragFEChipMask();
      
      for(unsigned int dmu=0;dmu<N_DMUS;dmu++) {
        
        m_arrays->m_ROD_DMUBCIDErr[rosL][drawer][dmu] = ((*itColl)->getFragBCID() >> dmu) & 1;
        m_arrays->m_ROD_DMUmemoryErr[rosL][drawer][dmu] = ((*itColl)->getFragMemoryPar() >> dmu) & 1;
        m_arrays->m_ROD_DMUSstrobeErr[rosL][drawer][dmu] = ((*itColl)->getFragSstrobe() >> dmu) & 1;
        m_arrays->m_ROD_DMUDstrobeErr[rosL][drawer][dmu]    = ((*itColl)->getFragDstrobe() >> dmu) & 1;
        m_arrays->m_ROD_DMUHeadformatErr[rosL][drawer][dmu] = ((*itColl)->getFragHeaderBit() >> dmu) & 1;
        m_arrays->m_ROD_DMUHeadparityErr[rosL][drawer][dmu] = ((*itColl)->getFragHeaderPar() >> dmu) & 1;
        m_arrays->m_ROD_DMUDataformatErr[rosL][drawer][dmu] = ((*itColl)->getFragSampleBit() >> dmu) & 1;
        m_arrays->m_ROD_DMUDataparityErr[rosL][drawer][dmu] = ((*itColl)->getFragSamplePar() >> dmu) & 1;
        m_arrays->m_ROD_DMUfeCRC[rosL][drawer][dmu] = ((*itColl)->getFragFEChipMask() >> dmu) & 1;
        m_arrays->m_ROD_DMUrodCRC[rosL][drawer][dmu] = ((*itColl)->getFragRODChipMask() >> dmu) & 1;
      }
    }
    
    if (m_compareMode) {
      m_l2Builder->SumE(ROS,drawer,m_rchUnit,E,gain,m_bad[rosL][drawer],sumE);
      if (dspCont) {
        m_sumEt_xx[m_l2Builder->idToIndex(fragId)] = sumE[0];
        m_sumEz_xx[m_l2Builder->idToIndex(fragId)] = sumE[1];
        m_sumE_xx[m_l2Builder->idToIndex(fragId)] = sumE[2];
      }
      else {
        m_sumEt_zz[m_l2Builder->idToIndex(fragId)] = sumE[0];
        m_sumEz_zz[m_l2Builder->idToIndex(fragId)] = sumE[1];
        m_sumE_zz[m_l2Builder->idToIndex(fragId)] = sumE[2];
      }
    }
  }
  
  if (m_compareMode && dspCont) {
    
    const TileL2Container* l2Cnt = SG::makeHandle(m_l2CntKey, ctx).get();
    
    TileL2Container::const_iterator it = l2Cnt->begin();
    TileL2Container::const_iterator end= l2Cnt->end();
    int i=0;
    for(; it != end; ++it) {
      m_sumEt_yy[i++] = (*it)->sumEt();
      m_sumEz_yy[i++] = (*it)->sumEz();
      m_sumE_yy[i++]  = (*it)->sumE();
    }
  }
  
  return StatusCode::SUCCESS;
}

StatusCode
TileAANtuple::storeMFRawChannels(const EventContext& ctx
                                 , const SG::ReadHandleKey<TileRawChannelContainer>& containerKey
                                 , float * ene
                                 , float * time
                                 , float chi2[N_ROS2][N_MODULES][N_CHANS]
                                 , float ped[N_ROS2][N_MODULES][N_CHANS]
                                 , bool fillAll)
{
  if (containerKey.empty()) {// empty name, nothing to do
    return StatusCode::FAILURE;
  }

  // get named container
  const TileRawChannelContainer* rcCnt = \
    SG::makeHandle (containerKey, ctx).get();
  
  TileRawChannelUnit::UNIT rChUnit = rcCnt->get_unit();
  ATH_MSG_VERBOSE( "RawChannel unit is " << rChUnit );
  
  bool dspCont = ( rChUnit >= TileRawChannelUnit::OnlineADCcounts );
  if (dspCont) { // this is container with DSP results
    m_dspUnit = rChUnit;
    m_dspFlags = rcCnt->get_bsflags() >> 16;
    ATH_MSG_VERBOSE( "DSP flag is 0x" << MSG::hex << m_dspFlags << MSG::dec
                    << " DSP unit is " << m_dspUnit);
    
  } else if ((m_useDspUnits || m_finalUnit >= TileRawChannelUnit::OnlineADCcounts)
             && rChUnit != TileRawChannelUnit::ADCcounts) {
    ATH_MSG_ERROR( "RawChannel units are not ADC counts, can't apply DSP-like calibration" );
    return StatusCode::FAILURE;
  }
  
  if (m_calibrateEnergy) {
    if (m_useDspUnits) { // calibrate a-la online
      m_rchUnit = m_dspUnit;
    } else { // convert to final units
      m_rchUnit = (TileRawChannelUnit::UNIT)m_finalUnit;
    }
  } else {
    m_rchUnit = rChUnit;
  }
  ATH_MSG_VERBOSE( "Final RawChannel unit is " << m_rchUnit );
  
  std::vector<float> sumE(3);
  float E[48];
  int gain[48];
  if (m_compareMode && dspCont) memset(m_bad, 0, sizeof(m_bad));
  
  // Get iterator for all TRCColl in TRCCont
  TileRawChannelContainer::const_iterator itColl = (*rcCnt).begin();
  TileRawChannelContainer::const_iterator itCollEnd = (*rcCnt).end();
  
  TileRawChannelCollection::const_iterator it, itEnd;
  
  // Go through all TileRawChannelCollections
  for (; itColl != itCollEnd; ++itColl) {
    int fragId = (*itColl)->identify();
    int drawerIdx = TileCalibUtils::getDrawerIdxFromFragId(fragId);
    int drawer = fragId & 0x3F;
    int ROS = (fragId >> 8);
    int rosI = ROS - 1;
    int rosL = rosI;
    int rosH = rosI + N_ROS;
    
    ATH_MSG_VERBOSE( "TRC ("<< containerKey.key()
                    <<") Event# "<< m_evtNr
                    << " Frag id 0x" << MSG::hex << fragId << MSG::dec
                    << " ROS " << ROS
                    << " drawer " << drawer );
    
    // go through all TileRawChannels in collection
    it = (*itColl)->begin();
    itEnd = (*itColl)->end();
    
    int cmpCounter = 0;
    if (m_compareMode) {
      memset(E, 0, sizeof(E));
      memset(gain, 0, sizeof(gain));
    }
    for(; it != itEnd; ++it) {
      const TileRawChannel* rch = (*it);
      
      HWIdentifier hwid=rch->adc_HWID();
      
      // determine channel
      int channel = m_tileHWID->channel(hwid);
      // convert channel number to PMT number if needed
      if (m_pmtOrder) channel = digiChannel2PMT(ROS,channel);
      
      // determine gain and set ros index accordingly
      int adc = m_tileHWID->adc(hwid);
      if (m_calibMode) {
        if (m_compareMode) {
          ++cmpCounter;
          if(cmpCounter>48) rosI = rosH;
          else              rosI = rosL;
        } else {
          if(adc == 1) rosI = rosH;
          else         rosI = rosL;
        }
      }

      /// final calibration
      float energy = 0.;
      for (int i = 0; i < 7; ++i) {
        energy = rch->amplitude(i);
        if (m_rchUnit != rChUnit) {
          if (m_rchUnit < TileRawChannelUnit::OnlineOffset)
            energy = m_tileToolEmscale->channelCalib(drawerIdx, channel, adc, energy, rChUnit, m_rchUnit);
          else
            energy = m_tileToolEmscale->channelCalibOnl(drawerIdx, channel, adc, energy, m_rchUnit);
        }
	      	
        ene[sample_ind(rosI,drawer,channel,i)] = energy;
        time[sample_ind(rosI,drawer,channel,i)] = rch->time(i);
      }
      chi2[rosI][drawer][channel] = rch->quality();
      ped[rosI][drawer][channel] = rch->pedestal();
      
      if (m_arrays->m_gain[rosI][drawer][channel] < 0)
        m_arrays->m_gain[rosI][drawer][channel] = adc;
      
      if (m_compareMode) { // filling array for SumEt calculations
        E[channel] = ene[sample_ind(rosI,drawer,channel,0)];
        gain[channel] = adc;
        if (dspCont) { // use bad flag from DSP container only
          m_bad[rosL][drawer][channel] = (rch->quality() > 15.99);
          //} else {
          //m_bad[rosL][drawer][channel] = m_tileBadChanTool->getAdcStatus(drawerIdx, channel, adc).isBad();
        }
      }
      
      if (msgLvl(MSG::VERBOSE)) {
        int index,pmt;
        rch->cell_ID_index(index,pmt);
        ATH_MSG_VERBOSE( "TRC ch " << channel
                         << " gain " << adc
                         << " type " << std::min(index,0)
                         << " ene=" << ene[sample_ind(rosI,drawer,channel,0)]
                         << " time=" << rch->time()
                         << " chi2=" << rch->quality()
                         << " ped=" << rch->pedestal()  );
      }
    }
    
    if (fillAll) {
      
      m_arrays->m_ROD_GlobalCRC[rosL][drawer] = (*itColl)->getFragGlobalCRC() & 1;
      m_arrays->m_ROD_BCID[rosL][drawer] = (*itColl)->getFragDSPBCID();
      m_arrays->m_ROD_DMUMask[rosL][drawer][0] = (*itColl)->getFragRODChipMask();
      m_arrays->m_ROD_DMUMask[rosL][drawer][1] = (*itColl)->getFragFEChipMask();
      
      for(unsigned int dmu=0;dmu<N_DMUS;dmu++) {
        
        m_arrays->m_ROD_DMUBCIDErr[rosL][drawer][dmu] = ((*itColl)->getFragBCID() >> dmu) & 1;
        m_arrays->m_ROD_DMUmemoryErr[rosL][drawer][dmu] = ((*itColl)->getFragMemoryPar() >> dmu) & 1;
        m_arrays->m_ROD_DMUSstrobeErr[rosL][drawer][dmu] = ((*itColl)->getFragSstrobe() >> dmu) & 1;
        m_arrays->m_ROD_DMUDstrobeErr[rosL][drawer][dmu]    = ((*itColl)->getFragDstrobe() >> dmu) & 1;
        m_arrays->m_ROD_DMUHeadformatErr[rosL][drawer][dmu] = ((*itColl)->getFragHeaderBit() >> dmu) & 1;
        m_arrays->m_ROD_DMUHeadparityErr[rosL][drawer][dmu] = ((*itColl)->getFragHeaderPar() >> dmu) & 1;
        m_arrays->m_ROD_DMUDataformatErr[rosL][drawer][dmu] = ((*itColl)->getFragSampleBit() >> dmu) & 1;
        m_arrays->m_ROD_DMUDataparityErr[rosL][drawer][dmu] = ((*itColl)->getFragSamplePar() >> dmu) & 1;
        m_arrays->m_ROD_DMUfeCRC[rosL][drawer][dmu] = ((*itColl)->getFragFEChipMask() >> dmu) & 1;
        m_arrays->m_ROD_DMUrodCRC[rosL][drawer][dmu] = ((*itColl)->getFragRODChipMask() >> dmu) & 1;
      }
    }
    
    if (m_compareMode) {
      m_l2Builder->SumE(ROS, drawer, m_rchUnit, E, gain, m_bad[rosL][drawer], sumE);
      if (dspCont) {
        m_sumEt_xx[m_l2Builder->idToIndex(fragId)] = sumE[0];
        m_sumEz_xx[m_l2Builder->idToIndex(fragId)] = sumE[1];
        m_sumE_xx[m_l2Builder->idToIndex(fragId)] = sumE[2];
      }
      else {
        m_sumEt_zz[m_l2Builder->idToIndex(fragId)] = sumE[0];
        m_sumEz_zz[m_l2Builder->idToIndex(fragId)] = sumE[1];
        m_sumE_zz[m_l2Builder->idToIndex(fragId)] = sumE[2];
      }
    }
  }
  
  if (m_compareMode && dspCont) {
    
    const TileL2Container* l2Cnt = SG::makeHandle(m_l2CntKey, ctx).get();
    
    TileL2Container::const_iterator it = l2Cnt->begin();
    TileL2Container::const_iterator end= l2Cnt->end();
    int i=0;
    for(; it != end; ++it) {
      m_sumEt_yy[i++] = (*it)->sumEt();
      m_sumEz_yy[i++] = (*it)->sumEz();
      m_sumE_yy[i++]  = (*it)->sumE();
    }
  }
  
  return StatusCode::SUCCESS;
}

/**
 /// Fill Ntuple with info from TileDigits
 /// Return true if the collection is empty
 */
StatusCode
TileAANtuple::storeDigits(const EventContext& ctx
                          , const SG::ReadHandleKey<TileDigitsContainer>& containerKey
                          , short *a_sample
                          , short a_gain[N_ROS2][N_MODULES][N_CHANS]
                          , bool fillAll)
{
  if (containerKey.empty()) // empty name, nothing to do
    return StatusCode::FAILURE;
  
  // Read Digits from TES
  const TileDigitsContainer* digitsCnt =
    SG::makeHandle (containerKey, ctx).get();
  
  bool emptyColl = true;
  
  // Get iterator for all TDColl in TDCont
  TileDigitsContainer::const_iterator itColl = (*digitsCnt).begin();
  TileDigitsContainer::const_iterator itCollEnd = (*digitsCnt).end();
  
  // Go through all TileDigitsCollections
  for(; itColl != itCollEnd; itColl++) {
    int fragId = (*itColl)->identify();
    int drawer = fragId & 0x3F;
    int ROS = (fragId>>8);
    int rosI = ROS-1;
    int rosL = rosI;
    int rosH = rosI + N_ROS;
    
    if (msgLvl(MSG::VERBOSE)) {
      ATH_MSG_VERBOSE( "Event# " << m_evtNr
                       << " Frag id 0x" << MSG::hex << fragId << MSG::dec
                       << " ROS " << ROS
                       << " drawer " << drawer  );
      
      if (fillAll) {
        ATH_MSG_VERBOSE( "       Size=" << (*itColl)->getFragSize()
                         << " BCID=" << (*itColl)->getFragBCID() << MSG::hex
                         << " CRC=0x" << ((*itColl)->getFragCRC()&0xffff)
                         << " DMUMask=0x" << ((*itColl)->getFragDMUMask()&0xffff) << MSG::dec  );
        
        ATH_MSG_VERBOSE( "       Lvl1ID=" << (*itColl)->getLvl1Id()
                         << " Lvl1Type=" << (*itColl)->getLvl1Type()
                         << " EvBCID=" << (*itColl)->getRODBCID()
                         << " EvType=" << (*itColl)->getDetEvType()  );
        
        ATH_MSG_VERBOSE("       Header=" << (*itColl)->getFragChipHeaderWords()  );
      }
    }
    
    /// Store ROD header info from collection
    /// (should be just one per ROD, i.e. 4 subsequent drawers give the same ROD number)
    int rod = (rosL*N_MODULES + drawer) >> 2;
    
    m_l1ID[rod] = (*itColl)->getLvl1Id();
    m_l1Type[rod] = (*itColl)->getLvl1Type();
    m_evType[rod] = (*itColl)->getDetEvType();
    m_evBCID[rod] = (*itColl)->getRODBCID();
    if (m_trigType == 0) m_trigType = -m_l1Type[rod]; // make negative to distinguish from TileCal internal trig types
    
    TileDigitsCollection::const_iterator it = (*itColl)->begin();
    TileDigitsCollection::const_iterator itEnd = (*itColl)->end();
    
    // non empty collection
    if(it != itEnd) {
      emptyColl = false;
      
      if (fillAll) {
        
        // store evtnr, bcid,crc, size
        // Same for lo and hi, because they come from the same fragment
        
        m_arrays->m_rodBCID[rosL][drawer] = (*itColl)->getRODBCID();
        m_arrays->m_fragSize[rosL][drawer]=(*itColl)->getFragSize();
        
        m_arrays->m_slinkCRC[rosL][drawer][0] = ((*itColl)->getFragCRC()>>16) & 0xffff;
        m_arrays->m_dmuMask[rosL][drawer][0] = ((*itColl)->getFragDMUMask()>>16) & 0xffff;
        m_arrays->m_slinkCRC[rosL][drawer][1] = (*itColl)->getFragCRC() & 0xffff;
        m_arrays->m_dmuMask[rosL][drawer][1] = (*itColl)->getFragDMUMask() & 0xffff;
        
        uint32_t CRCmask = (*itColl)->getFragDMUMask(); //mask of FE+ROD DMU crc check (16bit+16bit) 0xffffffff == All ok
        uint32_t fe_crc = CRCmask & 0xFFFF;
        uint32_t rod_crc = CRCmask >> 16;
        
        const std::vector<uint32_t> & headerVec = (*itColl)->getFragChipHeaderWords();
        unsigned int headsize = std::min(16U, static_cast<unsigned int>(headerVec.size()));
        
        for (unsigned int ih = 0; ih < headsize; ++ih) {
          
          m_arrays->m_DMUheader[rosL][drawer][ih] = headerVec[ih];   /// Full DMU header, stored for debugging
          m_arrays->m_DMUbcid[rosL][drawer][ih] = (headerVec[ih] & 0xFFF);             /// BCID in DMU header
          m_arrays->m_DMUformatErr[rosL][drawer][ih] = CheckDMUFormat(headerVec[ih]); /// bit_31==1 and bit_17==0
          m_arrays->m_DMUparityErr[rosL][drawer][ih] = CheckDMUParity(headerVec[ih]); /// parity must be an odd number
          m_arrays->m_DMUmemoryErr[rosL][drawer][ih] = (headerVec[ih] >> 25 & 0x1); /// memory parity error bit_25
          m_arrays->m_DMUSstrobeErr[rosL][drawer][ih] = (headerVec[ih] >> 24 & 0x1); /// single strobe error bit_24 (it is recovered)
          m_arrays->m_DMUDstrobeErr[rosL][drawer][ih] = (headerVec[ih] >> 23 & 0x1); /// double strobe error bit_23 (cannot be recovered)
          
          m_arrays->m_feCRC[rosL][drawer][ih] = (fe_crc >> ih & 0x1);
          m_arrays->m_rodCRC[rosL][drawer][ih] = (rod_crc >> ih & 0x1);
        }
        
        if (m_calibMode) {
          const std::vector<uint32_t> & headerVecHi = (*itColl)->getFragChipHeaderWordsHigh();
          unsigned int headsizehi = std::min(16U, static_cast<unsigned int>(headerVecHi.size()));
          
          for (unsigned int ih = 0; ih < headsizehi; ++ih) {
            
            m_arrays->m_DMUheader[rosH][drawer][ih] = headerVecHi[ih];                    /// Full DMU header, stored for debugging
            m_arrays->m_DMUbcid[rosH][drawer][ih] = (headerVecHi[ih] & 0xFFF);            /// BCID in DMU header
            m_arrays->m_DMUformatErr[rosH][drawer][ih] = CheckDMUFormat(headerVecHi[ih]); /// bit_31==1 and bit_17==0
            m_arrays->m_DMUparityErr[rosH][drawer][ih] = CheckDMUParity(headerVecHi[ih]); /// parity must be an odd number
            m_arrays->m_DMUmemoryErr[rosH][drawer][ih]  = (headerVecHi[ih] >> 25 & 0x1);  /// memory parity error bit_25
            m_arrays->m_DMUSstrobeErr[rosH][drawer][ih] = (headerVecHi[ih] >> 24 & 0x1);  /// single strobe error bit_24 (it is recovered)
            m_arrays->m_DMUDstrobeErr[rosH][drawer][ih] = (headerVecHi[ih] >> 23 & 0x1);  /// double strobe error bit_23 (cannot be recovered)
            
            m_arrays->m_feCRC[rosH][drawer][ih]  = -1 ; //Variables must be filled anyway, empty variables are not allowed
            m_arrays->m_rodCRC[rosH][drawer][ih] = -1;  //Variables must be filled anyway, empty variables are not allowed
          }
        }
      }
      
      
      int cmpCounter = 0;
      // go through all TileDigits in collection
      for (; it != itEnd; it++) {
        const TileDigits* digit = (*it);
        
        HWIdentifier hwid = digit->adc_HWID();
        
        // determine channel
        int channel = m_tileHWID->channel(hwid);
        // convert channel number to PMT number if needed
        if (m_pmtOrder) channel = digiChannel2PMT(ROS,channel);
        
        // determine gain and set ros index accordingly
        int gain = m_tileHWID->adc(hwid);
        if (m_calibMode) {
          if (m_compareMode) {
            ++cmpCounter;
            if (cmpCounter > 48) rosI = rosH;
            else rosI = rosL;
          } else {
            if (gain == 1) rosI = rosH;
            else rosI = rosL;
          }
        }
        
        a_gain[rosI][drawer][channel] = gain;
        
        // get digits
        const std::vector<float> & sampleVec = digit->samples();
        int siz = sampleVec.size();
        if (msgLvl(MSG::VERBOSE)) {
          int index,pmt;
          digit->cell_ID_index(index,pmt);
          msg(MSG::VERBOSE) << "TD ch " << channel
          << " gain "<< gain
          << " type " << std::min(index,0) << " {";
          for(int i=0;i<siz;i++) {
            msg(MSG::VERBOSE) <<(int)sampleVec[i] << " ";
          }
        }
	// changed N_SAMPLES to number of samples from tile configurator
        if (siz > m_nSamples) {
          siz = m_nSamples;
          if (msgLvl(MSG::VERBOSE))
            ATH_MSG_VERBOSE( "} ONLY " << siz << " digits saved to ntuple" );
        } else {
          if (msgLvl(MSG::VERBOSE))
            ATH_MSG_VERBOSE( "}"  );
        }
        for (int n = 0; n < siz; ++n) {
          a_sample[sample_ind(rosI,drawer,channel,n)] = (short) sampleVec[n];
        }
      }
    }
  }
  
  if (emptyColl) return StatusCode::FAILURE;
  else return StatusCode::SUCCESS;
}

StatusCode TileAANtuple::storeTMDBDecision(const EventContext& ctx) {

  const char * part[4] = {"LBA","LBC","EBA","EBC"};

  // Read Decision from TES
  //
  if (!m_tileMuRcvContainerKey.empty()){

    ATH_MSG_VERBOSE( "reading TMDB decision from " << m_tileMuRcvContainerKey.key() ); 

    const TileMuonReceiverContainer *decisionCnt =
      SG::makeHandle (m_tileMuRcvContainerKey, ctx).get();
  
    TileMuonReceiverContainer::const_iterator it = decisionCnt->begin();
    TileMuonReceiverContainer::const_iterator itLast = decisionCnt->end();
  
    // Go through all decisions
    for(; it != itLast; ++it) {

      const TileMuonReceiverObj * obj = (*it);

      const std::vector<bool> & decision = obj->GetDecision(); 
      int siz = decision.size();

      if (siz>0) {

        int fragId = (*it)->identify();
        int drawer = fragId & 0x3F;
        int ros    = (fragId>>8) - 1;
 
        if (siz > N_TMDBDECISIONS) {
          ATH_MSG_VERBOSE( "ONLY " << N_TMDBDECISIONS << " decisions saved to ntuple instead of " << siz);
          siz = N_TMDBDECISIONS;
        }

        for (int n = 0; n < siz; ++n) {
          m_arrays->m_decisionTMDB[ros][drawer][n] = (unsigned char) decision[n];
        }

        if (msgLvl(MSG::VERBOSE)) {
          std::stringstream ss;
          for (int n = 0; n < siz; ++n) {
            ss<<std::setw(5)<<(int)m_arrays->m_decisionTMDB[ros][drawer][n];
          }
          ATH_MSG_VERBOSE( "   0x" <<MSG::hex<< fragId <<MSG::dec<<" "<< part[ros] 
                           << std::setfill('0') << std::setw(2)
                           << drawer+1 << std::setfill(' ') 
                           << "      decision: " <<ss.str()  );
        }
      }
    }
  }

  return StatusCode::SUCCESS;
}
 
StatusCode TileAANtuple::storeTMDBDigits(const EventContext& ctx) {

  const char * part[4] = {"LBA","LBC","EBA","EBC"};

  // Read Digits from TES
  //
  if (!m_tileMuRcvDigitsContainerKey.empty()){

    ATH_MSG_VERBOSE( "reading TMDB digits from " << m_tileMuRcvDigitsContainerKey.key() ); 

    const TileDigitsContainer* digitsCnt =
      SG::makeHandle (m_tileMuRcvDigitsContainerKey, ctx).get();
  
    TileDigitsContainer::const_iterator itColl1 = (*digitsCnt).begin();
    TileDigitsContainer::const_iterator itCollEnd1 = (*digitsCnt).end();
  
    // Go through all TileDigitsCollections
    for(; itColl1 != itCollEnd1; ++itColl1) {

      TileDigitsCollection::const_iterator it1 = (*itColl1)->begin();
      TileDigitsCollection::const_iterator itEnd1 = (*itColl1)->end();

      if (it1!=itEnd1) {

        int fragId = (*itColl1)->identify();
        int drawer = fragId & 0x3F;
        int ros    = (fragId>>8) - 1;
        int ichannel = 0;
 
        ATH_MSG_VERBOSE( "   0x" <<MSG::hex<< fragId <<MSG::dec<<" "<< part[ros]
                         << std::setfill('0') << std::setw(2)
                          << drawer+1 << std::setfill(' ') ); 

        for (; it1 != itEnd1; ++it1) {

          if (ichannel>=N_TMDBCHANS) {
            ATH_MSG_WARNING("Too many channels in TMDB Digi container for frag 0x" <<MSG::hex<< fragId <<MSG::dec <<" keeping only first " << N_TMDBCHANS << " channels in ntuple ");
            break;
          }

          const TileDigits* digit = (*it1);
	  
          // get digits
          const std::vector<float> & sampleVec = digit->samples();
          int siz = sampleVec.size();

          if (siz > m_nSamples) {
            ATH_MSG_VERBOSE( "ONLY " << m_nSamples << " digits saved to ntuple instead of " << siz);
            siz = m_nSamples;
          }

          for (int n = 0; n < siz; ++n) {
            m_arrays->m_sampleTMDB[sample_ind_TMDB(ros,drawer,ichannel,n)] = (unsigned char) sampleVec[n];
          }

          if (msgLvl(MSG::VERBOSE)) {
            std::stringstream ss;
            for (int n = 0; n < siz; ++n) {
              ss<<std::setw(5)<<(int)m_arrays->m_sampleTMDB[sample_ind_TMDB(ros,drawer,ichannel,n)];
            }
            ATH_MSG_VERBOSE( "      dig: " <<ros+1<<"/"<<drawer<<"/"<<m_tileHWID->channel(digit->adc_HWID())<<": "<<ss.str()  );
          }
      
          ++ichannel;
        }
      }
    }
  }

  return StatusCode::SUCCESS;
}

StatusCode TileAANtuple::storeTMDBRawChannel(const EventContext& ctx) {

  const char * part[4] = {"LBA","LBC","EBA","EBC"};

  // Read Raw Channels from TDS
  //
  if (!m_tileMuRcvRawChannelContainerKey.empty()){

    ATH_MSG_VERBOSE( "reading TMDB energies from " << m_tileMuRcvRawChannelContainerKey.key() ); 

    const TileRawChannelContainer* rcCnt =
      SG::makeHandle (m_tileMuRcvRawChannelContainerKey, ctx).get();

    TileRawChannelContainer::const_iterator itColl2 = (*rcCnt).begin();
    TileRawChannelContainer::const_iterator itCollEnd2 = (*rcCnt).end();
  
    // Go through all TileDigitsCollections
    for(; itColl2 != itCollEnd2; ++itColl2) {

      TileRawChannelCollection::const_iterator it2 = (*itColl2)->begin();
      TileRawChannelCollection::const_iterator itEnd2 = (*itColl2)->end();

      if (it2!=itEnd2) {

        int fragId = (*itColl2)->identify();
        int drawer = fragId & 0x3F;
        int ros    = (fragId>>8) - 1;
        int ichannel = 0;
  
        ATH_MSG_VERBOSE( "   0x" <<MSG::hex<< fragId <<MSG::dec<<" "<< part[ros]
                         << std::setfill('0') << std::setw(2)
                         << drawer+1 << std::setfill(' ') ); 

        for (; it2 != itEnd2; ++it2) {

          if (ichannel>=N_TMDBCHANS) {
            ATH_MSG_WARNING("Too many channels in TMDB RCh container for frag 0x" <<MSG::hex<< fragId <<MSG::dec <<" keeping only first " << N_TMDBCHANS << " channels in ntuple ");
            break;
          }

          const TileRawChannel* rc = (*it2);
        
          m_arrays->m_eTMDB[ros][drawer][ichannel] =  rc -> amplitude();

          ATH_MSG_VERBOSE( "      rc: " <<ros+1<<"/"<<drawer<<"/"<<m_tileHWID->channel(rc->adc_HWID())<< ": " << m_arrays->m_eTMDB[ros][drawer][ichannel] );

          ++ichannel;
        }
      }
    }
  }

  return StatusCode::SUCCESS;
}

StatusCode
TileAANtuple::finalize() {
  
  if (m_arrays->m_sample) free(m_arrays->m_sample);
  if (m_arrays->m_sampleFlt) free(m_arrays->m_sampleFlt);
  if (m_arrays->m_sampleTMDB) free(m_arrays->m_sampleTMDB);
  if (m_arrays->m_eMF) free(m_arrays->m_eMF);
  if (m_arrays->m_tMF) free(m_arrays->m_tMF);

  ATH_MSG_INFO( "finalize() successfully" );
  return StatusCode::SUCCESS;
}

StatusCode
TileAANtuple::ntuple_clear() {
  
  TRIGGER_clearBranch();
  CISPAR_clearBranch();
  LASER_clearBranch();
  DIGI_clearBranch();
  TMDB_clearBranch();

  return StatusCode::SUCCESS;
}

StatusCode
TileAANtuple::initNTuple(const EventContext& ctx) {
  //Aux Ntuple creation
  
  if (m_compressSettings >= 0) {
    std::vector<std::string> files;
    m_fileMgr->getFiles(Io::ROOT, ( Io::WRITE | Io::CREATE ), files);
    for (std::string& file : files) {
      TFile* outFile = (TFile*) m_fileMgr->fptr(file);
      if (outFile) {
	ATH_MSG_INFO("Changing compressing settings to " << m_compressSettings << " for file: " << file);
	outFile->SetCompressionSettings(m_compressSettings);
      }
    }
  }

  if (m_ntupleID.size() > 0) {
    
    std::string ntupleID = m_ntupleID + "_map";
    
    TTree* ntuplePtr = new TTree(ntupleID.c_str(), "TileCal-CellMap");
    if(m_thistSvc->regTree("/" + m_streamName + "/" + ntupleID, ntuplePtr).isFailure()) {
      ATH_MSG_ERROR( "Problem registering TileRec CellMap Tree" );
    }
    
    fillCellMap(ntuplePtr);
    
    //Ntuple creation
    m_ntuplePtr = new TTree(m_ntupleID.c_str(), "TileCal-Ntuple");
    m_ntuplePtr->SetMaxTreeSize(m_treeSize);
    if (m_thistSvc->regTree("/" + m_streamName + "/" + m_ntupleID, m_ntuplePtr).isFailure()) {
      ATH_MSG_ERROR( "Problem registering TileRec Tree");
    }
    
    TRIGGER_addBranch();
    CISPAR_addBranch();
    if (!m_laserObjectKey.empty()) {
      const TileLaserObject* laserObj =
        SG::makeHandle(m_laserObjectKey, ctx).get();
      m_las_version = laserObj->getVersion();
      LASER_addBranch();
    }
    DIGI_addBranch();
    TMDB_addBranch();
  }
  
  //DCS Ntuple creation
  if (m_checkDCS) {
    std::string ntupleDCS = "Tile_DCS";
    m_DCSntuplePtr = new TTree(ntupleDCS.c_str(), "TileCal-DCS data");
    if (m_thistSvc->regTree("/" + m_streamName + "/" + ntupleDCS, m_DCSntuplePtr).isFailure()) {
      ATH_MSG_ERROR( "Problem registering TileRec DCS Tree");
    }
    DCS_addBranch();
  }
  
  return StatusCode::SUCCESS;
}


/**
 /////////////////////////////////////////////////////////////////////////////
 //
 ///    Variables Legenda
 ///
 ///  - C : a character string terminated by the 0 character
 ///  - B : an 8 bit signed integer
 ///  - b : an 8 bit unsigned integer                    2^8=256
 ///  - S : a 16 bit signed integer (i.e. a "short")
 ///  - s : a 16 bit unsigned integer                    2^16=65536
 ///  - I : a 32 bit signed integer (i.e an "int")
 ///  - i : a 32 bit unsigned integer                    2^32=4294967296
 ///  - F : a 32 bit floating point (i.e. a "float")
 ///  - D : a 64 bit floating point (i.e. a "double")
 ///  - L : a 64 bit signed integer
 ///  - l : a 64 bit unsigned integer
 ///  - O : a boolean
 //
 */

/**
 //////////////////////////////////////////////////////////////////////////////
 ///Create eta-phi map for all channels
 //
 //////////////////////////////////////////////////////////////////////////////
 */
void TileAANtuple::fillCellMap(TTree* ntuplePtr) {
  
  float eta[4][64][48];
  float phi[4][64][48];
  short tower[4][64][48];
  short sample[4][64][48];
  short ind[4][64][48];
  short pmt[4][64][48];
  short bad[4][64][48][2];
  
  CLEAR(eta);
  CLEAR(phi);
  CLEAR1(tower);
  CLEAR1(sample);
  CLEAR1(ind);
  CLEAR(pmt);
  CLEAR1(bad);
  
  ntuplePtr->Branch("eta", eta, "eta[4][64][48]/F");
  ntuplePtr->Branch("phi", phi, "phi[4][64][48]/F");
  ntuplePtr->Branch("tower", tower, "tower[4][64][48]/S");
  ntuplePtr->Branch("sample", sample, "sample[4][64][48]/S");
  ntuplePtr->Branch("ind", ind, "ind[4][64][48]/S");
  ntuplePtr->Branch("pmt", pmt, "pmt[4][64][48]/S");
  ntuplePtr->Branch("bad", bad, "bad[4][64][48][2]/S");
  
  TileDetDescrManager::calo_element_const_iterator itr = m_tileMgr->tile_cell_begin();
  TileDetDescrManager::calo_element_const_iterator end = m_tileMgr->tile_cell_end();
  
  for (; itr != end; ++itr) {
    const CaloDetDescrElement * caloDDE = (*itr);
    Identifier cell_id = caloDDE->identify();
    IdentifierHash hash[2] = { caloDDE->onl1(), caloDDE->onl2() };
    for (int i = 0; i < 2; ++i) {
      if (hash[i] != TileHWID::NOT_VALID_HASH) {
        HWIdentifier adc_id = m_tileHWID->adc_id(hash[i], 0);
        int ROS = m_tileHWID->ros(adc_id);
        int drawer = m_tileHWID->drawer(adc_id);
        int chan = m_tileHWID->channel(adc_id);
        int index, pm;
        m_cabling->h2s_cell_id_index(adc_id, index, pm); // index=-2 for MBTS or -1 for all disconnected channels
        if (index > -1) index = i; // just put back 0 or 1 for all connected channels
        pm = m_cabling->channel2hole(ROS, chan); // convert channel to PMT number, ignoring difference between drawers
        if ((ROS == 3 && drawer == 14) || (ROS == 4 && drawer == 17)) {
          if (pm < 0) pm = -pm; // crack scin in EBA15 EBC18 in non-standard place - recover positive pmt number
          if (chan == 2 || chan == 3) pm = -pm; // no D4 in special EBA15 EBC18 - put negative sign
        }
        if (m_pmtOrder) chan = digiChannel2PMT(ROS,chan); // convert channel to PMT-1
        int rosI = ROS-1; // make ros index from 0 to 3
        eta[rosI][drawer][chan] = caloDDE->eta();
        phi[rosI][drawer][chan] = caloDDE->phi();
        tower[rosI][drawer][chan] = m_tileID->tower(cell_id);
        sample[rosI][drawer][chan] = m_tileID->sample(cell_id);
        ind[rosI][drawer][chan] = index;
        pmt[rosI][drawer][chan] = pm;
      }
    }
  }
  
  const MbtsDetDescrManager* mbtsMgr; //!< Pointer to MbtsDetDescrManager
  if ( detStore()->retrieve(mbtsMgr).isFailure() ) {
    ATH_MSG_WARNING( "Unable to retrieve MbtsDetDescrManager from DetectorStore" );
    mbtsMgr = 0;
  }
  for (int ROS = 1; ROS < 5; ++ROS) {
    int rosI = ROS - 1;
    for (int drawer = 0; drawer < 64; ++drawer) {
      for (int chan = 0; chan < 48; ++chan) {
        for (int adc = 0; adc < 2; ++adc) {
          HWIdentifier adc_id = m_tileHWID->adc_id(ROS, drawer, chan, adc);
          bad[rosI][drawer][chan][adc] = (short) m_tileBadChanTool->encodeStatus(m_tileBadChanTool->getAdcStatus(adc_id));
          int index, pm;
          Identifier cell_id = m_cabling->h2s_cell_id_index(adc_id, index, pm);
          if (index == -2) { // MBTS
            ind[rosI][drawer][chan] = index;
            pmt[rosI][drawer][chan] = 1; // we know that all MBTS are in PMT 1
            if (mbtsMgr) {
              const CaloDetDescrElement * caloDDE = mbtsMgr->get_element(cell_id);
              if (caloDDE) {
                if (caloDDE->z() > 0.0)
                  eta[rosI][drawer][chan] = fabs(caloDDE->eta());
                else
                  eta[rosI][drawer][chan] = -fabs(caloDDE->eta());
                phi[rosI][drawer][chan] = caloDDE->phi();
              }
            }
          }
        }
      }
    }
  }
  
  ntuplePtr->Fill();
}


/**
 //////////////////////////////////////////////////////////////////////////////
 ///Add TRIGGER variables to the Tree
 //
 //////////////////////////////////////////////////////////////////////////////
 */
void TileAANtuple::TRIGGER_addBranch(void) {
  m_ntuplePtr->Branch("EvTime",&m_evTime,"EvTime/I");
  m_ntuplePtr->Branch("Run",&m_run,"Run/I");
  m_ntuplePtr->Branch("LumiBlock",&m_lumiBlock,"LumiBlock/I");
  m_ntuplePtr->Branch("Evt",&m_evt,"Evt/I");
  m_ntuplePtr->Branch("EvtNr",&m_evtNr,"EvtNr/I");
  m_ntuplePtr->Branch("Trig",&m_trigType,"Trig/I");
  m_ntuplePtr->Branch("DSPflags",&m_dspFlags,"DSPflags/i");
  m_ntuplePtr->Branch("DSPunits",&m_dspUnit,"DSPunits/S");
  m_ntuplePtr->Branch("OFLunits",&m_rchUnit,"OFLunits/S");
  
  if (m_bsInput) {
    m_ntuplePtr->Branch("L1ID",   m_l1ID,   "L1ID[65]/I");
    m_ntuplePtr->Branch("L1Type", m_l1Type, "L1Type[65]/I");
    m_ntuplePtr->Branch("EvType", m_evType, "EvType[65]/I");
    m_ntuplePtr->Branch("EvBCID", m_evBCID, "EvBCID[65]/I");
  }
}


/**
 //////////////////////////////////////////////////////////////////////////////
 //Clear Tree TRIGGER variables
 //
 //////////////////////////////////////////////////////////////////////////////
 */
void TileAANtuple::TRIGGER_clearBranch(void) {
  m_evTime=0;
  m_run=0;
  m_evt=0;
  m_trigType=0;
  m_dspFlags=0;
  
  CLEAR1(m_l1ID);
  CLEAR1(m_l1Type);
  CLEAR1(m_evType);
  CLEAR1(m_evBCID);
}


/**
 //////////////////////////////////////////////////////////////////////////////
 ///Add Tree CISPAR variables Tree
 //
 //////////////////////////////////////////////////////////////////////////////
 */
void TileAANtuple::CISPAR_addBranch(void) {
  if (!m_beamElemContainerKey.key().empty()) {
    m_ntuplePtr->Branch("cispar",m_cispar,"cispar[110]/i");
  }
}


/**
 //////////////////////////////////////////////////////////////////////////////
 //Clear Tree CISPAR variables
 //
 //////////////////////////////////////////////////////////////////////////////
 */
void TileAANtuple::CISPAR_clearBranch(void) {
  if (!m_beamElemContainerKey.key().empty()) {
    CLEAR(m_cispar);
  }
}


/**
 //////////////////////////////////////////////////////////////////////////////
 ///Add Tree LASER variables Tree
 //
 //////////////////////////////////////////////////////////////////////////////
 */
void TileAANtuple::LASER_addBranch(void) {
  
  if (!m_laserObjectKey.empty()) {
    
    const char* gainnames[2]  = {"LG","HG"};
    const char* channames[16] = {"Diode0","Diode1","Diode2","Diode3","Diode4","Diode5","Diode6","Diode7",
      "Diode8","Diode9","PMT1","ExtCIS0","IntCIS","DiodePhocal","PMT2","ExtCIS1"};
    
    m_ntuplePtr->Branch("LASER_BCID", &m_las_BCID, "LASER_BCID/I");
    
    m_ntuplePtr->Branch("LASER_FILTER", &m_las_Filt, "LASER_FILTER/I");
    m_ntuplePtr->Branch("LASER_REQAMP", &m_las_ReqAmp, "LASER_REQAMP/F");
    m_ntuplePtr->Branch("LASER_MEASAMP", &m_las_MeasAmp, "LASER_MEASAMP/F");
    
    if(m_las_version==2) {
      
      ATH_MSG_DEBUG("LASERII BRANCHING..");
      
      m_ntuplePtr->Branch(Form("LASER_QDCTIMEOUT"),&(m_qdctimeout),Form("LASER_QDCTIMEOUT/O"));
      m_ntuplePtr->Branch(Form("LASER_TDCTIMEOUT"),&(m_tdctimeout),Form("LASER_TDCTIMEOUT/O"));
      m_ntuplePtr->Branch(Form("LASER_DAQTYPE"),&(m_daqtype),Form("LASER_DAQTYPE/I"));
      for(int chan=0;chan<32;++chan){
        int ch=chan>>1;
        int gn=chan&1;
        m_ntuplePtr->Branch(Form("LASER_%s_%s_ADC",gainnames[gn],channames[ch]),&(m_arrays->m_chan[chan]),Form("LASER_%s_%s_ADC/I",gainnames[gn],channames[ch]));
        m_ntuplePtr->Branch(Form("LASER_%s_%s_Ped",gainnames[gn],channames[ch]),&(m_arrays->m_chan_Ped[chan]),Form("LASER_%s_%s_Ped/F",gainnames[gn],channames[ch]));
        m_ntuplePtr->Branch(Form("LASER_%s_%s_Led",gainnames[gn],channames[ch]),&(m_arrays->m_chan_Led[chan]),Form("LASER_%s_%s_Led/F",gainnames[gn],channames[ch]));
        m_ntuplePtr->Branch(Form("LASER_%s_%s_Ped1",gainnames[gn],channames[ch]),&(m_arrays->m_chan_Lin[chan]),Form("LASER_%s_%s_Ped1/F",gainnames[gn],channames[ch]));
        m_ntuplePtr->Branch(Form("LASER_%s_%s_Alpha",gainnames[gn],channames[ch]),&(m_arrays->m_chan_Alpha[chan]),Form("LASER_%s_%s_Alpha/F",gainnames[gn],channames[ch]));
        m_ntuplePtr->Branch(Form("LASER_%s_%s_Sigma_Ped",gainnames[gn],channames[ch]),&(m_arrays->m_chan_SPed[chan]),Form("LASER_%s_%s_Sigma_Ped/F",gainnames[gn],channames[ch]));
        m_ntuplePtr->Branch(Form("LASER_%s_%s_Sigma_Led",gainnames[gn],channames[ch]),&(m_arrays->m_chan_SLed[chan]),Form("LASER_%s_%s_Sigma_Led/F",gainnames[gn],channames[ch]));
        m_ntuplePtr->Branch(Form("LASER_%s_%s_Sigma_Ped1",gainnames[gn],channames[ch]),&(m_arrays->m_chan_SLin[chan]),Form("LASER_%s_%s_Sigma_Ped1/F",gainnames[gn],channames[ch]));
        m_ntuplePtr->Branch(Form("LASER_%s_%s_Sigma_Alpha",gainnames[gn],channames[ch]),&(m_arrays->m_chan_SAlpha[chan]),Form("LASER_%s_%s_Sigma_Alpha/F",gainnames[gn],channames[ch]));
      } // FOR
      
    } else {
      
      ATH_MSG_DEBUG("LASERI BRANCHING..");
      
      m_ntuplePtr->Branch("LASER_Diode_1_ADC", &(m_arrays->m_las_D_ADC[0][0]), "LASER_Diode_1_ADC/I");
      m_ntuplePtr->Branch("LASER_Diode_2_ADC", &(m_arrays->m_las_D_ADC[0][1]), "LASER_Diode_2_ADC/I");
      m_ntuplePtr->Branch("LASER_Diode_3_ADC", &(m_arrays->m_las_D_ADC[0][2]), "LASER_Diode_3_ADC/I");
      m_ntuplePtr->Branch("LASER_Diode_4_ADC", &(m_arrays->m_las_D_ADC[0][3]), "LASER_Diode_4_ADC/I");
      
      m_ntuplePtr->Branch("LASER_Diode_1_Ped", &(m_arrays->m_las_D_Ped[0][0]), "LASER_Diode_1_Ped/F");
      m_ntuplePtr->Branch("LASER_Diode_2_Ped", &(m_arrays->m_las_D_Ped[0][1]), "LASER_Diode_2_Ped/F");
      m_ntuplePtr->Branch("LASER_Diode_3_Ped", &(m_arrays->m_las_D_Ped[0][2]), "LASER_Diode_3_Ped/F");
      m_ntuplePtr->Branch("LASER_Diode_4_Ped", &(m_arrays->m_las_D_Ped[0][3]), "LASER_Diode_4_Ped/F");
      
      m_ntuplePtr->Branch("LASER_Diode_1_Ped_RMS", &(m_arrays->m_las_D_Ped_RMS[0][0]), "LASER_Diode_1_Ped_RMS/F");
      m_ntuplePtr->Branch("LASER_Diode_2_Ped_RMS", &(m_arrays->m_las_D_Ped_RMS[0][1]), "LASER_Diode_2_Ped_RMS/F");
      m_ntuplePtr->Branch("LASER_Diode_3_Ped_RMS", &(m_arrays->m_las_D_Ped_RMS[0][2]), "LASER_Diode_3_Ped_RMS/F");
      m_ntuplePtr->Branch("LASER_Diode_4_Ped_RMS", &(m_arrays->m_las_D_Ped_RMS[0][3]), "LASER_Diode_4_Ped_RMS/F");
      
      m_ntuplePtr->Branch("LASER_Diode_1_Alpha", &(m_arrays->m_las_D_Alpha[0][0]), "LASER_Diode_1_Alpha/F");
      m_ntuplePtr->Branch("LASER_Diode_2_Alpha", &(m_arrays->m_las_D_Alpha[0][1]), "LASER_Diode_2_Alpha/F");
      m_ntuplePtr->Branch("LASER_Diode_3_Alpha", &(m_arrays->m_las_D_Alpha[0][2]), "LASER_Diode_3_Alpha/F");
      m_ntuplePtr->Branch("LASER_Diode_4_Alpha", &(m_arrays->m_las_D_Alpha[0][3]), "LASER_Diode_4_Alpha/F");
      
      m_ntuplePtr->Branch("LASER_Diode_1_Alpha_RMS", &(m_arrays->m_las_D_Alpha_RMS[0][0]), "LASER_Diode_1_Alpha_RMS/F");
      m_ntuplePtr->Branch("LASER_Diode_2_Alpha_RMS", &(m_arrays->m_las_D_Alpha_RMS[0][1]), "LASER_Diode_2_Alpha_RMS/F");
      m_ntuplePtr->Branch("LASER_Diode_3_Alpha_RMS", &(m_arrays->m_las_D_Alpha_RMS[0][2]), "LASER_Diode_3_Alpha_RMS/F");
      m_ntuplePtr->Branch("LASER_Diode_4_Alpha_RMS", &(m_arrays->m_las_D_Alpha_RMS[0][3]), "LASER_Diode_4_Alpha_RMS/F");
      
      m_ntuplePtr->Branch("LASER_Diode_1_AlphaPed", &(m_arrays->m_las_D_AlphaPed[0][0]), "LASER_Diode_1_AlphaPed/F");
      m_ntuplePtr->Branch("LASER_Diode_2_AlphaPed", &(m_arrays->m_las_D_AlphaPed[0][1]), "LASER_Diode_2_AlphaPed/F");
      m_ntuplePtr->Branch("LASER_Diode_3_AlphaPed", &(m_arrays->m_las_D_AlphaPed[0][2]), "LASER_Diode_3_AlphaPed/F");
      m_ntuplePtr->Branch("LASER_Diode_4_AlphaPed", &(m_arrays->m_las_D_AlphaPed[0][3]), "LASER_Diode_4_AlphaPed/F");
      
      m_ntuplePtr->Branch("LASER_Diode_1_AlphaPed_RMS", &(m_arrays->m_las_D_AlphaPed_RMS[0][0]), "LASER_Diode_1_AlphaPed_RMS/F");
      m_ntuplePtr->Branch("LASER_Diode_2_AlphaPed_RMS", &(m_arrays->m_las_D_AlphaPed_RMS[0][1]), "LASER_Diode_2_AlphaPed_RMS/F");
      m_ntuplePtr->Branch("LASER_Diode_3_AlphaPed_RMS", &(m_arrays->m_las_D_AlphaPed_RMS[0][2]), "LASER_Diode_3_AlphaPed_RMS/F");
      m_ntuplePtr->Branch("LASER_Diode_4_AlphaPed_RMS", &(m_arrays->m_las_D_AlphaPed_RMS[0][3]), "LASER_Diode_4_AlphaPed_RMS/F");
      
      m_ntuplePtr->Branch("LASER_PMT_1_ADC", &(m_arrays->m_las_PMT_ADC[0][0]), "LASER_PMT_1_ADC/I");
      m_ntuplePtr->Branch("LASER_PMT_2_ADC", &(m_arrays->m_las_PMT_ADC[0][1]), "LASER_PMT_2_ADC/I");
      
      m_ntuplePtr->Branch("LASER_PMT_1_TDC", &(m_arrays->m_las_PMT_TDC[0][0]), "LASER_PMT_1_TDC/I");
      m_ntuplePtr->Branch("LASER_PMT_2_TDC", &(m_arrays->m_las_PMT_TDC[0][1]), "LASER_PMT_2_TDC/I");
      
      m_ntuplePtr->Branch("LASER_PMT_1_Ped", &(m_arrays->m_las_PMT_Ped[0][0]), "LASER_PMT_1_Ped/F");
      m_ntuplePtr->Branch("LASER_PMT_2_Ped", &(m_arrays->m_las_PMT_Ped[0][1]), "LASER_PMT_2_Ped/F");
      
      m_ntuplePtr->Branch("LASER_PMT_1_Ped_RMS", &(m_arrays->m_las_PMT_Ped_RMS[0][0]), "LASER_PMT_1_Ped_RMS/F");
      m_ntuplePtr->Branch("LASER_PMT_2_Ped_RMS", &(m_arrays->m_las_PMT_Ped_RMS[0][1]), "LASER_PMT_2_Ped_RMS/F");
      
    }
    
    m_ntuplePtr->Branch("LASER_HEAD_Temp", &m_las_Temperature, "LASER_HEAD_Temp/F");
  }
}


/**
 //////////////////////////////////////////////////////////////////////////////
 //Clear Tree LASER variables
 //
 //////////////////////////////////////////////////////////////////////////////
 */
void TileAANtuple::LASER_clearBranch(void) {
  
  if (!m_laserObjectKey.empty()) {
    
    m_las_BCID = 0;
    
    m_las_Filt = 0;
    m_las_ReqAmp = 0.0;
    m_las_MeasAmp = 0.0;
    m_las_Temperature = 0.0;
    
    // LASERII
    //    memset(m_arrays->m_chan, 0, sizeof(m_arrays->m_chan));
    //    memset(m_arrays->m_chan_Ped, 0, sizeof(m_arrays->m_chan_Ped));
    //    memset(m_arrays->m_chan_Led, 0, sizeof(m_arrays->m_chan_Led));
    //    memset(m_arrays->m_chan_Lin, 0, sizeof(m_arrays->m_chan_Lin));
    //    memset(m_arrays->m_chan_Alpha, 0, sizeof(m_arrays->m_chan_Alpha));
    //    memset(m_arrays->m_chan_SPed, 0, sizeof(m_arrays->m_chan_SPed));
    //    memset(m_arrays->m_chan_SLed, 0, sizeof(m_arrays->m_chan_SLed));
    //    memset(m_arrays->m_chan_SLin, 0, sizeof(m_arrays->m_chan_SLin));
    //    memset(m_arrays->m_chan_SAlpha, 0, sizeof(m_arrays->m_chan_SAlpha));
    
    // LASERI
    memset(m_arrays->m_las_D_ADC,          0, sizeof(m_arrays->m_las_D_ADC));
    memset(m_arrays->m_las_D_Ped,          0, sizeof(m_arrays->m_las_D_Ped));
    memset(m_arrays->m_las_D_Ped_RMS,      0, sizeof(m_arrays->m_las_D_Ped_RMS));
    memset(m_arrays->m_las_D_Alpha,        0, sizeof(m_arrays->m_las_D_Alpha));
    memset(m_arrays->m_las_D_Alpha_RMS,    0, sizeof(m_arrays->m_las_D_Alpha_RMS));
    memset(m_arrays->m_las_D_AlphaPed,     0, sizeof(m_arrays->m_las_D_AlphaPed));
    memset(m_arrays->m_las_D_AlphaPed_RMS, 0, sizeof(m_arrays->m_las_D_AlphaPed_RMS));
    
    memset(m_arrays->m_las_PMT_ADC,        0, sizeof(m_arrays->m_las_PMT_ADC));
    memset(m_arrays->m_las_PMT_TDC,        0, sizeof(m_arrays->m_las_PMT_TDC));
    memset(m_arrays->m_las_PMT_Ped,        0, sizeof(m_arrays->m_las_PMT_Ped));
    memset(m_arrays->m_las_PMT_Ped_RMS,    0, sizeof(m_arrays->m_las_PMT_Ped_RMS));
    
  }
}


/**
 //////////////////////////////////////////////////////////////////////////////
 ///Add Tree DIGI variables Tree
 //
 //////////////////////////////////////////////////////////////////////////////
 */
void TileAANtuple::DIGI_addBranch(void)
{
  std::string suf[3] = {"_lo","_hi",""};
  if (m_compareMode) {
    suf[0] = "_xx";
    suf[1] = "_yy";
    suf[2] = "_zz";
    m_ntuplePtr->Branch(NAME2("sumEt",suf[0]),       m_sumEt_xx,     NAME3("sumEt",    suf[0],"[4][64]/F")); // float
    m_ntuplePtr->Branch(NAME2("sumEz",suf[0]),       m_sumEz_xx,     NAME3("sumEz",    suf[0],"[4][64]/F")); // float
    m_ntuplePtr->Branch(NAME2("sumE", suf[0]),       m_sumE_xx,      NAME3("sumE",     suf[0],"[4][64]/F")); // float
    m_ntuplePtr->Branch(NAME2("sumEt",suf[1]),       m_sumEt_yy,     NAME3("sumEt",    suf[1],"[4][64]/F")); // float
    m_ntuplePtr->Branch(NAME2("sumEz",suf[1]),       m_sumEz_yy,     NAME3("sumEz",    suf[1],"[4][64]/F")); // float
    m_ntuplePtr->Branch(NAME2("sumE", suf[1]),       m_sumE_yy,      NAME3("sumE",     suf[1],"[4][64]/F")); // float
    m_ntuplePtr->Branch(NAME2("sumEt",suf[2]),       m_sumEt_zz,     NAME3("sumEt",    suf[2],"[4][64]/F")); // float
    m_ntuplePtr->Branch(NAME2("sumEz",suf[2]),       m_sumEz_zz,     NAME3("sumEz",    suf[2],"[4][64]/F")); // float
    m_ntuplePtr->Branch(NAME2("sumE", suf[2]),       m_sumE_zz,      NAME3("sumE",     suf[2],"[4][64]/F")); // float
  }

  int sample_size = N_ROS*N_MODULES*N_CHANS*m_nSamples;

  int imin = 2, imax = 3, ir = 0, is = 0;
  
  if (m_calibMode) {
    imin = 0;
    imax = 2;
  }
  
  for (int i = imin; i < imax; ++i) {
    
    std::string f_suf(suf[i]);
   
    if (m_fltDigitsContainerKey.empty() && m_digitsContainerKey.empty()
        && (!m_rawChannelContainerKey.empty()
            || !m_fitRawChannelContainerKey.empty()
            || !m_fitcRawChannelContainerKey.empty()
            || !m_optRawChannelContainerKey.empty()
            || !m_qieRawChannelContainerKey.empty()
            || !m_dspRawChannelContainerKey.empty()
            || !m_mfRawChannelContainerKey.empty()
            || !m_of1RawChannelContainerKey.empty()
            || !m_wienerRawChannelContainerKey.empty()
            || !m_bsInput) ) {
          
          m_ntuplePtr->Branch(NAME2("gain",f_suf),            m_arrays->m_gain[ir],          NAME3("gain",         f_suf,"[4][64][48]/S"));    // short
          
        } else {
          
          if (!m_fltDigitsContainerKey.empty()) {
            if (!m_digitsContainerKey.empty()) { // should use different names for two containers
              
              m_ntuplePtr->Branch(NAME2("sampleFlt",f_suf),   &(m_arrays->m_sample[is]),     NAME5("sampleFlt",    f_suf,"[4][64][48][",std::to_string(m_nSamples),"]/S")); // short 
              m_ntuplePtr->Branch(NAME2("gainFlt",f_suf),     m_arrays->m_gainFlt[ir],       NAME3("gainFlt",      f_suf,"[4][64][48]/S"));    // short
            } else {
              m_ntuplePtr->Branch(NAME2("sample",f_suf),      &(m_arrays->m_sampleFlt[is]),     NAME5("sampleFlt",    f_suf,"[4][64][48][",std::to_string(m_nSamples),"]/S")); // short 
              if (!m_rawChannelContainerKey.empty()
                  || !m_fitRawChannelContainerKey.empty()
                  || !m_fitcRawChannelContainerKey.empty()
                  || !m_optRawChannelContainerKey.empty()
                  || !m_qieRawChannelContainerKey.empty()
                  || !m_of1RawChannelContainerKey.empty()
                  || !m_dspRawChannelContainerKey.empty()
                  || !m_wienerRawChannelContainerKey.empty()
                  || m_bsInput) {
                
                m_ntuplePtr->Branch(NAME2("gain",f_suf),      m_arrays->m_gain[ir],          NAME3("gain",         f_suf,"[4][64][48]/S"));    // short
              } else {
                m_ntuplePtr->Branch(NAME2("gain",f_suf),      m_arrays->m_gainFlt[ir],       NAME3("gainFlt",      f_suf,"[4][64][48]/S"));    // short
              }
            }
          }
          
          if (!m_digitsContainerKey.empty()) {
            m_ntuplePtr->Branch(NAME2("sample",f_suf),          &(m_arrays->m_sample[is]),        NAME5("sample",       f_suf,"[4][64][48][",std::to_string(m_nSamples),"]/S")); // short 
            m_ntuplePtr->Branch(NAME2("gain",f_suf),            m_arrays->m_gain[ir],          NAME3("gain",         f_suf,"[4][64][48]/S"));    // short
            
            if (m_bsInput) {
              m_ntuplePtr->Branch(NAME2("DMUheader",f_suf),       m_arrays->m_DMUheader[ir],     NAME3("DMUheader",        f_suf,"[4][64][16]/i")); // uint32
              m_ntuplePtr->Branch(NAME2("DMUBCID",f_suf),         m_arrays->m_DMUbcid[ir],       NAME3("DMUBCID",          f_suf,"[4][64][16]/S")); // short
              m_ntuplePtr->Branch(NAME2("DMUmemoryErr",f_suf),    m_arrays->m_DMUmemoryErr[ir],  NAME3("DMUmemoryErr",     f_suf,"[4][64][16]/S")); // short
              m_ntuplePtr->Branch(NAME2("DMUSstrobeErr",f_suf),   m_arrays->m_DMUSstrobeErr[ir], NAME3("DMUSstrobeErr",    f_suf,"[4][64][16]/S")); // short
              m_ntuplePtr->Branch(NAME2("DMUDstrobeErr",f_suf),   m_arrays->m_DMUDstrobeErr[ir], NAME3("DMUDstrobeErr",    f_suf,"[4][64][16]/S")); // short
              m_ntuplePtr->Branch(NAME2("DMUheadformatErr",f_suf),m_arrays->m_DMUformatErr[ir],  NAME3("DMUheadformatErr", f_suf,"[4][64][16]/S")); // short
              m_ntuplePtr->Branch(NAME2("DMUheadparityErr",f_suf),m_arrays->m_DMUparityErr[ir],  NAME3("DMUheadparityErr", f_suf,"[4][64][16]/S")); // short
              
              m_ntuplePtr->Branch(NAME2("feCRC",f_suf),         m_arrays->m_feCRC[ir],         NAME3("feCRC",        f_suf,"[4][64][16]/S")); // short
              m_ntuplePtr->Branch(NAME2("rodCRC",f_suf),        m_arrays->m_rodCRC[ir],        NAME3("rodCRC",       f_suf,"[4][64][16]/S")); // short
              
              if (i == imin) { // common for low and high gain
                m_ntuplePtr->Branch("rodBCID",  m_arrays->m_rodBCID,  "rodBCID[4][64]/S");    // short
                m_ntuplePtr->Branch("fragSize", m_arrays->m_fragSize,"fragSize[4][64]/S");    // short
                m_ntuplePtr->Branch("DMUmask",  m_arrays->m_dmuMask,  "DMUmask[4][64][2]/s"); // unsigned short
                m_ntuplePtr->Branch("slinkCRC", m_arrays->m_slinkCRC,"slinkCRC[4][64][2]/s"); // unsigned short
              }
            }
          }
        }
    
    if (!m_rawChannelContainerKey.empty()) {
      m_ntuplePtr->Branch(NAME2("ene",f_suf),     m_arrays->m_ene[ir],          NAME3("ene",f_suf,"[4][64][48]/F")); // float
      m_ntuplePtr->Branch(NAME2("time",f_suf),    m_arrays->m_time[ir],        NAME3("time",f_suf,"[4][64][48]/F")); // float
      m_ntuplePtr->Branch(NAME2("ped",f_suf),     m_arrays->m_ped[ir],          NAME3("ped",f_suf,"[4][64][48]/F")); // float
      m_ntuplePtr->Branch(NAME2("chi2",f_suf),    m_arrays->m_chi2[ir],        NAME3("chi2",f_suf,"[4][64][48]/F")); // float
    }
    
    if (!m_fitRawChannelContainerKey.empty()) {
      m_ntuplePtr->Branch(NAME2("eFit",f_suf),    m_arrays->m_eFit[ir],        NAME3("eFit",f_suf,"[4][64][48]/F")); // float
      m_ntuplePtr->Branch(NAME2("tFit",f_suf),    m_arrays->m_tFit[ir],        NAME3("tFit",f_suf,"[4][64][48]/F")); // float
      m_ntuplePtr->Branch(NAME2("pedFit",f_suf),  m_arrays->m_pedFit[ir],    NAME3("pedFit",f_suf,"[4][64][48]/F")); // float
      m_ntuplePtr->Branch(NAME2("chi2Fit",f_suf), m_arrays->m_chi2Fit[ir],  NAME3("chi2Fit",f_suf,"[4][64][48]/F")); // float
    }
    
    if (!m_fitcRawChannelContainerKey.empty()) {
      m_ntuplePtr->Branch(NAME2("eFitc",f_suf),   m_arrays->m_eFitc[ir],      NAME3("eFitc",f_suf,"[4][64][48]/F")); // float
      m_ntuplePtr->Branch(NAME2("tFitc",f_suf),   m_arrays->m_tFitc[ir],      NAME3("tFitc",f_suf,"[4][64][48]/F")); // float
      m_ntuplePtr->Branch(NAME2("pedFitc",f_suf), m_arrays->m_pedFitc[ir],  NAME3("pedFitc",f_suf,"[4][64][48]/F")); // float
      m_ntuplePtr->Branch(NAME2("chi2Fitc",f_suf),m_arrays->m_chi2Fitc[ir],NAME3("chi2Fitc",f_suf,"[4][64][48]/F")); // float
    }
    
    if (!m_optRawChannelContainerKey.empty()) {
      m_ntuplePtr->Branch(NAME2("eOpt",f_suf),    m_arrays->m_eOpt[ir],        NAME3("eOpt",f_suf,"[4][64][48]/F")); // float
      m_ntuplePtr->Branch(NAME2("tOpt",f_suf),    m_arrays->m_tOpt[ir],        NAME3("tOpt",f_suf,"[4][64][48]/F")); // float
      m_ntuplePtr->Branch(NAME2("pedOpt",f_suf),  m_arrays->m_pedOpt[ir],    NAME3("pedOpt",f_suf,"[4][64][48]/F")); // float
      m_ntuplePtr->Branch(NAME2("chi2Opt",f_suf), m_arrays->m_chi2Opt[ir],  NAME3("chi2Opt",f_suf,"[4][64][48]/F")); // float
    }
    
    if (!m_qieRawChannelContainerKey.empty()) {
      m_ntuplePtr->Branch(NAME2("eQIE",f_suf),       m_arrays->m_eQIE[ir],             NAME3("eQIE",f_suf,"[4][64][48]/F")); // float
      m_ntuplePtr->Branch(NAME2("tQIE",f_suf),       m_arrays->m_tQIE[ir],             NAME3("tQIE",f_suf,"[4][64][48]/F")); // float
      m_ntuplePtr->Branch(NAME2("pedQIE",f_suf),     m_arrays->m_pedQIE[ir],         NAME3("pedQIE",f_suf,"[4][64][48]/F")); // float
      m_ntuplePtr->Branch(NAME2("chi2QIE",f_suf),    m_arrays->m_chi2QIE[ir],       NAME3("chi2QIE",f_suf,"[4][64][48]/F")); // float
    }

    if (!m_of1RawChannelContainerKey.empty()) {
      m_ntuplePtr->Branch(NAME2("eOF1",f_suf),       m_arrays->m_eOF1[ir],             NAME3("eOF1",f_suf,"[4][64][48]/F")); // float
      m_ntuplePtr->Branch(NAME2("tOF1",f_suf),       m_arrays->m_tOF1[ir],             NAME3("tOF1",f_suf,"[4][64][48]/F")); // float
      m_ntuplePtr->Branch(NAME2("pedOF1",f_suf),     m_arrays->m_pedOF1[ir],         NAME3("pedOF1",f_suf,"[4][64][48]/F")); // float
      m_ntuplePtr->Branch(NAME2("chi2OF1",f_suf),    m_arrays->m_chi2OF1[ir],       NAME3("chi2OF1",f_suf,"[4][64][48]/F")); // float
    }
    
    if (!m_dspRawChannelContainerKey.empty() && !m_reduced) {
      m_ntuplePtr->Branch(NAME2("eDsp",f_suf),       m_arrays->m_eDsp[ir],             NAME3("eDsp",f_suf,"[4][64][48]/F")); // float
      m_ntuplePtr->Branch(NAME2("tDsp",f_suf),       m_arrays->m_tDsp[ir],             NAME3("tDsp",f_suf,"[4][64][48]/F")); // float
      m_ntuplePtr->Branch(NAME2("pedDsp",f_suf),     m_arrays->m_pedDsp[ir],         NAME3("pedDsp",f_suf,"[4][64][48]/F")); // float
      m_ntuplePtr->Branch(NAME2("chi2Dsp",f_suf),    m_arrays->m_chi2Dsp[ir],       NAME3("chi2Dsp",f_suf,"[4][64][48]/F")); // float
    }

    if (!m_wienerRawChannelContainerKey.empty()) {
      m_ntuplePtr->Branch(NAME2("eWiener",f_suf),    m_arrays->m_eWiener[ir],       NAME3("eWiener",f_suf,"[4][64][48]/F")); // float
      m_ntuplePtr->Branch(NAME2("tWiener",f_suf),    m_arrays->m_tWiener[ir],       NAME3("tWiener",f_suf,"[4][64][48]/F")); // float
      m_ntuplePtr->Branch(NAME2("pedWiener",f_suf),  m_arrays->m_pedWiener[ir],   NAME3("pedWiener",f_suf,"[4][64][48]/F")); // float
      m_ntuplePtr->Branch(NAME2("chi2Wiener",f_suf), m_arrays->m_chi2Wiener[ir], NAME3("chi2Wiener",f_suf,"[4][64][48]/F")); // float
    }
    
    if (!m_mfRawChannelContainerKey.empty()) {
      m_ntuplePtr->Branch(NAME2("eMF",f_suf),        &(m_arrays->m_eMF[is]),               NAME3("eMF",f_suf,NAME3("[4][64][48][",std::to_string(m_nSamples),"]/F"))); // float 
      m_ntuplePtr->Branch(NAME2("tMF",f_suf),        &(m_arrays->m_tMF[is]),               NAME3("tMF",f_suf,NAME3("[4][64][48][",std::to_string(m_nSamples),"]/F"))); // float 
      m_ntuplePtr->Branch(NAME2("chi2MF",f_suf),     m_arrays->m_chi2MF[ir],         NAME3("chi2MF",f_suf,"[4][64][48]/F")); // float
      m_ntuplePtr->Branch(NAME2("pedMF",f_suf),      m_arrays->m_pedMF[ir],           NAME3("pedMF",f_suf,"[4][64][48]/F")); // float
    }
    
    if (m_bsInput) {
      if (i == imin) { // common for low and high gain
        m_ntuplePtr->Branch("ROD_GlobalCRC",       m_arrays->m_ROD_GlobalCRC,        "ROD_GlobalCRC[4][64]/S");            // short
        m_ntuplePtr->Branch("ROD_BCID",            m_arrays->m_ROD_BCID,             "ROD_BCID[4][64]/S");                 // short
        m_ntuplePtr->Branch("ROD_DMUBCIDErr",      m_arrays->m_ROD_DMUBCIDErr,       "ROD_DMUBCIDErr[4][64][16]/S");       // short
        m_ntuplePtr->Branch("ROD_DMUmemoryErr",    m_arrays->m_ROD_DMUmemoryErr,     "ROD_DMUmemoryErr[4][64][16]/S");     // short
        m_ntuplePtr->Branch("ROD_DMUSstrobeErr",   m_arrays->m_ROD_DMUSstrobeErr,    "ROD_DMUSstrobeErr[4][64][16]/S");    // short
        m_ntuplePtr->Branch("ROD_DMUDstrobeErr",   m_arrays->m_ROD_DMUDstrobeErr,    "ROD_DMUDstrobeErr[4][64][16]/S");    // short
        m_ntuplePtr->Branch("ROD_DMUheadformatErr",m_arrays->m_ROD_DMUHeadformatErr, "ROD_DMUheadformatErr[4][64][16]/S"); // short
        m_ntuplePtr->Branch("ROD_DMUheadparityErr",m_arrays->m_ROD_DMUHeadparityErr, "ROD_DMUheadparityErr[4][64][16]/S"); // short
        m_ntuplePtr->Branch("ROD_DMUdataformatErr",m_arrays->m_ROD_DMUDataformatErr, "ROD_DMUdataformatErr[4][64][16]/S"); // short
        m_ntuplePtr->Branch("ROD_DMUdataparityErr",m_arrays->m_ROD_DMUDataparityErr, "ROD_DMUdataparityErr[4][64][16]/S"); // short
        
        m_ntuplePtr->Branch("ROD_feCRC",           m_arrays->m_ROD_DMUfeCRC,         "ROD_feCRC[4][64][16]/S");            // short
        m_ntuplePtr->Branch("ROD_rodCRC",          m_arrays->m_ROD_DMUrodCRC,        "ROD_rodCRC[4][64][16]/S");           // short
        m_ntuplePtr->Branch("ROD_DMUmask",         m_arrays->m_ROD_DMUMask,          "ROD_DMUmask[4][64][2]/s");           // unsigned short
      }
    }
    ir += N_ROS;
    is += sample_size;
  }
}

/**
 //////////////////////////////////////////////////////////////////////////////
 ///Clear Tree DIGI variables
 //////////////////////////////////////////////////////////////////////////////
 */
void TileAANtuple::DIGI_clearBranch(void) {
  unsigned int size = (m_calibMode) ? 1 : 2;
  
  if (m_compareMode) {
    CLEAR(m_sumEt_xx);
    CLEAR(m_sumEt_yy);
    CLEAR(m_sumEt_zz);
    CLEAR(m_sumEz_xx);
    CLEAR(m_sumEz_yy);
    CLEAR(m_sumEz_zz);
    CLEAR(m_sumE_xx);
    CLEAR(m_sumE_yy);
    CLEAR(m_sumE_zz);
  }
  
  CLEAR3(m_arrays->m_gain, size);
  
  if (!m_fltDigitsContainerKey.empty()) {
    CLEAR5(m_arrays->m_sampleFlt, size); 
    CLEAR3(m_arrays->m_gainFlt, size);
  }
  
  if (!m_digitsContainerKey.empty()) {

    CLEAR5(m_arrays->m_sample,size);
    
    if (m_bsInput) {
      CLEAR2(m_arrays->m_DMUheader, size);
      CLEAR3(m_arrays->m_DMUbcid, size);
      CLEAR3(m_arrays->m_DMUformatErr, size);
      CLEAR3(m_arrays->m_DMUparityErr, size);
      CLEAR3(m_arrays->m_DMUmemoryErr, size);
      CLEAR3(m_arrays->m_DMUSstrobeErr, size);
      CLEAR3(m_arrays->m_DMUDstrobeErr, size);
      
      CLEAR3(m_arrays->m_feCRC, size);
      CLEAR3(m_arrays->m_rodCRC, size);
      
      CLEAR1(m_arrays->m_rodBCID);
      CLEAR1(m_arrays->m_fragSize);
      CLEAR(m_arrays->m_dmuMask);
      CLEAR(m_arrays->m_slinkCRC);
    }
    
  }
  
  if (!m_rawChannelContainerKey.empty()) {
    CLEAR2(m_arrays->m_ene, size);
    CLEAR2(m_arrays->m_time, size);
    CLEAR2(m_arrays->m_ped, size);
    CLEAR2(m_arrays->m_chi2, size);
  }
  
  if (!m_fitRawChannelContainerKey.empty()) {
    CLEAR2(m_arrays->m_eFit, size);
    CLEAR2(m_arrays->m_tFit, size);
    CLEAR2(m_arrays->m_pedFit, size);
    CLEAR2(m_arrays->m_chi2Fit, size);
  }
  
  if (!m_fitcRawChannelContainerKey.empty()) {
    CLEAR2(m_arrays->m_eFitc, size);
    CLEAR2(m_arrays->m_tFitc, size);
    CLEAR2(m_arrays->m_pedFitc, size);
    CLEAR2(m_arrays->m_chi2Fitc, size);
  }
  
  if (!m_optRawChannelContainerKey.empty()) {
    CLEAR2(m_arrays->m_eOpt, size);
    CLEAR2(m_arrays->m_tOpt, size);
    CLEAR2(m_arrays->m_pedOpt, size);
    CLEAR2(m_arrays->m_chi2Opt, size);
  }
  

  if (!m_qieRawChannelContainerKey.empty()) {
    CLEAR2(m_arrays->m_eQIE, size);
    CLEAR2(m_arrays->m_tQIE, size);
    CLEAR2(m_arrays->m_pedQIE, size);
    CLEAR2(m_arrays->m_chi2QIE, size);
  }

  if (!m_of1RawChannelContainerKey.empty()) {
    CLEAR2(m_arrays->m_eOF1, size);
    CLEAR2(m_arrays->m_tOF1, size);
    CLEAR2(m_arrays->m_pedOF1, size);
    CLEAR2(m_arrays->m_chi2OF1, size);
  }
  
  if (!m_dspRawChannelContainerKey.empty()) {
    CLEAR2(m_arrays->m_eDsp, size);
    CLEAR2(m_arrays->m_tDsp, size);
    CLEAR2(m_arrays->m_pedDsp, size);
    CLEAR2(m_arrays->m_chi2Dsp, size);
  }
  
  if (!m_mfRawChannelContainerKey.empty()) {
    CLEAR4(m_arrays->m_eMF, size);
    CLEAR4(m_arrays->m_tMF, size);
    CLEAR2(m_arrays->m_chi2MF, size);
    CLEAR2(m_arrays->m_pedMF, size);
  }

  if (!m_wienerRawChannelContainerKey.empty()) {
    CLEAR2(m_arrays->m_eWiener, size);
    CLEAR2(m_arrays->m_tWiener, size);
    CLEAR2(m_arrays->m_pedWiener, size);
    CLEAR2(m_arrays->m_chi2Wiener, size);
  }
  
  if (m_bsInput) {
    CLEAR1(m_arrays->m_ROD_GlobalCRC);
    CLEAR1(m_arrays->m_ROD_BCID);
    CLEAR1(m_arrays->m_ROD_DMUBCIDErr);
    CLEAR1(m_arrays->m_ROD_DMUmemoryErr);
    CLEAR1(m_arrays->m_ROD_DMUSstrobeErr);
    CLEAR1(m_arrays->m_ROD_DMUDstrobeErr);
    CLEAR1(m_arrays->m_ROD_DMUHeadformatErr);
    CLEAR1(m_arrays->m_ROD_DMUHeadparityErr);
    CLEAR1(m_arrays->m_ROD_DMUDataformatErr);
    CLEAR1(m_arrays->m_ROD_DMUDataparityErr);
    CLEAR1(m_arrays->m_ROD_DMUfeCRC);
    CLEAR1(m_arrays->m_ROD_DMUrodCRC);
    CLEAR(m_arrays->m_ROD_DMUMask);
  }
  
}

/*//////////////////////////////////////////////////////////////////////////////
 // TMDB variables
 //////////////////////////////////////////////////////////////////////////////
 */

void TileAANtuple::TMDB_addBranch(void)
{

  if (!m_tileMuRcvRawChannelContainerKey.empty()) {
    m_ntuplePtr->Branch("eTMDB", m_arrays->m_eTMDB, "eTMDB[4][64][8]/F");  // float m_arrays->m_eTMDB[N_ROS][N_MODULES][N_TMDBCHANS]
  }

  if (!m_tileMuRcvDigitsContainerKey.empty()) {
    m_ntuplePtr->Branch("sampleTMDB", &(m_arrays->m_sampleTMDB[0]), NAME3("sampleTMDB[4][64][8][",std::to_string(m_nSamples),"]/b")); // unsigned char m_arrays->m_sampleTMDB[N_ROS][N_MODULES][N_TMDBCHANS][N_SAMPLES] 
  }

  if (!m_tileMuRcvContainerKey.empty()) {
    m_ntuplePtr->Branch("decisionTMDB", m_arrays->m_decisionTMDB, "decisionTMDB[4][64][4]/b"); // unsigned char m_arrays->m_decisionTMDB[N_ROS][N_MODULES][N_TMDBDECISIONS] 
  }

}

void TileAANtuple::TMDB_clearBranch(void)
{
  if (!m_tileMuRcvRawChannelContainerKey.empty()) CLEAR(m_arrays->m_eTMDB);
  if (!m_tileMuRcvDigitsContainerKey.empty()) CLEAR6(m_arrays->m_sampleTMDB);
  if (!m_tileMuRcvContainerKey.empty()) CLEAR(m_arrays->m_decisionTMDB);
}

/*/////////////////////////////////////////////////////////////////////////////
 // DCS variables
 /////////////////////////////////////////////////////////////////////////////
 */

void TileAANtuple::DCS_addBranch() {
  bool br[9];
  int mask = m_DCSBranches;
  
  for (int i = 0; i < 9; ++i) {
    br[i] = (mask % 10);
    mask /= 10;
  }
  
  if (br[0]) {
    m_DCSntuplePtr->Branch("EvTime", &m_evTime,  "EvTime/I");
    m_DCSntuplePtr->Branch("Run",    &m_run,     "Run/I");
    m_DCSntuplePtr->Branch("LumiBlock",&m_lumiBlock,"LumiBlock/I");
    m_DCSntuplePtr->Branch("HHMMSS", &m_HHMMSS,  "HHMMSS/I");
    m_DCSntuplePtr->Branch("Evt",    &m_evt,     "Evt/I");
    m_DCSntuplePtr->Branch("EvtNr",  &m_evtNr,   "EvtNr/I");
  }
  
  if (br[1]) m_DCSntuplePtr->Branch("TEMP",    m_arrays->m_TEMP,    "TEMP[4][64][7]/F");
  if (br[2]) m_DCSntuplePtr->Branch("HV",      m_arrays->m_HV,      "HV[4][64][48]/F");
  if (br[3]) m_DCSntuplePtr->Branch("HVSET",   m_arrays->m_HVSET,   "HVSET[4][64][48]/F");
  if (br[4]) m_DCSntuplePtr->Branch("DRSTATES",m_arrays->m_DRSTATES,"DRSTATES[4][64]/I");
  if (br[5]) m_DCSntuplePtr->Branch("HVSTATUS",m_arrays->m_HVSTATUS,"HVSTATUS[4][64][48]/S");
  if (br[6]) m_DCSntuplePtr->Branch("DRSTATUS",m_arrays->m_DRSTATUS,"DRSTATUS[4][64]/S");
  if (br[7]) m_DCSntuplePtr->Branch("CHSTATUS",m_arrays->m_CHSTATUS,"CHSTATUS[4][64][48]/S");
  if (br[8]) {
    m_DCSntuplePtr->Branch("nBadDr",    &m_nBadDr,    "nBadDr/I");
    m_DCSntuplePtr->Branch("nBadHV",    &m_nBadHV,    "nBadHV/I");
    m_DCSntuplePtr->Branch("nBadDCS",   &m_nBadDCS,   "nBadDCS/I");
    m_DCSntuplePtr->Branch("nBadDB",    &m_nBadDB,    "nBadDB/I");
    m_DCSntuplePtr->Branch("nBadTotal", &m_nBadTotal, "nBadTotal/I");
  }
}

StatusCode TileAANtuple::storeDCS() {

  ATH_MSG_DEBUG( "Filling DCS ntuple:"
                 <<" evtCnt=" << m_evtNr
                 << " evt=" << m_evt
                 << " lumi=" << m_lumiBlock << "  " << m_dateTime );

  CLEAR(m_arrays->m_TEMP);
  CLEAR(m_arrays->m_HV);
  CLEAR(m_arrays->m_HVSET);
  CLEAR(m_arrays->m_DRSTATES);
  CLEAR(m_arrays->m_HVSTATUS);
  CLEAR(m_arrays->m_DRSTATUS);
  CLEAR(m_arrays->m_CHSTATUS);

  m_nBadDr = 0;
  m_nBadHV = 0;
  m_nBadDCS = 0;
  m_nBadDB  = 0;
  m_nBadTotal = 0;
  for (int ROS = 1; ROS < 5; ++ROS) {
    int rosI = ROS - 1;
      
    for (int drawer = 0; drawer < 64; ++drawer) {
      m_arrays->m_DRSTATES[rosI][drawer] = m_tileDCS->getDrawerStates(ROS, drawer);
      m_arrays->m_DRSTATUS[rosI][drawer] = m_tileDCS->getDCSStatus(ROS, drawer);
      bool drbad = m_tileDCS->isStatusBad(ROS, drawer);
        
      if (drbad) {
        ++m_nBadDr;
      }
        
      if (msgLvl(MSG::VERBOSE) || m_arrays->m_DRSTATUS[rosI][drawer] != TileDCSState::OK_DRAWER) {
        ATH_MSG_VERBOSE( "Module=" << TileCalibUtils::getDrawerString(ROS, drawer)
                         << " DRSTATES=" << m_arrays->m_DRSTATES[rosI][drawer]
                         << " DRSTATUS=" << m_arrays->m_DRSTATUS[rosI][drawer]
                         << " => " << ((drbad) ? "bad" : "good")  );
      }
        
      unsigned int drawerIdx = TileCalibUtils::getDrawerIdx(ROS,drawer);
      for (int channel=0; channel<48; ++channel){
        TileBchStatus chStat = m_tileBadChanTool->getChannelStatus(drawerIdx,channel);
        m_arrays->m_HV[rosI][drawer][channel]       = m_tileDCS->getChannelHV(ROS, drawer, channel);
        m_arrays->m_HVSET[rosI][drawer][channel]    = m_tileDCS->getChannelHVSet(ROS, drawer, channel);
        m_arrays->m_HVSTATUS[rosI][drawer][channel] = m_tileDCS->getDCSHVStatus(ROS, drawer, channel);
        m_arrays->m_CHSTATUS[rosI][drawer][channel] = m_tileDCS->getDCSStatus(ROS, drawer, channel)
          + 100 * m_tileBadChanTool->encodeStatus(chStat);
        bool chbad = m_tileDCS->isStatusBad(ROS, drawer, channel);
          
        if (chbad || chStat.isBad()) {
          ++m_nBadTotal;
          if (chbad) ++m_nBadDCS;
          if (chStat.isBad()) ++m_nBadDB;
        }
          
        if (m_tileDCS->isStatusHVBad(ROS, drawer, channel)) {
          ++m_nBadHV;
        }
          
        if (msgLvl(MSG::VERBOSE) || (chbad && !drbad)) {
          int pmt=abs(m_cabling->channel2hole(ROS,channel));
          ATH_MSG_VERBOSE( "Module=" << TileCalibUtils::getDrawerString(ROS, drawer)
                           << " channel=" << channel << " pmt=" << pmt
                           << " HV=" << m_arrays->m_HV[rosI][drawer][channel]
                           << " HVSET=" << m_arrays->m_HVSET[rosI][drawer][channel]
                           << " HVSTATUS=" << m_arrays->m_HVSTATUS[rosI][drawer][channel]
                           << " CHSTATUS=" << m_arrays->m_CHSTATUS[rosI][drawer][channel]
                           << " => " << ((chbad) ? "bad" : "good")  );
        }
      }
        
      for (int ind=0; ind<7; ++ind){
        m_arrays->m_TEMP[rosI][drawer][ind] = m_tileDCS->getChannelHV(ROS, drawer, ind+48);
        ATH_MSG_VERBOSE( "Module=" << TileCalibUtils::getDrawerString(ROS, drawer)
                         << " TEMP" << ind+1 << "=" << m_arrays->m_TEMP[rosI][drawer][ind] );
          
      }
    }
  }
    
  ATH_MSG_DEBUG( "BAD status in DCS: nBadDr=" << m_nBadDr
                 << " nBadHV=" << m_nBadHV
                 << " nBadDCS=" << m_nBadDCS
                 << " nBadDB="  << m_nBadDB
                 << " nBadTotal=" << m_nBadTotal );
    
  m_DCSntuplePtr->Fill();

  
  return StatusCode::SUCCESS;
}
