/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// ********************************************************************
//
// NAME:     TileRawChannelMonTool.cxx
// PACKAGE:  
//
// AUTHOR:   Lukas Plazak
//
//
// ********************************************************************

#include "TileMonitoring/TileRawChannelTimeMonTool.h"

#include "TileCalibBlobObjs/TileCalibUtils.h"
#include "TileConditions/ITileBadChanTool.h"
#include "TileEvent/TileRawChannelContainer.h"
#include "TileRecUtils/TileBeamInfoProvider.h"

#include "TProfile.h"
#include "TProfile2D.h"


#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>


/*---------------------------------------------------------*/
TileRawChannelTimeMonTool::TileRawChannelTimeMonTool(const std::string & type, const std::string & name, const IInterface* parent)
 : TileFatherMonTool(type, name, parent)
 , m_beamInfo("TileBeamInfoProvider")
 , m_tileBadChanTool("TileBadChanTool")
 , m_DQstatus(0)
 , m_doOnline(false)
 , m_old_lumiblock(-1)
 , m_delta_lumiblock(0)
 , m_nEvents(0)
 , m_bookProfHistOnce(5, false)


	 /*---------------------------------------------------------*/
{
  declareInterface<IMonitorToolBase>(this);

	// run type 1 - phys, 2 - las, 4 - ped, 8 - cis
  declareProperty("runType", m_runType = 0);
  declareProperty("LowGainThreshold", m_lowGainThreshold = 10.0);
  declareProperty("HiGainThreshold", m_hiGainThreshold = 40.0);
  declareProperty("bigain", m_bigain = true);
  declareProperty("doOnline"               , m_doOnline = false); //online mode
  declareProperty("TileRawChannelContainer", m_contName = "TileRawChannelOpt2"); //SG RC Container
  declareProperty("TileBadChanTool"        , m_tileBadChanTool);

  m_path = "/Tile/RawChannelTime"; 
}

/*---------------------------------------------------------*/
TileRawChannelTimeMonTool::~TileRawChannelTimeMonTool()
  /*---------------------------------------------------------*/
{
}

/*---------------------------------------------------------*/
StatusCode TileRawChannelTimeMonTool::initialize()
  /*---------------------------------------------------------*/
{

  ATH_MSG_INFO("in initialize()");

  CHECK(m_beamInfo.retrieve());
  CHECK(m_tileBadChanTool.retrieve());

  m_nEvents = 0;
  m_thresholds[0] = m_lowGainThreshold;
  m_thresholds[1] = m_hiGainThreshold;

  CHECK(TileFatherMonTool::initialize());

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode TileRawChannelTimeMonTool::bookHists()
  /*---------------------------------------------------------*/
{
  if (msgLvl(MSG::DEBUG)) {
    msg(MSG::DEBUG) << "in bookHists()" << endmsg;
    msg(MSG::DEBUG) << "Using base path " << m_path << endmsg;
  }

  for (unsigned int ros = 1; ros < TileCalibUtils::MAX_ROS; ++ros) {
    for (unsigned int drawer = 0; drawer < TileCalibUtils::MAX_DRAWER; ++drawer) {
      bookHists(ros, drawer);
    }
  }

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
void TileRawChannelTimeMonTool::bookHists(int ros, int drawer) {
  /*---------------------------------------------------------*/

  std::string partitionName[5] = { "AUX", "LBA", "LBC", "EBA", "EBC" };

  std::string moduleName = TileCalibUtils::getDrawerString(ros, drawer);
  std::string subDir = moduleName;
  std::string histName, histTitle;

  ATH_MSG_DEBUG("in bookHists() for module " << moduleName);

  for(int dig = 0; dig < 8; ++dig){

    histName = moduleName +  "_digitizer_" + std::to_string(dig + 1);
    histTitle = moduleName + " Time vs LumiBlock, Digi" + std::to_string(dig + 1);

    if(m_doOnline){
      profile_hist[ros][drawer][dig]=bookProfile(subDir, histName, histTitle, 100, -99.5, 0.5);
      profile_hist[ros][drawer][dig]->GetXaxis()->SetTitle("Last LumiBlocks");
    }else{
      profile_hist[ros][drawer][dig]=bookProfile(subDir, histName, histTitle, 1500, -0.5, 1499.5);
      profile_hist[ros][drawer][dig]->GetXaxis()->SetTitle("LumiBlock");
    }
  }

  if (!m_bookProfHistOnce[ros]) {

    histName = partitionName[ros];
    histName += "_avgTime";

    histTitle = partitionName[ros];
    histTitle += " Average Time ";

    subDir = "Summary";

    profile2d_hist[ros] = bookProfile2D(subDir, histName, histTitle, 64, 0.5, 64.5, 48, -0.5, 47.5, -50, 50);

    std::string module_name;
    std::string cell_name;
    std::string channel_name;

    for (unsigned int drw = 0; drw < TileCalibUtils::MAX_DRAWER; drw += 2) {
      module_name = TileCalibUtils::getDrawerString(ros, drw);
      profile2d_hist[ros]->GetXaxis()->SetBinLabel(drw + 1, module_name.c_str());
    }

    for (unsigned int channel = 0; channel < TileCalibUtils::MAX_CHAN; ++channel) {
      cell_name = getCellName(ros, channel);
      channel_name = cell_name + (cell_name.empty() ? "ch" : "_ch") + std::to_string(channel);
      profile2d_hist[ros]->GetYaxis()->SetBinLabel(channel + 1, channel_name.c_str());
    }

    m_bookProfHistOnce[ros] = true;
  }

}

/*---------------------------------------------------------*/
StatusCode TileRawChannelTimeMonTool::fillHists()
  /*---------------------------------------------------------*/
{

  static const int ch2digi[48] = {7, 7, 7,  7, 7, 7,
                                  6, 6, 6,  6, 6, 6,
                                  5, 5, 5,  5, 5, 5,
                                  4, 4, 4,  4, 4, 4,
                                  3, 3, 3,  3, 3, 3,
                                  2, 2, 2,  2, 2, 2,
                                  1, 1, 1,  1, 1, 1,
                                  0, 0, 0,  0, 0, 0};

  ATH_MSG_DEBUG("in fillHists()");
  fillEvtInfo();
  if (m_nEvents % 1000 == 0) ATH_MSG_INFO(m_nEvents<<" events processed so far");
  ++m_nEvents;

  m_DQstatus = m_beamInfo->getDQstatus();


  const TileRawChannelContainer* RawChannelCnt;
  CHECK(evtStore()->retrieve(RawChannelCnt, m_contName));

  std::vector<double> avgTimePerPart(TileCalibUtils::MAX_ROS, 0.0);
  std::vector<double> sumTimeCh(TileCalibUtils::MAX_ROS, 0.0);
  std::vector<int>  nCh(TileCalibUtils::MAX_ROS, 0);

  for (int k = 0; k < 2; k++) {  //k=0 - computing avg time, k=1 applying avg time

    TileRawChannelContainer::const_iterator collItr = RawChannelCnt->begin();
    TileRawChannelContainer::const_iterator lastColl = RawChannelCnt->end();

    for (; collItr != lastColl; ++collItr) {

      TileRawChannelCollection::const_iterator chItr = (*collItr)->begin();
      TileRawChannelCollection::const_iterator lastCh = (*collItr)->end();

      if (chItr != lastCh) {

        HWIdentifier adc_id = (*chItr)->adc_HWID();
        int ros = m_tileHWID->ros(adc_id);
        int drawer = m_tileHWID->drawer(adc_id);

        unsigned int drawerIdx = TileCalibUtils::getDrawerIdx(ros, drawer);

        for (; chItr != lastCh; ++chItr) {


          const TileRawChannel* rch = (*chItr);
          adc_id = rch->adc_HWID();
          unsigned int chan = m_tileHWID->channel(adc_id);
          
          if (isDisconnected(ros, drawer, chan)) continue;

          int adc = m_tileHWID->adc(adc_id);

          bool good  = m_DQstatus->isAdcDQgood(ros, drawer, chan, adc) && m_beamInfo->isChanDCSgood(ros, drawer, chan);

          if (good) {
            TileBchStatus status = m_tileBadChanTool->getAdcStatus(drawerIdx, chan, adc);
            good = ! ( status.isBad() ||  status.isBadTiming() );
          }

          if ( !good || (rch->amplitude()) < m_thresholds[adc] ) continue;

          double time = rch->time();
          double timeCorr = 0;

          if (k == 0) { 

            if ((ros == 3 || ros == 4)
                && (chan == 0 || chan == 1 || chan == 2 || chan == 3 || chan == 4 || chan == 5 || chan == 12 || chan == 13 || chan == 18
                  || chan == 19)) {
            } else {
              sumTimeCh[ros] += time;
              nCh[ros] += 1;
            }


          }  else if (k == 1) { //if k==0 
          
            timeCorr = time - avgTimePerPart[ros];      

            profile2d_hist[ros]->Fill(drawer + 1, chan, timeCorr);

            int32_t current_lumiblock = getLumiBlock();
            if(m_old_lumiblock == -1) {
              m_old_lumiblock = current_lumiblock;
            }

            if(m_doOnline) {
              m_delta_lumiblock = current_lumiblock - m_old_lumiblock;

              if(m_delta_lumiblock > 0) {//move bins

                for (unsigned int ros = 1; ros < TileCalibUtils::MAX_ROS; ++ros) {
                  for (unsigned int drawer = 0; drawer < TileCalibUtils::MAX_DRAWER; ++drawer) {
                    for (unsigned int digi = 0; digi < 8; ++digi) {
                      ShiftTprofile(profile_hist[ros][drawer][digi], m_delta_lumiblock);
                    }
                  }
                }

                m_old_lumiblock = current_lumiblock;
              }

              profile_hist[ros][drawer][ch2digi[chan]]->Fill(0., timeCorr);

            } else {// End of Online
              profile_hist[ros][drawer][ch2digi[chan]]->Fill(current_lumiblock, timeCorr);
            }
          }      //k==1 
        } // loop over channels
      }
    }

    if (k == 0) {
      for (unsigned int ros = 1; ros < TileCalibUtils::MAX_ROS; ++ros) {
        if (nCh[ros] != 0) {
          avgTimePerPart[ros] = sumTimeCh[ros] / nCh[ros];
        } else {
          avgTimePerPart[ros] = 0;
        }
        sumTimeCh[ros] = 0;
        nCh[ros] = 0;
      } //for
    } //if k==0
  } // loop over k 
  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode TileRawChannelTimeMonTool::finalHists()
  /*---------------------------------------------------------*/
{
  ATH_MSG_INFO("in finalHists()");

  return StatusCode::SUCCESS;
}
/*---------------------------------------------------------*/
StatusCode TileRawChannelTimeMonTool::checkHists(bool /* fromFinalize */)
  /*---------------------------------------------------------*/
{
  ATH_MSG_INFO("in checkHists()");

  return StatusCode::SUCCESS;
}




