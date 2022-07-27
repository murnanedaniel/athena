/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

//*****************************************************************************
//  Filename : TileHitToTTL1.cxx
//  Author   : F. Merritt
//  Created  : Feb, 2003
//
//  DESCRIPTION:
//      Created to simulate the Tile Level 1 Trigger Towers (TTL1), which are
//      hardware sums of the Tile channels, about 5 channels per tower.
//      The towers are read out in N time slices, with N=9 as a default.
//      Noise can be added to every sample and threshold can be applied
//
//  HISTORY:
//      August 2009: M. Dunford: Adding saturation of the MBTS pulses
//      March 2009: M. Dunford: Adding in bad channels. Now using the standard 
//                  pulse shape routines (instead of ascii files)   
//
//  BUGS:
//
//*****************************************************************************

// Tile includes
#include "TileSimAlgs/TileHitToTTL1.h"
#include "TileEvent/TileLogicalOrdering.h"
#include "TileIdentifier/TileHWID.h"
#include "TileConditions/TileInfo.h"
#include "TileConditions/TileCablingService.h"
#include "TileCalibBlobObjs/TileCalibUtils.h"

// Calo includes
#include "CaloIdentifier/TileID.h"
#include "CaloIdentifier/TileTBID.h"
#include "CaloIdentifier/CaloLVL1_ID.h"

// Atlas includes
#include "StoreGate/ReadHandle.h"
#include "StoreGate/WriteHandle.h"
#include "StoreGate/ReadCondHandle.h"
#include "AthenaKernel/errorcheck.h"
// For the Athena-based random numbers.
#include "AthenaKernel/IAthRNGSvc.h"
#include "AthenaKernel/RNGWrapper.h"

// Gaudi includes
#include "GaudiKernel/ISvcLocator.h"

//CLHEP includes
#include <CLHEP/Random/Randomize.h>

//C++ STL includes
#include <vector>
#include <algorithm>


using CLHEP::RandGaussQ;



//
// Alg standard initialize function
//
StatusCode TileHitToTTL1::initialize() {

  // retrieve CaloLVL1_ID, TileID, TileHWID helpers and TileIfno from det store

  ATH_CHECK( detStore()->retrieve(m_TT_ID) );

  ATH_CHECK( detStore()->retrieve(m_tileID) );

  ATH_CHECK( detStore()->retrieve(m_tileTBID) );

  ATH_CHECK( detStore()->retrieve(m_tileHWID) );

  //=== Get Tile Info
  ATH_CHECK( detStore()->retrieve(m_tileInfo, m_infoName) );

  //=== Get Tile Bad channel tool
  ATH_CHECK( m_tileBadChanTool.retrieve() );

  //=== get TileCondToolEmscale
  ATH_CHECK( m_tileToolEmscale.retrieve() );

  //=== Get rndm number service
  ATH_CHECK( m_rndmSvc.retrieve() );

  ATH_CHECK( m_samplingFractionKey.initialize() );

  ATH_CHECK( m_cablingSvc.retrieve() );
  m_cabling = m_cablingSvc->cablingService();

  // The 'cosmics' setting refers to running with the 
  // Tile/Chicago Cosmic boards instead of the L1Calo
  // This is a pre-2009 configuration
  // The 'standard' setting is the L1Calo configuration
  m_cosmicsType = (m_TileTTL1Type == "Cosmics");

  if (m_TileTTL1Type == "Cosmics") {
    ATH_MSG_INFO( "Cosmics TTL1 type selected" );

  } else if (m_TileTTL1Type == "Standard") {
    ATH_MSG_INFO( "Standard TTL1 type selected" );

  } else {
    ATH_MSG_INFO( "TileHitToTTL1 failed to recognize the TileTTL1Type: " << m_TileTTL1Type );
    return StatusCode::FAILURE;
  }

  if (m_mbtsTTL1ContainerKey.key().empty()) {
    ATH_MSG_INFO( "TileTTL1 from MBTS will not be produced" ); 
  } else {
    ATH_MSG_INFO( "Storing MBTS TileTTL1 in separate container " << m_mbtsTTL1ContainerKey.key() );
  }
  ATH_CHECK( m_mbtsTTL1ContainerKey.initialize(SG::AllowEmpty) );

  if (m_maskBadChannels) {
    ATH_MSG_INFO( "Bad Channel trigger status will be applied" );
  } else {
    ATH_MSG_INFO(  "Bad Channel trigger status will be ignored" );
  }

  /*............................................................................*/
  // Get all global parameters that will be needed for processing
  // For Cosmics, convolute sub-hits with full shape, producing a fine-grained
  // vector for each tower. After noise/pedestal, calculate peak value and 
  // put only that in TileTTL1 object, otherwise TileTTL1 will be too large.
  // For now, use the same shape for L1 and cosmics, but the latter should
  // be wider, will change soon.
  // For Standard L1, also convolute sub-hits with fine-grained shape, but 
  // produce sample-grained vector in TileTTL1.*/
  // Load the pulse shape (so as to not perform this on every execution)
  // The pulse shape for the MBTS, L1Calo (standard) config and the 
  // cosmics config is the same
  // Right now the phase is set to zero (the pulse is centered)
  m_nSamples = m_tileInfo->NdigitSamples(); // number of time slices for each chan
  m_iTrig = m_tileInfo->ItrigSample();   // index of the triggering time slice
  double phase = 0.0;
  std::vector<double> ttl1Shape(m_nSamples, 0.);
  m_tileInfo->ttl1Shape(m_nSamples, m_iTrig, phase, ttl1Shape);
  if (msgLvl(MSG::DEBUG)) {
    for (int jsamp = 0; jsamp < m_nSamples; ++jsamp) {
      msg(MSG::DEBUG) << "jsamp=" << jsamp << " ttl1shape=" << ttl1Shape[jsamp] << endmsg;
    } // end of pulse shape loading
  }

  // Data for MBTS
  // The MBTS have a special hardware configuration. The trigger uses the tower
  // output but the high gain from the 3&1 card
  m_MBTSiTrig = m_tileInfo->ItrigSample();
  m_MBTSnSamples = m_tileInfo->NdigitSamples();

  // Get TileNoise flag from TileInfo (true => generate noise in TileDigits) 
  m_tileNoise = m_tileInfo->TileNoise();
  // Get TileZeroSuppress flag from TileInfo (true => apply threshold to Digits) 
  m_tileThresh = m_tileInfo->TileZeroSuppress();

  if (msgLvl(MSG::DEBUG)) {
    msg(MSG::DEBUG) << " TTL1Shape[" << m_iTrig << "]=" << ttl1Shape[m_iTrig] << endmsg;

    msg(MSG::DEBUG) << " nSamples=" << m_nSamples
                    << ", iTrig=" << m_iTrig
                    << ", tileNoise=" << ((m_tileNoise) ? "true" : "false")
                    << ", tileThresh=" << ((m_tileThresh) ? "true" : "false") << endmsg;
  }

  ATH_CHECK( m_hitContainerKey.initialize() );
  ATH_CHECK( m_ttl1ContainerKey.initialize() );

  ATH_MSG_INFO( "TileHitToTTL1 initialization completed" );

  return StatusCode::SUCCESS;
}
/*==========================================================================*/
//
// Begin Execution Phase.
//
StatusCode TileHitToTTL1::execute(const EventContext &ctx) const {

  ATH_MSG_DEBUG( "Executing TileHitToTTL1" );

  // declare array for random number generation for noise in samples.
  double Rndm[16];      // Can't use variable size array

  // Prepare RNG Service
  ATHRNG::RNGWrapper* rngWrapper = m_rndmSvc->getEngine(this, m_randomStreamName);
  rngWrapper->setSeed( m_randomStreamName, ctx );

  /*........................................................................*/
  // Get hit container from TES and create TTL1 and MBTS container
  // Note that hit container has 256 collections (one for each drawer),
  // but TTL1 container has no collections and no structure
  /*........................................................................*/

  SG::ReadHandle<TileHitContainer> hitContainer(m_hitContainerKey, ctx);
  ATH_CHECK( hitContainer.isValid() );

  SG::WriteHandle<TileTTL1Container> ttl1Container(m_ttl1ContainerKey, ctx);
  // Register the TTL1 container in the TES
  ATH_CHECK( ttl1Container.record(std::make_unique<TileTTL1Container>()) );
  ATH_MSG_DEBUG( "TileTTL1Container registered successfully (" << m_ttl1ContainerKey.key() << ")" );

  std::unique_ptr<TileTTL1Container>  mbtsTTL1Container;
  if (!m_mbtsTTL1ContainerKey.key().empty()) {
    mbtsTTL1Container = std::make_unique<TileTTL1Container>();
  }

  SG::ReadCondHandle<TileSamplingFraction> samplingFraction(m_samplingFractionKey, ctx);
  ATH_CHECK( samplingFraction.isValid() );

  /*........................................................................*/
  // Create temporary arrays for processing signals.
  // Create array for all TT amplitudes in a single drawer
  /*........................................................................*/
  Identifier ttId[16]; // array of TT identifiers in a single drawer
  std::vector<double> ttAmp[16]; // array of all TT amplitudes in a single drawer
  std::vector<double> MBTSAmp; // MBTS amplitudes in a single drawer
  bool ttHit[16];      // array of TT occupancy in a single drawer
  bool MBTSHit;	       // MBTS occupancy in a single drawer
  int nTT;             // number of hit towers in this drawer.
  int nHit;            // number of hits in this drawer.
  int nIgnore;         // number of ignored hits in this drawer.
  int nTTTot = 0;        // total number of hit towers.
  int nHitTot = 0;       // total number of hits.
  int nIgnoreTot = 0;    // total number of ignored hits.
  double ttAmpTot = 0;   // total energy in good level-1 towers.
  double ttAmpTotIg = 0.;   // total energy in "ignored" level-1 towers.
  int minieta, maxieta, posneg;

  for (int ieta = 0; ieta < 16; ++ieta)
    ttAmp[ieta].resize(m_nSamples);
  MBTSAmp.resize(m_MBTSnSamples);

  // temporary array to make sub of all subhits in one channel 
  std::vector<double> hitSamples(m_nSamples);

  // Create array for the nSamples time-samples of a single tower.
  // If running the Cosmics configuration, put out only the peak value
  std::vector<float> ttL1samples;
  std::vector<float> MBTSsamples;
  if (m_cosmicsType) ttL1samples.resize(1);
  else ttL1samples.resize(m_nSamples);
  MBTSsamples.resize(m_MBTSnSamples);

  std::vector<double> ttl1Shape(m_nSamples, 0.);
  std::vector<double> ttl1MBTSShape(m_MBTSnSamples, 0.);
  /*........................................................................*/
  // Begin loop over all collections (collection = electronic drawer). 
  /*........................................................................*/

  for (const TileHitCollection* hitCollection : *hitContainer) {

    // get drawer and ros number
    HWIdentifier drawer_id = m_tileHWID->drawer_id(hitCollection->identify());
    int ros = m_tileHWID->ros(drawer_id);
    int drawer = m_tileHWID->drawer(drawer_id);
    int drawerIdx = TileCalibUtils::getDrawerIdx(ros, drawer);

    // check that this drawer is connected, if not skip it
    if (m_cosmicsType || m_cabling->connected(ros, drawer)) {
      ATH_MSG_VERBOSE( "ROS " << ros
                      << " drawer " << drawer
                      << " is connected");
    } else {
      continue;
    }

    // Find the partition for this drawer
    // Also set the min and max eta range for the towers in this partition
    switch (ros) {
      case TileHWID::BARREL_POS:
        posneg = +1;
        minieta = 0;
        maxieta = 8;
        break;
      case TileHWID::BARREL_NEG:
        posneg = -1;
        minieta = 0;
        maxieta = 8;
        break;
      case TileHWID::EXTBAR_POS:
        posneg = +1;
        minieta = 9;
        maxieta = m_lastTower;
        break;
      case TileHWID::EXTBAR_NEG:
        posneg = -1;
        minieta = 9;
        maxieta = m_lastTower;
        break;
      default:
        posneg = minieta = maxieta = 0;
        break;
    }

    // Zero temporary array of trigger tower amplitudes (TTL1amp) for this collection
    for (int ieta = 0; ieta < 16; ++ieta) {
      for (int js = 0; js < m_nSamples; ++js)
        ttAmp[ieta][js] = 0.0;
    }
    for (int js = 0; js < m_MBTSnSamples; ++js)
      MBTSAmp[js] = 0.0;

    MBTSHit = false;
    memset(ttHit, 0, sizeof(ttHit));
    nTT = nIgnore = nHit = 0;

    /*........................................................................*/
    // Iterate over all hits in this collection, summing amps for each tower.
    /*........................................................................*/
    for (const TileHit* tile_hit : *hitCollection) {

      // Get hit Identifier (= pmt_id)
      Identifier pmt_id = tile_hit->pmt_ID();
      // Get hit HWIdentifier (= channel_id)
      HWIdentifier pmt_HWid = tile_hit->pmt_HWID();
      // Get channel and ADC number
      int channel = m_tileHWID->channel(pmt_HWid);

      // conversion to hit energy after EMscale correction 
      double hit_calib = samplingFraction->getSamplingFraction(drawerIdx, channel);
      hit_calib = std::round(hit_calib * 1000) / 1000;
      // conversion to charge measured by digitizer
      // The trigger always uses the low gain
      double qfactor = hit_calib / m_tileToolEmscale->channelCalib(drawerIdx, channel
                                                                   , TileID::LOWGAIN, 1.
                                                                   , TileRawChannelUnit::PicoCoulombs
                                                                   , TileRawChannelUnit::MegaElectronVolts);

      // determine if the channel is good from the channel status DB
      bool is_good = true;
      if (m_maskBadChannels) {
        TileBchStatus status = m_tileBadChanTool->getAdcStatus(drawerIdx, channel, TileID::LOWGAIN);

        // if channel is bad, set qfactor to zero
        if (status.isNoGainL1()) {
          is_good = false;
          qfactor = 0.0;
        }

        // check if channel is half gain, scale the qfactor by 50%
        if (status.isHalfGainL1()) {
          is_good = false;
          qfactor *= 0.5;
        }

      } // end of bad channel status

      // Treat signal from MBTS differently
      if (m_tileTBID->is_tiletb(pmt_id)) {

        if (mbtsTTL1Container) {

          // Loop over the subhits for this channel.  For each one,
          // convolute with shaping function and add to digitSamples.
          for (int js = 0; js < m_MBTSnSamples; ++js)
            hitSamples[js] = 0.0;

          int n_hits = tile_hit->size();
          for (int ihit = 0; ihit < n_hits; ++ihit) {
            // Need to pass the negative of t_hit, this is because ttl1Shape returns the amplitude at 
            // a given phase, whereas the t_hit from t=0 when the hit took place
            double t_hit = -(tile_hit->time(ihit));
            m_tileInfo->ttl1Shape(m_MBTSnSamples, m_MBTSiTrig, t_hit, ttl1MBTSShape);

            double e_hit = tile_hit->energy(ihit);
            for (int js = 0; js < m_MBTSnSamples; ++js) {
              hitSamples[js] += e_hit * ttl1MBTSShape[js];
            } // end of loop over MBTS samples
          }  // end loop over sub-hits 

          if (MBTSHit) {
            for (int js = 0; js < m_MBTSnSamples; ++js)
              MBTSAmp[js] += qfactor * hitSamples[js];

          } else {
            MBTSHit = true;
            for (int js = 0; js < m_MBTSnSamples; ++js)
              MBTSAmp[js] = qfactor * hitSamples[js];
          }

          if (msgLvl(MSG::VERBOSE)) {

            // Diagnostic checks
            int side = m_tileTBID->type(pmt_id);
            int phi = m_tileTBID->module(pmt_id);
            int eta = m_tileTBID->channel(pmt_id);
            int channel = m_tileHWID->channel(pmt_HWid);

            msg(MSG::VERBOSE) << "New MBTS Hit:"
                              << " ros=" << ros
                              << ", drawer=" << drawer
                              << ", ch=" << channel
                              << ", side=" << side
                              << ", phi=" << phi
                              << ", eta=" << eta
                              << ", e0=" << hitSamples[m_MBTSiTrig] * hit_calib
                              << ", chan is good=" << is_good << endmsg;
          }
        }
        continue;
      } // end of MBTS loop

      // Get TT Identifier for this pmt 
      Identifier tt_id = tile_hit->tt_ID();

      // Get eta-phi indices of TTL1 for this channel.
      int ieta = m_TT_ID->eta(tt_id);
      int iphi = m_TT_ID->phi(tt_id); // (same as module number).
      if (iphi != drawer && m_tileID->sample(pmt_id) != TileID::SAMP_E)
        ATH_MSG_ERROR( "drawer=" << drawer << ", iphi=" << iphi );

      // Loop over the subhits for this channel.  For each one,
      // convolute with shaping function and add to digitSamples.
      for (int js = 0; js < m_nSamples; ++js)
        hitSamples[js] = 0.0;
      int n_hits = tile_hit->size();
      for (int ihit = 0; ihit < n_hits; ++ihit) {
        // Need to pass the negative of t_hit, this is because ttl1Shape returns the amplitude at 
        // a given phase, whereas the t_hit from t=0 when the hit took place
        double t_hit = -(tile_hit->time(ihit));
        m_tileInfo->ttl1Shape(m_nSamples, m_iTrig, t_hit, ttl1Shape);

        double e_hit = tile_hit->energy(ihit);
        for (int js = 0; js < m_nSamples; ++js) {
          hitSamples[js] += e_hit * ttl1Shape[js];
        } // end of loop over samples
      } // end loop over sub-hits

      // check if TT already exists, if so just add energy
      if (ttHit[ieta]) {
        for (int js = 0; js < m_nSamples; ++js)
          ttAmp[ieta][js] += qfactor * hitSamples[js];

        // if not create new TT
      } else {
        ttId[ieta] = tt_id;
        ttHit[ieta] = true;
        for (int js = 0; js < m_nSamples; ++js)
          ttAmp[ieta][js] = qfactor * hitSamples[js];

        if (ieta >= minieta && ieta <= maxieta)
          ++nTT; // count only valid TT
      }
      ++nHit;
      if (ieta < minieta || ieta > maxieta || !is_good)
        ++nIgnore;

      //Sum cell energy for comparison to other algos.
      if (ieta >= minieta && ieta <= maxieta && is_good) {
        ttAmpTot += hitSamples[m_iTrig] * hit_calib;
      } else {
        ttAmpTotIg += hitSamples[m_iTrig] * hit_calib;
      }

      if (msgLvl(MSG::VERBOSE)) {

        // Diagnostic checks:
        int side = m_tileID->side(pmt_id);
        int tower = m_tileID->tower(pmt_id);
        int sample = m_tileID->sample(pmt_id);
        int pmt = m_tileID->pmt(pmt_id);
        int channel = m_tileHWID->channel(pmt_HWid);

        msg(MSG::VERBOSE) << "New Hit:"
                          << " ros=" << ros
                          << ", drawer=" << drawer
                          << ", ch=" << channel
                          << ", side=" << side
                          << ", tower=" << tower
                          << ", sample=" << sample
                          << ", pmt=" << pmt
                          << ", e0=" << hitSamples[m_iTrig] * hit_calib
                          << ", ie=" << ieta
                          << ", ip=" << iphi
                          << ", chan is good= " << is_good;

        if (ieta >= minieta && ieta <= maxieta)
          msg(MSG::VERBOSE) << endmsg;
        else
          msg(MSG::VERBOSE) << " Outside limits" << endmsg;
      } // end of verbose printing 

    } // end loop over hits in this drawer.

    nTTTot += nTT;
    nHitTot += nHit;
    nIgnoreTot += nIgnore;

    ATH_MSG_VERBOSE( "      Statistics for"
                    << " ROS=" << ros
                    << ", drawer=" << drawer
                    << ";  posneg=" << posneg
                    << ", minieta=" << minieta
                    << ", maxieta=" << maxieta
                    << ";  nTT=" << nTT
                    << ", nHit=" << nHit
                    << ", nIgnore=" << nIgnore );


    /*........................................................................*/
    // We now have all the TTL1 amplitudes for this drawer.  
    // Loop over towers to produce the electronics signals (= time samples).
    // If tileNoise is requested, generate random numbers to give noise
    // For Cosmics configuration, only calculate peak value.
    // (no ADC on Chicago custom board, only discriminator)
    /*........................................................................*/

    if (mbtsTTL1Container) {
      Identifier MBTS_id = m_cabling->drawer2MBTS_id(drawer_id);
      if (MBTS_id.is_valid()) {
        bool Good = m_tileNoise || MBTSHit;
        if (Good) {
          double ttL1NoiseSigma = m_tileInfo->MBTSL1NoiseSigma(MBTS_id);
          double ttL1Thresh = m_tileInfo->MBTSL1Thresh(MBTS_id);
          double ttL1Ped = m_tileInfo->MBTSL1Ped(MBTS_id);
          double ttL1Calib = m_tileInfo->MBTSL1Calib(MBTS_id);
          double ttL1Max = m_tileInfo->MBTSL1Max(MBTS_id);

          if (m_tileNoise)
            RandGaussQ::shootArray(rngWrapper->getEngine(ctx), m_MBTSnSamples, Rndm);
          for (int jsamp = 0; jsamp < m_MBTSnSamples; ++jsamp) {
            MBTSAmp[jsamp] *= ttL1Calib; // convert pCb to mV
            MBTSsamples[jsamp] = MBTSAmp[jsamp] + ttL1Ped;
            if (m_tileNoise)
              MBTSsamples[jsamp] += ttL1NoiseSigma * Rndm[jsamp];

            // check if the voltage is above the saturation point,
            // if so, set the pulse value to the saturation point
            if (MBTSsamples[jsamp] > ttL1Max)
              MBTSsamples[jsamp] = ttL1Max;

          }  // end loop over samples

          if (m_tileThresh)
            if (MBTSsamples[m_MBTSiTrig] - ttL1Ped < ttL1Thresh)
              Good = false;

          if (Good) {
            std::unique_ptr<TileTTL1> mbtsTTL1 = std::make_unique<TileTTL1>(MBTS_id, MBTSsamples);
            mbtsTTL1Container->push_back(mbtsTTL1.release());
            ATH_MSG_DEBUG( "mbtsTTL1 saved. Is MBTS hit " << MBTSHit
                          << " Is noise " << m_tileNoise );
          }
        }
      }
    }

    for (int ieta = minieta; ieta <= maxieta; ++ieta) {
      int iphi = drawer;
      bool Good = m_tileNoise || ttHit[ieta];
      if (Good) {
        if (!ttHit[ieta])
          ttId[ieta] = m_TT_ID->tower_id(posneg, 1, 0, ieta, drawer);

        /* Include shaping fuction, pedestal, and noise. */

        // Chicago cosmic Trigger board
        if (m_cosmicsType) {
          double ttL1NoiseSigma = m_tileInfo->TTL1CosmicsNoiseSigma(ttId[ieta]);
          double ttL1Thresh = m_tileInfo->TTL1CosmicsThresh(ttId[ieta]);
          double ttL1Ped = m_tileInfo->TTL1CosmicsPed(ttId[ieta]);
          double ttL1Calib = m_tileInfo->TTL1CosmicsCalib(ttId[ieta]);

          double peakAmp = 0.0;
          int peakSamp = 0;
          for (int jsamp = 0; jsamp < m_nSamples; ++jsamp) {
            ttAmp[ieta][jsamp] *= ttL1Calib; // convert pCb to mV
            ttAmp[ieta][jsamp] += ttL1Ped;
            if (ttAmp[ieta][jsamp] > peakAmp) {
              peakAmp = ttAmp[ieta][jsamp];
              peakSamp = jsamp;
            }
          }  // end loop over samples

          if (m_tileNoise)
            peakAmp += ttL1NoiseSigma * RandGaussQ::shoot(rngWrapper->getEngine(ctx));
          ttL1samples[0] = peakAmp;
          if (m_tileThresh) {
            if (ttL1samples[0] - ttL1Ped < ttL1Thresh)
              Good = false;
          }
          if (msgLvl(MSG::DEBUG) && Good) {
            msg(MSG::DEBUG) << " TTL1:  "
                            << " ros=" << ros
                            << ", ieta=" << ieta
                            << ", iphi=" << iphi
                            << ", hitTrue=" << ttHit[ieta]
                            << ", Good=" << Good
                            << ", peak Amp=" << ttAmp[ieta][peakSamp]
                            << ", with noise=" << ttL1samples[0] << endmsg;
          }

          // ATLAS Level1 Calo Trigger
        } else {
          double ttL1NoiseSigma = m_tileInfo->TTL1NoiseSigma(ttId[ieta]);
          double ttL1Thresh = m_tileInfo->TTL1Thresh(ttId[ieta]);
          double ttL1Ped = m_tileInfo->TTL1Ped(ttId[ieta]);
          double ttL1Calib = m_tileInfo->TTL1Calib(ttId[ieta]);
          double ttL1Max = m_tileInfo->TTL1Max(ttId[ieta]);

          if (m_tileNoise)
            RandGaussQ::shootArray(rngWrapper->getEngine(ctx), m_nSamples, Rndm);
          for (int jsamp = 0; jsamp < m_nSamples; ++jsamp) {
            ttAmp[ieta][jsamp] *= ttL1Calib; // convert pCb to mV
            ttL1samples[jsamp] = ttAmp[ieta][jsamp] + ttL1Ped;
            if (m_tileNoise)
              ttL1samples[jsamp] += ttL1NoiseSigma * Rndm[jsamp];

            // check if the voltage is above the saturation point,
            // if so, set the pulse value to the saturation point
            if (ttL1samples[jsamp] > ttL1Max)
              ttL1samples[jsamp] = ttL1Max;

          }  // end loop over samples

          if (m_tileThresh) {
            if (ttL1samples[m_iTrig] - ttL1Ped < ttL1Thresh)
              Good = false;
          }
          if (msgLvl(MSG::DEBUG) && Good) {
            msg(MSG::DEBUG) << " TTL1:  "
                            << " ros=" << ros
                            << ", ieta=" << ieta
                            << ", iphi=" << iphi
                            << ", hitTrue=" << ttHit[ieta]
                            << ", Good=" << Good
                            << ", amp0=" << ttAmp[ieta][m_iTrig]
                            << ", digitIn=" << ttL1samples[m_iTrig] << endmsg;
          }
        }
      } // end first "Good" section.

      /* Create the new TTL1 object and store in TTL1Container. */
      if (Good) {
        std::unique_ptr<TileTTL1> ttl1 = std::make_unique<TileTTL1>(ttId[ieta], ttL1samples);
        ttl1Container->push_back(ttl1.release());
      }  // end second "Good" section.
    } // end loop over towers
  } // end loop over collections

  // sort all trigger towers according to identifier
  if (m_cosmicsType) {
    ATH_MSG_DEBUG( "Sorting container of size " << ttl1Container->size() );
    TileLogicalOrdering<TileTTL1> order;
    ttl1Container->sort(order);
  }

  if (mbtsTTL1Container) {
    SG::WriteHandle<TileTTL1Container> mbtsContainer(m_mbtsTTL1ContainerKey, ctx);
    ATH_CHECK( mbtsContainer.record(std::move(mbtsTTL1Container)));

    ATH_MSG_DEBUG( "MBTS TileTTL1Container registered successfully (" << m_mbtsTTL1ContainerKey.key() << ")" );
  }

  // Execution completed.
  if (msgLvl(MSG::DEBUG)) {
    msg(MSG::DEBUG) << "TileHitToTTL1 execution completed." << endmsg;
    msg(MSG::DEBUG) << " nTTTot=" << nTTTot
                    << " nHitTot=" << nHitTot
                    << " nIgnoreTot=" << nIgnoreTot
                    << "  ttAmpTot=" << ttAmpTot
                    << " ttAmpTotIg=" << ttAmpTotIg
                    << " =>eneTot=" << ttAmpTot + ttAmpTotIg << endmsg;
  }

  return StatusCode::SUCCESS;
}

StatusCode TileHitToTTL1::finalize() {

  ATH_MSG_INFO( "TileHitToTTL1::finalize() end" );

  return StatusCode::SUCCESS;
}
