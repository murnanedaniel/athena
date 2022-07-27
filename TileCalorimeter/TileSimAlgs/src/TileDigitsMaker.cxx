/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

//****************************************************************************
//  Filename : TileDigitsMaker.cxx
//  Author   : Zhifang
//  Created  : Feb 2006 from TileHitToDigits
//
//  DESCRIPTION:
//      Created to simulate the digits information (amplitudes of N time-slices,
//      with N about 9) received by the ROD's.)
//
//  HISTORY:
//
//  BUGS:
//
//*****************************************************************************

// Tile includes
#include "TileSimAlgs/TileDigitsMaker.h"
#include "TileEvent/TileMutableDigitsContainer.h"
#include "TileEvent/TileRawChannelContainer.h"
#include "TileIdentifier/TileHWID.h"
#include "TileConditions/TileInfo.h"
#include "TileCalibBlobObjs/TileCalibUtils.h"
#include "TileConditions/TileCablingService.h"
#include "TileConditions/TilePulseShapes.h"

// Calo includes
#include "CaloIdentifier/TileID.h"
#include "CaloIdentifier/TileTBID.h"

// Atlas include
#include "AthenaKernel/errorcheck.h"
// For the Athena-based random numbers.
#include "AthenaKernel/IAthRNGSvc.h"
#include "AthenaKernel/RNGWrapper.h"

#include "AthenaKernel/Units.h"
#include "StoreGate/ReadHandle.h"
#include "StoreGate/WriteHandle.h"
#include "StoreGate/ReadCondHandle.h"
#include "GaudiKernel/ThreadLocalContext.h"
// Pile up
#include "PileUpTools/PileUpMergeSvc.h"

// Gaudi includes

//CLHEP includes
#include <CLHEP/Random/Randomize.h>


#include "TMatrixF.h"
#include "TDecompChol.h"
#include "cmath"

#include <algorithm>

using CLHEP::RandGaussQ;
using CLHEP::RandFlat;
using Athena::Units::MeV;


//
// Alg standard initialize function
//
StatusCode TileDigitsMaker::initialize() {
  // retrieve TileID helper and TileInfo from det store

  ATH_CHECK( detStore()->retrieve(m_tileID) );

  ATH_CHECK( detStore()->retrieve(m_tileTBID) );

  ATH_CHECK( detStore()->retrieve(m_tileHWID) );

  ATH_CHECK( detStore()->retrieve(m_tileInfo, m_infoName) );

  //=== intialize TileEMScale condtion object
  ATH_CHECK( m_emScaleKey.initialize() );

  //=== initialize TileSampleNoise condition object
  ATH_CHECK( m_sampleNoiseKey.initialize() );

  ATH_CHECK( m_cablingSvc.retrieve() );
  m_cabling = m_cablingSvc->cablingService();

  /* Get needed parameters from tileInfo. */
  m_nSamples = m_tileInfo->NdigitSamples(); // number of time slices for each chan
  m_iTrig = m_tileInfo->ItrigSample();   // index of the triggering time slice
  m_i_ADCmax = m_tileInfo->ADCmax();     // adc saturation value used in assignment
  m_f_ADCmax = m_i_ADCmax;               // adc saturation value used in assignment
  m_f_ADCmaxHG = m_f_ADCmax - 0.5;       // value of switch from high to low gain
  m_ADCmaxMinusEps = m_f_ADCmax - 0.01;
  m_ADCmaxPlusEps = m_f_ADCmax + 0.01;
  m_f_ADCmaskValue = m_tileInfo->ADCmaskValue();               // indicates channels which were masked in background dataset
  m_tileNoise = m_tileInfo->TileNoise(); // (true => generate noise in TileDigits)
  m_tileCoherNoise = m_tileInfo->TileCoherNoise(); // (true => generate coherent noise in TileDigits)
  m_tileThresh = m_tileInfo->TileZeroSuppress(); // (true => apply threshold to Digits)
  m_tileThreshHi = m_tileInfo->ThresholdDigits(TileID::HIGHGAIN);
  m_tileThreshLo = m_tileInfo->ThresholdDigits(TileID::LOWGAIN);

  if (m_tileNoise || m_tileCoherNoise || m_rndmEvtOverlay) {
    ATH_CHECK( m_rndmSvc.retrieve());
  }

  ATH_MSG_DEBUG( "Event Overlay: " << ((m_rndmEvtOverlay)?"true":"false"));
  ATH_MSG_DEBUG( "Masking Channels: " << ((m_maskBadChannels)?"true":"false"));
  ATH_MSG_DEBUG( "Using pulse shapes from COOL: " << ((m_useCoolPulseShapes)?"true":"false"));

  ATH_CHECK( m_pulseShapeKey.initialize(m_useCoolPulseShapes) );
  ATH_CHECK( m_samplingFractionKey.initialize() );

  if (m_useCoolPulseShapes) {
    ATH_MSG_DEBUG( "Initializing pulse shape conditions object");
  } else {
    /* Get fine-grained shaping profile (0.5-ns bins) for both gains*/
    m_digitShapeHi = m_tileInfo->digitsFullShapeHi();
    m_digitShapeHi.push_back(0.0);
    m_digitShapeLo = m_tileInfo->digitsFullShapeLo();
    m_digitShapeLo.push_back(0.0);
  }

  //=== Initialize bad channels key
  ATH_CHECK( m_badChannelsKey.initialize(m_maskBadChannels || m_rndmEvtOverlay) );

  m_nShapeHi = m_tileInfo->digitsNBinsHi();
  m_nBinsPerXHi = m_tileInfo->digitsBinsPerXHi();
  m_binTime0Hi = m_tileInfo->digitsTime0BinHi();
  m_timeStepHi = 25.0 / m_nBinsPerXHi;

  m_nShapeLo = m_tileInfo->digitsNBinsLo();
  m_nBinsPerXLo = m_tileInfo->digitsBinsPerXLo();
  m_binTime0Lo = m_tileInfo->digitsTime0BinLo();
  m_timeStepLo = 25.0 / m_nBinsPerXLo;

  m_inputDigitContainerName = m_inputDigitContainerKey.key();
  ATH_CHECK( m_inputDigitContainerKey.initialize(!m_onlyUseContainerName && m_rndmEvtOverlay) );

  if (m_rndmEvtOverlay) {
    m_tileNoise = false;
    m_tileCoherNoise = false;
    m_tileThresh = false;
    m_calibRun = false;
    if (m_allChannels<0) m_allChannels = 2; // create all channels with noise in overlay by default

    ATH_MSG_INFO( "Pileup and/or noise added by overlaying digits of random events");

    // locate the PileUpMergeSvc and initialize our local ptr
    if (m_onlyUseContainerName) {
      ATH_CHECK( service("PileUpMergeSvc", m_mergeSvc) );

      ATH_MSG_INFO( "PileUpMergeSvc successfully initialized");
    }

    ATH_CHECK( m_DQstatusKey.initialize() );

  } else {
    ATH_CHECK( m_DQstatusKey.initialize(!m_DQstatusKey.empty()) );

    if (m_allChannels<0) m_allChannels = 0;                 // do not create all channels by default
    if (m_tileNoise || m_tileCoherNoise) m_allChannels = 2; // unless noise is set to True
    if (msgLvl(MSG::INFO)) {
      msg(MSG::INFO) << "Obtained info from TileInfo" << endmsg;
      msg(MSG::INFO) << "tileNoise=" << ((m_tileNoise) ? "true" : "false")
                     << ", tileCoherNoise=" << ((m_tileCoherNoise) ? "true" : "false")
                     << ", tileThresh=" << ((m_tileThresh) ? "true" : "false");
      if (m_tileThresh)
        msg(MSG::INFO) << ", thresh(hi,lo)=" << m_tileThreshHi << "," << m_tileThreshLo << endmsg;
      else
        msg(MSG::INFO) << endmsg;
    }
  }

  if (m_allChannels>1)
    ATH_MSG_INFO( "Create all channels with noise: true");
  else if (m_allChannels>0)
    ATH_MSG_INFO( "Create all channels without noise: true");
  else
    ATH_MSG_INFO( "Create all channels: false");

  if (m_calibRun) {
    m_filteredDigitsContainerKey = "";
  }

  if (!m_filteredDigitsContainerKey.key().empty()) {
    ATH_MSG_INFO( "Keep digits with hit energy above " << m_filterThreshold / MeV
                  << " MeV in " << m_filteredDigitsContainerKey.key() << " container");
    ATH_MSG_INFO( "Keep digits from MBTS with original G4 hit energy above "
                 << m_filterThresholdMBTS / MeV << " MeV ");

    ATH_CHECK( m_filteredDigitsContainerKey.initialize() );

  } else {
    m_filterThreshold = HUGE_VALL;
    m_filterThresholdMBTS = HUGE_VALL;
  }

  ATH_MSG_DEBUG( "nShapeHi=" << m_nShapeHi
                << " nBinsPerXHi=" << m_nBinsPerXHi
                << " timeStepHi=" << m_timeStepHi
                << " binTime0Hi=" << m_binTime0Hi);

  ATH_MSG_DEBUG( "nShapeLo=" << m_nShapeLo
                << " nBinsPerXLo=" << m_nBinsPerXLo
                << " timeStepLo=" << m_timeStepLo
                << " binTime0Lo=" << m_binTime0Lo);

  // decrease by 1, now they are indexes of last element in a vector
  --m_nShapeHi;
  --m_nShapeLo;

  /* ==================================*/
  // Store HWID's for all 12288 channels (48 channels in each of 64 drawers).
  IdContext drawer_context = m_tileHWID->drawer_context();
  int ndrawers = m_tileHWID->drawer_hash_max();
  const int nchMax = 48; // number of channels per drawer

  ATH_MSG_DEBUG( "ndrawers=" << ndrawers
                 << " nchMax=" << nchMax
                 << " HIGAIN=" << TileID::HIGHGAIN
                 << " LOWGAIN=" << TileID::LOWGAIN);

  /* Store all (12288) Identifiers for the calorimeter adc's for HIGHAIN */
  m_all_ids.reserve(ndrawers);
  for (int dr = 0; dr < ndrawers; ++dr) {
    HWIdentifier drawer_id;
    m_tileHWID->get_id(dr, drawer_id, &drawer_context);

    m_all_ids.push_back(std::make_unique<HWIdentifier[]>(nchMax));
    std::unique_ptr<HWIdentifier[]>& adc_ids = m_all_ids.back();

    int ros = m_tileHWID->ros(drawer_id);
    if (ros > 0) {
      int drawer = m_tileHWID->drawer(drawer_id);
      IdentifierHash idhash;
      m_tileHWID->get_hash(drawer_id, idhash, &drawer_context);
      int drawerIdx = TileCalibUtils::getDrawerIdx(ros, drawer);
      if (drawerIdx != (int)idhash) {
        ATH_MSG_ERROR("drawer " << m_tileHWID->to_string(drawer_id, -2)
                       << " hash " << idhash << " NOT EQUAL to idx " << drawerIdx);
      } else if (msgLvl(MSG::VERBOSE) && m_cabling->connected(ros, drawer)) {
        msg(MSG::VERBOSE) << "drawer " << m_tileHWID->to_string(drawer_id, -2)
                          << " hash " << idhash << endmsg;
      }

      for (int ch = 0; ch < nchMax; ++ch) {
        adc_ids[ch] = m_tileHWID->adc_id(drawer_id, ch, TileID::HIGHGAIN);
      }
    }
  }

  ATH_CHECK( m_hitContainer_DigiHSTruthKey.initialize(m_doDigiTruth) );
  ATH_CHECK( m_digitsContainer_DigiHSTruthKey.initialize(m_doDigiTruth) );

  ATH_CHECK( m_hitContainerKey.initialize() );
  ATH_CHECK( m_digitsContainerKey.initialize() );

  ATH_MSG_INFO( "TileDigitsMaker initialization completed");

  return StatusCode::SUCCESS;
}


StatusCode TileDigitsMaker::execute(const EventContext &ctx) const {
  ATH_MSG_DEBUG( "Executing TileDigitsMaker");

  // Prepare RNG service
  ATHRNG::RNGWrapper* rngWrapper = nullptr;
  CLHEP::HepRandomEngine* rngEngine = nullptr;
  if (m_tileNoise || m_tileCoherNoise || m_rndmEvtOverlay) {
    rngWrapper = m_rndmSvc->getEngine(this, m_randomStreamName);
    rngWrapper->setSeed( m_randomStreamName, ctx );
    rngEngine = rngWrapper->getEngine(ctx);
  }

  SG::ReadCondHandle<TileCalibDataFlt> sampleNoiseHandle(m_sampleNoiseKey, ctx);
  TileSampleNoise sampleNoise(*sampleNoiseHandle);

  bool first = (msgLvl(MSG::VERBOSE) && ctx.evt() == 0 && !m_rndmEvtOverlay );
  if (first) {
    ATH_MSG_VERBOSE( "Dumping 2G noise parameters");
    first = false;
    IdContext drawer_context = m_tileHWID->drawer_context();
    int ndrawers = m_tileHWID->drawer_hash_max();
    const int nchMax = 48; // number of channels per drawer
    for (int dr = 0; dr < ndrawers; ++dr) {
      HWIdentifier drawer_id;
      m_tileHWID->get_id(dr, drawer_id, &drawer_context);
      int ros = m_tileHWID->ros(drawer_id);
      int drawer = m_tileHWID->drawer(drawer_id);
      if (m_cabling->connected(ros, drawer)) {
        IdentifierHash idhash;
        m_tileHWID->get_hash(drawer_id, idhash, &drawer_context);
        for (int ch = 0; ch < nchMax; ++ch) {
          double pedSimHi = sampleNoise.getPed(idhash, ch, TileID::HIGHGAIN);
          double sigmaHi_Hfn1 = sampleNoise.getHfn1(idhash, ch, TileID::HIGHGAIN);
          double sigmaHi_Hfn2 = sampleNoise.getHfn2(idhash, ch, TileID::HIGHGAIN);
          double sigmaHi_Norm = sampleNoise.getHfnNorm(idhash, ch, TileID::HIGHGAIN);
          double pedSimLo = sampleNoise.getPed(idhash, ch, TileID::LOWGAIN);
          double sigmaLo_Hfn1 = sampleNoise.getHfn1(idhash, ch, TileID::LOWGAIN);
          double sigmaLo_Hfn2 = sampleNoise.getHfn2(idhash, ch, TileID::LOWGAIN);
          double sigmaLo_Norm = sampleNoise.getHfnNorm(idhash, ch, TileID::LOWGAIN);
          ATH_MSG_VERBOSE( "Channel " << m_tileHWID->to_string(drawer_id,-2) << "/" << ch
                           << " pedHi="<< pedSimHi
                           << " pedLo="<< pedSimLo
                           << " rmsHi="<< sigmaHi_Hfn1 << "," << sigmaHi_Hfn2 << "," << sigmaHi_Norm
                           << " rmsLo="<< sigmaLo_Hfn1 << "," << sigmaLo_Hfn2 << "," << sigmaLo_Norm);

        }
      }
    }
  }

  // declare array for random number generation
  double Rndm[16];       // Can't use variable size array,
  double RndmLo[16];     // Can't use variable size array,
  double Rndm_dG[1];    // uniform random number for the double gaussian
  double RndmLo_dG[1];  // uniform random number for the double gaussian

  // step1: Get hit container from TES 
  SG::ReadHandle<TileHitContainer> hitContainer(m_hitContainerKey, ctx);
  ATH_CHECK( hitContainer.isValid() );

  SG::ReadHandle<TileHitContainer> hitContainer_DigiHSTruth;
  if(m_doDigiTruth){
      hitContainer_DigiHSTruth = SG::ReadHandle<TileHitContainer> (m_hitContainer_DigiHSTruthKey, ctx);
      ATH_CHECK( hitContainer_DigiHSTruth.isValid() );
  }

  // Zero sums for monitoring.
  int nChSum = 0;
  int nChHiSum = 0;
  int nChLoSum = 0;
  int nChHiAcc = 0;
  int nChLoAcc = 0;
  int nChHiFlt = 0;
  int nChLoFlt = 0;
  int nChHiCut = 0;
  int nChLoCut = 0;
  double echtot_Acc = 0.;
  double echint_Acc = 0.;
  double echtot_Cut = 0.;
  double echint_Cut = 0.;
  double HitSum = 0.;
  double EneSum = 0.;
  double RChSum = 0.;

  /* step2: Set up  Digits container */

  auto digitsContainer = std::make_unique<TileMutableDigitsContainer>(true,
                                                                      TileFragHash::Digitizer,
                                                                      TileRawChannelUnit::ADCcounts);
  ATH_CHECK( digitsContainer->status() );

  std::unique_ptr<TileMutableDigitsContainer> digitsContainer_DigiHSTruth;
  if(m_doDigiTruth){
    digitsContainer_DigiHSTruth = std::make_unique<TileMutableDigitsContainer>(true,
                                                                               TileFragHash::Digitizer,
                                                                               TileRawChannelUnit::ADCcounts);
    ATH_CHECK( digitsContainer_DigiHSTruth->status() );
  }

  std::unique_ptr<TileMutableDigitsContainer> filteredContainer;
  if (!m_filteredDigitsContainerKey.key().empty()) {
    filteredContainer = std::make_unique<TileMutableDigitsContainer>(true,
                                                                     TileFragHash::Digitizer,
                                                                     TileRawChannelUnit::ADCcounts,
                                                                     SG::VIEW_ELEMENTS);
    ATH_CHECK( filteredContainer->status() );
  }

  /* Set up buffers for handling information in a single collection. */
  IdentifierHash idhash;
  IdContext drawer_context = m_tileHWID->drawer_context();
  const int nchMax = 48; // number of channels per drawer
  std::vector<int> igain(nchMax, -1);
  std::vector<int> ntot_ch(nchMax, 0);
  std::vector<double> ech_tot(nchMax, 0.0);
  std::vector<double> ech_int(nchMax, 0);
  std::vector<double> ech_int_DigiHSTruth(nchMax, 0);
  std::vector<int> over_gain(nchMax, -1);

  /* Make a vector of digits (to be filled at the end from m_drawerBuffer arrays) */
  std::vector<float> digitsBuffer(m_nSamples);
  std::vector<float> digitsBufferLo(m_nSamples); // for calib runs
  std::vector<float> digitsBuffer_DigiHSTruth(m_nSamples);
  std::vector<float> digitsBufferLo_DigiHSTruth(m_nSamples); // for calib runs

  std::vector<double> emptyBuffer;
  std::vector<std::vector<double>> drawerBufferHi(nchMax, std::vector<double>(m_nSamples));
  std::vector<std::vector<double>> drawerBufferLo(nchMax, std::vector<double>(m_nSamples));

  std::vector<std::vector<double>> drawerBufferHi_DigiHSTruth;
  std::vector<std::vector<double>> drawerBufferLo_DigiHSTruth;
  if (m_doDigiTruth) {
    drawerBufferHi_DigiHSTruth.resize(nchMax, std::vector<double>(m_nSamples));
    drawerBufferLo_DigiHSTruth.resize(nchMax, std::vector<double>(m_nSamples));
  }


  /* everything for calculation of coherent noise */
  // booleans for coherent noise
  Bool_t coherNoiseHi = false;
  Bool_t coherNoiseLo = false;
  TMatrixD CorrWeightHi;
  TMatrixD CorrWeightLo;
  std::vector<std::unique_ptr<double[]>> CorrRndmVec;
  std::vector<std::unique_ptr<double[]>> CorrRndmVecLo;
  if (m_tileCoherNoise) {
    for (int k = 0; k < m_nSamples; ++k) {
      CorrRndmVec.push_back(std::make_unique<double[]>(nchMax));
    }
    if (m_calibRun) {
      for (int k = 0; k < m_nSamples; ++k) {
        CorrRndmVecLo.push_back(std::make_unique<double[]>(nchMax));
      }
    }
  }

  TileMutableDigitsContainer::const_iterator collItrRndm;
  TileMutableDigitsContainer::const_iterator lastCollRndm;
  std::unique_ptr<TileMutableDigitsContainer> backgroundDigitContainer{};
  if (m_rndmEvtOverlay) {
    backgroundDigitContainer = std::make_unique<TileMutableDigitsContainer>(true,
                                                                            TileFragHash::Digitizer,
                                                                            TileRawChannelUnit::ADCcounts);
    ATH_CHECK( backgroundDigitContainer->status() );

    if (m_onlyUseContainerName) {
      typedef PileUpMergeSvc::TimedList<TileDigitsContainer>::type TimedDigitContList;
      TimedDigitContList digitContList;
      ATH_CHECK( m_mergeSvc->retrieveSubEvtsData(m_inputDigitContainerName, digitContList));
      ATH_MSG_DEBUG( "TileDigitsCnt successfully retrieved ");


      if (digitContList.size() == 0) {
        ATH_MSG_WARNING( "No overlay done ... ");
        return StatusCode::SUCCESS;
      }

      TimedDigitContList::iterator iTzeroDigitCont(digitContList.begin());
      for (const auto* digitCollection : *(iTzeroDigitCont->second)) {
        for (const auto* digit : *digitCollection) {
          auto pDigits = std::make_unique<TileDigits>(*digit);
          ATH_CHECK(backgroundDigitContainer->push_back(std::move(pDigits)));
        }
      }
    }
    else {
      SG::ReadHandle<TileDigitsContainer> tileDigitsContainerHandle(m_inputDigitContainerKey, ctx);
      if (tileDigitsContainerHandle.isValid()) {
        for (const auto* digitCollection : *tileDigitsContainerHandle) {
          for (const auto* digit : *digitCollection) {
            auto pDigits = std::make_unique<TileDigits>(*digit);
            ATH_CHECK(backgroundDigitContainer->push_back(std::move(pDigits)));
          }
        }
      }
      else {
        ATH_MSG_ERROR("ReadHandle to Background Digits is invalid.");
        return StatusCode::FAILURE;
      }
    }

    collItrRndm = backgroundDigitContainer->begin();
    lastCollRndm = backgroundDigitContainer->end();
  }

  SG::ReadCondHandle<TileEMScale> emScale(m_emScaleKey, ctx);
  ATH_CHECK( emScale.isValid() );

  const TileCalibDataFlt* pulseShape = nullptr;
  if (m_useCoolPulseShapes) {
    SG::ReadCondHandle<TileCalibDataFlt> pulseShapeHandle(m_pulseShapeKey, ctx);
    ATH_CHECK( pulseShapeHandle.isValid() );
    pulseShape = pulseShapeHandle.retrieve();
  }
  TilePulse pulse(pulseShape);

  SG::ReadCondHandle<TileSamplingFraction> samplingFraction(m_samplingFractionKey, ctx);
  ATH_CHECK( samplingFraction.isValid() );

  const TileDQstatus* dqStatus = nullptr;
  if (m_rndmEvtOverlay) {
    SG::ReadHandle<TileDQstatus> DQstatusHandle(m_DQstatusKey, ctx);
    ATH_CHECK( DQstatusHandle.isValid() );
    dqStatus = DQstatusHandle.get();
  }

  const TileBadChannels* badChannels = nullptr;
  if (m_maskBadChannels || m_rndmEvtOverlay) {
    SG::ReadCondHandle<TileBadChannels> badChannelsHandle(m_badChannelsKey,ctx);
    ATH_CHECK( badChannelsHandle.isValid() );
    badChannels = badChannelsHandle.retrieve();
  }

  // iterate over all collections in a container
  // Hit Container and signal hit container are the same size (1 entry per channel)
  TileHitContainer::const_iterator collItr_DigiHSTruth;
  if(m_doDigiTruth) collItr_DigiHSTruth = hitContainer_DigiHSTruth->begin();

  /* ----------------------------------------------------------------- */
  /* Begin loop over the Hit collections.  All collections are defined */
  /* (even if they have no hits), and all the digit information        */
  /* including pileup events are contained in the collection.          */
  /*-------------------------------------------------------------------*/
  for (const TileHitCollection* hitCollection : *hitContainer) {
    /* Get array of HWID's for this drawer (stored locally). */
    HWIdentifier drawer_id = m_tileHWID->drawer_id(hitCollection->identify());
    int ros = m_tileHWID->ros(drawer_id);
    int drawer = m_tileHWID->drawer(drawer_id);
    int drawerIdx = TileCalibUtils::getDrawerIdx(ros, drawer);
    if (m_cabling->connected(ros, drawer)) {
      ATH_MSG_VERBOSE( "ROS "<< ros << " drawer " << drawer << " is connected");
    } else {
      if (m_rndmEvtOverlay && collItrRndm != lastCollRndm) {
        ++collItrRndm; // skip also one drawer in digi overlay container
      }
      if (m_doDigiTruth) {
        ++collItr_DigiHSTruth;
      } // End DigiHSTruth stuff
      continue;
    }

    m_tileHWID->get_hash(drawer_id, idhash, &drawer_context);
    const std::unique_ptr<HWIdentifier[]>& adc_ids = m_all_ids[idhash];

    /* Initialize gain settings.  If noise is requested, all channels are */
    /* set to be active.  If not, set them all to be inactive (gain=-1).  */
    /* Only those which contain actual hits will be set active when the   */
    /* hits are read in.                                                  */
    int igainch = (m_allChannels) ? TileID::HIGHGAIN : -1;
    if (m_rndmEvtOverlay) {
      std::fill(over_gain.begin(), over_gain.end(), -1);
    } else if (m_tileNoise || m_tileCoherNoise) {
      igainch = TileID::HIGHGAIN;
    }

    std::fill(ech_tot.begin(), ech_tot.end(), 0.0);
    std::fill(ech_int.begin(), ech_int.end(), 0.0);
    std::fill(ntot_ch.begin(), ntot_ch.end(), 0);
    std::fill(igain.begin(), igain.end(), igainch);

    std::vector<std::reference_wrapper<std::vector<std::vector<double>>>> drawerBuffers{drawerBufferHi, drawerBufferLo};
    if (m_doDigiTruth) {
      drawerBuffers.push_back(drawerBufferHi_DigiHSTruth);
      drawerBuffers.push_back(drawerBufferLo_DigiHSTruth);
    }
    for (std::vector<std::vector<double>>& drawerBuffer : drawerBuffers) {
      for (std::vector<double>& digitsBuffer : drawerBuffer) {
        std::fill(digitsBuffer.begin(), digitsBuffer.end(), 0);
      }
    }

    if (m_rndmEvtOverlay && collItrRndm != lastCollRndm) {
      const TileDigitsCollection *bkgDigitCollection(*collItrRndm);
      ATH_CHECK(overlayBackgroundDigits(bkgDigitCollection, hitCollection, drawerBufferLo, drawerBufferHi,
                                        igain, ros, drawer, drawerIdx, over_gain, *emScale, sampleNoise, dqStatus, badChannels));
      ++collItrRndm; // skip to next digi collection
    }

    std::vector<bool> signal_in_channel(nchMax, false);
    std::vector<bool> signal_in_channel_DigiHSTruth(nchMax, false);
    ATH_CHECK(fillDigitCollection( hitCollection, drawerBufferLo, drawerBufferHi,
                                   igain, over_gain, ech_int, signal_in_channel, *emScale, *samplingFraction, pulse));
    if(m_doDigiTruth){
      ATH_CHECK(fillDigitCollection( *collItr_DigiHSTruth, drawerBufferLo_DigiHSTruth, drawerBufferHi_DigiHSTruth,
                                     igain, over_gain, ech_int_DigiHSTruth, signal_in_channel_DigiHSTruth, *emScale, *samplingFraction, pulse));
    } // End DigiHSTruth stuff

    /* Now all signals for this collection are stored in m_drawerBuffer, 
     accessed with digitSamplesHi and digitSampleLo. */
    if (msgLvl(MSG::VERBOSE)) {
      for (int ich = 0; ich < nchMax; ++ich) {
        if (igain[ich] > -1) {
          std::vector<double>& digitSamplesHi = drawerBufferHi[ich];
          std::vector<double>& digitSamplesLo = drawerBufferLo[ich];
          msg(MSG::VERBOSE) << "total:  ADC " << m_tileHWID->to_string(adc_ids[ich],-1) << "/" << igain[ich] 
                            << " nhit=" << ntot_ch[ich]
                            << " e_ch=" << ech_tot[ich]
                            << " AinTHi=" << digitSamplesHi[m_iTrig]
                            << " AinTLo=" << digitSamplesLo[m_iTrig] << endmsg;
        }
      }
    }

    /* ---------------------------------------------------------------      */
    /* Now all signals for this drawer are stored in m_drawerBuffer arrays, */
    /* and we are finished with TileHits for this collection.  Loop over    */
    /* channels to add noise and pedestal.  Check for saturation for        */
    /* each channel, and in case of saturation convert to low gain.         */
    /* -------------------------------------------------------------- */

    // =============CORRELATION MODIFICATION (F Spano)============== 
    //
    // Define CoVariance Matrix corresponding to noise.
    // TO UPDATE:: Such Matrix will have to be loaded from database in the future; 1 matrix per drawer.
    // NOW: 
    // a) define one covariance matrix
    // b) find Cholesky decomposition use for corrlation building
    // if this is set load the matrix
    if (m_tileCoherNoise) {
      ATH_MSG_VERBOSE( "Coherent noise for ROS " << ros
                       << " drawer " << drawer
                       << " with " << nchMax << " channels and "
                       << m_nSamples << "samples ");

      // get decomposed covariance matrix for hi gain
      coherNoiseHi = 1;
      if (coherNoiseHi) {
        CorrWeightHi.ResizeTo(*(m_tileInfo->DecoCovariance(ros, drawer, TileID::HIGHGAIN)));
        CorrWeightHi = *(m_tileInfo->DecoCovariance(ros, drawer, TileID::HIGHGAIN));
      }

      // get decomposed covariance matrix for low gain
      coherNoiseLo = 1;
      if (coherNoiseLo) {
        CorrWeightLo.ResizeTo(*(m_tileInfo->DecoCovariance(ros, drawer, TileID::LOWGAIN)));
        CorrWeightLo = *(m_tileInfo->DecoCovariance(ros, drawer, TileID::LOWGAIN));
      }

      //NOTE: ShootArray's inputs are : the engine, the size, the vector, the mean, the standard dev
      for (int k = 0; k < m_nSamples; ++k) {
        double* RndmVec = CorrRndmVec[k].get();
        RandGaussQ::shootArray(rngEngine, nchMax, RndmVec, 0.0, 1.0);
      }

      if (m_calibRun) {
        for (int k = 0; k < m_nSamples; ++k) {
          double * RndmVecLo = CorrRndmVecLo[k].get();
          RandGaussQ::shootArray(rngEngine, nchMax, RndmVecLo, 0.0, 1.0);
        }
      }
    }
    // =============CORRELATION MODIFICATION (F Spano)============== end

    // looping over channels
    for (int ich = 0; ich < nchMax; ++ich) {
      /* If igain<0, channel is inactive => skip it.                    */
      if (igain[ich] < 0)
        continue;

      /* Generate the nSamp Digits for high gain. Check each for saturation. */
      ++nChHiSum;
      HWIdentifier adc_id = adc_ids[ich];
      HWIdentifier adc_id_lo; // for calib runs
      Identifier pmt_id = m_cabling->h2s_pmt_id(adc_id);
      ATH_MSG_DEBUG( "Ch " << m_tileHWID->to_string(adc_id,-1)
                     << " PMT " << (pmt_id.is_valid() ? m_tileID->to_string(pmt_id,-1) : (signal_in_channel[ich] ? "fake gap" : "not connected"))
                     << " gain=" << igain[ich]);

      if (m_calibRun || m_maskBadChannels) {
        adc_id_lo = m_tileHWID->adc_id(drawer_id, ich, TileID::LOWGAIN);
        if (m_calibRun) {
          ++nChLoSum;
        }
      }

      bool chanLoIsBad = false;
      bool chanHiIsBad = false;
      if (m_maskBadChannels) {
        TileBchStatus statusLo = badChannels->getAdcStatus(adc_id_lo);
        TileBchStatus statusHi = badChannels->getAdcStatus(adc_id);
        chanLoIsBad = statusLo.isBad();
        chanHiIsBad = statusHi.isBad();
      }

      /* Get pedestal and noise values */
      double pedSimHi(0.), sigmaHi_Hfn1(0.), sigmaHi_Hfn2(0.), sigmaHi_Norm(0.), pedSimLo(0.),
          sigmaLo_Hfn1(0.), sigmaLo_Hfn2(0.), sigmaLo_Norm(0.);
      bool good_ch = (over_gain[ich]<9);
      bool overNoiseHG(over_gain[ich]!=TileID::HIGHGAIN && good_ch); // it's always true if no overlay
      bool overNoiseLG(over_gain[ich]!=TileID::LOWGAIN  && good_ch); // it's always true if no overlay
      bool tileNoiseHG(false),tileNoiseLG(false);

      if (overNoiseHG) {
        overNoiseHG &= (m_rndmEvtOverlay && m_allChannels>1); // set it to true only for overlay
        tileNoiseHG = m_tileNoise || overNoiseHG;

        pedSimHi = sampleNoise.getPed(idhash, ich, TileID::HIGHGAIN);
        // bug fix for wrong ped value in DB
        if (pedSimHi == 0.0 && (signal_in_channel[ich] || pmt_id.is_valid()))
          pedSimHi = 50.;

        sigmaHi_Hfn1 = sampleNoise.getHfn1(idhash, ich, TileID::HIGHGAIN);
        sigmaHi_Hfn2 = sampleNoise.getHfn2(idhash, ich, TileID::HIGHGAIN);
        if (sigmaHi_Hfn1>0 || sigmaHi_Hfn2) {
          sigmaHi_Norm = sigmaHi_Hfn1 / (sigmaHi_Hfn1
                       + sigmaHi_Hfn2 * sampleNoise.getHfnNorm(idhash, ich, TileID::HIGHGAIN));
        } else {
          sigmaHi_Hfn1 = sampleNoise.getHfn(idhash, ich, TileID::HIGHGAIN);
          sigmaHi_Norm = 1.;
        }
      }
      
      if (overNoiseLG) {
        overNoiseLG &= (m_rndmEvtOverlay && m_allChannels>1); // set it to true only for overlay
        tileNoiseLG = m_tileNoise || overNoiseLG;

        pedSimLo = sampleNoise.getPed(idhash, ich, TileID::LOWGAIN);
        // bug fix for wrong ped value in DB
        if (pedSimLo == 0.0 && (signal_in_channel[ich] || pmt_id.is_valid()))
          pedSimLo = 30.;

        sigmaLo_Hfn1 = sampleNoise.getHfn1(idhash, ich, TileID::LOWGAIN);
        sigmaLo_Hfn2 = sampleNoise.getHfn2(idhash, ich, TileID::LOWGAIN);
        if (sigmaLo_Hfn1 > 0 || sigmaLo_Hfn2) {
          sigmaLo_Norm = sigmaLo_Hfn1 / (sigmaLo_Hfn1
                       + sigmaLo_Hfn2 * sampleNoise.getHfnNorm(idhash, ich, TileID::LOWGAIN));
        } else {
          sigmaLo_Hfn1 = sampleNoise.getHfn(idhash, ich, TileID::LOWGAIN);
          sigmaLo_Norm = 1.;
        }
      }

      /* If tileNoise is requested, generate array of random numbers.   */
      if (tileNoiseLG) { // true if tileNoise is set or noise is needed for low gain in overlay
        RandGaussQ::shootArray(rngEngine, m_nSamples, Rndm, 0.0, 1.0);
        RandFlat::shootArray(rngEngine, 1, Rndm_dG, 0.0, 1.0);
        if (m_calibRun) {
          RandGaussQ::shootArray(rngEngine, m_nSamples, RndmLo, 0.0, 1.0);
          RandFlat::shootArray(rngEngine, 1, RndmLo_dG, 0.0, 1.0);
        }
      }

      std::vector<double>& digitSamplesHi = drawerBufferHi[ich];
      std::vector<double>& digitSamplesLo = drawerBufferLo[ich];
      std::vector<double>& digitSamplesHi_DigiHSTruth = (m_doDigiTruth) ? drawerBufferHi_DigiHSTruth[ich] : emptyBuffer;
      std::vector<double>& digitSamplesLo_DigiHSTruth = (m_doDigiTruth) ? drawerBufferLo_DigiHSTruth[ich] : emptyBuffer;

      ATH_MSG_DEBUG(" Channel " << ros << '/' << drawer << '/' << ich 
                     << " sampHi=" << digitSamplesHi[m_iTrig]
                     << " pedHi=" << pedSimHi
                     << " sampLo=" << digitSamplesLo[m_iTrig]
                     << " pedLo=" << pedSimLo);

      // looping over samples
      for (int js = 0; js < m_nSamples; ++js) {

        digitsBuffer[js] = digitSamplesHi[js] + pedSimHi;
        if(m_doDigiTruth) {
          digitsBuffer_DigiHSTruth[js] = digitSamplesHi_DigiHSTruth[js] + pedSimHi;
        }

        double noiseHi(0.0);
        // Full noise pattern, including coherent noise has priority over normal noise //F Spano'
        if (coherNoiseHi) {
          // get the js-th  correct random vector of 48 elements for the jsth sample k  //F Spano'
          std::unique_ptr<double[]>& CorVec = CorrRndmVec[js];
          // apply Y=C*Z where Z is the random vector of 48 normal indep variables, and C is the Cholesky decomposition //F Spano'
          for (int i = 0; i < nchMax; ++i) noiseHi += CorrWeightHi(i, ich) * CorVec[i];
        } else if (tileNoiseHG) {
          //using the same gaussian(sigma) for all samples in one channel in one event
          if (Rndm_dG[0] < sigmaHi_Norm) noiseHi = sigmaHi_Hfn1 * Rndm[js];
          else noiseHi = sigmaHi_Hfn2 * Rndm[js];
        }

        if (digitsBuffer[js] + noiseHi >= 0.0) {
          digitsBuffer[js] += noiseHi;
          if(m_doDigiTruth) digitsBuffer_DigiHSTruth[js] += noiseHi;
        } else {
          digitsBuffer[js] -= noiseHi;
          if(m_doDigiTruth) digitsBuffer_DigiHSTruth[js] -= noiseHi;
        }


        if (m_integerDigits) {
          digitsBuffer[js] = round(digitsBuffer[js]);
          if(m_doDigiTruth) digitsBuffer_DigiHSTruth[js] = round(digitsBuffer_DigiHSTruth[js]);
        }

        if (m_calibRun) { //Calculate also low gain
          digitsBufferLo[js] = digitSamplesLo[js] + pedSimLo;
          if(m_doDigiTruth) digitsBufferLo_DigiHSTruth[js] = digitSamplesLo_DigiHSTruth[js] + pedSimLo;
          double noiseLo(0.0);
          // Full noise pattern, including coherent noise has priority over normal noise //F Spano'
          if (coherNoiseLo) {
            // get the js-th  correct random vector of 48 elements for the jsth sample // F Spano'
            std::unique_ptr<double[]>& CorVecLo = CorrRndmVecLo[js];
            // apply Y=C*Z where Z is the random vector of 48 normal indep variables, and C is the Cholesky decomposition // F Spano'
            for (int i = 0; i < nchMax; ++i) noiseLo += CorrWeightLo(i, ich) * CorVecLo[i];
          } else if (tileNoiseLG) {
            //using the same gaussian (sigma) for all samples in one channel in one event
            if (RndmLo_dG[0] < sigmaLo_Norm) noiseLo = sigmaLo_Hfn1 * RndmLo[js];
            else noiseLo = sigmaLo_Hfn2 * RndmLo[js];
          }
          
          if (digitsBufferLo[js] + noiseLo >= 0.0) {
            digitsBufferLo[js] += noiseLo;
            if(m_doDigiTruth) digitsBufferLo_DigiHSTruth[js] += noiseLo;
          } else {
            digitsBufferLo[js] -= noiseLo;
            if(m_doDigiTruth) digitsBufferLo_DigiHSTruth[js] -= noiseLo;
          }
          
          if (m_integerDigits) {
            digitsBufferLo[js] = round(digitsBufferLo[js]);
            if(m_doDigiTruth) digitsBufferLo_DigiHSTruth[js] = round(digitsBufferLo_DigiHSTruth[js]);
          }


        } else if ((digitsBuffer[js] >= m_f_ADCmaxHG && good_ch) || igain[ich] == TileID::LOWGAIN) { // saturation of high gain in non-calib run
                                                                                                 // or low gain in digi overlay
          --nChHiSum;
          ++nChLoSum;
          igain[ich] = TileID::LOWGAIN;
          adc_id = m_tileHWID->adc_id(drawer_id, ich, TileID::LOWGAIN);

          // reset all samples in digitsBuffer[] to Low Gain values
          for (js = 0; js < m_nSamples; ++js) {
            digitsBuffer[js] = digitSamplesLo[js] + pedSimLo;
            if(m_doDigiTruth) digitsBuffer_DigiHSTruth[js] = digitSamplesLo_DigiHSTruth[js] + pedSimLo;
            double noiseLo(0.0);
            // Full noise pattern, including coherent noise has priority over normal noise //F Spano'
            if (coherNoiseLo) {
              // get the js-th  correct random vector of 48 elements for the jsth sample // F Spano'
              double* CorVec = CorrRndmVec[js].get(); // reuse the same rndm as for high gain
              // apply Y=C*Z where Z is the random vector of 48 normal indep variables, and C is the Cholesky decomposition // F Spano'
              for (int i = 0; i < nchMax; ++i) noiseLo += CorrWeightLo(i, ich) * CorVec[i];
            } else if (tileNoiseLG) {
              //using the same gaussian (sigma) for all samples in one channel in one event
	            // reuse the same rndm as for high gain
              if (Rndm_dG[0] < sigmaLo_Norm) noiseLo = sigmaLo_Hfn1 * Rndm[js];
              else noiseLo = sigmaLo_Hfn2 * Rndm[js];
            }

            if (digitsBuffer[js] + noiseLo >= 0.0) {
              digitsBuffer[js] += noiseLo;
              if(m_doDigiTruth) digitsBuffer_DigiHSTruth[js] += noiseLo;
            } else {
              digitsBuffer[js] -= noiseLo;
              if(m_doDigiTruth) digitsBuffer_DigiHSTruth[js] -= noiseLo;
            }

            if (digitsBuffer[js] > m_f_ADCmax && good_ch) {
              digitsBuffer[js] = m_f_ADCmax;
              if(m_doDigiTruth) digitsBuffer_DigiHSTruth[js] = m_f_ADCmax;
            }
            if (m_integerDigits) {
              digitsBuffer[js] = round(digitsBuffer[js]);
              if(m_doDigiTruth) digitsBuffer_DigiHSTruth[js] = round(digitsBuffer_DigiHSTruth[js]);
            }
          }

          overNoiseHG = false;

          if (msgLvl(MSG::VERBOSE)) {
            msg(MSG::VERBOSE) << "Channel " << ros << '/' << drawer << '/' << ich << "/" << igain[ich]
                              << "  Switch to low gain  Amp(lo)=" << digitsBuffer[m_iTrig] << endmsg;
            if (overNoiseLG) {
              if (sigmaLo_Norm<1.0) {
                msg(MSG::VERBOSE) << "LG Ped & noise from DB "
                                  << pedSimLo << " " << sigmaLo_Hfn1 << " " << sigmaLo_Hfn2 << " " << sigmaLo_Norm
                                  << ((Rndm_dG[0] < sigmaLo_Norm)?(" sig1 used"):(" sig2 used")) << endmsg;
              } else {
                msg(MSG::VERBOSE) << "LG Ped & noise from DB "
                                  << pedSimLo << " " << sigmaLo_Hfn1 << endmsg;
              }
            }
          }
          break;
        }
      }
      if (msgLvl(MSG::VERBOSE)) {
        if (overNoiseHG) {
          if (sigmaHi_Norm<1.0) {
            msg(MSG::VERBOSE) << "HG Ped & noise from DB " 
                              << pedSimHi << " " << sigmaHi_Hfn1 << " " << sigmaHi_Hfn2 << " " << sigmaHi_Norm
                              << ((Rndm_dG[0] < sigmaHi_Norm)?(" sig1 used"):(" sig2 used")) << endmsg;
          } else {
            msg(MSG::VERBOSE) << "HG Ped & noise from DB " 
                              << pedSimHi << " " << sigmaHi_Hfn1 << endmsg;
          }
        }
        if (m_calibRun && overNoiseLG) {
          if (sigmaLo_Norm<1.0) {
            msg(MSG::VERBOSE) << "LG Ped & noise from DB "
                              << pedSimLo << " " << sigmaLo_Hfn1 << " " << sigmaLo_Hfn2 << " " << sigmaLo_Norm
                              << ((RndmLo_dG[0] < sigmaLo_Norm)?(" sig1 used"):(" sig2 used")) << endmsg;
          } else {
            msg(MSG::VERBOSE) << "LG Ped & noise from DB "
                              << pedSimLo << " " << sigmaLo_Hfn1 << endmsg;
          }
        }
      }

      if (m_calibRun) { // calib run - keep both low and high gain

        if (chanHiIsBad) {
          std::fill(digitsBuffer.begin(), digitsBuffer.end(), m_f_ADCmaskValue);
          if (m_doDigiTruth) {
            std::fill(digitsBuffer_DigiHSTruth.begin(), digitsBuffer_DigiHSTruth.end(), m_f_ADCmaskValue);
          }
          ATH_MSG_DEBUG( "Masking Channel " << ros << '/' << drawer << '/' << ich << "/1 HG" );
        }

        auto pDigits = std::make_unique<TileDigits>(adc_id, digitsBuffer);
        ATH_CHECK( digitsContainer->push_back(std::move(pDigits)) );

        if(m_doDigiTruth && digitsContainer_DigiHSTruth){
          auto digits_DigiHSTruth = std::make_unique<TileDigits>(adc_id, digitsBuffer_DigiHSTruth);
          ATH_CHECK( digitsContainer_DigiHSTruth->push_back(std::move(digits_DigiHSTruth)) );
        }

        if (chanLoIsBad) {
          std::fill(digitsBufferLo.begin(), digitsBufferLo.end(), m_f_ADCmaskValue);
          if(m_doDigiTruth) {
            std::fill(digitsBufferLo_DigiHSTruth.begin(), digitsBufferLo_DigiHSTruth.end(), m_f_ADCmaskValue);
          }

          ATH_MSG_DEBUG( "Masking Channel " << ros << '/' << drawer << '/' << ich << "/0 LG");
        }

        auto pDigitsLo = std::make_unique<TileDigits>(adc_id_lo, digitsBufferLo);
        ATH_CHECK( digitsContainer->push_back(std::move(pDigitsLo)) );

        if(m_doDigiTruth && digitsContainer_DigiHSTruth){
          auto pDigitsLo_DigiHSTruth = std::make_unique<TileDigits>(adc_id_lo, digitsBufferLo_DigiHSTruth);
          ATH_CHECK( digitsContainer_DigiHSTruth->push_back(std::move(pDigitsLo_DigiHSTruth)) );
        }
      } else { //normal run

        bool hiGain = (igain[ich] == TileID::HIGHGAIN);

        // If tileThresh, apply threshold cut to the in-time Digits signal
        bool isChannelGood = true;
        if (m_tileThresh) {
          if (hiGain) { // make threshold only on high gain
            double ampInTime = digitsBuffer[m_iTrig] - pedSimHi;
            if (m_integerDigits)
              ampInTime = round(ampInTime);
            if (m_tileThreshHi < 0) {
              if (fabs(ampInTime) < fabs(m_tileThreshHi))
                isChannelGood = false;
            } else {
              if (ampInTime < m_tileThreshHi)
                isChannelGood = false;
            }
          }
        }
        // If channel is good, create TileDigits object and store in container.
        if (isChannelGood) {
          echtot_Acc += ech_tot[ich];
          echint_Acc += fabs(ech_int[ich]);
          if (hiGain) {
            ++nChHiAcc;
          } else {
            ++nChLoAcc;
          }

          if (hiGain) {
            //if (m_rndmEvtOverlay // not needed, because DQstatus have been already checked before
            //    && !(theDQstatus->isAdcDQgood(ros, drawer, ich, TileID::HIGHGAIN))) {
            //  chanHiIsBad = true;
            //  ATH_MSG_DEBUG( "BAD DQ Channel " << ros << '/' << drawer << '/' << ich << "/1 HG");
            //}
            if (chanHiIsBad) {
              if (pmt_id.is_valid()) {
                std::fill(digitsBuffer.begin(), digitsBuffer.end(), m_f_ADCmaskValue);
                if (m_doDigiTruth) {
                  std::fill(digitsBuffer_DigiHSTruth.begin(), digitsBuffer_DigiHSTruth.end(), m_f_ADCmaskValue);
                }
              } else if (good_ch) {
                ATH_MSG_DEBUG( "Disconnected Channel " << ros << '/' << drawer << '/' << ich);
                std::fill(digitsBuffer.begin(), digitsBuffer.end(), 0.);
                if (m_doDigiTruth) {
                  std::fill(digitsBuffer_DigiHSTruth.begin(), digitsBuffer_DigiHSTruth.end(), 0.);
                }
              }
              ATH_MSG_DEBUG( "Masking Channel " << ros << '/' << drawer << '/' << ich << "/1 HG");
            }
          } else {
            //if (m_rndmEvtOverlay // not needed, because DQstatus have been already checked before
            //    && !(theDQstatus->isAdcDQgood(ros, drawer, ich, TileID::LOWGAIN))) {
            //  chanLoIsBad = true;
            //  ATH_MSG_DEBUG( "BAD DQ Channel " << ros << '/' << drawer << '/' << ich << "/0 LG");
            //}
            if (chanLoIsBad) {
              if (pmt_id.is_valid()) {
                std::fill(digitsBuffer.begin(), digitsBuffer.end(), m_f_ADCmaskValue);
                if (m_doDigiTruth) {
                  std::fill(digitsBuffer_DigiHSTruth.begin(), digitsBuffer_DigiHSTruth.end(), m_f_ADCmaskValue);
                }
              } else if (good_ch) {
                ATH_MSG_DEBUG( "Disconnected Channel " << ros << '/' << drawer << '/' << ich);
                std::fill(digitsBuffer.begin(), digitsBuffer.end(), 0.);
                if (m_doDigiTruth) {
                  std::fill(digitsBuffer_DigiHSTruth.begin(), digitsBuffer_DigiHSTruth.end(), 0.);
                }
              }
              ATH_MSG_DEBUG( "Masking Channel " << ros << '/' << drawer << '/' << ich << "/0 LG");
            }
          }

          auto pDigits = std::make_unique<TileDigits>(adc_id, digitsBuffer);

          if (ech_int[ich] > m_filterThreshold || ech_int[ich] < -m_filterThresholdMBTS) {
            if (filteredContainer) ATH_CHECK( filteredContainer->push_back(pDigits.get()) );
            if (hiGain) {
              ++nChHiFlt;
            } else {
              ++nChLoFlt;
            }
          }

          ATH_CHECK( digitsContainer->push_back(std::move(pDigits)) );
          if(m_doDigiTruth && digitsContainer_DigiHSTruth){
            auto pDigits_DigiHSTruth = std::make_unique<TileDigits>(adc_id, digitsBuffer_DigiHSTruth);
            ATH_CHECK( digitsContainer_DigiHSTruth->push_back(std::move(pDigits_DigiHSTruth)) );
          }

          if (msgLvl(MSG::VERBOSE)) {
            double pedSim = ((hiGain) ? pedSimHi : pedSimLo);
            double ampInTime = digitsBuffer[m_iTrig] - pedSim;
            if (m_integerDigits)
              ampInTime = round(ampInTime);
            msg(MSG::VERBOSE) << ((ech_int[ich] > m_filterThreshold
                                  || ech_int[ich] < -m_filterThresholdMBTS) ? "AccFlt" : "Accept")
                              << " ADC " << m_tileHWID->to_string(adc_id)
                              << " AinT=" << ampInTime
                              << " ped=" << pedSim
                              << " Ech=" << ech_tot[ich]
                              << " EinT=" << ech_int[ich] << endmsg;
            msg(MSG::VERBOSE) << "digits";
            for (unsigned int i = 0; i < digitsBuffer.size(); ++i)
              msg(MSG::VERBOSE) << " " << digitsBuffer[i];
            msg(MSG::VERBOSE) << endmsg;
          }
        } else {
          echtot_Cut += ech_tot[ich];
          echint_Cut += ech_int[ich];
          if (hiGain) {
            ++nChHiCut;
          } else {
            ++nChLoCut;
          }

          if (msgLvl(MSG::VERBOSE)) {
            double pedSim = ((hiGain) ? pedSimHi : pedSimLo);
            double ampInTime = digitsBuffer[m_iTrig] - pedSim;
            if (m_integerDigits)
              ampInTime = round(ampInTime);
            msg(MSG::VERBOSE) << "Reject. ADC " << m_tileHWID->to_string(adc_id)
                              << " AinT=" << ampInTime
                              << " ped=" << pedSim
                              << " Ech=" << ech_tot[ich]
                              << " EinT=" << ech_int[ich] << endmsg;
          }
        }
      }
    }
    if(m_doDigiTruth) ++collItr_DigiHSTruth;
  }

  if (msgLvl(MSG::DEBUG)) {
    msg(MSG::DEBUG) << "TileDigitsMaker execution completed." << endmsg;
    msg(MSG::DEBUG) << " nCh=" << nChSum
                    << " nChH/L=" << nChHiSum << "/" << nChLoSum
                    << " nFltH/L=" << nChHiFlt << "/" << nChLoFlt
                    << " Hit=" << HitSum
                    << " Ene=" << EneSum
                    << " RChSum=" << RChSum << endmsg;
    if (m_tileThresh) {
      msg(MSG::DEBUG) << " Accepted:  nChLo/Hi=" << nChLoAcc << "/" << nChHiAcc
                      << " eTot=" << echtot_Acc
                      << " eInT=" << echint_Acc << endmsg;
      msg(MSG::DEBUG) << " Rejected: nChLo/Hi=" << nChLoCut << "/" << nChHiCut
                      << " eTot=" << echtot_Cut
                      << " eInT=" << echint_Cut << endmsg;
    }
  }


  // step3: register the Digit container in the TES
  SG::WriteHandle<TileDigitsContainer> digitsCnt(m_digitsContainerKey, ctx);
  ATH_CHECK( digitsCnt.record(std::move(digitsContainer)) );

  if(m_doDigiTruth && digitsContainer_DigiHSTruth){
    SG::WriteHandle<TileDigitsContainer> digits_DigiHSTruth(m_digitsContainer_DigiHSTruthKey, ctx);
    ATH_CHECK( digits_DigiHSTruth.record(std::move(digitsContainer_DigiHSTruth)) );
  }

  if (filteredContainer) {
    SG::WriteHandle<TileDigitsContainer> filteredDigitsContainer(m_filteredDigitsContainerKey, ctx);
    ATH_CHECK( filteredDigitsContainer.record(std::move(filteredContainer)) );
  }

  return StatusCode::SUCCESS;
}

StatusCode TileDigitsMaker::finalize() {
  ATH_MSG_INFO( "TileDigitsMaker finalized successfully");

  return StatusCode::SUCCESS;
}

StatusCode TileDigitsMaker::fillDigitCollection(const TileHitCollection* hitCollection,
                                                std::vector<std::vector<double>>& drawerBufferLo,
                                                std::vector<std::vector<double>>& drawerBufferHi,
                                                std::vector<int>& igain, std::vector<int>& over_gain, std::vector<double>& ech_int,
                                                std::vector<bool> &signal_in_channel, const TileEMScale* emScale,
                                                const TileSamplingFraction* samplingFraction, const TilePulse& pulse) const{
  
  constexpr int nchMax = 48; // number of channels per drawer
  std::array<int, nchMax> ntot_ch; ntot_ch.fill(0);
  std::array<double, nchMax> ech_tot; ech_tot.fill(0.0);
  //double ech_int[nchMax];

  IdContext drawer_context = m_tileHWID->drawer_context();

  /* Set up buffers for handling information in a single collection. */
  HWIdentifier drawer_id = m_tileHWID->drawer_id(hitCollection->identify());
  int ros = m_tileHWID->ros(drawer_id);
  int drawer = m_tileHWID->drawer(drawer_id);
  int drawerIdx = TileCalibUtils::getDrawerIdx(ros, drawer);


  // iterate over all hits in a collection
  for (const TileHit* tileHit : *hitCollection) {

    /* Get hit Identifier (= pmt_ID) and needed parameters for this channel */
    Identifier pmt_id = tileHit->pmt_ID();
    double mbts_extra_factor = (m_tileTBID->is_tiletb(pmt_id)) ? -1.0 : 1.0;
    HWIdentifier channel_id = tileHit->pmt_HWID();
    int ich = m_tileHWID->channel(channel_id);
    signal_in_channel[ich] = true;

    if (over_gain[ich] > 9) {
      if (msgLvl(MSG::DEBUG)) {
        int n_hits = tileHit->size();
        double e_hit(0.);
        for (int ihit = 0; ihit < n_hits; ++ihit) {
          e_hit += tileHit->energy(ihit);
        }
        e_hit *= std::round(samplingFraction->getSamplingFraction(drawerIdx, ich) * 1000) / 1000;
        ech_tot[ich] += e_hit;
        ntot_ch[ich] += n_hits;
        ATH_MSG_VERBOSE("BAD Overlay digits - skip hit in channel " << m_tileHWID->to_string(channel_id,-1));
      }
      continue;
    } else {
      ATH_MSG_VERBOSE("new hit in channel " << m_tileHWID->to_string(channel_id,-1));
    }

    /* Set gain=high and get digitSamples and calibration for this channel. */
    if (igain[ich] < 0)
      igain[ich] = TileID::HIGHGAIN;
    // conversion from scintillator energy to total cell energy (sampling fraction)  
    double hit_calib = samplingFraction->getSamplingFraction(drawerIdx, ich);
    hit_calib = std::round(hit_calib * 1000) / 1000;

    // conversion to ADC counts for high gain
    double efactorHi = hit_calib / emScale->calibrateChannel(drawerIdx, ich, TileID::HIGHGAIN, 1.
                                    , TileRawChannelUnit::ADCcounts, TileRawChannelUnit::MegaElectronVolts);
    // conversion to ADC counts for low gain
    double efactorLo = hit_calib / emScale->calibrateChannel(drawerIdx, ich, TileID::LOWGAIN, 1.
                                  , TileRawChannelUnit::ADCcounts, TileRawChannelUnit::MegaElectronVolts);

    std::vector<double>& digitSamplesHi = drawerBufferHi[ich];
    std::vector<double>& digitSamplesLo = drawerBufferLo[ich];
    /* Loop over the subhits for this channel.  For each one,
     convolute with shaping function and add to digitSamples.                  */
    int n_hits = tileHit->size();
    for (int ihit = 0; ihit < n_hits; ++ihit) {
      /* Get hit energy and convert to amplitude of high-gain and low-gain channel */
      double e_hit = tileHit->energy(ihit);
      double amp_ch = e_hit * efactorHi;
      double amp_ch_lo = e_hit * efactorLo;
      double ech_sub = e_hit * hit_calib;
      double t_hit = tileHit->time(ihit);

      ech_tot[ich] += ech_sub;
      if (fabs(t_hit) < 50.0) // ene within +/- 50 ns, used for filtered digits cut
        ech_int[ich] += ech_sub * mbts_extra_factor;
      ntot_ch[ich] += 1;

      // Assume time is in nanoseconds, use fine-grain shaping:
      int ishiftHi = (int) (t_hit / m_timeStepHi + 0.5);
      for (int js = 0; js < m_nSamples; ++js) {
        int k = m_binTime0Hi + (js - m_iTrig) * m_nBinsPerXHi - ishiftHi;
        if (k < 0)
          k = 0;
        else if (k > m_nShapeHi)
          k = m_nShapeHi;

        if (m_useCoolPulseShapes) {
          float phase = (k - m_binTime0Hi) * m_timeStepHi;
          float y, dy;
          pulse.getPulseShapeYDY(drawerIdx, ich, 1, phase, y, dy);
          double ampl = (double) y;
          digitSamplesHi[js] += amp_ch * ampl;
          ATH_MSG_VERBOSE( "Sample no.=" << js
                           << " Pulse index=" << k
                           << " Shape wt. =" << ampl
                           << " HIGAIN from COOL");

        } else {
          digitSamplesHi[js] += amp_ch * m_digitShapeHi[k];
          ATH_MSG_VERBOSE( "Sample no.=" << js
                           << " Pulse index=" << k
                           << " Shape wt. =" << m_digitShapeHi[k]
                           << " HIGAIN from TileInfo");
        }

      }
      int ishiftLo = (int) (t_hit / m_timeStepLo + 0.5);
      for (int js = 0; js < m_nSamples; ++js) {
        int k = m_binTime0Lo + (js - m_iTrig) * m_nBinsPerXLo - ishiftLo;
        if (k < 0)
          k = 0;
        else if (k > m_nShapeLo)
          k = m_nShapeLo;

        if (m_useCoolPulseShapes) {
          float phase = (k - m_binTime0Lo) * m_timeStepLo;
          float y, dy;
          pulse.getPulseShapeYDY(drawerIdx, ich, 0, phase, y, dy);
          double ampl = (double) y;
          digitSamplesLo[js] += amp_ch_lo * ampl;
          ATH_MSG_VERBOSE( "Sample no.=" << js
                           << " Pulse index=" << k
                           << " Shape wt. =" << ampl
                           << " LOGAIN from COOL");
        } else {
          digitSamplesLo[js] += amp_ch_lo * m_digitShapeLo[k];
          ATH_MSG_VERBOSE( "Sample no.=" << js
                           << " Pulse index=" << k
                           << " Shape wt. =" << m_digitShapeLo[k]
                           << " LOGAIN from TileInfo");
        }

      }

      if (msgLvl(MSG::VERBOSE)) {
        msg(MSG::VERBOSE) << "subHit:  ch=" << ich
                          << " e_hit=" << e_hit
                          << " t_hit=" << t_hit
                          << " SamplesHi[" << m_iTrig << "]=" << digitSamplesHi[m_iTrig]
                          << " SamplesLo[" << m_iTrig << "]=" << digitSamplesLo[m_iTrig] << endmsg;
      }
    } /* end loop over sub-hits */
  } /* end loop over hits for this collection. */


  return StatusCode::SUCCESS;

}

StatusCode TileDigitsMaker::overlayBackgroundDigits(const TileDigitsCollection *bkgDigitCollection,
                                                    const TileHitCollection* hitCollection,
                                                    std::vector<std::vector<double>>& drawerBufferLo,
                                                    std::vector<std::vector<double>>& drawerBufferHi,
                                                    std::vector<int>& igain, int ros, int drawer, int drawerIdx,
                                                    std::vector<int>& over_gain, const TileEMScale* emScale,
                                                    const TileSampleNoise& sampleNoise, const TileDQstatus* dqStatus,
                                                    const TileBadChannels* badChannels) const {

  if (hitCollection->identify() != bkgDigitCollection->identify()) {
    ATH_MSG_ERROR ( "Frag IDs for hit collection and digits overlay collection do not match "
                    << MSG::hex << hitCollection->identify() << " != " << bkgDigitCollection->identify()
                    << MSG::dec );
    return StatusCode::FAILURE;
  }


  // iterate over all digits in a collection
  for (const auto* bkgDigit : *bkgDigitCollection) {

    /* Get digit HWIdentifier (= adc_id) */
    HWIdentifier adcId = bkgDigit->adc_HWID();
    int channel = m_tileHWID->channel(adcId);
    int gain = m_tileHWID->adc(adcId);

    igain[channel] = gain;

    // get channel status
    bool good_dq = dqStatus->isAdcDQgood(ros, drawer, channel, gain);
    bool good_ch = (!badChannels->getAdcStatus(adcId).isBad());

    // get digits
    std::vector<float> digits = bkgDigit->samples();
    // get number of time samples & compare with nSamp
    int nSamp2 = digits.size();
    int goodDigits = nSamp2;
    float dig(m_f_ADCmaskValue),digmin(65536.),digmax(-65536.);
    if (goodDigits > 0) {
      auto minmax = std::minmax_element(digits.begin(), digits.end());
      digmin = *minmax.first;
      digmax = *minmax.second;
    }

    if (good_dq) {
      if (digmax > m_ADCmaxPlusEps) { // ignore everything in case of invalid digits
        dig = m_f_ADCmaskValue;
        goodDigits = 0;
      } else { // skip zeros or overflows
        float ADCmaxMinusEps = m_ADCmaxMinusEps;
        int badDigits = std::count_if(digits.begin(), digits.end(), [ADCmaxMinusEps](float dig){
                                                                      return (dig < 0.01) || (dig > ADCmaxMinusEps);});
        goodDigits -= badDigits;
        dig = digits.back();
      }
    } else if (goodDigits>0) {
      goodDigits = 0;
      dig = digits.back();
    }

    if (goodDigits>0) {
      over_gain[channel] = gain;
      if (nSamp2 != m_nSamples) {
        float lastDigit = digits.back();
        digits.resize(m_nSamples);
        // repeat last value in vector (nSamp-nSamp2) times
        std::fill(digits.begin() + nSamp2, digits.end(), lastDigit);
      }

      std::vector<double>& buffer = (gain == TileID::HIGHGAIN) ? drawerBufferHi[channel] : drawerBufferLo[channel];
      std::vector<double>& bufferLG = drawerBufferLo[channel];

      bool isFilledLG = false;
      if (gain == TileID::HIGHGAIN) {
        if (digmax - digmin > 5. && good_ch ) {// 5 ADC counts cut - to ignore pure noise in HG (less than 0.1 count effect in LG)
          float ratio = emScale->applyOnlineChargeCalibration(drawerIdx, channel, TileID::HIGHGAIN, 1.)
            / emScale->applyOnlineChargeCalibration(drawerIdx, channel, TileID::LOWGAIN, 1.); // ratio between low and high gain

          dig=std::min(digits[0],std::max(digmin, sampleNoise.getPed(drawerIdx, channel, TileID::HIGHGAIN)));

          std::transform(digits.begin(), digits.end(), bufferLG.begin(), [dig,ratio](float digit){return (digit - dig) * ratio;});
          isFilledLG = true;
        }
      }

      std::copy(digits.begin(), digits.end(), buffer.begin());

      if (msgLvl(MSG::VERBOSE)) {
        msg(MSG::VERBOSE) << "RNDM BG ADC " << m_tileHWID->to_string(adcId)
                          << " samples=";
        for (int js = 0; js < m_nSamples; ++js)
          msg(MSG::VERBOSE) << " " << buffer[js];
        if (!good_ch)
          msg(MSG::VERBOSE) << " BCH";
        if (!good_dq) {
          msg(MSG::VERBOSE) << " BDQ";
        } else if (isFilledLG) {
          msg(MSG::VERBOSE) << "  LG=";
          for (int js = 0; js < m_nSamples; ++js)
            msg(MSG::VERBOSE) << " " << int(bufferLG[js]*100)/100.;
        }
        msg(MSG::VERBOSE) << endmsg;
      }

    } else if (nSamp2 > 0) {
      over_gain[channel] = 10+gain; // flag problematic channel

      std::vector<double>& buffer = (gain == TileID::HIGHGAIN) ? drawerBufferHi[channel] : drawerBufferLo[channel];

      if (digmin != digmax || (dig!=0. && dig!=m_f_ADCmax)) {
        dig = m_f_ADCmaskValue; // keep only 0 or m_f_ADCmax as it is
      }
      std::fill(buffer.begin(), buffer.end(), dig);

      if (msgLvl(MSG::VERBOSE)) {
        msg(MSG::VERBOSE) << "BAD BG  ADC " << m_tileHWID->to_string(adcId)
                          << " samples=";
        for (int js = 0; js < nSamp2; ++js)
          msg(MSG::VERBOSE) << " " << digits[js];
        msg(MSG::VERBOSE) << ((good_ch)?"":" BCH") << ((good_dq)?"":" BDQ") << endmsg;
      }

    } else {
      ATH_MSG_VERBOSE( "NO BG   ADC " << m_tileHWID->to_string(adcId)
                       << " samples= 0 0 0 0 0 0 0"
                       << ((good_ch)?"":" BCH") << ((good_dq)?"":" BDQ") );
    }
  }
  return StatusCode::SUCCESS;
}
