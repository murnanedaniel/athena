/*
   Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
 */
/**
 * @file PixelAthHitMonAlg.cxx
 * @brief Reads Pixel RDO information and fills it into histograms 
 * @author Iskander Ibragimov
 **/
#include "PixelAthHitMonAlg.h"
#include <stdexcept>

PixelAthHitMonAlg::PixelAthHitMonAlg(const std::string& name, ISvcLocator* pSvcLocator) :
  AthMonitorAlgorithm(name, pSvcLocator)
{
  //jo flags
  declareProperty("doOnline", m_doOnline = false);
  declareProperty("doLumiBlock", m_doLumiBlock = false);
  declareProperty("doLowOccupancy", m_doLowOccupancy = false);
  declareProperty("doHighOccupancy", m_doHighOccupancy = false);
  declareProperty("doHeavyIonMon", m_doHeavyIonMon = false);
  declareProperty("doFEPlots", m_doFEPlots = false);
}

StatusCode PixelAthHitMonAlg::initialize() {
  ATH_CHECK( PixelAthMonitoringBase::initialize());
  ATH_CHECK(m_pixelRDOName.initialize());

  return StatusCode::SUCCESS;
}

StatusCode PixelAthHitMonAlg::fillHistograms(const EventContext& ctx) const {
  using namespace Monitored;

  int lb = GetEventInfo(ctx)->lumiBlock();
  auto lbval = Monitored::Scalar<int>("pixhitsmontool_lb", lb);
  auto hitGroup = getGroup("Hit");

  unsigned int bcid = GetEventInfo(ctx)->bcid();
  auto bcidval = Monitored::Scalar<unsigned int>("pixhitsmontool_bcid", bcid);
  auto rdocontainer = SG::makeHandle(m_pixelRDOName, ctx);
  if (!(rdocontainer.isValid())) {
    ATH_MSG_ERROR("Pixel Monitoring: Pixel RDO container " << m_pixelRDOName << " could not be found.");
    auto dataread_err = Monitored::Scalar<int>("hitdataread_err", DataReadErrors::ContainerInvalid);
    fill(hitGroup, dataread_err);
    return StatusCode::RECOVERABLE;
  } else {
    ATH_MSG_DEBUG("Pixel Monitoring: Pixel RDO container " << rdocontainer.name() << " is found.");
  }

  int nhits = 0;
  float nhits_layer[PixLayers::COUNT] = {
    0
  };

  int phiMod(-99);
  int etaMod(-99);
  bool copyFEval(false);

  AccumulatorArrays hitsPerEventArray = {{{0}}, {{0}}, {{0}}, {{0}}, {{0}}, {{0}}};
  SG::ReadHandle<InDet::SiDetectorElementStatus> pixel_status = getPixelDetElStatus(m_pixelDetElStatus, ctx);
  SG::ReadHandle<InDet::SiDetectorElementStatus> pixel_active = getPixelDetElStatus(m_pixelDetElStatusActiveOnly, ctx);
  // for inactive or bad modules init corresponding arrays entries to -1
  for (auto idIt = m_pixelid->wafer_begin(); idIt != m_pixelid->wafer_end(); ++idIt) {
    Identifier waferID = *idIt;
    IdentifierHash id_hash = m_pixelid->wafer_hash(waferID);
    int pixlayer = getPixLayersID(m_pixelid->barrel_ec(waferID), m_pixelid->layer_disk(waferID));
    if (pixlayer == 99) continue;
    if (isActive( !m_pixelDetElStatusActiveOnly.empty() ? pixel_active.cptr() : nullptr, id_hash) == false) {
      getPhiEtaMod(waferID, phiMod, etaMod, copyFEval);
      switch (pixlayer) {
      case PixLayers::kECA:
        hitsPerEventArray.DA[phiMod][etaMod] = -1;
        break;

      case PixLayers::kECC:
        hitsPerEventArray.DC[phiMod][etaMod] = -1;
        break;

      case PixLayers::kBLayer:
        hitsPerEventArray.B0[phiMod][etaMod] = -1;
        break;

      case PixLayers::kLayer1:
        hitsPerEventArray.B1[phiMod][etaMod] = -1;
        break;

      case PixLayers::kLayer2:
        hitsPerEventArray.B2[phiMod][etaMod] = -1;
        break;

      case PixLayers::kIBL:
        hitsPerEventArray.IBL[phiMod][etaMod] = -1;
        if (copyFEval) hitsPerEventArray.IBL[phiMod][++etaMod] = -1;
        break;
      }
    }
  }

  int nGood_layer[PixLayers::COUNT] = {
    0
  };
  int nActive_layer[PixLayers::COUNT] = {
    0
  };
  float avgocc_active_layer[PixLayers::COUNT] = {
    0
  };
  float avgocc_good_layer[PixLayers::COUNT] = {
    0
  };
  float avgocc_ratio_toIBL_layer[PixLayers::COUNT] = {
    0
  };

  PixelID::const_id_iterator idIt = m_pixelid->wafer_begin();
  PixelID::const_id_iterator idItEnd = m_pixelid->wafer_end();
  for (; idIt != idItEnd; ++idIt) {
    Identifier waferID = *idIt;
    IdentifierHash modHash = m_pixelid->wafer_hash(waferID);
    int pixlayer = getPixLayersID(m_pixelid->barrel_ec(waferID), m_pixelid->layer_disk(waferID));
    if (pixlayer == 99) continue;

    if (pixlayer == PixLayers::kIBL)
      {
	// normalization per FE for IBL
	//
	int nFE = getNumberOfFEs(pixlayer, m_pixelid->eta_module(waferID));
	int iblsublayer = (m_pixelid->eta_module(waferID) > -7 && m_pixelid->eta_module(waferID) < 6) ? PixLayers::kIBL2D : PixLayers::kIBL3D;
	for (int iFE=0; iFE<nFE; iFE++) {
	  Identifier pixID = m_pixelReadout->getPixelIdfromHash(modHash, iFE, 1, 1);
	  if (not pixID.is_valid()) continue;
	  auto [is_active,is_good] = isChipGood( !m_pixelDetElStatusActiveOnly.empty() ? pixel_active.cptr() : nullptr,
						 !m_pixelDetElStatus.empty() ? pixel_status.cptr() : nullptr,
						 modHash, iFE);
	  if (not is_active) continue;
	  nActive_layer[iblsublayer]++;
	  if (is_good) nGood_layer[iblsublayer]++;
	}
      }
    else
      {
	// normalization per module for the old pixel layers
	//
	bool is_active = isActive( !m_pixelDetElStatusActiveOnly.empty() ? pixel_active.cptr() :  nullptr, modHash);
	bool is_good   = isGood( !m_pixelDetElStatus.empty() ? pixel_status.cptr() :  nullptr, modHash);
	if (is_active) {
	  nActive_layer[pixlayer]++;
	  if (is_good) nGood_layer[pixlayer]++;
	}
      }
  }

  const int nChannels_mod[PixLayers::COUNT] = {
    46080, 46080, 46080, 46080, 46080, 26880, 26880
  };
  float nGoodChannels_layer[PixLayers::COUNT];
  float nActiveChannels_layer[PixLayers::COUNT];
  for (int i = 0; i < PixLayers::COUNT; i++) {
    nGoodChannels_layer[i] = 1.0 * nChannels_mod[i] * nGood_layer[i];
    nActiveChannels_layer[i] = 1.0 * nChannels_mod[i] * nActive_layer[i];
  }

  //*******************************************************************************
  //************************** Begin of filling Hit Histograms ********************
  //*******************************************************************************

  ATH_MSG_DEBUG("Filling Raw Hit (RDO) Monitoring Histograms");


  VecAccumulator2DMap HitMap(*this, "HitMap");
  VecAccumulator2DMap HitFEMap(*this, "HitFEMap");
  std::vector<int> hitLvl1a;
  std::unordered_map<int, std::vector<int> > hitLvl1aLayer;
  std::unordered_map<int, std::vector<int> > hitToTLayer;

  Identifier rdoID;

  for (auto colNext: *rdocontainer) {
    const InDetRawDataCollection<PixelRDORawData>* HitCollection(colNext);
    if (!HitCollection) {
      ATH_MSG_DEBUG("Pixel Monitoring: Pixel Hit container is empty.");
      auto dataread_err = Monitored::Scalar<int>("hitdataread_err", DataReadErrors::CollectionInvalid);
      fill(hitGroup, dataread_err);
      continue;
    }

    for (auto p_rdo: *HitCollection) {
      rdoID = p_rdo->identify();
      int pixlayer = getPixLayersID(m_pixelid->barrel_ec(rdoID), m_pixelid->layer_disk(rdoID));
      if (pixlayer == 99) continue;

      HitMap.add(pixlayer, rdoID, 1.0);
      if (m_doFEPlots) HitFEMap.add(pixlayer, rdoID, m_pixelReadout->getFE(rdoID, rdoID), 1.0);

      nhits++;
      hitLvl1a.push_back(p_rdo->getLVL1A());
      getPhiEtaMod(rdoID, phiMod, etaMod, copyFEval);
      switch (pixlayer) {
      case PixLayers::kECA:
        hitsPerEventArray.DA[phiMod][etaMod]++;
        break;

      case PixLayers::kECC:
        hitsPerEventArray.DC[phiMod][etaMod]++;
        break;

      case PixLayers::kBLayer:
        hitsPerEventArray.B0[phiMod][etaMod]++;
        break;

      case PixLayers::kLayer1:
        hitsPerEventArray.B1[phiMod][etaMod]++;
        break;

      case PixLayers::kLayer2:
        hitsPerEventArray.B2[phiMod][etaMod]++;
        break;

      case PixLayers::kIBL:
        hitsPerEventArray.IBL[phiMod][etaMod]++;
        break;
      }
      if (pixlayer == PixLayers::kIBL) 
	{
	  pixlayer = (m_pixelid->eta_module(rdoID) > -7 && m_pixelid->eta_module(rdoID) < 6) ? PixLayers::kIBL2D : PixLayers::kIBL3D;
	}
      nhits_layer[pixlayer]++;
      hitLvl1aLayer[pixlayer].push_back(p_rdo->getLVL1A());
      hitToTLayer[pixlayer].push_back(p_rdo->getToT());
    }
  }

  fill2DProfLayerAccum(HitMap);
  if (m_doFEPlots) fill2DProfLayerAccum(HitFEMap);

  auto vals = Monitored::Collection("Hit_LVL1A_pixel", hitLvl1a);
  fill(hitGroup, vals);
  for (const auto& itr : hitLvl1aLayer) {
    int layer = itr.first;
    try {
      //cppcheck-suppress containerOutOfBounds
      auto vals = Monitored::Collection("Hit_LVL1A_layer", hitLvl1aLayer.at(layer));
      fill(pixLayersLabel[layer], vals);
    } catch (std::out_of_range& e) {
      ATH_MSG_ERROR("Out of range access in PixelAthHitMonAlg::fillHistograms");
    }
  }
  for (const auto& itr : hitToTLayer) {
    int layer = itr.first;
    try {
      //cppcheck-suppress containerOutOfBounds
      auto vals = Monitored::Collection("HitToT_val", hitToTLayer.at(layer));
      fill(pixLayersLabel[layer], vals);
    } catch (std::out_of_range& e) {
      ATH_MSG_ERROR("Out of range access in PixelAthHitMonAlg::fillHistograms");
    }
  }


  auto nhitsval = Monitored::Scalar<int>("nhits_per_event", nhits);
  fill(hitGroup, lbval, nhitsval);
  fill1DProfLumiLayers("HitsPerLumi", lb, nhits_layer);

  fillFromArrays("HitOccupancyPP0", hitsPerEventArray, "OccupancyPerPixelEvent");

  for (int i = 0; i < PixLayers::COUNT; i++) {
    if (nGoodChannels_layer[i] > 0) avgocc_good_layer[i] = nhits_layer[i] / nGoodChannels_layer[i];
    auto val = Monitored::Scalar<float>("AvgOccPerBCID_val", avgocc_good_layer[i]);
    fill(pixLayersLabel[i], bcidval, val);
    if (nActiveChannels_layer[i] > 0) avgocc_active_layer[i] = nhits_layer[i] / nActiveChannels_layer[i];
  }
  fill1DProfLumiLayers("AvgOccActivePerLumi", lb, avgocc_active_layer);
  fill1DProfLumiLayers("AvgOccGoodPerLumi", lb, avgocc_good_layer);

  
  float nGoodChannels_ibl = nGoodChannels_layer[PixLayers::kIBL2D] + nGoodChannels_layer[PixLayers::kIBL3D];
  float avgocc_good_ibl(0);
  if (nGoodChannels_ibl>0) avgocc_good_ibl = (nhits_layer[PixLayers::kIBL2D] + nhits_layer[PixLayers::kIBL3D])/nGoodChannels_ibl;
  if (avgocc_good_ibl > 0) {
    for (int i = 0; i < PixLayers::COUNT; i++) {
      avgocc_ratio_toIBL_layer[i] = avgocc_good_layer[i] / avgocc_good_ibl;
    }
    fill1DProfLumiLayers("AvgOccRatioToIBLPerLumi", lb, avgocc_ratio_toIBL_layer);
  }

  if (nhits == 0) {
    auto dataread_err = Monitored::Scalar<int>("hitdataread_err", DataReadErrors::EmptyContainer);
    fill(hitGroup, dataread_err);
  }
  //*******************************************************************************
  //************************** End of filling Hit Histograms **********************
  //*******************************************************************************

  return StatusCode::SUCCESS;
}
