/*
   Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
 */
/**
 * @file PixelAthHitMonAlg.cxx
 * @brief Helper functions to fill different types of pixel histograms
 * @author Iskander Ibragimov
 **/
#include "PixelAthMonitoringBase.h"


//////////////////////////////////////////////
///
/// helper class to accumulate points to fill a 2D per-module plot with
///
void PixelAthMonitoringBase::VecAccumulator2DMap::add(const int layer, const Identifier& id,
                                                      float value) {
  m_pm[layer].push_back(m_host.m_pixelid->phi_module(id));
  int ld = m_host.m_pixelid->layer_disk(id);
  int em = ld;
  m_val[layer].push_back(value);

  bool copy = false;
  if (m_host.m_pixelid->barrel_ec(id) == 0) {
    em = m_host.m_pixelid->eta_module(id);
    if (ld == 0) {
      int feid = 0;
      int emf = 0;
      if (em < -6) {
        emf = em - 6;
      } else if (em > -7 && em < 6) {
        if (m_host.m_pixelid->eta_index(id) >= 80) feid = 1;
        emf = 2 * em + feid;
        copy = true;
      } else {
        emf = em + 6;
      }
      em = emf;
    }
  }
  m_em[layer].push_back(em);

  if (m_copy2DFEval && copy) {
    ++em;
    m_pm[layer].push_back(m_host.m_pixelid->phi_module(id));
    m_em[layer].push_back(em);
    m_val[layer].push_back(value);
  }
}

//////////////////////////////////////////////
///
/// helper class to accumulate points to fill a 2D per-FE plot with
///
void PixelAthMonitoringBase::VecAccumulator2DMap::add(const int layer, const Identifier& id,
                                                      int iFE, float value) {
  // value
  m_val[layer].push_back(value);

  // for old pixel see https://twiki.cern.ch/twiki/pub/Atlas/PixelCOOLoffline/pixel_modules_sketch.png
  //
  // phi (Y) coordinate
  if (layer == PixLayers::kIBL) m_pm[layer].push_back(m_host.m_pixelid->phi_module(id));
  else {
    if ((layer == PixLayers::kECA || layer == PixLayers::kECC) && (m_host.m_pixelid->phi_module(id) % 2 == 0)) {
      // even disk modules
      m_pm[layer].push_back(m_host.m_pixelid->phi_module(id) * 2 + iFE / 8);
    } else {
      m_pm[layer].push_back(m_host.m_pixelid->phi_module(id) * 2 + 1 - iFE / 8);
    }
  }
  // eta (X) coordinate
  int ld = m_host.m_pixelid->layer_disk(id);
  int em(0);
  // endcaps
  em = ld * 8;
  if (iFE < 8) em += (7 - iFE % 8);
  else em += iFE % 8;

  // barrel
  //
  if (m_host.m_pixelid->barrel_ec(id) == 0) {
    if (ld == 0) {  //ibl
      em = m_host.m_pixelid->eta_module(id);
      int emf;
      if (em < -6) {
        emf = em - 6;
      } else if (em > -7 && em < 6) {
        emf = 2 * em + iFE;
      } else {
        emf = em + 6;
      }
      em = emf;
    } else {
      em = m_host.m_pixelid->eta_module(id) * 8 - 4;
      if (iFE < 8) em += (7 - iFE % 8);
      else em += iFE % 8;
    }
  } // end barrel
  m_em[layer].push_back(em);
}

StatusCode PixelAthMonitoringBase::initialize() {
   ATH_CHECK( AthMonitorAlgorithm::initialize() );
   ATH_CHECK( m_pixelCondSummaryTool.retrieve( DisableTool{ !m_pixelDetElStatus.empty() && !m_pixelDetElStatusActiveOnly.empty() && !VALIDATE_STATUS_ARRAY_ACTIVATED} ) );
   if (m_pixelDetElStatus.empty() || m_pixelDetElStatusActiveOnly.empty() || VALIDATE_STATUS_ARRAY_ACTIVATED) {
      ATH_CHECK( m_pixelReadout.retrieve() );
   }
   ATH_CHECK(detStore()->retrieve(m_pixelid, "PixelID"));
   ATH_CHECK( m_pixelDetElStatus.initialize( !m_pixelDetElStatus.empty()) );
   ATH_CHECK( m_pixelDetElStatusActiveOnly.initialize( !m_pixelDetElStatusActiveOnly.empty()) );
   return StatusCode::SUCCESS;
}

//////////////////////////////////////////////
///
/// take VecAccumulator2DMap and fill the corresponding group
///
void PixelAthMonitoringBase::fill2DProfLayerAccum(const VecAccumulator2DMap& accumulator) const {
  // iterate over all actually filled layers
  for (const auto& itr : accumulator.m_pm) {
    // Define the monitored variables
    int layer = itr.first;
    auto pm = Monitored::Collection(accumulator.m_prof2Dname + "_pm", accumulator.m_pm.at(layer));
    auto val = Monitored::Collection(accumulator.m_prof2Dname + "_val", accumulator.m_val.at(layer));
    auto em = Monitored::Collection(accumulator.m_prof2Dname + "_em", accumulator.m_em.at(layer));
    fill(pixBaseLayersLabel[layer], pm, em, val);
  }
}

///
/// filling 1DProf per-lumi per-layer histograms ["ECA","ECC","BLayer","Layer1","Layer2","IBL","IBL2D","IBL3D"]
///
void PixelAthMonitoringBase::fill1DProfLumiLayers(const std::string& prof1Dname, int lumiblock, float* values,
                                                  int nlayers) const {
  ATH_MSG_VERBOSE("in fill1DProfLumiLayers()");

  // Define the monitored variables
  auto lb = Monitored::Scalar<int>(prof1Dname + "_lb", lumiblock);
  auto val = Monitored::Scalar<float>(prof1Dname + "_val", 1.0);

  int i_start = 0;
  int i_end = PixLayers::COUNT;
  if (nlayers == PixLayers::NFEI3LAYERS) i_end = nlayers;
  if (nlayers == PixLayers::COUNT - PixLayers::NFEI3LAYERS) i_start = PixLayers::NFEI3LAYERS;
  for (int i = i_start; i < i_end; i++) {
    val = values[i];
    fill(pixLayersLabel[i], lb, val);
  }
}

//////////////////////////////////////////////

///
/// filling 2DProf per-lumi per-layer histograms ["ECA","ECC","BLayer","Layer1","Layer2","IBL2D","IBL3D"]
///
void PixelAthMonitoringBase::fill2DProfLumiLayers(const std::string& prof2Dname, int lumiblock,
                                                  float(*values)[PixLayers::COUNT], const int* nCategories) const {
  ATH_MSG_VERBOSE("in fill2DProfLumiLayers()");

  // Define the monitored variables
  auto lb = Monitored::Scalar<int>(prof2Dname + "_lb", lumiblock);
  auto val = Monitored::Scalar<float>(prof2Dname + "_val", 1.0);
  auto cat = Monitored::Scalar<int>(prof2Dname + "_cat");

  for (int i = 0; i < PixLayers::COUNT; i++) {
    for (cat = 0; cat < nCategories[i]; cat++) {
      val = values[cat][i];
      fill(pixLayersLabel[i], lb, cat, val);
    }
  }
}

//////////////////////////////////////////////


///
/// filling 1DProfile per-pp0(ROD) histograms for ["ECA","ECC","BLayer","Layer1","Layer2","IBLA","IBLC"]
///
void PixelAthMonitoringBase::fillFromArrays(const std::string& namePP0, AccumulatorArrays& pixarrays,
                                            const std::string& name2DMap) const {
  ATH_MSG_VERBOSE("in fillFromArrays()");

  const float weightPix = 1.0 / 46080.0;
  const float weightIBL = 1.0 / 26880.0;

  bool fillPP0only(name2DMap == "");
  std::string pospp0varx = namePP0 + "_pospp0x";
  std::string valvarp = namePP0 + "_val";
  std::string posvarx = name2DMap + "_em";
  std::string posvary = name2DMap + "_pm";
  std::string valvarm = name2DMap + "_val";

  for (unsigned int a = 0; a < PixMon::kNumModulesDisk; ++a) {
    auto posy = Monitored::Scalar<int>(posvary, a);
    for (unsigned int b = 0; b < PixMon::kNumLayersDisk; ++b) {
      // to find out (and fill together into one PP0-histogram bin)
      // array content of the modules belonging to the same sector (or PP0)
      // the translation (a-1)/6 is used
      // to show PP0 values from other disks of the same endcap
      // in the same plot
      // the shift (b-1)*8 applies per disk counter b
      // (there are in total 8 sectors/disk)
      auto pospp0x = Monitored::Scalar<int>(pospp0varx, a / 6 + b * 8);
      auto posx = Monitored::Scalar<int>(posvarx, b);
      auto valp = Monitored::Scalar<float>(valvarp, pixarrays.DA[a][b]);
      auto valm = Monitored::Scalar<float>(valvarm, pixarrays.DA[a][b] * weightPix);
      if (pixarrays.DA[a][b] > -1) {
        fill("ECA", pospp0x, valp);
        if (!fillPP0only) fill("ECA", posx, posy, valm);
      }
      valp = pixarrays.DC[a][b];
      valm = pixarrays.DC[a][b] * weightPix;
      if (pixarrays.DC[a][b] > -1) {
        fill("ECC", pospp0x, valp);
        if (!fillPP0only) fill("ECC", posx, posy, valm);
      }
    }
  }

  for (unsigned int b = 0; b < PixMon::kNumModulesBarrel; ++b) {
    // translating array index into old Pixel module eta on a stave
    // i.e. 0..12 into -6..6 so that standard per-layer histograms
    // declared by define2DProfHist method can be filled
    auto posx = Monitored::Scalar<int>(posvarx, b - 6);

    for (unsigned int a = 0; a < PixMon::kNumStavesL0; ++a) {
      auto posy = Monitored::Scalar<int>(posvary, a);
      auto pospp0x = Monitored::Scalar<int>(pospp0varx, a);
      auto valp = Monitored::Scalar<float>(valvarp, pixarrays.B0[a][b]);
      auto valm = Monitored::Scalar<float>(valvarm, pixarrays.B0[a][b] * weightPix);
      if (pixarrays.B0[a][b] > -1) {
        fill("BLayer", pospp0x, valp);
        if (!fillPP0only) fill("BLayer", posx, posy, valm);
      }
    }
    for (unsigned int a = 0; a < PixMon::kNumStavesL1; ++a) {
      auto posy = Monitored::Scalar<int>(posvary, a);
      auto pospp0x = Monitored::Scalar<int>(pospp0varx, a);
      auto valp = Monitored::Scalar<float>(valvarp, pixarrays.B1[a][b]);
      auto valm = Monitored::Scalar<float>(valvarm, pixarrays.B1[a][b] * weightPix);
      if (pixarrays.B1[a][b] > -1) {
        fill("Layer1", pospp0x, valp);
        if (!fillPP0only) fill("Layer1", posx, posy, valm);
      }
    }
    for (unsigned int a = 0; a < PixMon::kNumStavesL2; ++a) {
      auto posy = Monitored::Scalar<int>(posvary, a);
      auto pospp0x = Monitored::Scalar<int>(pospp0varx, a);
      auto valp = Monitored::Scalar<float>(valvarp, pixarrays.B2[a][b]);
      auto valm = Monitored::Scalar<float>(valvarm, pixarrays.B2[a][b] * weightPix);
      if (pixarrays.B2[a][b] > -1) {
        fill("Layer2", pospp0x, valp);
        if (!fillPP0only) fill("Layer2", posx, posy, valm);
      }
    }
  }
  unsigned int nbina = PixMon::kNumStavesIBL;
  unsigned int nbinb = PixMon::kNumFEsIBL;
  for (unsigned int a = 0; a < nbina; ++a) {
    auto posy = Monitored::Scalar<int>(posvary, a);
    auto pospp0x = Monitored::Scalar<int>(pospp0varx, a);
    for (unsigned int b = 0; b < nbinb; ++b) {
      // translating array index into IBL frontend eta on a stave
      // i.e. 0..31 into -16..15 so that standard per-layer histograms
      // declared by define2DProfHist method can be filled
      auto posx = Monitored::Scalar<int>(posvarx, b - 16);
      auto valp = Monitored::Scalar<float>(valvarp, pixarrays.IBL[a][b]);
      auto valm = Monitored::Scalar<float>(valvarm, pixarrays.IBL[a][b] * weightIBL);
      if (pixarrays.IBL[a][b] > -1) {
        if (b > (0.5 * nbinb - 1)) {
          fill("IBLA", pospp0x, valp);
        } else {
          fill("IBLC", pospp0x, valp);
        }
        if (!fillPP0only) fill("IBL", posx, posy, valm);
      }
    }
  }
}

//////////////////////////////////////////////


///
/// helper function to get layers ID
///
int PixelAthMonitoringBase::getPixLayersID(int ec, int ld) const {
  int layer = 99;

  if (ec == 2) {
    layer = PixLayers::kECA;
  } else if (ec == -2) {
    layer = PixLayers::kECC;
  } else if (ec == 0) {
    if (ld == 0) layer = PixLayers::kIBL;
    if (ld == 1) layer = PixLayers::kBLayer;
    if (ld == 2) layer = PixLayers::kLayer1;
    if (ld == 3) layer = PixLayers::kLayer2;
  }
  return layer;
}

///
/// helper function to check if module is IBL planar based on pixel hash ID
///
bool PixelAthMonitoringBase::isIBL2D(int hashID) const {
  bool result(false);
  if ( hashID>=156 && hashID<=435 ) // IBL
    { 
      int module = (hashID-156) % 20;
      if (module>3 && module<16)
	{ 
	  result = true;
	}
    }
  return result;
}

//////////////////////////////////////////////

///
/// helper function to check if module is IBL 3D based on pixel hash ID
///
bool PixelAthMonitoringBase::isIBL3D(int hashID) const {
  bool result(false);
  if ( hashID>=156 && hashID<=435 ) // IBL
    { 
      int module = (hashID-156) % 20;
      if (module<4 || module>15)
	{ 
	  result = true;
	}
    }
  return result;
}

//////////////////////////////////////////////

///
/// helper function to get number of FEs per module
///
int PixelAthMonitoringBase::getNumberOfFEs(int pixlayer, int etaMod) const {
  int nFE(16);

  if (pixlayer == PixLayers::kIBL) {
    nFE = 1; // IBL 3D
    if (etaMod > -7 && etaMod < 6) nFE = 2; // IBL Planar
  }
  return nFE;
}

//////////////////////////////////////////////


///
/// helper function to get eta phi coordinates of per-layer arrays
///
void PixelAthMonitoringBase::getPhiEtaMod(Identifier& id, int& phiMod, int& etaMod,
                                          bool& copyFE) const {
  phiMod = m_pixelid->phi_module(id);

  int layerDisk = m_pixelid->layer_disk(id);
  etaMod = layerDisk;
  copyFE = false;
  if (m_pixelid->barrel_ec(id) == 0) {
    etaMod = m_pixelid->eta_module(id);
    if (layerDisk == 0) {
      if (etaMod < -6) {
        etaMod = etaMod - 6;
      } else if (etaMod > -7 && etaMod < 6) {
        int feid = 0;
        if (m_pixelid->eta_index(id) >= 80) feid = 1;
        etaMod = 2 * etaMod + feid;
        copyFE = true;
      } else {
        etaMod = etaMod + 6;
      }
      etaMod = etaMod + 16;
    } else etaMod = etaMod + 6;
  }
}

//////////////////////////////////////////////
///
/// checks if hit is on track
///

bool PixelAthMonitoringBase::isHitOnTrack(Identifier id, std::vector<Identifier> const& RDOIDs) const {
  return binary_search(RDOIDs.begin(), RDOIDs.end(), id);

  ;
}

//////////////////////////////////////////////


///
/// checks if cluster is on track
///

bool PixelAthMonitoringBase::isClusterOnTrack(Identifier id, std::vector<std::pair<Identifier,
                                                                                   double> > const& ClusterIDs) const {
  bool onTrack = false;

  std::pair<Identifier, double> searchVal = std::make_pair(id, -1.0);
  onTrack = std::binary_search(ClusterIDs.begin(), ClusterIDs.end(), searchVal,
                               [](std::pair<Identifier, double> l, std::pair<Identifier, double> r) -> bool {
    return l.first < r.first;
  });
  return onTrack;
}

//////////////////////////////////////////////

///
/// checks if cluster is on track and returns its cosalpha
///
bool PixelAthMonitoringBase::isClusterOnTrack(Identifier id, std::vector<std::pair<Identifier,
                                                                                   double> > const& ClusterIDs,
                                              double& cosalpha) const {
  bool onTrack(false);

  std::pair<Identifier, double> searchVal = std::make_pair(id, -1.0);
  auto it = std::lower_bound(ClusterIDs.begin(), ClusterIDs.end(), searchVal,
                             [](std::pair<Identifier, double> l, std::pair<Identifier, double> r) -> bool {
    return l.first < r.first;
  });

  if (it != ClusterIDs.end() && !(id < (*it).first)) {
    onTrack = true;
    cosalpha = (*it).second;
  }
  return onTrack;
}

//////////////////////////////////////////////
