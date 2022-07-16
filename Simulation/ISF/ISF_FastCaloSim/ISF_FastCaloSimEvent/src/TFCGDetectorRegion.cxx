 /*
  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
// TFCGDetectorRegion.cxx, (c) ATLAS Detector software             //
///////////////////////////////////////////////////////////////////

// class header include
#include "ISF_FastCaloSimEvent/TFCGDetectorRegion.h"

TFCGDetectorRegion::DetectorRegion(int id, int pid, int etaMin, int etaMax){
  m_id = id;
  m_pid = pid;
  m_etaMin = etaMin;
  m_etaMax = etaMax;  
}

std::vector<int> TFCGDetectorRegion::GetRelevantLayers() { 
  return m_relevantlayers;
}

TFCGDetectorRegion::Binning TFCGDetectorRegion::GetBinning() { 
  return m_binning;
}

bool TFCGDetectorRegion::IsSymmetrisedAlpha() { 
  return m_symmetrisedAlpha;
}

void TFCGDetectorRegion::SetRelevantLayers(std::vector<int> relevantlayers){
  m_relevantlayers = relevantlayers;
}

void TFCGDetectorRegion::SetBinning(TFCGDetectorRegion::Binning binning){
  m_binning = binning;
}

void TFCGDetectorRegion::SetSymmetrisedAlpha(bool symmetrisedAlpha){
  m_symmetrisedAlpha = symmetrisedAlpha;
}

    