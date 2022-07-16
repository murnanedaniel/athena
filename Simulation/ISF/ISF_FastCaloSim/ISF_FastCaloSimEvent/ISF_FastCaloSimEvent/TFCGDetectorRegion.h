/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/
 
//////////////////////////////////////////////////////////////////
// TFCGDetectorRegion.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef ISF_TFCGDETECTORREGION_H
#define ISF_TFCGDETECTORREGION_H 1

#include <map>
#include <vector>

#include "TH2D.h"

class TFCGDetectorRegion 
{
  public:
  typedef std::map<int, TH2D> Binning;

  TFCGDetectorRegion(int id, int pid, int etaMin, int etaMax);
  
  std::vector<int> GetRelevantLayers();
  Binning GetBinning();
  
  void SetRelevantLayers(std::vector<int> relevantlayers);
  void SetBinning(Binning binning);
  void SetSymmetrisedAlpha(bool symmetrisedAlpha);
  
  bool IsSymmetrisedAlpha();
  int GetID();
  
  private:
  int m_id;
  int m_pid;
  int m_etaMin;
  int m_etaMax;
  bool m_symmetrisedAlpha;

  Binning m_binning;
  std::vector<int> m_relevantlayers;
};

#endif //> !ISF_TFCGDETECTORREGION_H