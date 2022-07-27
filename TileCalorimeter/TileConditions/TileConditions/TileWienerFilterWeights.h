/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TILECONDITIONS_TILEWIENERFILTERWEIGHTS_H
#define TILECONDITIONS_TILEWIENERFILTERWEIGHTS_H

#include "GaudiKernel/MsgStream.h"

typedef struct
{
  int luminosity;
  double generalWeights[4][48][8]; // one set of weights for each TileCall cell
  double optimalWeights[7][8]; // weights for E3 and E4 only
} TileWienerFilterWeightsStruct;

class  TileWienerFilterWeights {

  friend class TileInfoLoader;

 public:

  TileWienerFilterWeights();
  ~TileWienerFilterWeights();

  void loadWeights(MsgStream &log);
  const TileWienerFilterWeightsStruct * TileWFstruct() const { return m_weights; }
  bool loaded()    { return m_loaded; }

 private:

  //variables
  int  m_Luminosity;
  int  m_NSamples_Phys = 0;
  bool m_loaded;
  TileWienerFilterWeightsStruct * m_weights;

};

#endif // TILECONDITIONS_TILEWIENERFILTERWEIGHTS_H
