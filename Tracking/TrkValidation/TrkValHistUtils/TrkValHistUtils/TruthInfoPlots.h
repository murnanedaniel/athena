/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TRKVALHISTUTILS_TRUTHINFOPLOTS_H
#define TRKVALHISTUTILS_TRUTHINFOPLOTS_H

#include "PlotBase.h"
#include "xAODTruth/TruthParticle.h"

namespace Trk{

class TruthInfoPlots: public PlotBase {
  public:
    TruthInfoPlots(PlotBase *pParent, const std::string& sDir):PlotBase(pParent, sDir){ init();}
  void fill(const xAOD::TruthParticle& truthprt, float weight=1.0);
 		
    TH1* truthType;
    TH1* origin;

  private:
    void init();
    void initializePlots();
			
};

}

#endif

