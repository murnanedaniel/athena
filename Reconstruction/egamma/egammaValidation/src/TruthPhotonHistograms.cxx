/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
*/

#include "TruthPhotonHistograms.h"

using namespace egammaMonitoring;

#include "GaudiKernel/ServiceHandle.h"
#include "AsgTools/AnaToolHandle.h"

StatusCode TruthPhotonHistograms::initializePlots() {

  ATH_CHECK(ParticleHistograms::initializePlots());

  histoMap["convRadius"] = new TH1D(Form("%s_%s",m_name.c_str(),"convRadius"), ";Conversion Radius [mm]; Conversion Radius Events", 14, m_cR_bins);

  ATH_CHECK(m_rootHistSvc->regHist(m_folder+"convRadius", histoMap["convRadius"]));

  return StatusCode::SUCCESS;

}

void TruthPhotonHistograms::fill(const xAOD::IParticle& phrec) {

  float trueR = -999;

  ParticleHistograms::fill(phrec);

  const xAOD::TruthParticle *tmp  = xAOD::TruthHelpers::getTruthParticle(phrec);

  if (tmp) {
    if (tmp->pdgId() == 22 && tmp->hasDecayVtx()) {

      float x = tmp->decayVtx()->x();
      float y = tmp->decayVtx()->y();
      trueR = sqrt( x*x + y*y );

    }
  }

  histoMap["convRadius"]->Fill(trueR);


} // fill
