/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

/*
 *
   author: Max Baugh
   15/9/15
 *
 */

#include "InDet_BadMatchRate.h"
#include <utility>
#include "TProfile.h"

using namespace TMath;

InDet_BadMatchRate::InDet_BadMatchRate(InDetPlotBase *pParent, std::string sDir) :
  InDetPlotBase(pParent, sDir),
  m_BadMatchRate{},
  m_BMR_vs_logpt{},
  m_ReallyFakeRate{},
  m_trackinjet_badmatchrate_vs_dr_gr_j100{} {
  //
}

void
InDet_BadMatchRate::initializePlots() {
  // Bad Match Rate plots, Truth Matching Probability < 50.1%
  book(m_BadMatchRate, "BadMatchRate");
  book(m_BMR_vs_logpt, "BadMatchRate_vs_logpt");

  // Really Fake Rate plots, Truth Matching Probability < 50.0%
  book(m_ReallyFakeRate, "ReallyFakeRate");

  // TrackinJet Bad Match Rate plots
  book(m_trackinjet_badmatchrate_vs_dr_gr_j100, "trackinjet_badmatchrate_vs_dr_gr_j100");
}

void
InDet_BadMatchRate::fillBMR(const xAOD::TrackParticle &particle, float weight) {
  float trketa = particle.eta();
  float trkpt = particle.pt();
  float logpt = Log10(trkpt) - 3.0; // -3 converts from MeV to GeV

  m_BadMatchRate->Fill(trketa, weight);
  m_BMR_vs_logpt->Fill(logpt, weight);
}

void
InDet_BadMatchRate::fillRF(const xAOD::TrackParticle &particle, float weight) {
  float trketa = particle.eta();

  m_ReallyFakeRate->Fill(trketa, weight);
}

void
InDet_BadMatchRate::jetBMR(const xAOD::TrackParticle &track, const xAOD::Jet &jet, float weight) {
  float jet_et = jet.pt() / 1000.; // divide by 1000 to convert to GeV
  float dR = jet.p4().DeltaR(track.p4());

  if (jet_et > 100) {
    m_trackinjet_badmatchrate_vs_dr_gr_j100->Fill(dR, weight);
  }
}
