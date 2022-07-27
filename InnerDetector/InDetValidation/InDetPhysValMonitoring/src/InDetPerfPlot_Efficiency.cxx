/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#include "InDetPerfPlot_Efficiency.h"
// #include "TrkValHistUtils/EfficiencyPurityCalculator.h"
#include "xAODTruth/TruthParticle.h"
#include "xAODTruth/TruthVertex.h"
#include "InDetPhysValMonitoringUtilities.h"
#include <cmath>
using namespace IDPVM;

InDetPerfPlot_Efficiency::InDetPerfPlot_Efficiency(InDetPlotBase* pParent, const std::string& sDir) :
  InDetPlotBase(pParent, sDir),
  m_efficiency_vs_pteta{},
  m_efficiency_vs_eta{},
  m_efficiency_vs_pt{},
  m_efficiency_vs_pt_low{},
  m_efficiency_vs_pt_high{},
  m_efficiency_vs_pt_log{},
  m_efficiency_vs_phi{},
  m_efficiency_vs_d0{},
  m_efficiency_vs_z0{},
  m_efficiency_vs_R{},
  m_efficiency_vs_Z{},
  m_efficiency_vs_mu{},
  m_extended_efficiency_vs_d0{},
  m_extended_efficiency_vs_z0{},
  m_efficiency_vs_prodR{},
  m_efficiency_vs_prodR_extended{},
  m_efficiency_vs_prodZ{},
  m_efficiency_vs_prodZ_extended{} {
  // nop
}

void
InDetPerfPlot_Efficiency::initializePlots() {

  book(m_efficiency_vs_pteta, "efficiency_vs_pteta");

  book(m_efficiency_vs_eta, "efficiency_vs_eta");
  book(m_efficiency_vs_pt, "efficiency_vs_pt");
  book(m_efficiency_vs_pt_low, "efficiency_vs_pt_low");
  book(m_efficiency_vs_pt_high, "efficiency_vs_pt_high");
  book(m_efficiency_vs_phi, "efficiency_vs_phi");
  book(m_efficiency_vs_d0, "efficiency_vs_d0");
  book(m_efficiency_vs_z0, "efficiency_vs_z0");
  book(m_efficiency_vs_R, "efficiency_vs_R");
  book(m_efficiency_vs_Z, "efficiency_vs_Z");
  book(m_efficiency_vs_mu, "efficiency_vs_mu");

  book(m_extended_efficiency_vs_d0, "extended_efficiency_vs_d0");
  book(m_extended_efficiency_vs_z0, "extended_efficiency_vs_z0");
  book(m_efficiency_vs_prodR, "efficiency_vs_prodR");
  book(m_efficiency_vs_prodR_extended, "efficiency_vs_prodR_extended");
  book(m_efficiency_vs_prodZ, "efficiency_vs_prodZ");
  book(m_efficiency_vs_prodZ_extended, "efficiency_vs_prodZ_extended");

  book(m_efficiency_vs_pt_log, "efficiency_vs_pt_log");
  const TH1* h = m_efficiency_vs_pt_log->GetTotalHistogram();
  int nbins = h->GetNbinsX();
  std::vector<double> logptbins = IDPVM::logLinearBinning(nbins, h->GetBinLowEdge(1), h->GetBinLowEdge(nbins + 1), false);
  m_efficiency_vs_pt_log->SetBins(nbins, logptbins.data());
}

void
InDetPerfPlot_Efficiency::fill(const xAOD::TruthParticle& truth, const bool isGood, float weight, float mu) {
  double eta = truth.eta();
  double pt = truth.pt() / Gaudi::Units::GeV; // convert MeV to GeV
  double phi = truth.phi();

  fillHisto(m_efficiency_vs_pteta, pt, eta, isGood, weight);

  fillHisto(m_efficiency_vs_eta, eta, isGood, weight);
  fillHisto(m_efficiency_vs_pt, pt, isGood, weight);
  fillHisto(m_efficiency_vs_pt_low, pt, isGood, weight);
  fillHisto(m_efficiency_vs_pt_high, pt, isGood, weight);
  fillHisto(m_efficiency_vs_phi, phi, isGood, weight);
  fillHisto(m_efficiency_vs_pt_log, pt, isGood, weight);

  double d0 = truth.auxdata<float>("d0");
  double z0 = truth.auxdata<float>("z0");
  double R = truth.auxdata<float>("prodR");
  double Z = truth.auxdata<float>("prodZ");
  fillHisto(m_efficiency_vs_d0, d0, isGood, weight);
  fillHisto(m_efficiency_vs_z0, z0, isGood, weight);
  fillHisto(m_efficiency_vs_R, R, isGood, weight);
  fillHisto(m_efficiency_vs_Z, Z, isGood, weight);

  fillHisto(m_extended_efficiency_vs_d0, d0, isGood, weight);
  fillHisto(m_extended_efficiency_vs_z0, z0, isGood, weight);
  fillHisto(m_efficiency_vs_mu, mu, isGood, weight);

  if (truth.hasProdVtx()) {
    const xAOD::TruthVertex* vtx = truth.prodVtx();
    double prod_rad = vtx->perp();
    double prod_z = vtx->z();
    fillHisto(m_efficiency_vs_prodR, prod_rad, isGood, weight);
    fillHisto(m_efficiency_vs_prodR_extended, prod_rad, isGood, weight);
    fillHisto(m_efficiency_vs_prodZ, prod_z, isGood, weight);
    fillHisto(m_efficiency_vs_prodZ_extended, prod_z, isGood, weight);
  }
}

void
InDetPerfPlot_Efficiency::finalizePlots() {
}
