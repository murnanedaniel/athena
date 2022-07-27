/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef PFONEUTRALVALIDATIONPLOTS_H
#define PFONEUTRALVALIDATIONPLOTS_H

#include "TrkValHistUtils/PlotBase.h"
#include "PFOHistUtils/PFOPlots.h"
#include "PFOHistUtils/PFOClusterMomentPlots.h"
#include "PFOHistUtils/PFOCalibHitClusterMomentPlots.h"
#include "PFOHistUtils/PFOAttributePlots.h"
#include "PFOHistUtils/FlowElement_LinkerPlots.h"
#include "xAODPFlow/FlowElement.h"
#include "xAODEventInfo/EventInfo.h"

class PFONeutralValidationPlots : public PlotBase {

 public:

  /** Standard Constructor */
  PFONeutralValidationPlots(PlotBase* pParent, const std::string& sDir, const std::string& sFEContainerName);

  /** fill the histograms up */  
  void fill(const xAOD::FlowElement& theFE, const xAOD::EventInfo& eventInfo);

 private:
  /** 4-vector histograms */
  PFO::PFOPlots m_FEPlots;
  /** Cluster Moment histograms */
  PFO::PFOClusterMomentPlots m_FEClusterMomentPlots;
  /** CalibHit Cluster Moment histograms */
  /**  PFO::PFOCalibHitClusterMomentPlots m_FECalibHitClusterMomentPlots; // MC doesn't generally have the relevant calibhits saved. To add at a later date if needed */
  /** FE attributes */
  PFO::PFOAttributePlots m_FEAttributePlots;  
  
  /** Flow element linkers to leptons/photons */
  PFO::FlowElement_LinkerPlots m_FELinkerPlots;
  
};
#endif
