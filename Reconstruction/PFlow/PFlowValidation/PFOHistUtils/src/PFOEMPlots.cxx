/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#include "PFOHistUtils/PFOEMPlots.h"

namespace PFO {

  PFOEMPlots::PFOEMPlots(PlotBase* pParent, std::string sDir, std::string sPFOContainerName, std::string sFEContainerName) : PlotBase(pParent, sDir), m_sPFOContainerName(sPFOContainerName), m_sFEContainerName(sFEContainerName) {
    m_PFO_ptEM = nullptr;
    m_PFO_etaEM = nullptr;
    m_PFO_phiEM = nullptr;
    m_PFO_mEM = nullptr;
    m_FE_ptEM = nullptr;
    m_FE_etaEM = nullptr;
    m_FE_phiEM = nullptr;
    m_FE_mEM = nullptr;
  }

  void PFOEMPlots::initializePlots(){
    if(!m_sPFOContainerName.empty()){
      m_PFO_ptEM = Book1D("_PtEM",m_sPFOContainerName + "_PtEM (Entries/1 GeV)",30,-10.0,20.0);
      m_PFO_etaEM = Book1D("_EtaEM",m_sPFOContainerName + "_EtaEM (Entries/0.1)",100,-5.0,5.0);
      m_PFO_phiEM = Book1D("_PhiEM",m_sPFOContainerName + "_PhiEM (Entries/0.1)",64,-3.2,3.2);
      m_PFO_mEM = Book1D("_mEM",m_sPFOContainerName + "_mEM (Entries/100 MeV)",10,0.0,0.5);
    }
    if(!m_sFEContainerName.empty()){
      m_FE_ptEM = Book1D("_PtEM",m_sFEContainerName + "_PtEM (Entries/1 GeV)",30,-10.0,20.0);
      m_FE_etaEM = Book1D("_EtaEM",m_sFEContainerName + "_EtaEM (Entries/0.1)",100,-5.0,5.0);
      m_FE_phiEM = Book1D("_PhiEM",m_sFEContainerName + "_PhiEM (Entries/0.1)",64,-3.2,3.2);
      m_FE_mEM = Book1D("_mEM",m_sFEContainerName + "_mEM (Entries/100 MeV)",10,0.0,0.5);
    }
  }

  void PFOEMPlots::fill(const xAOD::PFO& PFO, const xAOD::EventInfo& eventInfo){
    m_PFO_ptEM->Fill(PFO.ptEM()/1000.0,eventInfo.beamSpotWeight());
    m_PFO_etaEM->Fill(PFO.etaEM(),eventInfo.beamSpotWeight());
    m_PFO_phiEM->Fill(PFO.phiEM(),eventInfo.beamSpotWeight());
    m_PFO_mEM->Fill(PFO.mEM()/1000.0,eventInfo.beamSpotWeight());
  }  
}
