/*
 * Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
 */

#ifndef IDPERFMON_ELECTRONSELECTOR_H
#define IDPERFMON_ELECTRONSELECTOR_H

//==============================================================================
// Include files...
//==============================================================================
#include "InDetPerformanceMonitoring/EventAnalysis.h"

#include <map>
#include "TH1.h"

#include "xAODMuon/Muon.h"
#include "xAODTracking/Vertex.h"
#include "xAODTracking/VertexContainer.h"

#include "GaudiKernel/ServiceHandle.h"
#include "AsgTools/ToolHandle.h" 
#include "xAODEgamma/Electron.h"
#include "xAODEgamma/ElectronContainer.h"
#include "ElectronPhotonSelectorTools/AsgElectronLikelihoodTool.h"
//==============================================================================
// Forward class declarations...
//==============================================================================
class ElectronSelector : public EventAnalysis
{
 public:
  ElectronSelector();
  ~ElectronSelector();

  void setDebug (bool debug) {m_doDebug = debug;}
 
  // Override functions from EventAnalysis
  void       Init();
  void       PrepareElectronList (const xAOD::ElectronContainer* pxElecContainer);
  bool       RecordElectron (const xAOD::Electron *);
  void       OrderElectronList ();
  // virtual bool Reco();

 protected:
  // virtual void BookHistograms();

 private:
  typedef EventAnalysis PARENT;

  static unsigned int s_uNumInstances;

  // functions 
  void   Clear();

  // message stream
  MsgStream * m_msgStream;

  // Class variables
  const xAOD::Muon*           m_pxElectron;
  std::vector<const xAOD::TrackParticle*>  m_pxElTrackList; 

  // 
  bool m_doDebug;
  // 
  float m_ptCut;

  // Electron likelihood tool:
  AsgElectronLikelihoodTool* m_LHTool2015; //!

  // 
  int m_elecneg1;
  int m_elecneg2;
  int m_elecpos1;
  int m_elecpos2;
};

#endif

