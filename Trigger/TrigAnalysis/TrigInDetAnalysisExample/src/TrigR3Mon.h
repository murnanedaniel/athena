/// emacs: this is -*- c++ -*- 
/**
 **     @file    TrigR3Mon.h
 **
 **     @author  mark sutton
 **     @date    Tue  8 Feb 2022 09:08:26 GMT
 **
 **     Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
 **/


#ifndef TIDAEXAMPLE_TRIGR3MON_H
#define TIDAEXAMPLE_TRIGR3MON_H

#include "GaudiKernel/ToolHandle.h"
#include "AthenaMonitoringKernel/GenericMonitoringTool.h"
#include "AthenaMonitoring/AthMonitorAlgorithm.h"
#include "AthenaMonitoringKernel/Monitored.h"

#include <string>

#include "TrigInDetAnalysis/TrackFilter.h"
#include "TrigInDetAnalysis/TIDARoiDescriptor.h"
#include "TrigInDetAnalysis/TIDDirectory.h"
#include "TrigInDetAnalysis/Efficiency.h"

#include "TrigInDetAnalysisUtils/Filter_Track.h"
#include "TrigInDetAnalysisUtils/Filter_RoiSelector.h"
#include "TrigInDetAnalysisUtils/T_AnalysisConfig.h"
#include "TrigInDetAnalysisUtils/Associator_BestMatch.h"
#include "TrigInDetAnalysisUtils/TrackMatchDeltaR.h"
#include "TrigInDetAnalysisUtils/TrackMatchDeltaRCosmic.h"

#include "TrigInDetAnalysisExample/AnalysisConfig_Tier0.h"


class TrigR3Mon : public AthMonitorAlgorithm {


public:

  TrigR3Mon( const std::string & name, ISvcLocator* pSvcLocator);

  virtual ~TrigR3Mon();

  virtual StatusCode initialize() override;
  virtual StatusCode fillHistograms(const EventContext &context) const override;
  virtual StatusCode finalize() override;

  virtual StatusCode bookHistograms();

  void addMonGroupFromBase( const std::string&  ) { }

protected:

  // track selector cuts

  // test tracks
  double  m_pTCut;
  double  m_etaCut;
  double  m_d0Cut;
  double  m_z0Cut;

  int  m_siHits; // total number of si hits
  int m_pixHits; // pixel hits
  int m_sctHits; // sct hits

  int   m_trtHits; // high threshold hits
  int m_strawHits; // total number of straws

  // reference tracks
  double  m_tauEtCutOffline;
  double  m_doTauThreeProng;
  double  m_pTCutOffline;
  double  m_etaCutOffline;
  double  m_d0CutOffline;
  double  m_z0CutOffline;

  int m_siHitsOffline; // total number of si hits
  int m_pixHitsOffline; // pixel hits
  int m_sctHitsOffline; // sct hits
  int m_blayerHitsOffline;

  int m_pixHolesOffline; // pixel holes
  int m_sctHolesOffline;  // sct holes
  int m_siHolesOffline;   // total pix+sct holes

  int   m_trtHitsOffline; // high threshold hits
  int m_strawHitsOffline; // total number of straws

  // matching parameters
  double m_matchR;   // for DeltaR matcher
  double m_matchPhi; // for DeltaPhi matcher

  ToolHandle<Trig::TrigDecisionTool> m_tdt;

  std::vector<TrackFilter*>  m_filters;
  std::vector<TrackAssociator*>                 m_associators;

  /// do we need this ??? why not the base class ???
  std::vector<T_AnalysisConfig<AthReentrantAlgorithm>*>   m_sequences;

  std::vector<std::string> m_chainNames;
  std::vector<std::string> m_ntupleChainNames;
  std::string              m_releaseMetaData;

  bool m_buildNtuple;
  bool m_mcTruthIn;

  std::string m_analysis_config;
  std::string m_outputFileName;

  bool        m_genericFlag;

  bool        m_initialisePerRun;
  mutable bool        m_firstRun;

  //pdgId
  int m_selectTruthPdgId;

  int m_selectParentTruthPdgId;

  /// kepp events even if they fail the requested trigger chains
  bool  m_keepAllEvents;

  /// if an ntple file open?
  bool m_fileopen;

  /// is this the first event
  mutable bool m_first; 

  /// use only the highest pt tracks
  bool m_useHighestPT;

  /// if performing the vertex analysis, the index of the 
  /// offline vertex to look for
  int m_vtxIndex;

  /// also run purity analyses
  bool m_runPurity;

  /// determine whether this should be treated as a shifter chain
  bool m_shifter;

  /// max number of shifter chains to use - must be < 2 at the moment
  int m_shifterChains;

  /// additional string for the histogram directory
  std::string  m_sliceTag;

  /// do we want basic, or rigorous roi track containment
  bool         m_containTracks;

  bool         m_legacy;

  /// ntuple building variables

  double       m_fiducial_radius;

  bool         m_requireDecision;

  bool         m_filter_on_roi;

  ToolHandleArray<GenericMonitoringTool> m_monTools { this, "MonTools", {} }; // insane configuration paradigm ?

};



#endif //  TIDAEXAMPLE_TRIGR3MON_H
