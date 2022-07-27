/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////////////////////////////
// Package : MMRawDataMonitoring
// Author:  M. Biglietti, E. Rossi (Roma Tre)
//
// DESCRIPTION:
// Subject: MM-->Offline Muon Data Quality
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MMRawDataMonAlg_H
#define MMRawDataMonAlg_H

//Core Include
#include "AthenaMonitoring/AthMonitorAlgorithm.h"
#include "AthenaMonitoringKernel/Monitored.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h" 
//Helper Includes

#include "MuonAnalysisInterfaces/IMuonSelectionTool.h"
#include "MuonPrepRawData/MuonPrepDataContainer.h"
#include "MuonPrepRawData/MMPrepDataCollection.h"
#include "MuonIdHelpers/IMuonIdHelperSvc.h"
#include "MuonReadoutGeometry/MuonDetectorManager.h"
#include "AthenaMonitoring/DQAtlasReadyFilterTool.h"
#include "xAODTrigger/MuonRoIContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "MuonPrepRawData/MMPrepDataContainer.h"
#include "MuonPrepRawData/MMPrepData.h"
#include "StoreGate/ReadHandleKey.h"


namespace Muon {
  class MMPrepData;
  }

namespace {
  struct MMOverviewHistogramStruct;
  struct MMSummaryHistogramStruct;
  struct   MMByPhiStruct;
}

//stl includes                                                                                              
#include <string>

class MMRawDataMonAlg: public AthMonitorAlgorithm {
 public:

  MMRawDataMonAlg( const std::string& name, ISvcLocator* pSvcLocator );

  //  virtual ~MMRawDataMonAlg();
  virtual ~MMRawDataMonAlg()=default;
  virtual StatusCode initialize() override;
  virtual StatusCode fillHistograms(const EventContext& ctx) const override;
  
 private:

  ServiceHandle<Muon::IMuonIdHelperSvc> m_idHelperSvc {this, "MuonIdHelperSvc", "Muon::MuonIdHelperSvc/MuonIdHelperSvc"};

  ToolHandle<CP::IMuonSelectionTool> m_muonSelectionTool{this,"MuonSelectionTool","CP::MuonSelectionTool/MuonSelectionTool"};
  SG::ReadCondHandleKey<MuonGM::MuonDetectorManager> m_DetectorManagerKey {this, "DetectorManagerKey", "MuonDetectorManager","Key of input MuonDetectorManager condition data"};
  SG::ReadHandleKey<Trk::SegmentCollection> m_segm_type{this,"Eff_segm_type","TrackMuonSegments","muon segments"};
  SG::ReadHandleKey<Muon::MMPrepDataContainer> m_MMContainerKey{this,"MMPrepDataContainerName","MM_Measurements"};
  SG::ReadHandleKey<xAOD::MuonContainer> m_muonKey{this,"MuonKey","Muons","muons"};
  SG::ReadHandleKey<xAOD::TrackParticleContainer> m_meTrkKey{this, "METrkContainer", "ExtrapolatedMuonTrackParticles"};

  virtual StatusCode  fillMMOverviewVects(const Muon::MMPrepData*, MMOverviewHistogramStruct& vects, MMByPhiStruct (&occupancyPlots)[16][2]) const;
  virtual void  fillMMOverviewHistograms(const MMOverviewHistogramStruct& vects, MMByPhiStruct (&occupancyPlots)[16][2], const int lb) const;
  virtual StatusCode  fillMMSummaryVects( const Muon::MMPrepData*, MMSummaryHistogramStruct (&vects)[2][16][2][2][4]) const; //[side][stationPhi][stationEta][multiplet][gas_gap]
  virtual StatusCode  fillMMHistograms( const Muon::MMPrepData* ) const;                                      
  virtual StatusCode  fillMMSummaryHistograms( const MMSummaryHistogramStruct (&vects)[2][16][2][2][4]) const;

  void clusterFromTrack(const xAOD::TrackParticleContainer*,const int lb) const;
  void clusterFromSegments(const Trk::SegmentCollection*, const int lb) const;
  
  int get_PCB_from_channel(const int channel) const;
  int get_sectorPhi_from_stationPhi_stName(const int stationPhi, const std::string& stName) const;
  int get_sectorEta_from_stationEta(const int stationEta) const;

  int get_bin_for_occ_CSide_hist(const int stationEta, const int multiplet, const int gas_gap) const;
  int get_bin_for_occ_ASide_hist(const int stationEta, const int multiplet, const int gas_gap) const;
  int get_bin_for_occ_CSide_pcb_eta2_hist(const int stationEta, const int multiplet, const int gas_gap, const int PCB) const;
  int get_bin_for_occ_CSide_pcb_eta1_hist(const int stationEta, const int multiplet, const int gas_gap, const int PCB) const;
  int get_bin_for_occ_ASide_pcb_eta2_hist(const int stationEta, const int multiplet, const int gas_gap, const int PCB) const;
  int get_bin_for_occ_ASide_pcb_eta1_hist(const int stationEta, const int multiplet, const int gas_gap, const int PCB) const;
  int get_bin_for_occ_lb_CSide_pcb_eta2_hist(const int stationEta, const int multiplet, const int gas_gap, const int PCB,const int isector) const;
  int get_bin_for_occ_lb_CSide_pcb_eta1_hist(const int stationEta, const int multiplet, const int gas_gap, const int PCB,int isector) const;
  int get_bin_for_occ_lb_ASide_pcb_eta1_hist(const int stationEta, const int multiplet, const int gas_gap, const int PCB,int isector) const;
  int get_bin_for_occ_lb_ASide_pcb_eta2_hist(const int stationEta, const int multiplet, const int gas_gap, const int PCB, const int isector) const;
  int get_bin_for_occ_lb_pcb_hist(const int multiplet, const int gas_gap, const int PCB) const;

  void MMEfficiency(const xAOD::TrackParticleContainer*) const;


  Gaudi::Property<bool> m_doMMESD{this,"DoMMESD",true};
  Gaudi::Property<bool> m_do_mm_overview{this,"do_mm_overview",true};
  Gaudi::Property<bool> m_do_stereoCorrection{this,"do_stereoCorrection",false};
  Gaudi::Property<float> m_cut_pt{this,"cut_pt",15000};

};    
#endif
