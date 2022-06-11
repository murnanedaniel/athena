/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Package : sTGCRawDataMonAlg
// Author: Sebastian Fuenzalida Garrido
// Local supervisor: Edson Carquin Lopez
// Technical supervisor: Gerardo Vasquez
//
// DESCRIPTION:
// Subject: sTGC --> sTGC raw data monitoring
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef StgcRawDataMonAlg_H
#define StgcRawDataMonAlg_H

//Core Include
#include "AthenaMonitoring/AthMonitorAlgorithm.h"
#include "AthenaMonitoringKernel/Monitored.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h" 

//Helper Includes
#include "MuonAnalysisInterfaces/IMuonSelectionTool.h"
#include "MuonPrepRawData/MuonPrepDataContainer.h"
#include "MuonPrepRawData/sTgcPrepDataCollection.h"
#include "MuonIdHelpers/IMuonIdHelperSvc.h"
#include "MuonReadoutGeometry/MuonDetectorManager.h"
#include "AthenaMonitoring/DQAtlasReadyFilterTool.h"
#include "xAODTrigger/MuonRoIContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "MuonPrepRawData/sTgcPrepDataContainer.h"
#include "MuonPrepRawData/sTgcPrepData.h"
#include "StoreGate/ReadHandleKey.h"

//stl includes                                                                                 
#include <string>

namespace Muon 
{
  class sTgcPrepData;
}

namespace GeometricSectors
{
  static const std::vector<std::string> sTGC_Side = {"CSide", "ASide"}; 
}

namespace Histograms
{
  struct sTGCSummaryHistogramStruct 
  {
    std::vector<int> strip_charges_vec;
    
    std::vector<int> stationEta_perPhi_vec;
    std::vector<short unsigned int> strip_numbers_perPhi_vec;
    std::vector<int> charge_perPhi_vec;

    std::vector<int> charge_vec;
    std::vector<int> stationPhi_vec;
    std::vector<int> stationEta_vec;
  };
}


class StgcRawDataMonAlg: public AthMonitorAlgorithm 
{
 public:
  
  StgcRawDataMonAlg( const std::string& name, ISvcLocator* pSvcLocator );
  
  virtual ~StgcRawDataMonAlg()=default;
  virtual StatusCode initialize() override;
  virtual StatusCode fillHistograms(const EventContext& ctx) const override;
  
 private:  

  ServiceHandle<Muon::IMuonIdHelperSvc> m_idHelperSvc {this, "MuonIdHelperSvc", "Muon::MuonIdHelperSvc/MuonIdHelperSvc"};

  ToolHandle<CP::IMuonSelectionTool> m_muonSelectionTool{this,"MuonSelectionTool","CP::MuonSelectionTool/MuonSelectionTool"};
  SG::ReadCondHandleKey<MuonGM::MuonDetectorManager> m_DetectorManagerKey {this, "DetectorManagerKey",
      "MuonDetectorManager","Key of input MuonDetectorManager condition data"};
 
  virtual void fillsTGCOverviewHistograms(const std::vector<const Muon::sTgcPrepData*> &prd, const Muon::sTgcPrepData*) const;
  
  virtual StatusCode fillsTGCHistograms(const Muon::sTgcPrepData*) const;                                      
  virtual void fillsTGCSummaryHistograms(const Muon::sTgcPrepData*, Histograms::sTGCSummaryHistogramStruct (&vects)[2][2][4]) const; //[side][multiplet][gas_gap]
  
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

  SG::ReadHandleKey<Muon::sTgcPrepDataContainer> m_sTGCContainerKey{this,"sTGCPrepDataContainerName","STGC_Measurements"};
  SG::ReadHandleKey<xAOD::MuonContainer> m_muonKey{this,"MuonKey","Muons","muons"};

  Gaudi::Property<bool> m_doSTGCESD{this,"DoSTGCESD",true};
  Gaudi::Property<bool> m_do_sTgc_overview{this,"do_sTgc_overview",true};
  Gaudi::Property<bool> m_do_stereoCorrection{this,"do_stereoCorrection",false};
};    
#endif
