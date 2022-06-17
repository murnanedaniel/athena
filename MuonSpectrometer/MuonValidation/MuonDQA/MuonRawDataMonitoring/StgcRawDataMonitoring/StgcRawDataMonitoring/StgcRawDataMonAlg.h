/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Package : sTgcRawDataMonAlg
// Author: Sebastian Fuenzalida Garrido
// Local supervisor: Edson Carquin Lopez
// Technical supervisor: Gerardo Vasquez
//
// DESCRIPTION:
// Subject: sTgc --> sTgc raw data monitoring
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef sTgcRawDataMonAlg_H
#define sTgcRawDataMonAlg_H

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
  static const std::array<std::string, 2> sTgc_Side   = {"CSide", "ASide"};
}

namespace Histograms
{
  struct sTgcSummaryHistogramStruct 
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


class sTgcRawDataMonAlg: public AthMonitorAlgorithm 
{
 public:
  
  sTgcRawDataMonAlg(const std::string& name, ISvcLocator* pSvcLocator);
  
  virtual ~sTgcRawDataMonAlg()=default;
  virtual StatusCode initialize() override;
  virtual StatusCode fillHistograms(const EventContext& ctx) const override;
  
 private:  

  ServiceHandle<Muon::IMuonIdHelperSvc> m_idHelperSvc {this, "MuonIdHelperSvc", "Muon::MuonIdHelperSvc/MuonIdHelperSvc"};

  SG::ReadCondHandleKey<MuonGM::MuonDetectorManager> m_DetectorManagerKey {this, "DetectorManagerKey", "MuonDetectorManager","Key of input MuonDetectorManager condition data"};
 
  void fillsTgcOverviewHistograms(const std::vector<const Muon::sTgcPrepData*> &prd, const Muon::sTgcPrepData*) const;

  void fillsTgcSummaryHistograms(const Muon::sTgcPrepData*, Histograms::sTgcSummaryHistogramStruct (&vects)[2][2][4]) const; //[side][multiplet][gas_gap]

  int get_sectorPhi_from_stationPhi_stName(const int stationPhi, const std::string& stName) const;

  SG::ReadHandleKey<xAOD::TrackParticleContainer> m_meTrkKey{this, "METrkContainer", "ExtrapolatedMuonTrackParticles"};
  SG::ReadHandleKey<Muon::sTgcPrepDataContainer> m_sTgcContainerKey{this,"sTGCPrepDataContainerName", "STGC_Measurements"};

  SG::ReadHandleKey<xAOD::MuonContainer> m_muonKey{this, "MuonKey", "Muons", "muons"};

  Gaudi::Property<bool> m_dosTgcESD{this,"dosTgcESD",true};
  Gaudi::Property<bool> m_dosTgcOverview{this,"dosTgcOverview",true};
};    
#endif
