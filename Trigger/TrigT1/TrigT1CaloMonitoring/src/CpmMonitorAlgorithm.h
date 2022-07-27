/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/
#ifndef TRIGT1CALOMONITORING_CPMMONITORALGORITHM_H
#define TRIGT1CALOMONITORING_CPMMONITORALGORITHM_H

#include "AthenaMonitoring/AthMonitorAlgorithm.h"
#include "AthenaMonitoringKernel/Monitored.h"
#include "StoreGate/ReadHandleKey.h"
#include "StoreGate/WriteHandleKey.h"
#include "xAODTrigL1Calo/CPMTowerContainer.h" 
#include "xAODTrigL1Calo/CPMTobRoIContainer.h"
#include "xAODTrigL1Calo/CMXCPTobContainer.h"
#include "xAODTrigL1Calo/CMXCPHitsContainer.h"
#include "TrigT1Interfaces/TrigT1CaloDefs.h"


class CpmMonitorAlgorithm : public AthMonitorAlgorithm {
 public:CpmMonitorAlgorithm( const std::string& name, ISvcLocator* pSvcLocator );
  virtual ~CpmMonitorAlgorithm()=default;
  virtual StatusCode initialize() override;
  virtual StatusCode fillHistograms( const EventContext& ctx ) const override;

  // monitoring trigger tower structs for various towers
  struct MonitorTT{
    const xAOD::TriggerTower_v2* ttower{};
    double phi_scaled{}; // rescaled for 2D plots
  };

  struct MonitorCpmTT{
    const xAOD::CPMTower_v2* ttower{};
    // some modified/derived information 
    double phi_scaled{}; // rescaled for 2D plots
    int slice{}; // crate * m_maxSlices + emEnergyVec()).size() - 1;
    // errors
    bool emParityError{};
    bool emLinkDownError{};
    bool hadParityError{};
    bool hadLinkDownError{};   
  };

  struct MonitorTobRoI{
    const xAOD::CPMTobRoI_v1* tobroi;
    // some modified/derived information 
    double etaMod{}; 
    double phiMod{}; 
  };

  struct MonitorCmxCpTob{
    const xAOD::CMXCPTob_v1* tob;
    // some modified/derived information 
    uint8_t x{}; // crate * m_modules + cpm - 1
    uint8_t y{}; // chip * 4 + location
    int ybase{}; // cmx * 5
    // errors required to be used as masks
    bool parityError{};
    int ybaseError{}; 
  };

  struct MonitorCmxCpHits{
    const xAOD::CMXCPHits_v1* hit{};
    // some modified/derived information 
    uint8_t crateSlices{}; // crate * m_maxSlices + slices - 1
    uint8_t crateCmx{}; // crate * 2 + cmx
    // source flag
    bool srcTopoCheckSum{};
  };


private:

  // Phi scale for trigger tower eta/phi plots
  double m_phiScaleTT{};

  StringProperty m_packageName{this,"PackageName","CpmMonitor","group name for histograming"};

  Gaudi::Property<int> m_crates{this,"s_crates", 4,  "Number of CPM crates"};
  Gaudi::Property<int> m_modules{this,"s_modules", 14, "Number of modules per crate (modules numbered 1-14)"};
  Gaudi::Property<int> m_maxSlices{this,"s_maxSlices", 5,  "Maximum number of slices"};
  Gaudi::Property<int> m_tobsPerCPM{this,"s_tobsPerCPM", 5,  "Maximum number of TOBs per CPM sent to CMX"};
  Gaudi::Property<int> m_isolBits{this,"s_isolBits", 5,  "Number of bits for encoded isolation"};
  Gaudi::Property<int> m_threshBits{this,"s_threshBits", 3,  "Number of bits per threshold for hit sums"};
  Gaudi::Property<int> m_thresholds{this,"s_thresholds", 16, "Number of EM/Tau threshold bits"};
  Gaudi::Property<int> m_maxTobsPerCmx{this,"MaxTOBsPerCMX", 70,  "Maximum number of TOBs per CMX plotted"};

  // Error vector StoreGate key
  SG::WriteHandleKey<std::vector<int>> m_errorLocation{this,"ErrorLocation","L1CaloCPMErrorVector","Error vector name"};

  // Error summary plot bins
  enum SummaryErrors { EMParity, EMLink, HadParity, HadLink, CPMStatus,
                       TOBParity, SumParity, CMXStatus, NumberOfSummaryBins };


  // container keys including steering parameter and description
  SG::ReadHandleKey<xAOD::TriggerTowerContainer> m_xAODTriggerTowerContainerName{this, "BS_xAODTriggerTowerContainer",LVL1::TrigT1CaloDefs::xAODTriggerTowerLocation,"Trigger Tower Container"};
  SG::ReadHandleKey<xAOD::CPMTowerContainer> m_cpmTowerLocation{this, "CPMTowerLocation", LVL1::TrigT1CaloDefs::CPMTowerLocation, "CPM container"};
  SG::ReadHandleKey<xAOD::CPMTowerContainer> m_cpmTowerLocationOverlap{this, "CPMTowerLocationOverlap",LVL1::TrigT1CaloDefs::CPMTowerLocation + "Overlap", "CPM Overlap container"};
  SG::ReadHandleKey<xAOD::CPMTobRoIContainer> m_cpmTobRoiLocation{this, "CPMTobRoILocation", LVL1::TrigT1CaloDefs::CPMTobRoILocation, "CPMTobRoI container"};
  SG::ReadHandleKey<xAOD::CMXCPTobContainer> m_cmxCpTobLocation{this, "CMXCPTobLocation", LVL1::TrigT1CaloDefs::CMXCPTobLocation, "CMXCPTob container"};
  SG::ReadHandleKey<xAOD::CMXCPHitsContainer> m_cmxCpHitsLocation{this, "CMXCPHitsLocation", LVL1::TrigT1CaloDefs::CMXCPHitsLocation, "CMXCPHits container"};

  std::vector<bool> getIsolationBits(int val, int nThresh, int nBits) const;


  StatusCode fillCpmTowerVectors(SG::ReadHandle<xAOD::CPMTowerContainer> &cpmTower,
				 std::vector<MonitorCpmTT> &monCpmTTs_em, std::vector<MonitorCpmTT> &monCpmTTs_had,
				 std::vector<int> &errorsCPM, 
				 bool core,
				 Monitored::Scalar<int> &cpmLoc,
				 Monitored::Scalar<int> &GLinkParityError
				 ) const;


};
#endif
