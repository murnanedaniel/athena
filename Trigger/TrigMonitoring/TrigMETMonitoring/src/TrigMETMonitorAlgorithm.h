/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TRIGMETMONITORING_TRIGMETMONITORALGORITHM_H
#define TRIGMETMONITORING_TRIGMETMONITORALGORITHM_H

#include "AthenaMonitoring/AthMonitorAlgorithm.h"
#include "AthenaMonitoringKernel/Monitored.h"

#include "StoreGate/ReadHandleKey.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODTrigger/EnergySumRoI.h" 
#include <xAODTrigger/jFexMETRoI.h>
#include <xAODTrigger/jFexMETRoIContainer.h>
#include <xAODTrigger/jFexSumETRoI.h>
#include <xAODTrigger/jFexSumETRoIContainer.h>
#include <xAODTrigger/gFexGlobalRoI.h>
#include <xAODTrigger/gFexGlobalRoIContainer.h>
#include "xAODTrigMissingET/TrigMissingETContainer.h" 
#include "xAODTrigMissingET/TrigMissingETAuxContainer.h" 
#include "xAODMuon/MuonContainer.h"
#include "xAODMuon/Muon.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/Electron.h"
#include "xAODCaloEvent/CaloClusterContainer.h"
#include "xAODCaloEvent/CaloCluster.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTracking/TrackParticle.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/Vertex.h"

#include "xAODEventInfo/EventInfo.h"

#include "TrigDecisionInterface/ITrigDecisionTool.h"

class TrigMETMonitorAlgorithm : public AthMonitorAlgorithm {
 public:
  TrigMETMonitorAlgorithm( const std::string& name, ISvcLocator* pSvcLocator );
  virtual ~TrigMETMonitorAlgorithm();
  virtual StatusCode initialize() override;
  virtual StatusCode fillHistograms( const EventContext& ctx ) const override;

 private:
  double signed_log(double e, double epsilon) const;

  SG::ReadHandleKey<xAOD::MissingETContainer> m_offline_met_key{this, "offline_met_key", "MET_Reference_AntiKt4EMPFlow", "Offline met container name"};
  
  SG::ReadHandleKey<xAOD::ElectronContainer> m_hlt_electron_key{this, "hlt_electron_key", "HLT_egamma_Electrons_GSF", "HLT electron container name"};
  SG::ReadHandleKey<xAOD::MuonContainer> m_hlt_muon_key{this, "hlt_muon_key", "HLT_MuonsCB_FS", "HLT muon container name"};

  SG::ReadHandleKey<xAOD::CaloClusterContainer> m_topoclusters_key{this, "topoclusters_key", "HLT_TopoCaloClustersLCFS", "HLT topoclusters container name"};
  SG::ReadHandleKey<xAOD::TrackParticleContainer> m_tracks_key{this, "tracks_key", "HLT_IDTrack_FS_FTF", "HLT tracks container name"};
  SG::ReadHandleKey<xAOD::VertexContainer> m_vertex_key{this, "vertex_key", "HLT_IDVertex_FS", "HLT vertex container name"};
  SG::ReadHandleKey<xAOD::VertexContainer> m_offline_vertex_key{this, "offline_vertex_key", "PrimaryVertices", "Offline vertex container name"};

  SG::ReadHandleKey<xAOD::EnergySumRoI> m_lvl1_roi_key{this, "l1_roi_key", "LVL1EnergySumRoI", "L1 energy sum ROI container name"};
  SG::ReadHandleKey<xAOD::EnergySumRoI> m_lvl1_jnc_key{this, "l1_jnc_key", "jNOISECUT_MET", "jNOISE Cut MET container name"};
  SG::ReadHandleKey<xAOD::EnergySumRoI> m_lvl1_jrho_key{this, "l1_jrho_key", "jXERHO_MET", "jXERHO MET container name"};
  SG::ReadHandleKey<xAOD::EnergySumRoI> m_lvl1_gnc_key{this, "l1_gnc_key", "gXENOISECUT_MET", "gXENOISECUT MET container name"};
  SG::ReadHandleKey<xAOD::EnergySumRoI> m_lvl1_grho_key{this, "l1_grho_key", "gXERHO_MET", "gXERHO MET container name"};
  SG::ReadHandleKey<xAOD::EnergySumRoI> m_lvl1_gjwoj_key{this, "l1_gjwoj_key", "gXEJWOJ_MET", "gXERJWOW MET container name"};
  SG::ReadHandleKey<xAOD::EnergySumRoI> m_lvl1_gpufit_key{this, "l1_gpufit_key", "gXEPUFIT_MET", "gXEPUFIT MET container name"};

  SG::ReadHandleKey<xAOD::jFexMETRoIContainer> m_l1_jFexMet_key{this, "l1_jFexMet_key", "L1_jFexMETRoI", "L1 jFex MET container name"};
  SG::ReadHandleKey<xAOD::jFexSumETRoIContainer> m_l1_jFexSumEt_key{this, "l1_jFexSumEt_key", "L1_jFexSumETRoI", "L1 jFex sumEt container name"};
  SG::ReadHandleKey<xAOD::gFexGlobalRoIContainer> m_l1_gFexJwojScalar_key{this, "l1_gFexJwojScalar_key", "L1_gScalarEJwoj", "L1 gFex JWOJ Et and sumEt container name"};
  SG::ReadHandleKey<xAOD::gFexGlobalRoIContainer> m_l1_gFexJwojMETComponents_key{this, "l1_gFexJwojMETComponents_key", "L1_gMETComponentsJwoj", "L1 gFex JWOJ Ex and Ey container name"};
  SG::ReadHandleKey<xAOD::gFexGlobalRoIContainer> m_l1_gFexJwojMHTComponents_key{this, "l1_gFexJwojMHTComponents_key", "L1_gMHTComponentsJwoj", "L1 gFex JWOJ Hard Term Ex and Ey container name"};
  SG::ReadHandleKey<xAOD::gFexGlobalRoIContainer> m_l1_gFexJwojMSTComponents_key{this, "l1_gFexJwojMSTComponents_key", "L1_gMSTComponentsJwoj", "L1 gFex JWOJ Soft Term Ex and Ey container name"};
  SG::ReadHandleKey<xAOD::gFexGlobalRoIContainer> m_l1_gFexNCMETScalar_key{this, "l1_gFexNCMETScalar_key", "L1_gScalarENoiseCut", "L1 gFex NC Et and sumEt container name"};
  SG::ReadHandleKey<xAOD::gFexGlobalRoIContainer> m_l1_gFexNCMETComponents_key{this, "l1_gFexNCMETComponents_key", "L1_gMETComponentsNoiseCut", "L1 gFex NC Ex and Ey container name"};
  SG::ReadHandleKey<xAOD::gFexGlobalRoIContainer> m_l1_gFexRhoMETScalar_key{this, "l1_gFexRhoMETScalar_key", "L1_gScalarERms", "L1 gFex Rho Et and sumEt container name"};
  SG::ReadHandleKey<xAOD::gFexGlobalRoIContainer> m_l1_gFexRhoMETComponents_key{this, "l1_gFexRhoMETComponents_key", "L1_gMETComponentsRms", "L1 gFex Rho Ex and Ey container name"};

  SG::ReadHandleKey<xAOD::TrigMissingETContainer> m_hlt_cell_met_key{this, "hlt_cell_key", "HLT_MET_cell", "HLT Cell MET container name"};
  SG::ReadHandleKey<xAOD::TrigMissingETContainer> m_hlt_mht_met_key{this, "hlt_mht_key", "HLT_MET_mht", "HLT MHT MET container name"};
  SG::ReadHandleKey<xAOD::TrigMissingETContainer> m_hlt_tc_met_key{this, "hlt_tc_key", "HLT_MET_tc", "HLT TC MET container name"};
  SG::ReadHandleKey<xAOD::TrigMissingETContainer> m_hlt_tc_em_met_key{this, "hlt_tc_em_key", "HLT_MET_tc_em", "HLT TC EM MET container name"};
  SG::ReadHandleKey<xAOD::TrigMissingETContainer> m_hlt_tcpufit_met_key{this, "hlt_tcpufit_key", "HLT_MET_tcpufit", "HLT TCPufit MET container name"};
  SG::ReadHandleKey<xAOD::TrigMissingETContainer> m_hlt_tcpufit_sig30_met_key{this, "hlt_tcpufit_sig30_key", "HLT_MET_tcpufit_sig30", "HLT TCPufit sig30 MET container name"};
  SG::ReadHandleKey<xAOD::TrigMissingETContainer> m_hlt_trkmht_met_key{this, "hlt_trkmht_key", "HLT_MET_trkmht", "HLT TrkMht MET container name"};
  SG::ReadHandleKey<xAOD::TrigMissingETContainer> m_hlt_pfsum_met_key{this, "hlt_pfsum_key", "HLT_MET_pfsum", "HLT Pfsum MET container name"};
  SG::ReadHandleKey<xAOD::TrigMissingETContainer> m_hlt_pfsum_cssk_met_key{this, "hlt_pfsum_cssk_key", "HLT_MET_pfsum_cssk", "HLT Pfsum CSSK MET container name"};
  SG::ReadHandleKey<xAOD::TrigMissingETContainer> m_hlt_pfsum_vssk_met_key{this, "hlt_pfsum_vssk_key", "HLT_MET_pfsum_vssk", "HLT Pfsum VSSK MET container name"};
  SG::ReadHandleKey<xAOD::TrigMissingETContainer> m_hlt_pfopufit_met_key{this, "hlt_pfopufit_key", "HLT_MET_pfopufit", "HLT PfoPufit MET container name"};
  SG::ReadHandleKey<xAOD::TrigMissingETContainer> m_hlt_pfopufit_sig30_met_key{this, "hlt_pfopufit_sig30_key", "HLT_MET_pfopufit_sig30", "HLT PfoPufit sig30 MET container name"};
  SG::ReadHandleKey<xAOD::TrigMissingETContainer> m_hlt_cvfpufit_met_key{this, "hlt_cvfpufit_key", "HLT_MET_cvfpufit", "HLT CvfPufit MET container name"};
  SG::ReadHandleKey<xAOD::TrigMissingETContainer> m_hlt_mhtpufit_pf_met_key{this, "hlt_mhtpufit_pf_key", "HLT_MET_mhtpufit_pf", "HLT MhtPufitPf MET container name"};
  SG::ReadHandleKey<xAOD::TrigMissingETContainer> m_hlt_mhtpufit_em_met_key{this, "hlt_mhtpufit_em_key", "HLT_MET_mhtpufit_em", "HLT MhtPufitEm MET container name"};
  SG::ReadHandleKey<xAOD::TrigMissingETContainer> m_hlt_met_nn_key{this, "hlt_met_nn_key", "HLT_MET_nn", "HLT MET NN container name"};

  Gaudi::Property<std::vector<std::string>> m_l1Chains{this, "L1Chains", {}, "The L1 chains to monitor"};
  Gaudi::Property<std::vector<std::string>> m_hltChains{this, "HLTChains", {}, "The HLT shifter chains to monitor"};
  Gaudi::Property<std::vector<std::string>> m_hltChainsVal{this, "HLTChainsVal", {}, "The HLT val chains to monitor"};
  Gaudi::Property<std::vector<std::string>> m_hltChainsT0{this, "HLTChainsT0", {}, "The HLT t0 chains to monitor"};
  Gaudi::Property<std::vector<std::string>> m_algsL1{this, "algsL1", {}, "L1 algorithms to monitor"};
  Gaudi::Property<std::vector<std::string>> m_algsHLT{this, "algsHLT", {}, "HLT algorithms to monitor"};
  Gaudi::Property<std::vector<std::string>> m_algsHLT2d{this, "algsHLT2d", {}, "HLT algorithms for 2d eta-phi plots"};
  Gaudi::Property<std::vector<std::string>> m_algsL1Expert{this, "algsL1Expert", {}, "L1 algorithms for Expert"};
  Gaudi::Property<std::vector<std::string>> m_algsHLTExpert{this, "algsHLTExpert", {}, "HLT algorithms for Expert"};
  Gaudi::Property<std::vector<std::string>> m_algsMET2d_tcpufit{this, "algsMET2d_tcpufit", {}, "HLT algorithms for 2D MET wrt tcpufit"};
  Gaudi::Property<std::vector<std::string>> m_compNames{this, "compNames", {}, "Calorimeter component names"};
  Gaudi::Property<std::vector<std::string>> m_bitNames{this, "bitNames", {}, "Status bit names"};

  Gaudi::Property<double> m_electronPtCut{this, "electronPtCut", 0.0, "Electron pt cut for leading electron"};
  Gaudi::Property<double> m_electronEtaCut{this, "electronEtaCut", 0.0, "Electron eta cut for leading electron"};
  Gaudi::Property<double> m_muonPtCut{this, "muonPtCut", 0.0, "Muon pt cut for leading muon"};
  Gaudi::Property<double> m_muonEtaCut{this, "muonEtaCut", 0.0, "Muon eta cut for leading muon"};
  Gaudi::Property<std::vector<std::string>> m_signalLepAlgs{this, "signalLepAlgs", {}, "Signal lepton MET histograms"};
};
#endif
