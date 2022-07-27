/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef  TRIGL2MUONSA_MUFASTSTEERING_H
#define  TRIGL2MUONSA_MUFASTSTEERING_H

#include "AthenaBaseComps/AthReentrantAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/ServiceHandle.h"

#include "TrigSteeringEvent/TrigRoiDescriptor.h"

#include "MuFastDataPreparator.h"
#include "MuFastPatternFinder.h"
#include "MuFastStationFitter.h"
#include "MuFastTrackFitter.h"
#include "MuFastTrackExtrapolator.h"
#include "RecMuonRoIUtils.h"
#include "MuCalStreamerTool.h"
#include "CscSegmentMaker.h"
#include "FtfRoadDefiner.h"

//adding a part of DataHandle for AthenaMT
#include "AthenaBaseComps/AthReentrantAlgorithm.h"
#include "StoreGate/ReadHandleKey.h"
#include "StoreGate/WriteHandleKey.h"

#include "xAODEventInfo/EventInfo.h"
#include "xAODTrigMuon/L2StandAloneMuonContainer.h"
#include "xAODTrigMuon/L2CombinedMuonContainer.h"
#include "xAODTrigger/TrigCompositeAuxContainer.h"
#include "xAODTrigger/TrigCompositeContainer.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTrigger/MuonRoIContainer.h"
#include "AthenaMonitoringKernel/GenericMonitoringTool.h"

#include "GaudiKernel/IIncidentListener.h"

class IIncidentSvc;

enum ECRegions{ Bulk, WeakBFieldA, WeakBFieldB };

class MuFastSteering : public AthReentrantAlgorithm , public IIncidentListener
{
 public:
  enum {
    ITIMER_DATA_PREPARATOR=0,
    ITIMER_PATTERN_FINDER,
    ITIMER_STATION_FITTER,
    ITIMER_TRACK_FITTER,
    ITIMER_TRACK_EXTRAPOLATOR,
    ITIMER_CALIBRATION_STREAMER,
    ITIMER_TOTAL_PROCESSING
  };


 public:

  /** Constructor */
  MuFastSteering(const std::string& name, ISvcLocator* svc);

  virtual StatusCode initialize() override;
  virtual StatusCode stop() override;

  /** execute(), main code of the algorithm for AthenaMT*/
  virtual StatusCode execute(const EventContext& ctx) const override;

  /** findMuonSignature(), includes reconstract algorithms **/
  /** this function can be called from both execute() **/
  StatusCode findMuonSignature(const std::vector<const TrigRoiDescriptor*>&	roi,
			       const std::vector<const LVL1::RecMuonRoI*>& 	muonRoIs,
                               DataVector<xAOD::L2StandAloneMuon>& 		outputTracks,
			       TrigRoiDescriptorCollection&	 		outputID,
			       TrigRoiDescriptorCollection&	 		outputMS,
			       const bool                                       dynamicDeltaRpc,
			       const EventContext&                              ctx ) const;

  StatusCode findMuonSignature(const std::vector<const TrigRoiDescriptor*>&	roi,
			       const std::vector<const xAOD::MuonRoI*>& 	muonRoIs,
                               DataVector<xAOD::L2StandAloneMuon>& 		outputTracks,
			       TrigRoiDescriptorCollection&	 		outputID,
			       TrigRoiDescriptorCollection&	 		outputMS,
			       const bool                                       dynamicDeltaRpc,
			       const EventContext&                              ctx ) const;

  /** findMuonSignatureIO(), includes reconstract algorithms for inside-out mode **/
  StatusCode findMuonSignatureIO(const xAOD::TrackParticleContainer&            idtracks,
				 const std::vector<const TrigRoiDescriptor*>&    roids,
				 const std::vector<const LVL1::RecMuonRoI*>&     muonRoIs,
				 DataVector<xAOD::L2CombinedMuon>&              outputCBs,
				 DataVector<xAOD::L2StandAloneMuon>&            outputSAs,
				 const bool                                     dynamicDeltaRpc,
				 const EventContext&                            ctx ) const;

  StatusCode findMuonSignatureIO(const xAOD::TrackParticleContainer&            idtracks,
				 const std::vector<const TrigRoiDescriptor*>&    roids,
				 const std::vector<const xAOD::MuonRoI*>&        muonRoIs,
				 DataVector<xAOD::L2CombinedMuon>&              outputCBs,
				 DataVector<xAOD::L2StandAloneMuon>&            outputSAs,
				 const bool                                     dynamicDeltaRpc,
				 const EventContext&                            ctx ) const;

  /** findMultiTrackSignature(), includes reconstract algorithms for multi-track mode **/
  StatusCode findMultiTrackSignature(const std::vector<const TrigRoiDescriptor*>&	roi,
			             const std::vector<const LVL1::RecMuonRoI*>& 	muonRoIs,
                                     DataVector<xAOD::L2StandAloneMuon>& 		outputTracks,
                                     const bool                                         dynamicDeltaRpc,
                                     const EventContext&                                ctx) const;

  StatusCode findMultiTrackSignature(const std::vector<const TrigRoiDescriptor*>&	roi,
			             const std::vector<const xAOD::MuonRoI*>& 	        muonRoIs,
                                     DataVector<xAOD::L2StandAloneMuon>& 		outputTracks,
                                     const bool                                         dynamicDeltaRpc,
                                     const EventContext&                                ctx) const;

  int L2MuonAlgoMap(const std::string& name) const;

  virtual void handle(const Incident& incident) override;

 protected:

  /**
     Called at the end of the algorithm processing to set the steering
     navigation properly
  */
  bool updateOutputObjects(const LVL1::RecMuonRoI*                        roi,
                           const TrigRoiDescriptor*                       roids,
                           const TrigL2MuonSA::MuonRoad&                  muonRoad,
                           const TrigL2MuonSA::MdtRegion&                 mdtRegion,
                           const TrigL2MuonSA::RpcHits&                   rpcHits,
                           const TrigL2MuonSA::TgcHits&                   tgcHits,
                           const TrigL2MuonSA::RpcFitResult&              rpcFitResult,
                           const TrigL2MuonSA::TgcFitResult&              tgcFitResult,
                           const TrigL2MuonSA::MdtHits&                   mdtHits,
                           const TrigL2MuonSA::CscHits&                   cscHits,
			   const TrigL2MuonSA::StgcHits&                  stgcHits,
			   const TrigL2MuonSA::MmHits&                    mmHits,
                           const std::vector<TrigL2MuonSA::TrackPattern>& trackPatterns,
			   DataVector<xAOD::L2StandAloneMuon>&	          outputTracks,
			   TrigRoiDescriptorCollection&  	          outputID,
			   TrigRoiDescriptorCollection&   	          outputMS,
			   const EventContext&                            ctx) const;

  bool updateOutputObjects(const xAOD::MuonRoI*                           roi,
                           const TrigRoiDescriptor*                       roids,
                           const TrigL2MuonSA::MuonRoad&                  muonRoad,
                           const TrigL2MuonSA::MdtRegion&                 mdtRegion,
                           const TrigL2MuonSA::RpcHits&                   rpcHits,
                           const TrigL2MuonSA::TgcHits&                   tgcHits,
                           const TrigL2MuonSA::RpcFitResult&              rpcFitResult,
                           const TrigL2MuonSA::TgcFitResult&              tgcFitResult,
                           const TrigL2MuonSA::MdtHits&                   mdtHits,
                           const TrigL2MuonSA::CscHits&                   cscHits,
			   const TrigL2MuonSA::StgcHits&                  stgcHits,
			   const TrigL2MuonSA::MmHits&                    mmHits,
                           const std::vector<TrigL2MuonSA::TrackPattern>& trackPatterns,
			   DataVector<xAOD::L2StandAloneMuon>&	          outputTracks,
			   TrigRoiDescriptorCollection&  	          outputID,
			   TrigRoiDescriptorCollection&   	          outputMS,
			   const EventContext&                            ctx) const;

  bool storeMuonSA(const LVL1::RecMuonRoI*             roi,
                   const TrigRoiDescriptor*            roids,
               	   const TrigL2MuonSA::MuonRoad&       muonRoad,
               	   const TrigL2MuonSA::MdtRegion&      mdtRegion,
               	   const TrigL2MuonSA::RpcHits&        rpcHits,
               	   const TrigL2MuonSA::TgcHits&        tgcHits,
               	   const TrigL2MuonSA::RpcFitResult&   rpcFitResult,
               	   const TrigL2MuonSA::TgcFitResult&   tgcFitResult,
               	   const TrigL2MuonSA::MdtHits&        mdtHits,
               	   const TrigL2MuonSA::CscHits&        cscHits,
		   const TrigL2MuonSA::StgcHits&       stgcHits,
		   const TrigL2MuonSA::MmHits&         mmHits,
               	   const TrigL2MuonSA::TrackPattern&   pattern,
                   DataVector<xAOD::L2StandAloneMuon>& outputTracks,
                   const EventContext&                 ctx) const;

  bool storeMuonSA(const xAOD::MuonRoI*                roi,
                   const TrigRoiDescriptor*            roids,
               	   const TrigL2MuonSA::MuonRoad&       muonRoad,
               	   const TrigL2MuonSA::MdtRegion&      mdtRegion,
               	   const TrigL2MuonSA::RpcHits&        rpcHits,
               	   const TrigL2MuonSA::TgcHits&        tgcHits,
               	   const TrigL2MuonSA::RpcFitResult&   rpcFitResult,
               	   const TrigL2MuonSA::TgcFitResult&   tgcFitResult,
               	   const TrigL2MuonSA::MdtHits&        mdtHits,
               	   const TrigL2MuonSA::CscHits&        cscHits,
		   const TrigL2MuonSA::StgcHits&       stgcHits,
		   const TrigL2MuonSA::MmHits&         mmHits,
               	   const TrigL2MuonSA::TrackPattern&   pattern,
                   DataVector<xAOD::L2StandAloneMuon>& outputTracks,
                   const EventContext&                 ctx) const;

  bool storeMSRoiDescriptor(const TrigRoiDescriptor*                  roids,
		            const TrigL2MuonSA::TrackPattern&         pattern,
                            const DataVector<xAOD::L2StandAloneMuon>& outputTracks,
		            TrigRoiDescriptorCollection&	      outputMS) const;


  bool storeIDRoiDescriptor(const TrigRoiDescriptor*                  roids,
		            const TrigL2MuonSA::TrackPattern&         pattern,
                            const DataVector<xAOD::L2StandAloneMuon>& outputTracks,
		            TrigRoiDescriptorCollection&	      outputID) const;

  /**
     Update monitoring variables
  */
  StatusCode updateMonitor(const LVL1::RecMuonRoI*                  roi,
			   const TrigL2MuonSA::MdtHits&             mdtHits,
                           std::vector<TrigL2MuonSA::TrackPattern>& trackPatterns ) const;
  StatusCode updateMonitor(const xAOD::MuonRoI*                     roi,
			   const TrigL2MuonSA::MdtHits&             mdtHits,
                           std::vector<TrigL2MuonSA::TrackPattern>& trackPatterns ) const;



 protected:

  // Tools
  ToolHandle<TrigL2MuonSA::MuFastDataPreparator>     m_dataPreparator {
	this, "DataPreparator", "TrigL2MuonSA::MuFastDataPreparator", "data preparator" };
  ToolHandle<TrigL2MuonSA::MuFastPatternFinder>      m_patternFinder {
	this, "PatternFinder", "TrigL2MuonSA::MuFastPatternFinder", "pattern finder" };
  ToolHandle<TrigL2MuonSA::MuFastStationFitter>      m_stationFitter {
	this, "StationFitter", "TrigL2MuonSA::MuFastStationFitter", "station fitter" };
  ToolHandle<TrigL2MuonSA::MuFastTrackFitter>        m_trackFitter {
	this, "TrackFitter", "TrigL2MuonSA::MuFastTrackFitter", "track fitter" };
  ToolHandle<TrigL2MuonSA::MuFastTrackExtrapolator>  m_trackExtrapolator {
	this, "TrackExtrapolator", "TrigL2MuonSA::MuFastTrackExtrapolator", "track extrapolator" };
  ToolHandle<TrigL2MuonSA::FtfRoadDefiner>           m_ftfRoadDefiner {
        this, "FtfRoadDefiner", "TrigL2MuonSA::FtfRoadDefiner", "ftf road definer" };

  /** Handle to MuonBackExtrapolator tool */
  ToolHandle<ITrigMuonBackExtrapolator> m_backExtrapolatorTool {
	this, "BackExtrapolator", "TrigMuonBackExtrapolator", "public tool for back extrapolating the muon tracks to the IV" };

  // calibration streamer tool
  ToolHandle<TrigL2MuonSA::MuCalStreamerTool> m_calStreamer {
  	this, "CalibrationStreamer", "TrigL2MuonSA::MuCalStreamerTool", "calibration stream" };

  // Utils
  TrigL2MuonSA::RecMuonRoIUtils  m_recMuonRoIUtils;

  //Tools for CSC
  ToolHandle<TrigL2MuonSA::CscSegmentMaker> m_cscsegmaker {
	this, "CscSegmentMaker", "TrigL2MuonSA::CscSegmentMaker", "" };


 private:

  ServiceHandle< IIncidentSvc > m_incidentSvc{this, "IncidentSvc", "IncidentSvc"};      
  ServiceHandle<Gaudi::Interfaces::IOptionsSvc> m_jobOptionsSvc {this, "JobOptionsSvc", "JobOptionsSvc",    "Job options service to retrieve DataFlowConfig"   };
  // Property
  Gaudi::Property< float > m_scaleRoadBarrelInner { this, "Scale_Road_BarrelInner", 1 };
  Gaudi::Property< float > m_scaleRoadBarrelMiddle { this, "Scale_Road_BarrelMiddle", 1 };
  Gaudi::Property< float > m_scaleRoadBarrelOuter { this, "Scale_Road_BarrelOuter", 1 };

  Gaudi::Property< bool > m_use_mcLUT { this, "UseLUTForMC", true};
  Gaudi::Property< bool > m_use_new_segmentfit { this, "USE_NEW_SEGMENTFIT", true};
  Gaudi::Property< bool > m_use_rpc { this, "USE_RPC", true};
  Gaudi::Property< bool > m_use_stgc { this, "USE_STGC", true};
  Gaudi::Property< bool > m_use_mm { this, "USE_MM", true};
  Gaudi::Property< bool > m_use_RoIBasedDataAccess_MDT  { this, "USE_ROIBASEDACCESS_MDT",  true};
  Gaudi::Property< bool > m_use_RoIBasedDataAccess_RPC  { this, "USE_ROIBASEDACCESS_RPC",  true};
  Gaudi::Property< bool > m_use_RoIBasedDataAccess_TGC  { this, "USE_ROIBASEDACCESS_TGC",  true};
  Gaudi::Property< bool > m_use_RoIBasedDataAccess_CSC  { this, "USE_ROIBASEDACCESS_CSC",  true};
  Gaudi::Property< bool > m_use_RoIBasedDataAccess_STGC { this, "USE_ROIBASEDACCESS_STGC", true};
  Gaudi::Property< bool > m_use_RoIBasedDataAccess_MM   { this, "USE_ROIBASEDACCESS_MM",   true};
  Gaudi::Property< bool > m_doCalStream { this, "DoCalibrationStream", true};
  Gaudi::Property< bool > m_calDataScouting { this, "MuonCalDataScouting", false};
  Gaudi::Property< bool > m_rpcErrToDebugStream { this, "RpcErrToDebugStream", false};
  Gaudi::Property< bool > m_use_endcapInnerFromBarrel { this, "UseEndcapInnerFromBarrel", false};


  Gaudi::Property< int > m_esd_rpc_size { this, "ESD_RPC_size", 100 };
  Gaudi::Property< int > m_esd_tgc_size { this, "ESD_TGC_size", 50 };
  Gaudi::Property< int > m_esd_mdt_size { this, "ESD_MDT_size", 100 };
  Gaudi::Property< int > m_esd_csc_size { this, "ESD_CSC_size", 100 };
  Gaudi::Property< int > m_esd_stgc_size { this, "ESD_STGC_size", 100 };
  Gaudi::Property< int > m_esd_mm_size { this, "ESD_MM_size", 100 };

  Gaudi::Property< double > m_rWidth_RPC_Failed { this, "R_WIDTH_RPC_FAILED", 400 };
  Gaudi::Property< double > m_rWidth_TGC_Failed { this, "R_WIDTH_TGC_FAILED", 200 };

  Gaudi::Property< double > m_winPt { this, "WinPt", 4.0 };

  Gaudi::Property< bool > m_insideOut { this, "InsideOutMode", false, "" };
  Gaudi::Property< bool > m_multiTrack { this, "multitrackMode", false, "" };
  Gaudi::Property< bool > m_doEndcapForl2mt { this, "doEndcapForl2mt", false, "" };
  Gaudi::Property< float > m_ftfminPt { this, "FTFminPt", 3500, "pT [MeV] threshold to FTF tracks for L2Muon Inside-out mode" };
  Gaudi::Property< bool > m_topoRoad { this, "topoRoad", false, "create road in barrel not to highly overlap surrounding L1 RoIs" };
  Gaudi::Property< float > m_dPhisurrRoI { this, "dPhisurrRoI", 99, "phi range to find surrounding L1 RoIs" };
  Gaudi::Property< float > m_dEtasurrRoI { this, "dEtasurrRoI", 99, "eta range to find surrounding L1 RoIs" };

  float getRoiSizeForID(bool isEta, const xAOD::L2StandAloneMuon* muonSA) const;

  Gaudi::Property< bool > m_allowOksConfig { this, "AllowOksConfig", true};
  Gaudi::Property< std::string > m_calBufferName { this, "MuonCalBufferName", "/tmp/testOutput"};
  Gaudi::Property< int > m_calBufferSize { this, "MuonCalBufferSize", 1024*1024};

  // Enable to fill FS RoI for ID (cosmic run)
  Gaudi::Property< bool > m_fill_FSIDRoI { this, "FILL_FSIDRoI", false, "Fill FS RoI for ID (will be used in cosmic run)"};

  Gaudi::Property< bool > m_useRun3Config { this, "UseRun3Config", false, "use Run3 L1Muon EDM; xAOD::MuonRoI"};

  //adding a part of DataHandle for AthenaMT
  //ReadHandle xAOD::EventInfo
  SG::ReadHandleKey<xAOD::EventInfo> m_eventInfoKey{
	this, "EventInfo", "EventInfo", "Name of the xAOD::EventInfo object"};

  //ReadHandle MURoIs
  SG::ReadHandleKey<TrigRoiDescriptorCollection> m_roiCollectionKey{
	this, "MuRoIs", "HLT_MURoIs", "Name of the input data from HLTSeeding"};

  //ReadHandle RecMuonRoIs
  SG::ReadHandleKey<DataVector<LVL1::RecMuonRoI>> m_run2recRoiCollectionKey{
	this, "Run2RecMuonRoI", "HLT_RecMURoIs", "Name of the input data on LVL1::RecMuonRoI produced by HLTSeeding"};
  SG::ReadHandleKey<xAOD::MuonRoIContainer> m_recRoiCollectionKey{
	this, "RecMuonRoI", "LVL1MuonRoIs", "Name of the input data on xAOD::MuonRoI"};

  //ReadHandle FTF Tracks
  SG::ReadHandleKey<xAOD::TrackParticleContainer> m_FTFtrackKey{
	this, "TrackParticlesContainerName", "HLT_xAODTracks_Muon", "Name of the input data on xAOD::TrackParticleContainer produced by FTF for Inside-Out mode"};

  //WriteHandle <xAOD::L2StandAloneMuonContainer>
  SG::WriteHandleKey<xAOD::L2StandAloneMuonContainer> m_muFastContainerKey{
	this, "MuonL2SAInfo", "MuonL2SAInfo", "Name of the output data on xAOD::L2StandAloneMuonContainer"};

  //WriteHandle <xAOD::L2StandAloneMuonContainer>
  SG::WriteHandleKey<xAOD::TrigCompositeContainer> m_muCompositeContainerKey{
	this, "MuonCalibrationStream", "MuonCalibrationStream", "Name of the decisions object attached by MuFastSteering"};

  //WriteHandle <TrigRoiDescriptor> for ID
  SG::WriteHandleKey<TrigRoiDescriptorCollection> m_muIdContainerKey{
	this, "forID", "forID", "Name of the output data for Inner Detector"};

  //WriteHandle <TrigRoiDescriptor> for MS
  SG::WriteHandleKey<TrigRoiDescriptorCollection> m_muMsContainerKey{
	this, "forMS", "forMS", "Name of the output data for MS"};

  //WriteHandle <xAOD::L2CombinedMuonContainer> for Inside-out Mode
  SG::WriteHandleKey<xAOD::L2CombinedMuonContainer> m_outputCBmuonCollKey{
	this, "L2IOCB", "MuonL2CBInfo", "output CB Muon container name"};

  // Monitor system
  ToolHandle< GenericMonitoringTool > m_monTool { this, "MonTool", "", "Monitoring tool" };



};

#endif // MUFASTSTEERING_H
