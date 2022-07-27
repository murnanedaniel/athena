/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/
//////////////////////////////////////////////////////////////////////////////////////////////
// Package : RpcRawDataMonitoring
// Author:   N. Benekos(Illinois) - M. Bianco(INFN-Lecce) - G. Chiodini(INFN-Lecce)
// Sept. 2007
//
// DESCRIPTION:
// Subject: RPCLV1-->Offline Muon Data Quality
/////////////////////////////////////////////////////////////////////////////////////////////

#ifndef RpcLv1RawDataValAlg_H
#define RpcLv1RawDataValAlg_H

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/StatusCode.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/NTuple.h"

#include "AthenaMonitoring/AthenaMonManager.h"
#include "AthenaMonitoring/ManagedMonitorToolBase.h"
#include "MuonDQAUtils/MuonDQAHistMap.h"

#include "xAODEventInfo/EventInfo.h"

#include "GaudiKernel/ServiceHandle.h"
#include "MuonIdHelpers/IMuonIdHelperSvc.h"
#include "RPC_CondCabling/RpcCablingCondData.h"
#include "StoreGate/ReadCondHandleKey.h"

#include "MuonReadoutGeometry/MuonDetectorManager.h"
#include "MuonReadoutGeometry/RpcReadoutElement.h"
#include "MuonReadoutGeometry/MuonReadoutElement.h"
#include "MuonReadoutGeometry/RpcReadoutSet.h"
#include "MuonGeoModel/MYSQL.h"

#include "MuonPrepRawData/MuonPrepDataContainer.h"

#include "MuonRDO/RpcFiredChannel.h"
#include "MuonRDO/RpcCoinMatrix.h"
#include "MuonRDO/RpcPad.h"
#include "MuonRDO/RpcPadContainer.h"
#include "MuonRDO/RpcSectorLogicContainer.h"

#include "RpcGlobalUtilities.h"

#include "StoreGate/ReadHandleKey.h"

#include <sstream>
#include <string.h>
#include <vector>
#include <map>
 
class TFile;
template <class ConcreteAlgorithm> class AlgFactory;
/////////////////////////////////////////////////////////////////////////////

class RpcLv1RawDataValAlg: public ManagedMonitorToolBase {
 
 public:

  RpcLv1RawDataValAlg ( const std::string & type, const std::string & name, const IInterface* parent );
  virtual ~RpcLv1RawDataValAlg()=default;
  StatusCode initialize(); 
 
  virtual StatusCode bookHistogramsRecurrent();
  virtual StatusCode fillHistograms();
  virtual StatusCode procHistograms();

 private:

  //Function for histogram booking and parameterd for fitting
  StatusCode bookRPCLV1cmatimevschHistograms(const std::string& m_sectorlogic_name, const std::string& m_tower_name, const std::string& m_cma_name);    
  StatusCode bookRPCLV1TriggerRoadHistograms(const std::string& m_sectorlogic_name, const std::string& m_tower_name, const std::string& m_cma_name, const std::string& m_thr_name);
  StatusCode bookRPCLV1ProfilesHistograms(int m_i_sector, const std::string& m_sectorlogic_name, const std::string& m_cma_name,  int m_i_ijk, const std::string& m_ijk_name) ;  

  MuonDQAHistMap m_stationHists;

  StatusCode StoreTriggerType();
  int GetTriggerType() { return m_trigtype; }
  int m_trigtype = 0;
  
  int m_sector;
  int m_side;

  SG::ReadHandleKey<xAOD::EventInfo> m_eventInfo{this,"EventInfo","EventInfo","event info"};
  SG::ReadHandleKey<RpcPadContainer> m_rpcRdoKey{this,"RpcRdo","RPCPAD","RPC RDO"};

  std::vector<std::string> m_sectorlogicTowerCma_name_list   ;
  std::vector<std::string> m_sectorlogicTowerCma_name_list2  ;
  std::vector<std::string> m_profile_list                    ;
  
  bool m_doClusters         ;
  bool m_checkCabling       ;
  bool m_rpclv1file         ;  
  bool m_rpclv1hist         ; 
  bool m_rpclv1prof         ;
  int  m_rpclv1reducenbins  ;
     
  // MuonDetectorManager from the conditions store
  SG::ReadCondHandleKey<MuonGM::MuonDetectorManager> m_DetectorManagerKey {this, "DetectorManagerKey", 
      "MuonDetectorManager", 
      "Key of input MuonDetectorManager condition data"};    

  ServiceHandle<Muon::IMuonIdHelperSvc> m_idHelperSvc {this, "MuonIdHelperSvc", "Muon::MuonIdHelperSvc/MuonIdHelperSvc"};
  
  SG::ReadCondHandleKey<RpcCablingCondData> m_readKey{this, "ReadKey", "RpcCablingCondData", "Key of RpcCablingCondData"};

  void  bookRPCCoolHistograms(std::vector<std::string>::const_iterator &m_iter, int, int, const std::string& m_layer);
  bool m_doCoolDB;
  
  std::string m_chamberName;
  std::string m_StationSize;
  int m_StationEta;
  int m_StationPhi;
  int m_lastEvent;
  int m_cosmicStation;
};


#endif 
