/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#include "TrigT1RPCRecRoiTool.h"
#include "MuonIdHelpers/IMuonIdHelperSvc.h"
#include "MuonReadoutGeometry/MuonDetectorManager.h"
#include "MuonReadoutGeometry/RpcReadoutElement.h"
#include "RPC_CondCabling/RpcCablingCondData.h"
#include "RPCcablingInterface/IRPCcablingServerSvc.h"
#include <fstream>

namespace LVL1 {
  TrigT1RPCRecRoiTool::TrigT1RPCRecRoiTool( const std::string& type, const std::string& name, const IInterface* parent) : base_class(type,name,parent)
  {
  }
  
  TrigT1RPCRecRoiTool::~TrigT1RPCRecRoiTool() {
  }
  
  StatusCode TrigT1RPCRecRoiTool::initialize() {
    ATH_CHECK(m_DetectorManagerKey.initialize());
    ATH_CHECK( detStore()->retrieve(m_muonMgr) );
    ATH_CHECK( m_idHelperSvc.retrieve() );
    ServiceHandle<IRPCcablingServerSvc> rpcCabGet ("RPCcablingServerSvc", name());
    ATH_CHECK( rpcCabGet.retrieve() );
    ATH_CHECK( rpcCabGet->giveCabling(m_cabling) );
    if(m_useRun3Config){
      ATH_MSG_INFO("update to Run 3 bit mask");
      updateBitMask( Run3 );
    }
    return StatusCode::SUCCESS;
  }
 
  TrigT1MuonRecRoiData TrigT1RPCRecRoiTool::roiData( const unsigned int & roiWord )const
  {
    TrigT1MuonRecRoiData data;

    // subsystem ID and  sector ID
    data.set_side( getBitMaskValue(&roiWord,SubSysIDMask()) );
    data.set_sector( getBitMaskValue(&roiWord,BarrelSectorIDMask()) );
    // RoI
    data.set_roi( getBitMaskValue(&roiWord,BarrelRoIMask()) );
    
    // Get the strips delimiting the RoIs from rPCcablingSvc
    Identifier EtaLowBorder_id;
    Identifier EtaHighBorder_id;
    Identifier PhiLowBorder_id;
    Identifier PhiHighBorder_id;
    
    Amg::Vector3D EtaLowBorder_pos(0.,0.,0.);
    Amg::Vector3D EtaHighBorder_pos(0.,0.,0.);
    Amg::Vector3D PhiLowBorder_pos(0.,0.,0.);
    Amg::Vector3D PhiHighBorder_pos(0.,0.,0.);
    
    if(m_cabling->give_RoI_borders_id(data.side(), data.sector(), data.roi(),
				      EtaLowBorder_id, EtaHighBorder_id,
				      PhiLowBorder_id, PhiHighBorder_id)) 
      {

	const MuonGM::MuonDetectorManager* muonMgr = m_muonMgr;
	if(m_useConditionData){
	  SG::ReadCondHandle<MuonGM::MuonDetectorManager> DetectorManagerHandle{m_DetectorManagerKey};
	  muonMgr = DetectorManagerHandle.cptr(); 
	  if(muonMgr==nullptr){
	    ATH_MSG_ERROR("Null pointer to the read MuonDetectorManager conditions object. Use the one from DetectorStore");
	    muonMgr = m_muonMgr;
	  }
	}
  
	const MuonGM::RpcReadoutElement* EtaLowBorder_descriptor =
	  muonMgr->getRpcReadoutElement(EtaLowBorder_id);
	EtaLowBorder_pos = EtaLowBorder_descriptor->stripPos(EtaLowBorder_id);
	
	const MuonGM::RpcReadoutElement* EtaHighBorder_descriptor =
	  muonMgr->getRpcReadoutElement(EtaHighBorder_id);
	EtaHighBorder_pos = EtaHighBorder_descriptor->stripPos(EtaHighBorder_id);
	
	const MuonGM::RpcReadoutElement* PhiLowBorder_descriptor =
	  muonMgr->getRpcReadoutElement(PhiLowBorder_id);
	PhiLowBorder_pos = PhiLowBorder_descriptor->stripPos(PhiLowBorder_id);
	
	const MuonGM::RpcReadoutElement* PhiHighBorder_descriptor =
	  muonMgr->getRpcReadoutElement(PhiHighBorder_id);
	PhiHighBorder_pos =   PhiHighBorder_descriptor->stripPos(PhiHighBorder_id);
	
	data.set_etaMin( EtaLowBorder_pos.eta() );
	data.set_etaMax( EtaHighBorder_pos.eta() );
	data.set_eta( (data.etaMin()+data.etaMax())/2. );
	
	data.set_phiMin( PhiLowBorder_pos.phi() );
	data.set_phiMax( PhiHighBorder_pos.phi() );
	data.set_phi( (data.phiMin()+data.phiMax())/2. );

	data.update(); 
	
      }
    return data;
  }
  void TrigT1RPCRecRoiTool::RoIsize(const unsigned int & roiWord,
				    double & etaMin_LowHigh, double & etaMax_LowHigh,
				    double & phiMin_LowHigh, double & phiMax_LowHigh) const
  {
    double etaMin_Low=0;
    double etaMin_High=0;
    double etaMax_Low=0;
    double etaMax_High=0;
    
    SG::ReadCondHandle<RpcCablingCondData> readHandle{m_readKey};
    std::unique_ptr<const RpcCablingCondData> readCdo(*readHandle);
    
    auto data = roiData(roiWord);
    phiMin_LowHigh=data.phiMin();
    phiMax_LowHigh=data.phiMax();
    
    bool low  = etaDimLow(data,etaMin_Low,etaMax_Low, std::move(readCdo));
    bool high = etaDimHigh(data,etaMin_High,etaMax_High, std::move(readCdo));
    
    if (low&&high) {
      etaMin_LowHigh=std::min(etaMin_Low,etaMin_High);
      etaMax_LowHigh=std::max(etaMax_Low,etaMax_High);
    } else if (low) {
      etaMin_LowHigh =std::min(etaMin_Low,data.etaMin());
      etaMax_LowHigh =std::max(etaMax_Low,data.etaMax());
    } else if (high) {
      etaMin_LowHigh =std::min(etaMin_High,data.etaMin());
      etaMax_LowHigh =std::max(etaMax_High,data.etaMax());
    } else  {
      etaMin_LowHigh = data.etaMin();
      etaMax_LowHigh = data.etaMax();
    }
    return;
  }

  bool TrigT1RPCRecRoiTool::dumpRoiMap(const std::string& filename) const
  {
    const unsigned int maxSubsystem = 2;
    const unsigned int maxLogicSector = 32;
    const unsigned int maxRoI = 32;
    
    SG::ReadCondHandle<RpcCablingCondData> readHandle{m_readKey};
    std::unique_ptr<const RpcCablingCondData> readCdo(*readHandle);
    
    std::ofstream roi_map;
    roi_map.open(filename.c_str(), std::ios::out );
    if(!roi_map){
      ATH_MSG_WARNING( "Unable to open ROI_Mapping file!" );
    } else {
      roi_map <<"# side     sector   roi      etaMin       etaMax       phiMin       phiMax     etaMinLow    etaMaxLow    etaMinHigh   etaMaxHigh "<< std::endl;
      roi_map <<"# ------------------------------------------------------------------------------------------------------------------------------ "<< std::endl;
      for(unsigned int side=0;side < maxSubsystem; side++){
	for(unsigned int sector=0;sector < maxLogicSector; sector++){
	  for (unsigned int roi=0; roi<maxRoI; roi++){
	    unsigned long int roIWord = (m_useRun3Config) ? (roi+(side<<21)+(sector<<22)) : ((roi<<2)+(side<<14)+(sector<<15));
	    auto data = roiData(roIWord);
	    double etaMinLow(0),etaMaxLow(0),etaMinHigh(0),etaMaxHigh(0);
	    etaDimLow (data,etaMinLow, etaMaxLow, std::move(readCdo));
	    etaDimHigh(data,etaMinHigh,etaMaxHigh, std::move(readCdo));
	    roi_map << std::setw(8)  << side     << " "
		    << std::setw(8)  << sector   << " "
		    << std::setw(8)  << roi      << " "
		    << std::setw(12) << data.etaMin() << " "
		    << std::setw(12) << data.etaMax() << " "
		    << std::setw(12) << data.phiMin() << " "
		    << std::setw(12) << data.phiMax() << " "
		    << std::setw(8) << etaMinLow << " "
		    << std::setw(8) << etaMaxLow << " "
		    << std::setw(8) << etaMinHigh << " "
		    << std::setw(8) << etaMaxHigh <<  std::endl;
	  } 
	}
      }    
      roi_map.close();
    }
    return true;
  }
  
  
  bool TrigT1RPCRecRoiTool::etaDimLow (const TrigT1MuonRecRoiData& data,
				       double& etaMin, double& etaMax,
				       std::unique_ptr<const RpcCablingCondData> readCdo) const
  {
    // Get the strips delimiting the RoIs from rPCcablingSvc
    Identifier EtaLowBorder_id;
    Identifier EtaHighBorder_id;
    Identifier PhiLowBorder_id;
    Identifier PhiHighBorder_id;
    Amg::Vector3D EtaLowBorder_pos(0.,0.,0.);
    Amg::Vector3D EtaHighBorder_pos(0.,0.,0.);
    
    if(!readCdo->give_LowPt_borders_id(data.side(), data.sector(), data.roi(),
                                       EtaLowBorder_id, EtaHighBorder_id,
                                       PhiLowBorder_id, PhiHighBorder_id,
				       &m_idHelperSvc->rpcIdHelper())) return false;
    

    const MuonGM::MuonDetectorManager* muonMgr = m_muonMgr;
    if(m_useConditionData){
      SG::ReadCondHandle<MuonGM::MuonDetectorManager> DetectorManagerHandle{m_DetectorManagerKey};
      muonMgr = DetectorManagerHandle.cptr(); 
      if(muonMgr==nullptr){
	ATH_MSG_ERROR("Null pointer to the read MuonDetectorManager conditions object. Use the one from DetectorStore");
	muonMgr = m_muonMgr;
      }
    }
  
    const MuonGM::RpcReadoutElement* EtaLowBorder_descriptor =
      muonMgr->getRpcReadoutElement(EtaLowBorder_id);
    EtaLowBorder_pos = EtaLowBorder_descriptor->stripPos(EtaLowBorder_id);
    
    const MuonGM::RpcReadoutElement* EtaHighBorder_descriptor =
      muonMgr->getRpcReadoutElement(EtaHighBorder_id);
    EtaHighBorder_pos = EtaHighBorder_descriptor->stripPos(EtaHighBorder_id);
    
    etaMin=EtaLowBorder_pos.eta();
    etaMax=EtaHighBorder_pos.eta();
    if (etaMin>etaMax){
      double tmp=etaMin;
      etaMin=etaMax;
      etaMax=tmp;
    }
    return true;
  }
  
  bool TrigT1RPCRecRoiTool::etaDimHigh (const TrigT1MuonRecRoiData& data,
					double& etaMin, double& etaMax,
					std::unique_ptr<const RpcCablingCondData> readCdo) const
  {
    // Get the strips delimiting the RoIs from rPCcablingSvc
    Identifier EtaLowBorder_id;
    Identifier EtaHighBorder_id;
    Identifier PhiLowBorder_id;
    Identifier PhiHighBorder_id;
    Amg::Vector3D EtaLowBorder_pos(0.,0.,0.);
    Amg::Vector3D EtaHighBorder_pos(0.,0.,0.);
    
    if(!readCdo->give_HighPt_borders_id(data.side(), data.sector(), data.roi(),
					EtaLowBorder_id, EtaHighBorder_id,
					PhiLowBorder_id, PhiHighBorder_id,
					&m_idHelperSvc->rpcIdHelper())) return false;
    
    const MuonGM::MuonDetectorManager* muonMgr = m_muonMgr;
    if(m_useConditionData){
      SG::ReadCondHandle<MuonGM::MuonDetectorManager> DetectorManagerHandle{m_DetectorManagerKey};
      muonMgr = DetectorManagerHandle.cptr(); 
      if(muonMgr==nullptr){
	ATH_MSG_ERROR("Null pointer to the read MuonDetectorManager conditions object. Use the one from DetectorStore");
	muonMgr = m_muonMgr;
      }
    }

    const MuonGM::RpcReadoutElement* EtaLowBorder_descriptor =
      muonMgr->getRpcReadoutElement(EtaLowBorder_id);
    EtaLowBorder_pos = EtaLowBorder_descriptor->stripPos(EtaLowBorder_id);
    
    const MuonGM::RpcReadoutElement* EtaHighBorder_descriptor =
      muonMgr->getRpcReadoutElement(EtaHighBorder_id);
    EtaHighBorder_pos = EtaHighBorder_descriptor->stripPos(EtaHighBorder_id);
    
    etaMin=EtaLowBorder_pos.eta();
    etaMax=EtaHighBorder_pos.eta();
    if (etaMin>etaMax){
      double tmp=etaMin;
      etaMin=etaMax;
      etaMax=tmp;
    }
    return true;
  }
  
  
}
