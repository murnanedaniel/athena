/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#ifndef JIVEXML_TRIGRPCDATARETRIEVER_H
#define JIVEXML_TRIGRPCDATARETRIEVER_H

#include <string>

#include "JiveXML/IDataRetriever.h"

#include "AthenaBaseComps/AthAlgTool.h"

#include "MuonReadoutGeometry/MuonDetectorManager.h"
#include "MuonReadoutGeometry/RpcReadoutElement.h"
#include "MuonRPC_CnvTools/IRPC_RDO_Decoder.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"

#include "GaudiKernel/ServiceHandle.h"
#include "MuonIdHelpers/IMuonIdHelperSvc.h"

#include "RPC_CondCabling/RpcCablingCondData.h"
#include "StoreGate/ReadCondHandleKey.h"

class IRPCcablingSvc;

namespace JiveXML {

  class TrigRpcDataRetriever : virtual public IDataRetriever, public AthAlgTool {
    
  public:

    /// Standard Constructor
    TrigRpcDataRetriever(const std::string& type,const std::string& name,const IInterface* parent);

      /// Retrieve all the data
    virtual StatusCode retrieve(ToolHandle<IFormatTool> &FormatTool);

    /// Return the name of the data type
    virtual std::string dataTypeName() const { return m_typeName; };

    ///Default AthAlgTool methods
    StatusCode initialize();
   
  private:
    
   ///The data type that is generated by this retriever
    const std::string m_typeName;

    ///The storegate key for the CSC collection
    std::string m_sgKey;

    ServiceHandle<Muon::IMuonIdHelperSvc> m_idHelperSvc {this, "MuonIdHelperSvc", "Muon::MuonIdHelperSvc/MuonIdHelperSvc"};

    const IRPCcablingSvc *m_rpcCabling;

    SG::ReadCondHandleKey<RpcCablingCondData> m_readKey{this, "ReadKey", "RpcCablingCondData", "Key of RpcCablingCondData"};

    /// RPC decoder
    ToolHandle<Muon::IRPC_RDO_Decoder> m_rpcDecoder{this,"RpcRDO_Decoder","Muon::RpcRDO_Decoder"}; 
    
    ///Geo Model
    SG::ReadCondHandleKey<MuonGM::MuonDetectorManager> m_DetectorManagerKey {this, "DetectorManagerKey", 
	"MuonDetectorManager", 
	"Key of input MuonDetectorManager condition data"};    
  };
  
}
#endif
