/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#ifndef JIVEXML_RPCPREPDATARETRIEVER_H
#define JIVEXML_RPCPREPDATARETRIEVER_H

#include "JiveXML/IDataRetriever.h"

#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"

#include "MuonIdHelpers/IMuonIdHelperSvc.h"
#include "MuonPrepRawData/RpcPrepDataContainer.h"

#include <string>

namespace JiveXML {
  
  class RpcPrepDataRetriever : virtual public IDataRetriever, public AthAlgTool {
    
  public:

    /// Standard Constructor
    RpcPrepDataRetriever(const std::string& type, const std::string& name, const IInterface* parent);

     /// Retrieve all the data
    virtual StatusCode retrieve(ToolHandle<IFormatTool> &FormatTool);

    /// Return the name of the data type
    virtual std::string dataTypeName() const { return m_typeName; };

    ///Default AthAlgTool methods
    StatusCode initialize();

  private:
    
   ///The data type that is generated by this retriever
    const std::string m_typeName;

    ///The storegate key for the RPC collection
    SG::ReadHandleKey<Muon::RpcPrepDataContainer> m_sgKey{this, "StoreGateKey", "RPC_Measurements", "Name of the RpcPrepDataContainer"};

    ServiceHandle<Muon::IMuonIdHelperSvc> m_idHelperSvc {this, "MuonIdHelperSvc", "Muon::MuonIdHelperSvc/MuonIdHelperSvc"};

  };
  
}
#endif
