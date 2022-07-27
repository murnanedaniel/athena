/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef  TRIGL2MUONSA_RPCCLUSTERPREPARATOR_H
#define  TRIGL2MUONSA_RPCCLUSTERPREPARATOR_H


#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"

#include "TrigT1Interfaces/RecMuonRoI.h"
#include "TrigT1Interfaces/ITrigT1MuonRecRoiTool.h"
#include "TrigSteeringEvent/TrigRoiDescriptor.h"

#include "MuonIdHelpers/IMuonIdHelperSvc.h"
#include "MuonPrepRawData/MuonPrepDataContainer.h"
#include "MuonPrepRawData/RpcPrepDataContainer.h"
#include "MuonPrepRawData/RpcPrepDataCollection.h"

#include "ClusterPatFinder.h"


// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------

namespace TrigL2MuonSA {

class RpcClusterPreparator: public AthAlgTool
{
  public:

    RpcClusterPreparator(const std::string& type,
          		 const std::string& name,
			 const IInterface*  parent);

    virtual StatusCode initialize() override;

    StatusCode clusteringRPCs(const bool doMultiMuon,
                              const std::vector<const Muon::RpcPrepDataCollection*>& rpcCols,
                              const TrigRoiDescriptor*          p_roids,
                              const ToolHandle<ClusterPatFinder>*     clusterPatFinder,
                              TrigL2MuonSA::RpcLayerClusters&   rpcLayerClusters) const;

    typedef std::map<int, const Muon::RpcPrepData*, std::less<int> > pattern;
    int buildPatterns(const bool doMultiMuon,
                      const TrigRoiDescriptor*        p_roids,
                      const Muon::RpcPrepDataCollection* rpcCollection,
                      std::map<Identifier, pattern>&     digits) const;

    void buildClusters(const ToolHandle<ClusterPatFinder>*   clusterPatFinder,
                       std::map<Identifier, pattern>&  digits,
                       TrigL2MuonSA::RpcLayerClusters& rpcLayerClusters) const;
 private:

    ServiceHandle<Muon::IMuonIdHelperSvc> m_idHelperSvc {this, "MuonIdHelperSvc", "Muon::MuonIdHelperSvc/MuonIdHelperSvc"};
    ToolHandle<LVL1::ITrigT1MuonRecRoiTool> m_recRPCRoiTool{ this, "TrigT1RPCRecRoiTool", "LVL1::TrigT1RPCRecRoiTool/TrigT1RPCRecRoiTool"};

};

} // namespace TrigL2MuonSA

#endif  // TRIGL2MUONSA_RPCCLUSTERPREPARATOR_H
