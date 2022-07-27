/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef  TRIGL2MUONSA_MUCALSTREAMERTOOL_H
#define  TRIGL2MUONSA_MUCALSTREAMERTOOL_H

#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/ServiceHandle.h"

#include "MdtData.h"
#include "RpcData.h"
#include "TgcData.h"
#include "TrackData.h"
#include "MuCalCircClient.h"
#include "MuonCnvToolInterfaces/IMuonRawDataProviderTool.h"
#include "RPC_CondCabling/RpcCablingCondData.h"
#include "StoreGate/ReadCondHandleKey.h"
#include "IRegionSelector/IRegSelTool.h"
#include "MuonRDO/TgcRdoContainer.h"
#include "StoreGate/ReadHandleKey.h"
#include "xAODTrigger/MuonRoI.h"
#include "xAODEventInfo/EventInfo.h"
#include "MuonRDO/RpcPadContainer.h"

#include <fstream>
#include <string>
#include <vector>

// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------

namespace LVL1 {
  class RecMuonRoI;
}

namespace LVL2_MUON_CALIBRATION {
  class MdtCalibFragment;
  class RpcCalibFragment;
  class TgcCalibFragment;
}

class RpcPad;
class ROBDataProviderSvc;

inline unsigned bcBitmap(uint16_t bcTag) { 
  return (bcTag == 1 ? 4 : 
          ((bcTag == 0 || bcTag == 2) ? 2 :
           (bcTag == 3 ? 1 : 0)));
}
 
namespace TrigL2MuonSA {
  
  class MuCalStreamerTool: public AthAlgTool
  {
  public:
    
    MuCalStreamerTool(const std::string& type, 
		      const std::string& name,
		      const IInterface*  parent);

    virtual StatusCode initialize() override;

    std::string instanceName() const {return m_algInstanceName;}
    void setInstanceName(std::string name) { m_algInstanceName = name; }

    // set the properties
    void setBufferName(std::string buffName) {m_calBufferName=buffName;}

    std::vector<int>* getLocalBuffer()   {return &m_localBuffer;}
    int getLocalBufferSize() const {return m_localBuffer.size();}
    /* void clearLocalBuffer(); */

    //
    // initialize the stream
    StatusCode openStream(int calBufferSize);

    //
    // finalize and close the stream
    StatusCode closeStream();

    bool isStreamOpen();
    
    //
    // create the fragment corresponding to an roi
    StatusCode createRoiFragment(const xAOD::MuonRoI* roi, 
				 TrigL2MuonSA::TrackPattern& trackPattern, 
				 TrigL2MuonSA::MdtHits& mdtHits, 
				 TrigL2MuonSA::RpcHits& rpcHits, 
				 TrigL2MuonSA::TgcHits& tgcHits, 
				 //int calBufferSize, 
				 bool doDataScouting,
				 bool& updateTriggerElement, const EventContext&                        ctx) const; 

  private:

    Gaudi::Property< bool > m_writeToFile { this, "WriteToFile", false, "" };
    SG::ReadHandleKey<TgcRdoContainer> m_tgcRdoKey{this, "TgcRdoContainer", "TGCRDO", "Tgc RDO Input"};
    SG::ReadHandleKey<xAOD::EventInfo> m_eventInfoKey{ this, "EventInfoKey", "EventInfo", "" };
    SG::ReadHandleKey<RpcPadContainer> m_rpcPadKey{this, "RpcPadContainerKey", "RPCPAD", "RpcPad container fro MuCalStreamerTool"};

    // name of the calibration buffer or of the 
    // output file
    std::string m_calBufferName;
    std::string m_algInstanceName;

    // output file 
    std::ofstream m_outputFile;



    // the region selector

    ToolHandle<IRegSelTool> m_regSel_MDT{this, "RegSel_MDT", "RegSelTool/RegSelTool_MDT", "MDT Region Selector Tool"};
    ToolHandle<IRegSelTool> m_regSel_TGC{this, "RegSel_TGC", "RegSelTool/RegSelTool_TGC", "TGC Region Selector Tool"};
    


    SG::ReadCondHandleKey<RpcCablingCondData> m_readKey{this, "ReadKey", "RpcCablingCondData", "Key of RpcCablingCondData"};



    // local buffer for the TrigComposite object
    int m_localBufferSize = 0;
    std::vector<int> m_localBuffer;

    //
    // create the MDT fragment
    StatusCode createMdtFragment(TrigL2MuonSA::MdtHits& mdtHits,  
				 LVL2_MUON_CALIBRATION::MdtCalibFragment& mdtFragment, float phi) const;

    //
    // create the RPC fragment
    StatusCode createRpcFragment(const xAOD::MuonRoI* roi,
				 LVL2_MUON_CALIBRATION::RpcCalibFragment& rpcFragment, const EventContext&                        ctx) const;

    //
    // create the TGC fragment
    StatusCode createTgcFragment(std::vector<uint32_t>& tgcRobIdList,
				 LVL2_MUON_CALIBRATION::TgcCalibFragment& tgcFragment) const;



    TrigL2MuonSA::MuCalCircClient *m_circ = nullptr;

  };
  
} // namespace TrigL2MuonSA

#endif  // 
