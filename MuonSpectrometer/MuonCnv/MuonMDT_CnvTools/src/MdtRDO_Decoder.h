/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
*/

#ifndef MUONBYTESTREAMMDTRDODECODER_H
#define MUONBYTESTREAMMDTRDODECODER_H

#include "AthenaBaseComps/AthAlgTool.h"

#include "MuonRDO/MdtAmtHit.h"
#include "MuonDigitContainer/MdtDigit.h"
#include "MuonIdHelpers/MdtIdHelper.h"
#include "MuonCablingData/MuonMDT_CablingMap.h"
#include "StoreGate/ReadCondHandleKey.h"

#include "MuonMDT_CnvTools/IMDT_RDO_Decoder.h"

#include <inttypes.h>
#include <vector>


namespace Muon {
  // Decoder class for conversion from MDT RDOs to MDT digits
  // Stefano Rosati
  // CERN Jan 2004
  
  class MdtRDO_Decoder: public AthAlgTool, virtual public IMDT_RDO_Decoder {
    
  public:
    
    MdtRDO_Decoder( const std::string& type, const std::string& name,
        const IInterface* parent ) ;

    virtual StatusCode initialize();
    virtual StatusCode finalize();    
    
    MdtDigit * getDigit(const MdtAmtHit * amtHit, uint16_t& subdetId, 
			uint16_t& mrodId, uint16_t& csmId) const;
    
    Identifier getOfflineData(const MdtAmtHit * amtHit, uint16_t& subdetId, 
			      uint16_t& mrodId, uint16_t& csmId, int& tdc, int& width) const;
    
  private:
    
    const MdtIdHelper * m_mdtIdHelper;
    SG::ReadCondHandleKey<MuonMDT_CablingMap> m_readKey{this, "ReadKey", "MuonMDT_CablingMap", "Key of MuonMDT_CablingMap"};
    
  };
  
}

inline MdtDigit* Muon::MdtRDO_Decoder::getDigit(const MdtAmtHit* amtHit, uint16_t& subdetId, 
						uint16_t& mrodId, uint16_t& csmId) const
{
  int stationEta = 0;
  int stationPhi = 0;
  int stationName = 0;
  int multiLayer = 0;
  int tubeLayer  = 0;
  int tube = 0;
  
  uint16_t tdc     = amtHit->tdcId();
  uint16_t chan    = amtHit->channelId();
  uint16_t coarse  = amtHit->coarse();
  uint16_t fine    = amtHit->fine();
  int width   = (int) amtHit->width();
  
  SG::ReadCondHandle<MuonMDT_CablingMap> readHandle{m_readKey};
  const MuonMDT_CablingMap* readCdo{*readHandle};
  if(readCdo==nullptr){
    ATH_MSG_ERROR("Null pointer to the read conditions object");
    return NULL;
  }
  bool cab = readCdo->getOfflineId((uint8_t) subdetId, (uint8_t) mrodId, 
				   (uint8_t) csmId, (uint8_t) tdc, 
				   (uint8_t) chan,
				   stationName, stationEta, stationPhi,
				   multiLayer, tubeLayer, tube);
    
  if (!cab) return NULL;
  
  Identifier chanId = m_mdtIdHelper->channelID(stationName, stationEta, stationPhi, 
					       multiLayer, tubeLayer, tube);
  
  int tdcCounts = coarse*32+fine;
  
  MdtDigit* mdtDigit = new MdtDigit(chanId,tdcCounts,width,amtHit->isMasked());
  
  return mdtDigit;
}


inline Identifier Muon::MdtRDO_Decoder::getOfflineData(const MdtAmtHit* amtHit, uint16_t& subdetId, 
						       uint16_t& mrodId, uint16_t& csmId, int& tdcCounts, 
						       int& width) const
{
  int stationEta = 0;
  int stationPhi = 0;
  int stationName = 0;
  int multiLayer = 0;
  int tubeLayer  = 0;
  int tube = 0;
  
  uint16_t tdc     = amtHit->tdcId();
  uint16_t chan    = amtHit->channelId();
  uint16_t coarse  = amtHit->coarse();
  uint16_t fine    = amtHit->fine();
  width   = (int) amtHit->width();

  Identifier chanIdDefault;
  SG::ReadCondHandle<MuonMDT_CablingMap> readHandle{m_readKey};
  const MuonMDT_CablingMap* readCdo{*readHandle};
  if(readCdo==nullptr){
    ATH_MSG_ERROR("Null pointer to the read conditions object");
    return chanIdDefault;
  }
  bool cab = readCdo->getOfflineId(subdetId, mrodId, csmId, tdc, chan,
					stationName, stationEta, stationPhi,
					multiLayer, tubeLayer, tube);
  
  if (!cab) {
    
  }
  
  Identifier chanId = m_mdtIdHelper->channelID(stationName, stationEta, stationPhi, 
					       multiLayer, tubeLayer, tube);
  
  tdcCounts = coarse*32+fine;
  
  return chanId;
}

#endif
