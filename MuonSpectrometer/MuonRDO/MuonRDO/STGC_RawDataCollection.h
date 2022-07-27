/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef MUONRDO_STGC_RawDataCOLLECTION_H
#define MUONRDO_STGC_RawDataCOLLECTION_H

#include "MuonRDO/STGC_RawData.h"
#include "AthContainers/DataVector.h"
#include "Identifier/IdentifierHash.h"

namespace Muon {
  class STGC_RawDataContainerCnv_p1;
  class STGC_RawDataContainerCnv_p2;
  class STGC_RawDataContainerCnv_p3;
    
  class STGC_RawDataCollection : public DataVector<STGC_RawData>
  {
  public:
    friend class Muon::STGC_RawDataContainerCnv_p1;
    friend class Muon::STGC_RawDataContainerCnv_p2;
    friend class Muon::STGC_RawDataContainerCnv_p3;
    
    STGC_RawDataCollection(IdentifierHash hash) : m_idHash(hash) {}

    const IdentifierHash& identifyHash() const { return m_idHash; }
  private:
  
    /** Offline IdentifierHash for this collection*/
    IdentifierHash m_idHash;
  };
}

#endif
