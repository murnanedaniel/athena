/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#ifndef MUONCSC_CNVTOOLS_CSCRDOTOCSCPREPDATATOOLMT_H
#define MUONCSC_CNVTOOLS_CSCRDOTOCSCPREPDATATOOLMT_H 

#include "CscRdoToCscPrepDataToolCore.h"

#include <string>

////////////////////////////////////////////////////////////////////////////////////////
/// Author: Ketevi A. Assamagan
/// BNL, April 03, 2005
///
/// algorithm to decode RDO into CscStripPrepData
/// get the RDO container from Storegate
/// loop over the RDO
/// Decode RDO into PrepRawData
/// loop over the PrepRawData and build the PrepRawData container
/// store the PrepRawData container in StoreGate
////////////////////////////////////////////////////////////////////////////////////////

namespace Muon {
  class CscRdoToCscPrepDataToolMT : public CscRdoToCscPrepDataToolCore {

  public:
    
    CscRdoToCscPrepDataToolMT(const std::string& type, const std::string& name,
                            const IInterface* parent);
    
    /** destructor 
     */ 
    virtual ~CscRdoToCscPrepDataToolMT();
    
    virtual StatusCode initialize() override;
    
    virtual StatusCode finalize() override;

    using CscRdoToCscPrepDataToolCore::decode;
    
    virtual StatusCode decode(std::vector<IdentifierHash>& givenIdhs, std::vector<IdentifierHash>& decodedIdhs) override;
  };
}
#endif /// MUONCSC_CNVTOOL_CSCRDOTOCSCPREPDATA_H


