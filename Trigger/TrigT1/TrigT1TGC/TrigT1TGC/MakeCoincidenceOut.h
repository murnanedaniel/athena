/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// MakeCoincidenceOut.h
#ifndef __MAKECOINCIDENCEOUT_H__
#define __MAKECOINCIDENCEOUT_H__

// Gaudi includes
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/Property.h"
#include "StoreGate/StoreGateSvc.h"
#include "GaudiKernel/NTuple.h"

class TgcIdHelper;

namespace LVL1TGCTrigger {

  /** 
   *
   *   @short CWN implementation for TGC digits and coincidence(@ LowPt and HightPt) !Obsolete!
   *   
   *   @author Naoko Kanaya <Naoko.Kanaya@cern.ch>
   *   @author Hisaya Kurashige <Hisaya.Kurashige@cern.ch>
   *
   */

  class MakeCoincidenceOut : public Algorithm
  {
    
  public:
    
    MakeCoincidenceOut( const std::string& name, ISvcLocator* pSvcLocator ) ;
    
    ~MakeCoincidenceOut();
    
    StatusCode initialize() ;
    StatusCode execute() ;
    StatusCode finalize() ;
    StatusCode bookHistos() ;
    
  private:
    ServiceHandle<StoreGateSvc> m_sgSvc;
    StringProperty  m_key ;
    BooleanProperty m_WriteMCtruth;
    const TgcIdHelper* m_tgcIdHelper;
    NTuple::Tuple* m_ntuplePtr;
    
    NTuple::Item<long>   m_runNumber, m_eventNumber;
    NTuple::Item<long>   m_ncol, m_ndig;
    NTuple::Array<long>  m_stationName, m_stationEta, m_stationPhi, m_gasGap, m_isStrip, m_channel;
    NTuple::Item<long>   m_nhpt;
    NTuple::Array<long>  m_hbid, m_hsec, m_hmod, m_hreg, m_hssc, m_hssr, m_hssphi, m_hssid, m_hptr, m_hdr, m_hptphi, m_hdphi, m_hptmax;
    // MC
    NTuple::Item<long>   m_nmuMC;
    NTuple::Array<long>  m_idMC;
    NTuple::Array<double>  m_ptMC, m_etaMC, m_phiMC;
    
  }; // class MakeCoincidenceOut
 
} // namespace LVL1TGCTrigger

#endif // __MAKECOINCIDENCEOUT_H__
