/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#include "LArCalibTools/LArNoise2Ntuple.h"
#include "LArRawConditions/LArNoiseComplete.h"
#include "LArRawConditions/LArNoiseMC.h"
#include "CaloIdentifier/CaloGain.h"
/*
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/SmartDataPtr.h"
*/

//#include <fstream>

LArNoise2Ntuple::LArNoise2Ntuple(const std::string& name, ISvcLocator* pSvcLocator): 
  LArCond2NtupleBase(name, pSvcLocator) { 
  declareProperty("ContainerKey",m_contKey);
  //declareProperty("IsMC",m_isMC = false);

  m_ntTitle="Noise";
  m_ntpath="/NTUPLES/FILE1/NOISE";

}

LArNoise2Ntuple::~LArNoise2Ntuple() 
{}

StatusCode LArNoise2Ntuple::stop() {
  //const LArNoiseComplete* larNoiseComplete = NULL;
  //const LArNoiseMC* larNoiseMC = NULL;
  const ILArNoise* larNoise = NULL;
  StatusCode sc;
  sc=m_detStore->retrieve(larNoise,m_contKey);
  if (sc!=StatusCode::SUCCESS) {
     (*m_log)  << MSG::ERROR << "Unable to retrieve ILArNoise with key " 
               << m_contKey << " from DetectorStore" << endreq;
     return StatusCode::FAILURE;
  }

 NTuple::Item<long> cellIndex,gain;
 NTuple::Item<float> noise;

 sc=m_nt->addItem("icell",cellIndex,0,2000);
 if (sc!=StatusCode::SUCCESS)
   {(*m_log)  << MSG::ERROR << "addItem 'Cell Index' failed" << endreq;
    return StatusCode::FAILURE;
   }

 sc=m_nt->addItem("gain",gain,0,3);
 if (sc!=StatusCode::SUCCESS)
   {(*m_log) << MSG::ERROR << "addItem 'gain' failed" << endreq;
    return StatusCode::FAILURE;
   }


 sc=m_nt->addItem("noise",noise);
 if (sc!=StatusCode::SUCCESS)
   {(*m_log)  << MSG::ERROR << "addItem 'noise' failed" << endreq;
    return StatusCode::FAILURE;
   }


 unsigned cellCounter=0;
 for(long igain=CaloGain::LARHIGHGAIN; igain<CaloGain::LARNGAIN; igain++) {
   std::vector<HWIdentifier>::const_iterator itOnId = m_onlineId->channel_begin();
   std::vector<HWIdentifier>::const_iterator itOnIdEnd = m_onlineId->channel_end();
   for(; itOnId!=itOnIdEnd;++itOnId){
     const HWIdentifier hwid = *itOnId;
     if ( m_larCablingSvc->isOnlineConnected(hwid)) {
	 fillFromIdentifier(hwid);       
	 cellIndex = cellCounter;
         gain=igain;     
	 noise = larNoise->noise(hwid,igain);
	 sc=ntupleSvc()->writeRecord(m_nt);
	 if (sc!=StatusCode::SUCCESS) {
	   (*m_log)  << MSG::ERROR << "writeRecord failed" << endreq;
	   return StatusCode::FAILURE;
	 }
     }//end if isConnected
     cellCounter++;
  }//end loop over online ID
 } // ovr gains

 (*m_log)  << MSG::INFO << "LArNoise2Ntuple has finished." << endreq;
 return StatusCode::SUCCESS;
}// end finalize-method.
   
