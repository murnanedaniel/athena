/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#include "LArCalibTools/LArPhysWaves2Ntuple.h"
#include "CaloIdentifier/CaloGain.h"


LArPhysWaves2Ntuple::LArPhysWaves2Ntuple(const std::string& name, ISvcLocator* pSvcLocator): LArWaves2Ntuple(name, pSvcLocator)
{ 
  m_keylist.clear() ;
  declareProperty("KeyList",        m_keylist);
  declareProperty("NtupleName",     m_ntName); 
  declareProperty("AddCorrUndo",   m_addCorrUndo=true);
  declareProperty("ApplyCorrection",m_applyCorr=false);
}

LArPhysWaves2Ntuple::~LArPhysWaves2Ntuple()
{}


StatusCode LArPhysWaves2Ntuple::initialize() {
  m_ntTitle="Calibration Wave";
  m_ntpath=std::string("/NTUPLES/FILE1/")+m_ntName;
  return LArWaves2Ntuple::initialize();
}


StatusCode LArPhysWaves2Ntuple::stop() 
{ 
  
  StatusCode sc;

  sc=m_nt->addItem("gain",m_gain,0,3);
  if (sc!=StatusCode::SUCCESS) {
    (*m_log) << MSG::ERROR << "addItem 'gain' failed" << endreq;
    return StatusCode::FAILURE;
  }

  sc=m_nt->addItem("timeOffset",m_timeOffset,0.,100.);
  if (sc!=StatusCode::SUCCESS) {
    (*m_log) << MSG::ERROR << "addItem 'timeOffset' failed" << endreq;
    return StatusCode::FAILURE;
  }


  if (m_addCorrUndo) {
    sc=m_nt->addItem("corrUndo",m_corrUndo,0,1);
    if (sc!=StatusCode::SUCCESS) {
      (*m_log) << MSG::ERROR << "addItem 'corrUndo' failed" << endreq;
      return StatusCode::FAILURE;
    }
  }

  for ( unsigned k=0 ; k<m_keylist.size() ; k++ ) {
    const std::string& key = m_keylist[k] ;


    (*m_log) << MSG::INFO << "Processing WaveContainer from StoreGate! key=" << m_keylist[k] << endreq; 
    const LArPhysWaveContainer* physWaveContainer;
    LArPhysWaveContainer* physWaveContainer_nc;
    StatusCode sc;
    if (m_applyCorr) {
      sc = m_detStore->retrieve(physWaveContainer_nc,key);	    
      if (sc.isFailure()) {
	(*m_log) << MSG::ERROR << "Cannot read non-const LArPhysWaveContainer from StoreGate! key=" << key << endreq;
	(*m_log) << MSG::ERROR << "The ApplyCorrection option works only for non-const containers!" << endreq;
	return sc;
      }
      else
	(*m_log) << MSG::INFO << "Read non-const LArPhysWaveContainer from StoreGate! key= "  << key << endreq;

      if (!physWaveContainer_nc->correctionsApplied()) {
	sc=physWaveContainer_nc->applyCorrections();
	if (sc.isFailure()) {
	  (*m_log) << MSG::ERROR << "Failed to apply corrections to LArPhysWaveContainer!" << endreq;
	}
      }
      physWaveContainer=physWaveContainer_nc;
    }
    else {
      sc = m_detStore->retrieve(physWaveContainer,key);	    
      if (sc.isFailure()) {
	(*m_log) << MSG::ERROR << "Cannot read LArPhysWaveContainer from StoreGate! key=" << key << endreq;
	return StatusCode::FAILURE;
      } else 
	(*m_log) << MSG::INFO << "Read const LArPhysWaveContainer from StoreGate! key= "  << key << endreq;
    }
    
    for (unsigned igain=CaloGain::LARHIGHGAIN;igain<CaloGain::LARNGAIN ;++igain){
      std::vector<HWIdentifier>::const_iterator itOnId = m_onlineId->channel_begin();
      std::vector<HWIdentifier>::const_iterator itOnIdEnd = m_onlineId->channel_end();
      for(; itOnId!=itOnIdEnd;++itOnId){
	const HWIdentifier chid = *itOnId;
	m_gain=(long)igain;
	const LArPhysWave& wave=physWaveContainer->get(chid,igain);
	if (wave.isEmpty()) continue;
	if (m_addCorrUndo) m_corrUndo=0;
	fillWave(chid,wave); //Fill method from base-class
	m_timeOffset = wave.getTimeOffset();
	sc=ntupleSvc()->writeRecord(m_nt);
	if (sc!=StatusCode::SUCCESS) {
	  (*m_log) << MSG::ERROR << "writeRecord failed" << endreq;
	  return sc;
	}
      }//end loop over identifiers
    }//end loop over gains

    if (m_addCorrUndo) {
      for (unsigned igain=CaloGain::LARHIGHGAIN;igain<CaloGain::LARNGAIN ;++igain) {
	LArPhysWaveContainer::ConstCorrectionIt itUndo=physWaveContainer->undoCorrBegin(igain);
	LArPhysWaveContainer::ConstCorrectionIt itUndo_e=physWaveContainer->undoCorrEnd(igain);
	for(;itUndo!=itUndo_e;itUndo++) {
	  const HWIdentifier chid(itUndo->first);
	  const LArPhysWave& wave = itUndo->second;
	  m_gain  = (long)igain;
	  m_corrUndo = 1;
	  fillWave(chid,wave); //Fill method from base-class
	  sc=ntupleSvc()->writeRecord(m_nt);
	  if (sc!=StatusCode::SUCCESS) {
	    (*m_log) << MSG::ERROR << "writeRecord failed" << endreq;
	    return sc;
	  }
	}//end loop over corrections
      }//end loop over gain
    }//end if addCorrUndo
  }//end loop over container keys
  (*m_log) << MSG::INFO << "LArWave2Ntuple has finished." << endreq;
  return StatusCode::SUCCESS;
} // end finalize-method.

