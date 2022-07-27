/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "LArEventTest/LArDigitsToNtuple.h"
#include "CaloIdentifier/CaloGain.h"
#include "CaloIdentifier/CaloCell_ID.h"
#include "TBEvent/TBPhase.h"
#include "TBEvent/TBTriggerPatternUnit.h"
#include "TBEvent/TBScintillatorCont.h"
#include <fstream>


LArDigitsToNtuple::LArDigitsToNtuple(const std::string& name, ISvcLocator* pSvcLocator)
  : AthAlgorithm(name, pSvcLocator),
    m_onlineHelper(0),
    m_emId(0),
    m_hecId(0),
    m_fcalId(0),
    m_nt(0),
    m_nsamples(0),
    m_accept(0),
    m_reject(0),
    m_phase(0),
    m_scint(0),
    m_trigger(0),
    m_sca(0),
    m_ped(0)
{ 
  declareProperty("NSamples",m_nsamples=5);
  // By default to get only pedestal events, we need to set m_reject=1
  // But to reject non flaged events (physics events occuring out of
  // the beam signal...) we need also to set m_accept=14.
  declareProperty("accept", m_accept=14);
  declareProperty("reject", m_reject=1);
  declareProperty("ReadPhase", m_phase=1);
  declareProperty("ReadScint", m_scint=1);
  declareProperty("ReadTrigger", m_trigger=1);
  declareProperty("ReadSCA", m_sca=1);
  declareProperty("ReadPedestal", m_ped=1);
}

LArDigitsToNtuple::~LArDigitsToNtuple() 
{}

StatusCode LArDigitsToNtuple::initialize()
{
  const CaloCell_ID* idHelper = nullptr;
  ATH_CHECK( detStore()->retrieve (idHelper, "CaloCell_ID") );
  m_emId=idHelper->em_idHelper();
  m_fcalId=idHelper->fcal_idHelper();
  m_hecId=idHelper->hec_idHelper();

  if (!m_emId) {
    ATH_MSG_ERROR ( "Could not access lar EM ID helper" );
    return StatusCode::FAILURE;
  }
  if (!m_fcalId) {
    ATH_MSG_ERROR ( "Could not access lar FCAL ID helper" );
    return StatusCode::FAILURE;
  }
  if (!m_hecId) {
    ATH_MSG_ERROR ( "Could not access lar HEC ID helper" );
    return StatusCode::FAILURE;
  }

  ATH_CHECK( m_cablingKey.initialize() );

  ATH_CHECK( m_pedestalKey.initialize(m_ped) );
  ATH_CHECK( m_contKey.initialize() );
  ATH_CHECK( m_hdrContKey.initialize(m_sca) );

  ATH_CHECK( detStore()->retrieve(m_onlineHelper, "LArOnlineID") );
  ATH_MSG_DEBUG ( " Found the LArOnlineID helper. " );

  // Book ntuple
  NTupleFilePtr file1(ntupleSvc(),"/NTUPLES/FILE1");
  if (!file1) {
    ATH_MSG_ERROR ( "Booking of NTuple failed" );
    return StatusCode::FAILURE;
  }
  NTuplePtr nt(ntupleSvc(),"/NTUPLES/FILE1/DIGITS");
  if (!nt) {
    nt=ntupleSvc()->book("/NTUPLES/FILE1/DIGITS",CLID_ColumnWiseTuple,"Digits");
  }
  if (!nt) {
    ATH_MSG_ERROR ( "Booking of NTuple failed" );
    return StatusCode::FAILURE;
  }

  ATH_CHECK( nt->addItem("event",m_nt_event) );
  ATH_CHECK( nt->addItem("layer",m_nt_layer,0,4) );
  ATH_CHECK( nt->addItem("ieta",m_nt_eta,0,510) );
  ATH_CHECK( nt->addItem("iphi",m_nt_phi,0,1023) );
  ATH_CHECK( nt->addItem("region",m_nt_region,0,1) );
  ATH_CHECK( nt->addItem("barrel_ec",m_nt_barrel_ec,0,1) );
  ATH_CHECK( nt->addItem("pos_neg",m_nt_pos_neg,0,1) );
  ATH_CHECK( nt->addItem("detector",m_nt_detector,0,2) );
  ATH_CHECK( nt->addItem("FT",m_nt_FT,0,32) );
  ATH_CHECK( nt->addItem("slot",m_nt_slot,1,15) );
  ATH_CHECK( nt->addItem("channel",m_nt_channel,0,127) );
  ATH_CHECK( nt->addItem("gain",m_nt_gain,0,3) );
  ATH_CHECK( nt->addItem("samples",m_nsamples,m_nt_samples) );
  ATH_CHECK( nt->addItem("onlId",m_nt_onlId) );

  if(m_ped)     ATH_CHECK( nt->addItem("ped",m_nt_ped) );
  if(m_sca)     ATH_CHECK( nt->addItem("sca",m_nsamples,m_nt_sca) );
  if(m_phase)   ATH_CHECK( nt->addItem("tdc",m_nt_tdc) );
  if(m_trigger) ATH_CHECK( nt->addItem("trigger",m_nt_trigger) );
  if(m_scint)   ATH_CHECK( nt->addItem("S1",m_nt_S1) );

  m_nt=nt;

  return StatusCode::SUCCESS;
}  

StatusCode LArDigitsToNtuple::execute()
{
  int triggerword;
  double S1Adc,tdcphase;
  const EventContext& ctx = Gaudi::Hive::currentContext();
  const int eventnumber = ctx.eventID().event_number();

  

  // Retrieve the TBScintillators
  if(m_scint) {
    S1Adc=-999.0;
    const TBScintillatorCont * theTBScint;
    StatusCode sc =evtStore()->retrieve(theTBScint,"ScintillatorCont");
    if (sc.isFailure()) 
      {
	ATH_MSG_ERROR ( " Cannot read TBScintillatorCont from StoreGate! " );
	//return StatusCode::FAILURE;
      }
    else {
      for (const TBScintillator* scint : *theTBScint) {
	const std::string name = scint->getDetectorName();
	if (name=="S1") {
	  S1Adc = scint->getSignal();
	  break;
	}
      } //end loop over scintillator-container
    }
  }
  else S1Adc=0; 

  //Retrieve the TBPhase
  if(m_phase) {
    const TBPhase* theTBPhase=nullptr;
    StatusCode sc = evtStore()->retrieve(theTBPhase, "TBPhase");
    
    if (sc.isFailure()) {
      tdcphase = -999.0;
      ATH_MSG_ERROR( "cannot allocate TBPhase " );
      //return StatusCode::FAILURE;
    } else {
      tdcphase = theTBPhase->getPhase();
    }
  }
  else tdcphase = 0;

  //Retrieve the TBTriggerPatternUnit
  if(m_trigger) {
    const TBTriggerPatternUnit* theTBTriggerPatternUnit;
    StatusCode sc = evtStore()->retrieve(theTBTriggerPatternUnit, "TBTrigPat");

    if (sc.isFailure()) {
      triggerword = 0;
      ATH_MSG_ERROR( "cannot allocate TBTriggerPatternUnit" );
      //return StatusCode::FAILURE;
    } else {
      triggerword = theTBTriggerPatternUnit->getTriggerWord();
    }

    // As a trigger can be both physics and pedestal (Indeed!)
    // we need to make two checks to cover all possibilities
    // Check whether we should reject this trigger
    if(m_reject>0) if(triggerword&m_reject) return StatusCode::SUCCESS;
    // Check whether we should accept this trigger
    if(m_accept>0) if((triggerword&m_accept)==0) return StatusCode::SUCCESS;
  }
  else triggerword=0;

  // Retrieve the pedestal
  //Pointer to conditions data objects 
  const ILArPedestal* larPedestal=nullptr;
  if (m_ped) {
    SG::ReadCondHandle<ILArPedestal> pedHdl(m_pedestalKey,ctx);
    larPedestal= *pedHdl;
  }

  // Retrieve LArDigitContainer
  SG::ReadHandle<LArDigitContainer> digit_cont{m_contKey,ctx};
  ATH_MSG_INFO ( "Retrieved LArDigitContainer from StoreGate! key=" << m_contKey.key() );

  // Retrieve LArFebHeaderContainer
  const LArFebHeaderContainer *larFebHeaderContainer=nullptr;
  if(m_sca) {
    SG::ReadHandle<LArFebHeaderContainer> febHdrHdl{m_hdrContKey,ctx};
    larFebHeaderContainer=&(*febHdrHdl);
  }
  SG::ReadCondHandle<LArOnOffIdMapping> larCablingHdl{m_cablingKey,ctx};
  const LArOnOffIdMapping* cabling=*larCablingHdl;
  if(!cabling) {
     ATH_MSG_ERROR("Could not get LArOnOffIdMapping !!");
     return StatusCode::FAILURE;
  }

  // Fill ntuple
  LArDigitContainer::const_iterator it = digit_cont->begin(); 
  LArDigitContainer::const_iterator it_e = digit_cont->end(); 
  if(it==it_e) {
    ATH_MSG_DEBUG ( "LArDigitContainer is empty..." );
  }
  for(; it!=it_e; ++it){
    const HWIdentifier hwid=(*it)->channelID();//hardwareID();
    // Fill detector geometry information
    m_nt_event   = eventnumber;
    m_nt_onlId   = hwid.get_identifier32().get_compact();
    if(m_phase)   m_nt_tdc     = tdcphase;
    if(m_trigger) m_nt_trigger = triggerword;
    if(m_scint)   m_nt_S1      = S1Adc;
    try {
      Identifier id=cabling->cnvToIdentifier(hwid);
      if (m_emId->is_lar_em(id)) {
	m_nt_eta       = m_emId->eta(id);
	m_nt_phi       = m_emId->phi(id);
	m_nt_layer     = m_emId->sampling(id);
	m_nt_region    = m_emId->region(id);
	m_nt_detector  = 0;
      }
      else if (m_hecId->is_lar_hec(id)) {
	m_nt_eta       = m_hecId->eta(id);
	m_nt_phi       = m_hecId->phi(id);
	m_nt_layer     = m_hecId->sampling(id);
	m_nt_region    = m_hecId->region(id);
	m_nt_detector  = 1;
      }
      else if (m_fcalId->is_lar_fcal(id)) {
	m_nt_eta       = m_fcalId->eta(id);
	m_nt_phi       = m_fcalId->phi(id);
	m_nt_layer     = m_fcalId->module(id);
	m_nt_region    = 0;
	m_nt_detector  = 2;
      }
    }
    catch (LArID_Exception & except) {
      m_nt_eta       = -1;
      m_nt_phi       = -1;
      m_nt_layer     = -1;
      m_nt_region    = -1;
      m_nt_detector  = -1;
    }
    // Fill hardware information
    m_nt_barrel_ec = m_onlineHelper->barrel_ec(hwid);
    m_nt_pos_neg   = m_onlineHelper->pos_neg(hwid);
    m_nt_FT        = m_onlineHelper->feedthrough(hwid);
    m_nt_slot      = m_onlineHelper->slot(hwid);
    m_nt_channel   = m_onlineHelper->channel(hwid);

    // Fill pedestal information
    float thePedestal=-1;    
    if (larPedestal) {
      float  DBpedestal=larPedestal->pedestal(hwid,(*it)->gain());
      if (DBpedestal >= (1.0+LArElecCalib::ERRORCODE))
	thePedestal=DBpedestal;
      if (thePedestal<0) {
	thePedestal = -999;
	ATH_MSG_DEBUG ( "No pedestal found for this cell. Use default value " << thePedestal );
      }
      m_nt_ped = thePedestal;
    }
  
    // Fill raw data samples and gain
    for(unsigned int i=0;i<(*it)->samples().size();i++) {
      if((int)i>=m_nsamples) break;
      m_nt_samples[i]=(*it)->samples()[i];
    }
    m_nt_gain=(*it)->gain();
  
    // Fill SCA numbers
    if(m_sca) {
      const HWIdentifier febid=m_onlineHelper->feb_Id(hwid);
      for (const LArFebHeader* feb : *larFebHeaderContainer) {
	const HWIdentifier this_febid=feb->FEBId();
	
	if(this_febid!=febid) continue;
	for(unsigned int i=0;i<feb->SCA().size();i++) {
	  if((int)i>=m_nsamples) break;
	  m_nt_sca[i]=feb->SCA()[i];
	}
	break;
      } // End FebHeader loop
    }
    ATH_CHECK(ntupleSvc()->writeRecord(m_nt));
 
  }

  return StatusCode::SUCCESS;
}

StatusCode LArDigitsToNtuple::finalize()
{
  ATH_MSG_INFO ( "LArDigitsToNtuple has finished." );
  return StatusCode::SUCCESS;
}// end finalize-method.
