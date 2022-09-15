/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

//==================================================================================
// Include files...
//==================================================================================

// This files header
#include "InDetPerformanceMonitoring/ZmumuEvent.h"

// Standard headers

// Package Headers
#include "InDetPerformanceMonitoring/PerfMonServices.h"

// ATLAS headers
#include "StoreGate/StoreGateSvc.h"

//#include "muonEvent/MuonParamDefs.h"

#include "CLHEP/Random/RandFlat.h"

#include "xAODMuon/Muon.h"
#include "xAODMuon/MuonContainer.h"

#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Units/PhysicalConstants.h"
using CLHEP::GeV;

//==================================================================================
// Public Methods
//==================================================================================

ZmumuEvent::ZmumuEvent()
{
  m_xSampleName = "ZMM";

  m_container = PerfMonServices::MUID_COLLECTION;

  m_doDebug = false;

  // Setup the muon tags
  m_uMuonTags   = 2;
  m_uTrackMatch = 0;
  m_bLooseMatch = true;  // will use combined fit otherwise.
  m_etaCut      = 1.05;
  m_LeadingMuonPtCut = 20.;   
  m_SecondMuonPtCut = 15.;
  m_MassWindowLow = 60.0;
  m_MassWindowHigh = 120.0;
  m_OpeningAngleCut = 0.2; // in radians
  m_Z0GapCut = 5.0; // in mm
  m_SelectMuonByIso = true;
  m_SelectMuonByIP = true;
  m_analyzedEventCount = 0;
  m_eventsWithoutEnoughMuonsCount = 0;
  m_eventsWithEnoughMuonsCount = 0;
  m_acceptedEventCount = 0;
  m_testedMuonCount = 0;
  m_acceptedMuonCount = 0;
  m_eventselectioncount_toofewmuons = 0;
  m_eventselectioncount_notallmuonsfilled = 0;
  m_eventselectioncount_morethantwomuons = 0;
  m_eventselectioncount_ptofleadingmuon = 0;
  m_eventselectioncount_ptofsecondmuon = 0;
  m_eventselectioncount_masswindow = 0;
  m_eventselectioncount_openingangle = 0;
  m_eventselectioncount_dimuoncharge = 0;
  m_skipMScheck = false;
}

//==================================================================================
ZmumuEvent::~ZmumuEvent()
{
}

//==================================================================================
void ZmumuEvent::Init()
{
  

  m_xMuonID.Init();
  
  PARENT::Init();
}



//==================================================================================
const std::string ZmumuEvent::getRegion() const{

  const double eta1 = std::abs(m_pxRecMuon[m_muon1]->eta());
  const double eta2 = std::abs(m_pxRecMuon[m_muon2]->eta());

  if ( eta1 < m_etaCut && eta2 < m_etaCut )
    return "BB";

  else if( (eta1 < m_etaCut && eta2 > m_etaCut) || (eta1 > m_etaCut && eta2 < m_etaCut) )
    return "BE";

  else return "EE";
}


//==================================================================================
bool ZmumuEvent::Reco()
{
  if (m_doDebug) { std::cout << " * ZmumuEvent * ZmumuEvent::Reco() starting " << std::endl; }
  m_analyzedEventCount++;

  // Clear out the previous events record.
  this->Clear();

  const xAOD::MuonContainer* pxMuonContainer = PerfMonServices::getContainer<xAOD::MuonContainer>( m_container );

  // START patch by Anthony to avoid crash when MuonSpetrometer::Pt was not defined mainly for data16 
  // WARNING thus is necessary for data16 !!
  if (false) {
    for( auto muon :  *pxMuonContainer ){
      const xAOD::TrackParticle* idtrk(nullptr);
      const xAOD::TrackParticle* metrk(nullptr);
      idtrk = muon->trackParticle(xAOD::Muon::InnerDetectorTrackParticle);
      metrk = muon->trackParticle(xAOD::Muon::ExtrapolatedMuonSpectrometerTrackParticle);
      if (idtrk && metrk) {
	muon->auxdecor<float>("InnerDetectorPt") = idtrk->pt();      
	muon->auxdecor<float>("MuonSpectrometerPt") = metrk->pt();
      }
    }
  }
  // END patch 

  if (pxMuonContainer != nullptr) {
    if (m_doDebug) {std::cout << " * ZmumuEvent * track list has "<< pxMuonContainer->size() << " muon in xAOD::MuonContainer " << m_container <<std::endl; }
    xAOD::MuonContainer::const_iterator xMuonItr  = pxMuonContainer->begin();
    xAOD::MuonContainer::const_iterator xMuonItrE  = pxMuonContainer->end();

    int attemptedMuonCount = 0;
    int acceptedMuonCount = 0;

    while ( xMuonItr != xMuonItrE ){ // this loops on the muons in the pxMuonContainer
      const xAOD::Muon* pxCMuon = *xMuonItr;
      attemptedMuonCount++;
      m_testedMuonCount++;
      if(m_doDebug){std::cout << " * ZmumuEvent * Reco() ** attempt on xMuonItr number "<< attemptedMuonCount << " (pointer: "<< *xMuonItr <<")" << std::endl; }
      // Apply muon cuts
      if ( m_xMuonID.passSelection( pxCMuon) ) {
	if (RecordMuon( pxCMuon )) {
	  acceptedMuonCount++;
	  m_acceptedMuonCount++;
	  if (m_doDebug) {std::cout << "                          This muon is accepeted !! this is muon number " << acceptedMuonCount << " & full pass" << m_numberOfFullPassMuons << std::endl; }
	}
      }
      ++xMuonItr;
    } // end loop on muons
    if (m_doDebug) {std::cout << " * ZmumuEvent * accepted " << acceptedMuonCount << " muons from the input list of "<< pxMuonContainer->size() <<std::endl; }

    if (acceptedMuonCount < 2) m_eventsWithoutEnoughMuonsCount++;
    if (acceptedMuonCount >= 2) m_eventsWithEnoughMuonsCount++;
  } // muon container exist
  else {
    std::cout << " * ZmumuEvent * Can't retrieve combined muon collection (container: " << m_container <<") " << std::endl;
    return false;
  }
  
  // Ordering of muons
  this->OrderMuonList();
  
  // Reconstruct the invariant mass ( based on mu-sys pt ).
  ReconstructKinematics();
  
  m_passedSelectionCuts = EventSelection(ID);
  m_DiMuonPairInvMass =  m_fInvariantMass[ID];
  
  if (m_passedSelectionCuts) m_acceptedEventCount++;
  
  if(m_doDebug) {
    if ( m_passedSelectionCuts) std::cout << " * ZmumuEvent * Reco() * result of analyzed event " << m_analyzedEventCount << " --> Selected event :) " << std::endl;
    if (!m_passedSelectionCuts) std::cout << " * ZmumuEvent * Reco() * result of analyzed event " << m_analyzedEventCount << " --> Rejected event :( " << std::endl;
  }
  if (m_doDebug) {
    std::cout << " * ZmumuEvent::Reco * COMPLETED * Event has " << m_numberOfFullPassMuons 
	      << " muons. " << m_acceptedEventCount 
	      << " events accpeted out of " <<  m_analyzedEventCount 
	      << " tested ";
    if (m_passedSelectionCuts) std::cout << " This m= " << m_DiMuonPairInvMass; 
    std::cout << " * return " << m_passedSelectionCuts << std::endl; 
  }
  return m_passedSelectionCuts;
}


//==================================================================================
// Protected Methods
//==================================================================================
void ZmumuEvent::BookHistograms()
{
}

//==================================================================================
// Private Methods
//==================================================================================
bool ZmumuEvent::EventSelection(ZTYPE eType)
{
  if(m_doDebug) {std::cout <<" * ZmumuEvent * Event selection ** START ** for type: " << eType << " m_NumberOfFullPassMuons: " << m_numberOfFullPassMuons << std::endl;}

  // First require two muon-id's with cuts pre-applied.
  if ( m_numberOfFullPassMuons < 2 ) {
    if (m_doDebug) {std::cout <<" * ZmumuEvent * Failing number of good muons == 2 :( " <<  m_numberOfFullPassMuons << std::endl;}
    return false;
  }
  m_eventselectioncount_toofewmuons++;
  
  // crosscheck all muons have been properly filled
  bool allMuonsGood = true;
  for (unsigned int muonid=0; muonid < m_numberOfFullPassMuons; muonid++) {
    if (!m_pxMSTrack[muonid]) {
      allMuonsGood = false;
    } 
  }
  if (!allMuonsGood){
    return false;
  } 
  m_eventselectioncount_notallmuonsfilled++;
   
  if ( m_numberOfFullPassMuons > 2 ) {
    if (m_doDebug) {std::cout <<" * ZmumuEvent * Failing number of good muons == 2 :( " <<  m_numberOfFullPassMuons << std::endl;}
    return false;
  }
  m_eventselectioncount_morethantwomuons++;
  
  // momentum of the muons
  double leadingMuonPt, secondMuonPt;
  switch ( eType )
    {
    case MS :
      {
	leadingMuonPt = m_pxMSTrack[m_muon1]->pt();
	secondMuonPt = m_pxMSTrack[m_muon2]->pt();
	break;
      }
    case ME:
      {
	leadingMuonPt = m_pxMETrack[m_muon1]->pt();
	secondMuonPt = m_pxMETrack[m_muon2]->pt();
	break;
      }
    case CB:
      {
	leadingMuonPt = m_pxRecMuon[m_muon1]->pt();
	secondMuonPt = m_pxRecMuon[m_muon2]->pt();
	break;
      }
    case ID:
      {
	leadingMuonPt = m_pxIDTrack[m_muon1]->pt();
	secondMuonPt = m_pxIDTrack[m_muon2]->pt();
	break;
      }
          
    default:
      leadingMuonPt = m_pxRecMuon[m_muon1]->pt();
      secondMuonPt = m_pxRecMuon[m_muon2]->pt();
    }
  // up to here the leading and second pt are not really in the right order.
  // order the muon pt:
  if (secondMuonPt > leadingMuonPt) {
    double tempPt = leadingMuonPt;
    leadingMuonPt = secondMuonPt;
    secondMuonPt = tempPt;
  }
 
  // muon pt cut
  // if ( !(leadingMuonPt > m_LeadingMuonPtCut*CLHEP::GeV &&  secondMuonPt >  m_SecondMuonPtCut*CLHEP::GeV ) ) {
  if ( leadingMuonPt < m_LeadingMuonPtCut*CLHEP::GeV ) {
    if(m_doDebug){  std::cout <<" * ZmumuEvent * Failing 1st muon pt cut * Reco Pt:  " << leadingMuonPt << " < " << m_LeadingMuonPtCut*CLHEP::GeV << std::endl;}
    return false;
  }
  m_eventselectioncount_ptofleadingmuon++;
  
  if ( secondMuonPt < m_SecondMuonPtCut*CLHEP::GeV ) {
    if(m_doDebug){  std::cout <<" * ZmumuEvent * Failing 2nd muon pt cut * Reco Pt:  " << secondMuonPt << " < " <<  m_SecondMuonPtCut*CLHEP::GeV << std::endl;}
    return false;
  }
  m_eventselectioncount_ptofsecondmuon++;

  // Invariant mass window
  if ( m_fInvariantMass[eType]  < m_MassWindowLow  ) {
    if(m_doDebug) {std::cout <<" * ZmumuEvent * Failing mass window low cut:  reco m= " << m_fInvariantMass[eType] << " > " <<  m_MassWindowLow << std::endl;}
    return false;
  }
  if ( m_fInvariantMass[eType]  > m_MassWindowHigh ) {
    if(m_doDebug) {std::cout <<" * ZmumuEvent * Failing mass window high cut:  reco m= " << m_fInvariantMass[eType] << " > " <<  m_MassWindowHigh << std::endl;}
    return false;
  }
  m_eventselectioncount_masswindow++;

  // opening angle
  if ( m_fMuonDispersion[eType] <  m_OpeningAngleCut  ) {        
    if(m_doDebug) {std::cout <<" * ZmumuEvent * Failing opening angle cut. Opening angle " << m_fMuonDispersion[eType] << " < " <<  m_OpeningAngleCut << std::endl;}
    return false;
  }
  m_eventselectioncount_openingangle++;

  // opposite charge
  if ( getZCharge(eType) != 0  ) {
    if(m_doDebug) {
      std::cout <<" * ZmumuEvent * Failing get ZCharge != 0 cut * Reco q1= " << m_pxRecMuon[m_muon1]->charge()*m_pxRecMuon[m_muon1]->pt() <<std::endl;
	//<< "  q2= " <<  m_pxRecMuon[m_muon2]->charge()*m_pxRecMuon[m_muon2]->pt() << std::endl; //This might not exist!
      std::cout <<"                                             * ID   q1= " << m_pxIDTrack[m_muon1]->charge()*m_pxIDTrack[m_muon1]->pt() << std::endl;
	//	<< "  q2= " <<  m_pxIDTrack[m_muon2]->charge()*m_pxIDTrack[m_muon2]->pt() << std::endl; //This might not exist!
    }
    return false;
  }
  m_eventselectioncount_dimuoncharge++;

  //
  // both muons should come from the same vertex
  // if the vertex information is used, that is already guaranteed, but if not, one has to check the z0
  if (eType == ID) {
    double z0_muon1 = m_pxIDTrack[m_muon1]->vz() +  m_pxIDTrack[m_muon1]->z0();
    double z0_muon2 = m_pxIDTrack[m_muon2]->vz() +  m_pxIDTrack[m_muon2]->z0();
    if(m_doDebug) {
      std::cout << " * ZmumuEvent *  z0_muon1= " << z0_muon1 << "  z0_muon2= " << z0_muon2 << "  delta= " << z0_muon1-z0_muon2 << std::endl;
    }
    if ( std::abs(z0_muon1 - z0_muon2) > m_Z0GapCut) {
      if(m_doDebug) {
	std::cout << " * ZmumuEvent * Failing common vertex cut. z.vtx1= " << m_pxIDTrack[m_muon1]->vz() << "  z.vtx2= " << m_pxIDTrack[m_muon2]->vz() << std::endl;
	std::cout << " * ZmumuEvent * Failing common vertex cut. IDTrk.z0_1= " << m_pxIDTrack[m_muon1]->z0() << "  IDTrk.z0_2= " << m_pxIDTrack[m_muon2]->z0() << std::endl;
	std::cout << " * ZmumuEvent * z0_muon1= " << z0_muon1 << "  z0_muon2= " << z0_muon2 << "  delta= " << z0_muon1-z0_muon2 << " > " << m_Z0GapCut << " (cut)" << std::endl;
      } 
      return false;
    }
  }
  

  if(m_doDebug) {
    std::cout <<" * ZmumuEvent * Good muon pair: pt= " <<  leadingMuonPt/1000 
	      << " & " << secondMuonPt/1000 
	      << " GeV   dimuon invariant mass = " << m_fInvariantMass[eType] << " GeV " << std::endl;
  }
  return true;
}

//==================================================================================
void ZmumuEvent::Clear()
{
  m_numberOfFullPassMuons = 0;
  m_passedSelectionCuts   = false;
  m_DiMuonPairInvMass = -1.; // flag as no reconstructed inv mass yet
  m_muon1 = MUON1; // point to the first two
  m_muon2 = MUON2;

  for ( unsigned int u = 0; u < NUM_MUONS; ++u ) {
      m_pxRecMuon[u] = nullptr;
      m_pxMSTrack[u] = nullptr;
      m_pxMETrack[u] = nullptr;
      m_pxIDTrack[u] = nullptr;
  }
  for ( unsigned int v = 0; v < NUM_TYPES; ++v ) {
    m_fZPt[v]            = -999.9f;
    m_fZEtaDir[v]        = -999.9f;
    m_fZPhiDir[v]        = -999.9f;
    m_fInvariantMass[v]  = -999.9f;
    m_fMuonDispersion[v] = -999.9f;
  }
  return;
}

//==================================================================================
bool ZmumuEvent::RecordMuon( const xAOD::Muon* pxMuon )
{
  if(m_doDebug) { std::cout <<" * ZmumuEvent * RecordMuon * START ** muons recorded so far "<< m_numberOfFullPassMuons << " up to a maximum of " << NUM_MUONS << std::endl;}

  // This shouldn't really ever happen but just in case.
  if ( !pxMuon ) {
    if(m_doDebug) { std::cout <<" * ZmumuEvent * RecordMuon * bad pxMuon --> EXIT "<< std::endl;}
    return false;
  }

  if ( m_numberOfFullPassMuons < NUM_MUONS ) {
      // The main Muon
      m_pxRecMuon[m_numberOfFullPassMuons] = pxMuon;
      // Tracking Muon Spectrometer ( raw )
      const xAOD::TrackParticle* pxMSTrack   = pxMuon->trackParticle(xAOD::Muon::MuonSpectrometerTrackParticle);
      if (!pxMSTrack) {
	if(m_doDebug){  std::cout <<" * ZmumuEvent * RecordMuon * bad pxMSmuon --> EXIT "<< std::endl;}
	return false;
      } 
      m_pxMSTrack[m_numberOfFullPassMuons] = pxMSTrack;

      // Tracking ID ( fix later to include loose match track conditions )
      const xAOD::TrackParticle*  pxIDTrack  = pxMuon->trackParticle(xAOD::Muon::InnerDetectorTrackParticle);
      if (!pxIDTrack) {
	return false;
      }      
      m_pxIDTrack[m_numberOfFullPassMuons] = pxIDTrack;
      //
      if(m_doDebug){  std::cout <<"                m_pxRecMuon[" << m_numberOfFullPassMuons <<"]" 
				<< "  pt= " << m_pxRecMuon[m_numberOfFullPassMuons]->pt() << "  q= "<< m_pxRecMuon[m_numberOfFullPassMuons]->charge() << std::endl;}
      if(m_doDebug){  std::cout <<"                m_pxMSTrack[" << m_numberOfFullPassMuons <<"]" 
				<<"   pt= " << m_pxMSTrack[m_numberOfFullPassMuons]->pt() << "  q= " << m_pxMSTrack[m_numberOfFullPassMuons]->charge() << std::endl;}
      if(m_doDebug){  std::cout <<"                m_pxIDTrack[" << m_numberOfFullPassMuons <<"]"
				<<"   pt= " << m_pxIDTrack[m_numberOfFullPassMuons]->pt() << "  q= " << m_pxIDTrack[m_numberOfFullPassMuons]->charge() << std::endl;}
      // update count
      ++m_numberOfFullPassMuons;
  }
  if(m_doDebug){  std::cout <<" * ZmumuEvent * RecordMuon * return with a total of " << m_numberOfFullPassMuons << std::endl;}
  return true;
}


//==================================================================================
void ZmumuEvent::ReconstructKinematics()
{
  // Three ways. No checks here. Thus make sure the pointers are ok before this.
  if ( m_numberOfFullPassMuons == 2 )
    {
      // Note that all the util. functions will check the pointers & return -999.9f on failure.
      m_fInvariantMass[MS]      = EvalDiMuInvMass( m_pxMSTrack[m_muon1], m_pxMSTrack[m_muon2] );
      m_fMuonDispersion[MS]     = EvaluateAngle(   m_pxMSTrack[m_muon1], m_pxMSTrack[m_muon2] );
      m_fZPt[MS]                = EvalPt(          m_pxMSTrack[m_muon1], m_pxMSTrack[m_muon2] );
      m_fZEtaDir[MS]            = EvalEta(         m_pxMSTrack[m_muon1], m_pxMSTrack[m_muon2] );
      m_fZPhiDir[MS]            = EvalPhi(         m_pxMSTrack[m_muon1], m_pxMSTrack[m_muon2] );

      m_fInvariantMass[CB]      = EvalDiMuInvMass( m_pxRecMuon[m_muon1], m_pxRecMuon[m_muon2] );
      m_fMuonDispersion[CB]     = EvaluateAngle(   m_pxRecMuon[m_muon1], m_pxRecMuon[m_muon2] );
      m_fZPt[CB]                = EvalPt(          m_pxRecMuon[m_muon1], m_pxRecMuon[m_muon2] );
      m_fZEtaDir[CB]            = EvalEta(         m_pxRecMuon[m_muon1], m_pxRecMuon[m_muon2] );
      m_fZPhiDir[CB]            = EvalPhi(         m_pxRecMuon[m_muon1], m_pxRecMuon[m_muon2] );

      m_fInvariantMass[ID]      = EvalDiMuInvMass( m_pxIDTrack[m_muon1], m_pxIDTrack[m_muon2]);
      m_fMuonDispersion[ID]     = EvaluateAngle(   m_pxIDTrack[m_muon1], m_pxIDTrack[m_muon2] );
      m_fZPt[ID]                = EvalPt(          m_pxIDTrack[m_muon1], m_pxIDTrack[m_muon2] );
      m_fZEtaDir[ID]            = EvalEta(         m_pxIDTrack[m_muon1], m_pxIDTrack[m_muon2] );
      m_fZPhiDir[ID]            = EvalPhi(         m_pxIDTrack[m_muon1], m_pxIDTrack[m_muon2] );
    }
}

//==================================================================================
float ZmumuEvent::getPtImbalance( ZTYPE eType )
{
  // First determine what's positive
  if ( m_numberOfFullPassMuons == 2 )
    {
      switch ( eType )
	{
	case MS :
	  {
	    return EvalPtDiff( m_pxMSTrack[m_muon1], m_pxMSTrack[m_muon2] );
	  }
	case ME:
	  {
	    return EvalPtDiff( m_pxMETrack[m_muon1], m_pxMETrack[m_muon2] );
	  }
	case CB:
	  {
	    return EvalPtDiff( m_pxRecMuon[m_muon1], m_pxRecMuon[m_muon2] );
	  }
	case ID:
	  {
	    return EvalPtDiff( m_pxIDTrack[m_muon1], m_pxIDTrack[m_muon2] );
	  }
	default:
	  return -999.0;
	}
    }
  else
    {
      return -999.0;
    }
}

//==================================================================================
int ZmumuEvent::getZCharge( ZTYPE eType )
{
  switch ( eType )
    {
    case MS :
      {
	return ( static_cast<int>( EvalCharge( m_pxMSTrack[m_muon1], m_pxMSTrack[m_muon2] ) ) );
      }
    case ME:
      {
	return ( static_cast<int>( EvalCharge( m_pxMETrack[m_muon1], m_pxMETrack[m_muon2] ) ) );
      }
    case CB:
      {
	return ( static_cast<int>( EvalCharge( m_pxRecMuon[m_muon1], m_pxRecMuon[m_muon2] ) ) );
      }
    case ID:
      {
	return ( static_cast<int>( EvalCharge( m_pxIDTrack[m_muon1], m_pxIDTrack[m_muon2] ) ) );
      }
    default:
      return -999;
    }
}

//==================================================================================
unsigned int ZmumuEvent::getPosMuon( ZTYPE eType )
{
  if ( getNumberOfTaggedMuons() != 2 ) return 999;
  if ( getZCharge(eType) != 0        ) return 999;

  switch ( eType )
    {
    case MS :
      {
	if ( !m_pxMSTrack[m_muon1] || !m_pxMSTrack[m_muon2] ) return 999;
	return ( static_cast<int>( m_pxMSTrack[m_muon1]->charge() ) == 1  ? m_muon1 : m_muon2 );
      }
    case ME:
      {
	if ( !m_pxMETrack[m_muon1] || !m_pxMETrack[m_muon2] ) return 999;
	return ( static_cast<int>( m_pxMETrack[m_muon1]->charge() ) == 1  ? m_muon1 : m_muon2 );
      }
    case CB:
      {
	if ( !m_pxRecMuon[m_muon1] || !m_pxRecMuon[m_muon2] ) return 999;
	return ( static_cast<int>( m_pxRecMuon[m_muon1]->charge() ) == 1  ? m_muon1 : m_muon2 );
      }
    case ID:
      {
	if ( !m_pxIDTrack[m_muon1] || !m_pxIDTrack[m_muon2] ) return 999;
	return ( static_cast<int>( m_pxIDTrack[m_muon1]->charge() ) == 1 ? m_muon1 : m_muon2 );
      }
    default:
      return 999;
    }
}

//==================================================================================
unsigned int ZmumuEvent::getNegMuon( ZTYPE eType )
{
  int uTmp = getPosMuon( eType );
  if ( uTmp == 999 )
    {
      return 999;
    }
  else
    {
      return ( ( uTmp == m_muon1 ) ? m_muon2 : m_muon1 );
    }
}

//==================================================================================
const xAOD::TrackParticle*  ZmumuEvent::getLooseIDTk( unsigned int /*uPart*/ )
{
  const xAOD::TrackParticleContainer*  pxTrackContainer =
    PerfMonServices::getContainer<xAOD::TrackParticleContainer>( PerfMonServices::TRK_COLLECTION );

  if ( pxTrackContainer )
    {
      xAOD::TrackParticleContainer::const_iterator xTrkItr  = pxTrackContainer->begin();
      xAOD::TrackParticleContainer::const_iterator xTrkItrE  = pxTrackContainer->end();
      while ( xTrkItr != xTrkItrE )
	{
	  const xAOD::TrackParticle* pxTrack = *xTrkItr;
	  if ( !pxTrack ) continue;
	  const Trk::Track* pxTrkTrack = pxTrack->track();
	  if(!pxTrkTrack) continue;	  
	  const Trk::Perigee* pxPerigee = pxTrkTrack->perigeeParameters() ;
	  if ( !pxPerigee ) continue;
	  const float fTrkPhi   = pxPerigee->parameters()[Trk::phi];
	  const float fTrkEta   = pxPerigee->eta();

	  float fDPhi = std::abs( fTrkPhi -  m_pxMETrack[m_muon1]->phi() );
	  float fDEta = std::abs( fTrkEta -  m_pxMETrack[m_muon2]->eta() );
	  float fDR = sqrt( fDPhi*fDPhi + fDEta*fDEta );

	  if ( fDR < 0.3f )
	    {
	      return pxTrack;
	    }

	  ++xTrkItr;
	}
    }
  // if ()
  return nullptr;
}

//==================================================================================
void ZmumuEvent::SetLeadingMuonPtCut (double newvalue)
{
  // first set the new pt cut value
  m_LeadingMuonPtCut = newvalue;
  

  // the second muon pt cut can not be higher than the leading muon pt cut:
  if (m_LeadingMuonPtCut < m_SecondMuonPtCut) this->SetSecondMuonPtCut(m_LeadingMuonPtCut);

  // this has to be translated to the MuonSelector
  // but there one has to use the minimum momentum --> second muon
  //this->SetMuonPtCut(m_SecondMuonPtCut);
  if(m_doDebug && false){std::cout <<" * ZmumuEvent * SetLeadingMuonPtCut * new Pt cuts:  " << m_LeadingMuonPtCut << " & " << m_SecondMuonPtCut << "  MuonSelector: " << m_xMuonID.GetPtCut() << std::endl;}
  return;
}

//==================================================================================
void ZmumuEvent::SetSecondMuonPtCut (double newvalue) 
{
  m_SecondMuonPtCut = newvalue;

  // second muon pt shouldn't be higher than the leading muon pt
  if (m_LeadingMuonPtCut < m_SecondMuonPtCut) this->SetLeadingMuonPtCut(m_LeadingMuonPtCut);

  // this has to be translated to the MuonSelector
  this->SetMuonPtCut(m_SecondMuonPtCut);

  if(m_doDebug && false){std::cout <<" * ZmumuEvent * SetSecondMuonPtCut * new Pt cuts:  " << m_LeadingMuonPtCut << " & " << m_SecondMuonPtCut << "  MuonSelector: " << m_xMuonID.GetPtCut() << std::endl;}

  return;
}

//==================================================================================
void ZmumuEvent::OrderMuonList()
{
  int muPlusId = -9;
  int muMinusId = -9;
  double muPlusPt = 0.;
  double muMinusPt = 0.;

  if (m_doDebug) {std::cout <<" * ZmumuEvent * OrderMuonList * START *  input number of muons: " << m_numberOfFullPassMuons << std::endl;}

  if (m_numberOfFullPassMuons >= 2) {
    for (int imuon=0; imuon < (int) m_numberOfFullPassMuons; imuon++) {
      if (m_pxRecMuon[imuon] != nullptr) {
	
	if (m_pxRecMuon[imuon]->charge()== 1 && m_pxRecMuon[imuon]->pt()> muPlusPt) {
	  muPlusPt = m_pxRecMuon[imuon]->pt();
	  muPlusId = imuon;
	} 
	if (m_pxRecMuon[imuon]->charge()==-1 && m_pxRecMuon[imuon]->pt()> muMinusPt) {
	  muMinusPt = m_pxRecMuon[imuon]->pt();
	  muMinusId = imuon;
	} 
      } // muon exist
    } // for (int imuon
  } // if (m_numberOfFullPassMuons >= 2)
  if (muPlusId>=0 && muMinusId>=0) {
    m_muon1 = muPlusId;
    m_muon2 = muMinusId;
    m_numberOfFullPassMuons = 2; // the two muons have been selected. Let's pretend we have only two muons then.
  }
  if (m_doDebug) {std::cout <<" * ZmumuEvent * OrderMuonList * COMPLETED ** numberOfFullPassMuons: " << m_numberOfFullPassMuons 
			    << "  mu+: " << muPlusId << " mu-: " << muMinusId << std::endl;}
  return;
}


//==================================================================================
void ZmumuEvent::finalize()
{
  m_xMuonID.finalize();
  
  std::cout << " ** ZmumuEvent ** -- STATS -- " << std::endl
	    << "    Analyzed events           : " << m_analyzedEventCount << std::endl
	    << "    Tested muons              : " << m_testedMuonCount << std::endl
	    << "    Accepted muons            : " << m_acceptedMuonCount << std::endl
	    << "    Events without enough muon: " << m_eventsWithoutEnoughMuonsCount << std::endl
	    << "    Events with enough muons  : " << m_eventsWithEnoughMuonsCount << std::endl
            << "    pass few muons            : " << m_eventselectioncount_toofewmuons << std::endl
	    << "    pass not filled muons     : " << m_eventselectioncount_notallmuonsfilled << std::endl
	    << "    pass more than 2 muons    : " << m_eventselectioncount_morethantwomuons << std::endl
	    << "    pass pt lead              : " << m_eventselectioncount_ptofleadingmuon << std::endl
	    << "    pass pt 2nd               : " << m_eventselectioncount_ptofsecondmuon << std::endl
	    << "    pass mass window          : " << m_eventselectioncount_masswindow << std::endl	
	    << "    pass opening angle        : " << m_eventselectioncount_openingangle << std::endl
	    << "    pass dimuon charge        : " << m_eventselectioncount_dimuoncharge << std::endl
	    << "    Accepted events           : " << m_acceptedEventCount << std::endl
	    << std::endl;
  return;
}

