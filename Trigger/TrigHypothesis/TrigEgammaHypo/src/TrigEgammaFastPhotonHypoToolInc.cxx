
/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/


#include "TrigCompositeUtils/HLTIdentifier.h"
#include "TrigCompositeUtils/Combinators.h"
#include "AthenaMonitoringKernel/Monitored.h"
#include "TrigEgammaFastPhotonHypoToolInc.h"

namespace TCU = TrigCompositeUtils;

TrigEgammaFastPhotonHypoToolInc::TrigEgammaFastPhotonHypoToolInc( const std::string& type, 
            const std::string& name, 
            const IInterface* parent ) 
  : AthAlgTool( type, name, parent ),
    m_decisionId( HLT::Identifier::fromToolName( name ) ) {}



StatusCode TrigEgammaFastPhotonHypoToolInc::initialize()  {
  
  if ( !m_monTool.empty() ) 
    ATH_CHECK( m_monTool.retrieve() );

  ATH_MSG_DEBUG( "Initialization completed successfully:" );
  ATH_MSG_DEBUG( "AcceptAll      = " << ( m_acceptAll==true ? "True" : "False" ) );
  ATH_MSG_DEBUG( "EtaBins        = " << m_etabin      );
  ATH_MSG_DEBUG( "ETthr          = " << m_eTthr    << "(lo)/" << m_eT2thr    << "(hi)"     );
  ATH_MSG_DEBUG( "HADETthr       = " << m_hadeTthr << "(lo)/" << m_hadeT2thr << "(hi)"     );
  ATH_MSG_DEBUG( "CARCOREthr     = " << m_carcorethr  );
  ATH_MSG_DEBUG( "CAERATIOthr    = " << m_caeratiothr );
  ATH_MSG_DEBUG( "dPHICLUSTERthr = " << m_dphicluster );
  ATH_MSG_DEBUG( "dETACLUSTERthr = " << m_detacluster );


  std::vector<size_t> sizes( {m_eTthr.size(), m_eT2thr.size(), 
                              m_hadeTthr.size(), m_hadeT2thr.size(),
                              m_carcorethr.size(), m_caeratiothr.size() } );   

  if ( *std::min_element( sizes.begin(), sizes.end() ) != *std::max_element( sizes.begin(), sizes.end() )  ) {     
    ATH_MSG_ERROR( "Missconfiguration, cut properties listed above ( when DEBUG ) have different dimensions shortest: " 
                    <<  *std::min_element( sizes.begin(), sizes.end() ) << " longest " 
                    << *std::max_element( sizes.begin(), sizes.end() ) );     
    return StatusCode::FAILURE;   
  }

  return StatusCode::SUCCESS;
}

//==================================================================

StatusCode TrigEgammaFastPhotonHypoToolInc::decide( std::vector<PhotonInfo>& input)  const {
  for ( auto& i: input ) {
    if ( TCU::passed ( m_decisionId.numeric(), i.previousDecisionIDs ) ) {
      if ( decide( i.photon ) ) {
        TCU::addDecisionID( m_decisionId, i.decision );
      }
    }
  }
  return StatusCode::SUCCESS;
}

//==================================================================

bool TrigEgammaFastPhotonHypoToolInc::decide( const xAOD::TrigPhoton* photon ) const {

  auto cutCounter = Monitored::Scalar<int>( "CutCounter", -1 );
  auto PhEt =  Monitored::Scalar( "PhEt", -99. );
  auto PhEta = Monitored::Scalar( "PhEta", -99. );
  auto PhPhi = Monitored::Scalar( "PhPhi", -99. );
  auto dEta = Monitored::Scalar( "dEta", -99. );
  auto dPhi = Monitored::Scalar( "dPhi", -99. );
  auto PhRcore = Monitored::Scalar( "PhRcore", -99. );
  auto PhEratio = Monitored::Scalar( "PhRcore", -99. );
  auto PhHadEt = Monitored::Scalar( "PhHadEt", -99. );
  auto PhF1 = Monitored::Scalar( "PhF1", -99. );
  auto monitorIt  = Monitored::Group( m_monTool,
                                      cutCounter, 
                                      PhEt,
                                      PhEta, PhPhi,
                                      dEta, dPhi,
                                      PhRcore, PhEratio,
                                      PhHadEt, PhF1 );

  float EmET         = -99.0;
  float HadEmRatio   = -99.0;
  float Reta         = -99.0;
  float Eratio       = -99.0;
  float f1           = -99.0;
  float HadET        = -99.0;

  if(!photon) return false;

  cutCounter++;
  
  // Determine which eta bin to apply the cuts                                                                                                                                                              
  float absEta = std::abs( photon->eta() );

  int etaBin = findCutIndex(absEta);

  // getting photon variable           
  Eratio     = photon->eratio();
  Reta       = photon->rcore();
  EmET       = photon->pt();
  HadET      = photon->etHad();
  f1         = photon->f1();

  if(m_acceptAll) {
    ATH_MSG_DEBUG ( "Accept all property is set: TrigPhoton: ET_em=" << EmET << " cut in etaBin " 
                    << etaBin << " is ET_em >= " << m_eTthr[0] );
    return true;
  }
  
  // eta range                                                                                                                                                                                              
  if ( etaBin==-1 ) {
    ATH_MSG_INFO( "Photon eta: " << absEta << " outside eta range " << m_etabin[m_etabin.size()-1] );
    return true;
  } else {
    ATH_MSG_INFO( "eta bin used for cuts " << etaBin );
  }
  cutCounter++; // passed eta cut       


  // Reta (was previously called Rcore) 
  if ( Reta > m_carcorethr[etaBin] ){
    ATH_MSG_INFO( "TrigPhoton Reta=" << Reta << " cut in etaBin " 
                      << etaBin << " is Reta >= "  << m_carcorethr[etaBin]  );
    return  false;
  }
  cutCounter++;


  //  // Eratio           
  bool inCrack = ( absEta > 2.37 || ( absEta > 1.37 && absEta < 1.52) );
  if ( inCrack || f1<m_F1thr[etaBin] ) {
    ATH_MSG_INFO(  "TrigPhoton: InCrack= " << inCrack << " F1=" << f1
                     << " Eratio cut not being applied" );
  } else {
    if ( Eratio > m_caeratiothr[etaBin] ) {
      return false;
    }
  }
  cutCounter++;
  if(inCrack)  Eratio  = -1; //Set default value in crack for monitoring.                                                                                                                                   

  // ET_em                                                                                                                                                                                                  
  if ( EmET < m_eTthr[etaBin]) {
    ATH_MSG_INFO( "TrigPhoton: ET_em=" << EmET
		     << " not in etaBin " << etaBin << " is ET_em < " << m_eTthr[etaBin] );
    return false;
  }
  cutCounter++;

  // ET_had                                                                                                                                                                                                 
  // find which ET_had to apply : this depends on the ET_em and the eta bin                                                                                                                                 
  float hadET_cut=-1;

  if ( EmET >  m_eT2thr[etaBin] ) {
    hadET_cut = m_hadeT2thr[etaBin] ;
    ATH_MSG_INFO( "ET_em>" << m_eT2thr[etaBin] );
  } else {
    hadET_cut = m_hadeTthr[etaBin];
    ATH_MSG_INFO( "ET_em<" << m_eT2thr[etaBin] );
  }

  HadEmRatio = (EmET!=0) ? HadET/EmET : -1.0;

  if ( HadEmRatio < hadET_cut ){
    ATH_MSG_INFO( "TrigPhoton: ET_had=" <<  HadEmRatio
		     << "  not in etaBin " << etaBin );
    return false;

  }
  cutCounter++;
  
  return true;
  
}

//==================================================================

int TrigEgammaFastPhotonHypoToolInc::findCutIndex( float eta ) const {
  const float absEta = std::abs(eta);
  auto binIterator = std::adjacent_find( m_etabin.begin(), m_etabin.end(), [=](float left, float right){ return left < absEta and absEta < right; }  );
  if ( binIterator == m_etabin.end() ) {
    return -1;
  }
  return  binIterator - m_etabin.begin();
}


