///////////////////////// -*- C++ -*- /////////////////////////////

/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// CaloEnergyCnv_p2.cxx 
// Implementation file for class CaloEnergyCnv_p2
// Author: S.Binet<binet@cern.ch>
/////////////////////////////////////////////////////////////////// 


// STL includes

// AthenaPoolCnvSvc includes
#include "AthenaPoolCnvSvc/T_AthenaPoolTPConverter.h"

// muonEvent includes
#define private public
#define protected public
#include "muonEvent/DepositInCalo.h"
#include "muonEvent/CaloEnergy.h"
#undef private
#undef protected

// RecTPCnv includes
#include "RecTPCnv/DepositInCaloCnv_p2.h"
#include "RecTPCnv/CaloEnergyCnv_p2.h"
#include "TrkEventTPCnv/TrkMaterialOnTrack/EnergyLossCnv_p1.h"

// For converter.
#include "RecTPCnv/DepositInCalo_p1.h"
#include "RecTPCnv/DepositInCalo_p2.h"
#include "RootConversions/VectorConverter.h"
#include "RootConversions/TConverterRegistry.h"

// pre-allocate converters
static DepositInCaloCnv_p2 depositCnv;
static EnergyLossCnv_p1 energyLossCnv;

CaloEnergyCnv_p2::CaloEnergyCnv_p2()
{
}

/////////////////////////////////////////////////////////////////// 
// Non-Const methods: 
///////////////////////////////////////////////////////////////////

void CaloEnergyCnv_p2::persToTrans( const CaloEnergy_p2* pers,
				    CaloEnergy* trans, 
				    MsgStream& msg ) 
{
//   msg << MSG::DEBUG << "Loading CaloEnergy from persistent state..."
//       << endmsg;

      /// energy loss in calorimeter
  energyLossCnv.persToTrans( &pers->m_energyLoss,
			     trans,
			     msg );

  trans->m_energyLossType     = static_cast<CaloEnergy::EnergyLossType>(pers->m_energyLossType);
  trans->m_caloLRLikelihood   = pers->m_caloLRLikelihood;
  trans->m_caloMuonIdTag      = pers->m_caloMuonIdTag;
  trans->m_fsrCandidateEnergy = pers->m_fsrCandidateEnergy;

  // reserve enough space so no re-alloc occurs
  trans->m_deposits.reserve( pers->m_deposits.size() );

  typedef std::vector<DepositInCalo_p2> Deposits_t;
  for ( Deposits_t::const_iterator 
	  itr  = pers->m_deposits.begin(),
	  iEnd = pers->m_deposits.end();
	itr != iEnd;
	++itr ) {
    trans->m_deposits.push_back( DepositInCalo() );
    depositCnv.persToTrans( &*itr, &trans->m_deposits.back(), msg );
  }

  trans->m_etCore = pers->m_etCore;

//   msg << MSG::DEBUG << "Loaded CaloEnergy from persistent state [OK]"
//       << endmsg;

  return;
}

void CaloEnergyCnv_p2::transToPers( const CaloEnergy* trans, 
				    CaloEnergy_p2* pers, 
				    MsgStream& msg ) 
{
//   msg << MSG::DEBUG << "Creating persistent state of CaloEnergy..."
//       << endmsg;

  energyLossCnv.transToPers ( trans, &pers->m_energyLoss,  msg );

  pers->m_energyLossType     = trans->m_energyLossType;
  pers->m_caloLRLikelihood   = trans->m_caloLRLikelihood;
  pers->m_caloMuonIdTag      = trans->m_caloMuonIdTag;
  pers->m_fsrCandidateEnergy = trans->m_fsrCandidateEnergy;
  
  // reserve enough space so no re-alloc occurs
  pers->m_deposits.reserve( trans->m_deposits.size() );

  typedef std::vector<DepositInCalo> Deposits_t;
  for ( Deposits_t::const_iterator 
	  itr  = trans->m_deposits.begin(),
	  iEnd = trans->m_deposits.end();
	itr != iEnd;
	++itr ) {
    pers->m_deposits.push_back( DepositInCalo_p2() );
    depositCnv.transToPers( &*itr, &pers->m_deposits.back(), msg );
  }

  pers->m_etCore = trans->m_etCore;

//   msg << MSG::DEBUG << "Created persistent state of CaloEnergy [OK]"
//       << endmsg;
  return;
}


/**
 * @brief Register a streamer converter for backwards compatibility
 *        for the vector<DepositInCalo_p1> -> vector<DepositInCalo_p2>
 *        change.
 */
void CaloEnergyCnv_p2::registerStreamerConverter()
{
  // m_deposits was changed from vector<DepositInCalo_p1> to
  // vector<DepositInCalo_p2>.  Install a converter so that we can read
  // the old data.
  TConverterRegistry::Instance()->AddStreamerConverter
    ("vector<DepositInCalo_p1>", "vector<DepositInCalo_p2>",
     new RootConversions::VectorConverter<DepositInCalo_p1,DepositInCalo_p2>
     ("DepositInCalo_p1"));
}

/**
 * @brief register a C function to be executed at library loading time
 *        this is b/c there is no more the ability run the above static method
 *        `CaloEnergyCnv_p2::registerStreamerConverter` from pyroot (persistent
 *        classes' dictionaries are generated by reflex with --dataonly)
 */
extern "C" {
void caloenergy_cnv_p2_register_streamer()
{
  CaloEnergyCnv_p2::registerStreamerConverter();
}

}
