///////////////////////// -*- C++ -*- /////////////////////////////

/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// PhotonCnv_p5.cxx 
// Implementation file for class PhotonCnv_p5
/////////////////////////////////////////////////////////////////// 


// STL includes

// egammaEvent includes
#include "egammaEvent/Photon.h"

// egammaEventTPCnv includes
#include "egammaEventTPCnv/PhotonCnv_p5.h"

/////////////////////////////////////////////////////////////////// 
// Public methods: 
/////////////////////////////////////////////////////////////////// 

// Constructors
////////////////

// Destructor
///////////////

/////////////////////////////////////////////////////////////////// 
// Const methods: 
///////////////////////////////////////////////////////////////////

void PhotonCnv_p5::persToTrans( const Photon_p5* pers,
				Analysis::Photon* trans, 
				MsgStream& msg ) 
{
//   msg << MSG::DEBUG << "Loading Photon from persistent state..."
//       << endreq;

  // base class
  m_egammaCnv.persToTrans( &pers->m_egamma, trans, msg );

//   msg << MSG::DEBUG << "Loaded Photon from persistent state [OK]"
//       << endreq;

  return;
}

void PhotonCnv_p5::transToPers( const Analysis::Photon* trans, 
				Photon_p5* pers, 
				MsgStream& msg ) 
{
//   msg << MSG::DEBUG << "Creating persistent state of Photon..."
//       << endreq;

  // base class
  m_egammaCnv.transToPers( trans, &pers->m_egamma, msg );

//   msg << MSG::DEBUG << "Created persistent state of Photon [OK]"
//       << endreq;
  return;
}
