///////////////////////// -*- C++ -*- /////////////////////////////

/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// Photon_p2.h 
// Header file for class Photon_p2
/////////////////////////////////////////////////////////////////// 
#ifndef EGAMMAEVENTTPCNV_PHOTON_P2_H 
#define EGAMMAEVENTTPCNV_PHOTON_P2_H 1

// egammaEventTPCnv includes
#include "egammaEventTPCnv/egamma_p2.h"

// forward declarations
class PhotonCnv_p2;

class Photon_p2
{
  /////////////////////////////////////////////////////////////////// 
  // Friend classes
  /////////////////////////////////////////////////////////////////// 

  // Make the AthenaPoolCnv class our friend
  friend class PhotonCnv_p2;

  /////////////////////////////////////////////////////////////////// 
  // Public methods: 
  /////////////////////////////////////////////////////////////////// 
public: 

  /** Default constructor: 
   */
  Photon_p2();

  /** Destructor: 
   */
  ~Photon_p2();

  /////////////////////////////////////////////////////////////////// 
  // Private data: 
  /////////////////////////////////////////////////////////////////// 
private: 

  /// the egamma part 
  egamma_p2 m_egamma;
  
}; 

/////////////////////////////////////////////////////////////////// 
// Inline methods: 
/////////////////////////////////////////////////////////////////// 

inline Photon_p2::Photon_p2() :
  m_egamma()
{}

#endif //> EGAMMAEVENTTPCNV_PHOTON_P2_H
