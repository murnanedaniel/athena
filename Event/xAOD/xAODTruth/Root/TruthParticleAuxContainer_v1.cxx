/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// $Id: TruthParticleAuxContainer_v1.cxx 574253 2013-12-06 15:21:49Z krasznaa $

// Local include(s):
#include "xAODTruth/versions/TruthParticleAuxContainer_v1.h"

namespace xAOD {

   TruthParticleAuxContainer_v1::TruthParticleAuxContainer_v1()
   : AuxContainerBase() {

      AUX_VARIABLE( pdgId );
      AUX_VARIABLE( barcode );
      AUX_VARIABLE( status );
      AUX_VARIABLE( prodVtxLink );
      AUX_VARIABLE( decayVtxLink );
      AUX_VARIABLE( px );
      AUX_VARIABLE( py );
      AUX_VARIABLE( pz );
      AUX_VARIABLE( e );
      AUX_VARIABLE( m );
   }

} // namespace xAOD
