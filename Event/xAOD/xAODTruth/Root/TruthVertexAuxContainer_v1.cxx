/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// $Id: TruthVertexAuxContainer_v1.cxx 595827 2014-05-07 13:18:18Z jmitrevs $

// Local include(s):
#include "xAODTruth/versions/TruthVertexAuxContainer_v1.h"

namespace xAOD {

   TruthVertexAuxContainer_v1::TruthVertexAuxContainer_v1()
   : AuxContainerBase() {

      AUX_VARIABLE( id );
      AUX_VARIABLE( barcode );
      AUX_VARIABLE( incomingParticleLinks );
      AUX_VARIABLE( outgoingParticleLinks );
      AUX_VARIABLE( x );
      AUX_VARIABLE( y );
      AUX_VARIABLE( z );
      AUX_VARIABLE( t );
      AUX_VARIABLE( weights );
   }

} // namespace xAOD
