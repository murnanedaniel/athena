/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// System include(s):
#include <stdexcept>
#include <vector>

// Gaudi/Athena include(s):
#include "GaudiKernel/MsgStream.h"

//EDM include(s):
#include "AthContainers/tools/copyAuxStoreThinned.h"
#include "xAODTrigRinger/versions/TrigRNNOutputContainer_v1.h"
#include "xAODTrigRinger/TrigRNNOutputContainer.h"


// Local include(s):
#include "xAODTrigRNNOutputAuxContainerCnv_v1.h"

/// Convenience macro for setting the level of output messages
#define MSGLVL MSG::INFO

/// Another convenience macro for printing messages in the converter
#define ATH_MSG( MSG )                          \
   do {                                         \
      if( log.level() <= MSGLVL ) {             \
         log << MSGLVL << MSG << endmsg;        \
      }                                         \
   } while( 0 )

xAODTrigRNNOutputAuxContainerCnv_v1::xAODTrigRNNOutputAuxContainerCnv_v1(): T_AthenaPoolTPCnvBase< xAOD::TrigRNNOutputAuxContainer,
                                                                                                   xAOD::TrigRNNOutputAuxContainer_v1 >() {

}

void xAODTrigRNNOutputAuxContainerCnv_v1::persToTrans( const xAOD::TrigRNNOutputAuxContainer_v1* oldObj,
                                                             xAOD::TrigRNNOutputAuxContainer* newObj,
                                                             MsgStream& /*log*/ ) {

   // Greet the user:
   //ATH_MSG( "Converting xAOD::TrigRNNOutputAuxContainer_v1 to current version..." );

   // Remove this line once the converter is "ready":
   //ATH_MSG( "WARNING Converter is not complete yet!" );

   // Clear the transient object:
   newObj->resize( 0 );

   SG::copyAuxStoreThinned( *oldObj, *newObj, 0 );

   //The old  uses v1
   xAOD::TrigRNNOutputContainer_v1 oldInt;
   for( size_t i = 0; i < oldObj->size(); ++i ) {
     oldInt.push_back( new xAOD::TrigRNNOutput_v1() );
   }
   oldInt.setStore( oldObj );

   // The new uses v2
   xAOD::TrigRNNOutputContainer newInt;
   for (size_t i = 0; i < newObj->size(); ++i) {
      newInt.push_back( new xAOD::TrigRNNOutput() );
   }

   newInt.setStore( newObj);

   for (size_t i = 0; i < newObj->size(); ++i) {

     std::vector<float> dec = oldInt[i]->decision();
     newInt[i]->setRnnDecision( dec );
     newInt[i]->auxdata<uint32_t>("RoIword") = (uint32_t)oldInt[i]->RoIword();
     newInt[i]->auxdata<float>("et") = oldInt[i]->et();

   }

   return;
}

/// This function should never be called, as we are not supposed to convert
/// object before writing.
///
void xAODTrigRNNOutputAuxContainerCnv_v1::transToPers( const xAOD::TrigRNNOutputAuxContainer*,
                                                             xAOD::TrigRNNOutputAuxContainer_v1*,
                                                             MsgStream& log ) 
{

   log << MSG::ERROR << "Somebody called xAODTrigRNNOutputAuxContainerCnv_v1::transToPers" << endmsg;
   throw std::runtime_error( "Somebody called xAODTrigRNNOutputAuxContainerCnv_v1::"
                             "transToPers" );

   return;
}
