/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// $Id$
	
// System include(s):
#include <stdexcept>
	
// Gaudi/Athena include(s):
#include "GaudiKernel/MsgStream.h"
	
// Core EDM include(s):
#include "AthContainers/tools/copyAuxStoreThinned.h"
	
// Local include(s):
#include "xAODTrackParticleAuxContainerCnv_v3.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTracking/versions/TrackParticleContainer_v1.h"
	
/// Convenience macro for setting the level of output messages
#define MSGLVL MSG::DEBUG
	
/// Another convenience macro for printing messages in the converter
#define ATH_MSG( MSG )                          \
do {                                         \
  if( log.level() <= MSGLVL ) {             \
    log << MSGLVL << MSG << endmsg;        \
  }                                         \
} while( 0 )
	
xAODTrackParticleAuxContainerCnv_v3::xAODTrackParticleAuxContainerCnv_v3()
    : T_AthenaPoolTPCnvBase< xAOD::TrackParticleAuxContainer, xAOD::TrackParticleAuxContainer_v3 >() 
{
}
	
void xAODTrackParticleAuxContainerCnv_v3::
persToTrans(  const xAOD::TrackParticleAuxContainer_v3* oldObj, 
              xAOD::TrackParticleAuxContainer* newObj,
              MsgStream& log ) {
	
  // Greet the user:
  ATH_MSG( "Converting xAOD::TrackParticleAuxContainer_v3 to current version..." );
	
  // Clear the transient object:
  newObj->resize( 0 );
	
  // Copy the payload of the v3 object into the latest one by misusing
  // the thinning code a bit...
  SG::copyAuxStoreThinned( *oldObj, *newObj, 0 );

  
  // FIXME - what do we do about the identifier?
	
  // Print what happened:
  ATH_MSG( "Converting xAOD::TrackParticleAuxContainer_v3 to current version "
    "[OK]" );
	
  return;
}
	
/// This function should never be called, as we are not supposed to convert
/// object before writing.
///
void xAODTrackParticleAuxContainerCnv_v3::transToPers( const xAOD::TrackParticleAuxContainer*,
                                                       xAOD::TrackParticleAuxContainer_v3*,
                                                       MsgStream& log ) {
	
  log << MSG::ERROR
      << "Somebody called xAODTrackParticleAuxContainerCnv_v3::transToPers"
      << endmsg;
  throw std::runtime_error( "Somebody called xAODTrackParticleAuxContainerCnv_v3::"
    "transToPers" );
  return;
}
