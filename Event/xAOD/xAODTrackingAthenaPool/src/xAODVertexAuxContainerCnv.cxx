/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// $Id: xAODVertexAuxContainerCnv.cxx 635799 2014-12-12 22:08:09Z ssnyder $

// System include(s):
#include <exception>

// Local include(s):
#include "xAODVertexAuxContainerCnv.h"
#include "AthContainers/tools/copyThinned.h"
#include "AthenaKernel/IThinningSvc.h"

xAODVertexAuxContainerCnv::
xAODVertexAuxContainerCnv( ISvcLocator* svcLoc )
   : xAODVertexAuxContainerCnvBase( svcLoc ) {

}

xAOD::VertexAuxContainer*
xAODVertexAuxContainerCnv::
createPersistent( xAOD::VertexAuxContainer* trans ) {

   // Create a copy of the container:
   return SG::copyThinned (*trans, IThinningSvc::instance());
}

xAOD::VertexAuxContainer*
xAODVertexAuxContainerCnv::createTransient() {

   // The known ID(s) for this container:
   static const pool::Guid v1_guid( "B1F73A82-9B4E-4508-8EB0-EF7D6E05BA57" );

   // Check which version of the container we're reading:
   if( compareClassGuid( v1_guid ) ) {
      // It's the latest version, read it directly:
      return poolReadObject< xAOD::VertexAuxContainer >();
   }

   // If we didn't recognise the ID:
   throw std::runtime_error( "Unsupported version of "
                             "xAOD::VertexAuxContainer found" );
   return 0;
}
