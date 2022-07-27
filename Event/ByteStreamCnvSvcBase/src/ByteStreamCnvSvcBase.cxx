/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "ByteStreamCnvSvcBase/ByteStreamCnvSvcBase.h"
#include "ByteStreamCnvSvcBase/ByteStreamAddress.h"

#include "GaudiKernel/IOpaqueAddress.h"
#include "GaudiKernel/GenericAddress.h"
#include "GaudiKernel/IConverter.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/IClassIDSvc.h"

//______________________________________________________________________________
ByteStreamCnvSvcBase::ByteStreamCnvSvcBase(const std::string& name, ISvcLocator* pSvcLocator) :
   ::AthCnvSvc(name, pSvcLocator, ByteStreamAddress::storageType())
{
   declareProperty("InitCnvs", m_initCnvs); 
}
//______________________________________________________________________________
/// Standard Destructor
ByteStreamCnvSvcBase::~ByteStreamCnvSvcBase()   {
}
//______________________________________________________________________________
/// Initialize the service.
StatusCode ByteStreamCnvSvcBase::initialize()     {
   if (!::AthCnvSvc::initialize().isSuccess()) {
      ATH_MSG_FATAL("Cannot initialize AthCnvSvc base class.");
      return(StatusCode::FAILURE);
   }

   ServiceHandle<IIncidentSvc> incsvc("IncidentSvc", this->name());
   if (!incsvc.retrieve().isSuccess()) {
      ATH_MSG_FATAL("Cannot get IncidentSvc.");
      return(StatusCode::FAILURE);
   }
   incsvc->addListener(this, "BeginRun", 0, false, true); // true for singleshot
   return(StatusCode::SUCCESS);
}
//_______________________________________________________________________
StatusCode ByteStreamCnvSvcBase::queryInterface(const InterfaceID& riid, void** ppvInterface) {
   if (IByteStreamEventAccess::interfaceID().versionMatch(riid)) {
      *ppvInterface = dynamic_cast<IByteStreamEventAccess*>(this);
   } else {
      // Interface is not directly available: try out a base class
      return(::AthCnvSvc::queryInterface(riid, ppvInterface));
   }
   addRef();
   return(StatusCode::SUCCESS);
}
//______________________________________________________________________________
StatusCode ByteStreamCnvSvcBase::updateServiceState(IOpaqueAddress* pAddress) {
   if (pAddress != 0) {
      GenericAddress* pAddr = dynamic_cast<GenericAddress*>(pAddress);
      if (pAddr != 0) {
         return(StatusCode::SUCCESS);
      }
   }
   return(StatusCode::FAILURE);
}
//______________________________________________________________________________
void ByteStreamCnvSvcBase::handle(const Incident& /*incident*/) {
   ServiceHandle<IClassIDSvc> clidSvc("ClassIDSvc", name());
   if (!clidSvc.retrieve().isSuccess()) {
      ATH_MSG_ERROR("Cannot get ClassIDSvc.");
      return;
   }
   // Initialize the converters
   for (const std::string& cnv : m_initCnvs) {
      ATH_MSG_DEBUG("Accessing Converter for " << cnv);
      CLID id;
      if (!clidSvc->getIDOfTypeName(cnv, id).isSuccess()) {
         ATH_MSG_WARNING("Cannot get CLID for " << cnv);
      } else {
         IConverter* cnv = converter(id);
         if (cnv == 0) {
	    ATH_MSG_WARNING("Cannot get converter for " << cnv);
         } 
      }
   } 
   return;
}
