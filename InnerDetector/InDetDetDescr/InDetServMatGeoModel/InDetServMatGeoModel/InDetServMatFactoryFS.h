/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef INDETSERVMATGEOMODEL_INDETSERVMATFACTORYFS_H
#define INDETSERVMATGEOMODEL_INDETSERVMATFACTORYFS_H

#include "AthenaBaseComps/AthMessaging.h"
#include "GeoModelKernel/GeoVDetectorFactory.h"
#include "InDetServMatGeoModel/InDetServMatManager.h"

#include "GaudiKernel/ServiceHandle.h"

class StoreGateSvc;
class IRDBAccessSvc;

// InDet service material factory for Frozen Showers (FS)
class InDetServMatFactoryFS : public GeoVDetectorFactory, public AthMessaging  {

 public:
  
  // Constructor:
  InDetServMatFactoryFS(StoreGateSvc  *pDetStore,
			ServiceHandle<IRDBAccessSvc> pRDBAccess);
  
  // Destructor:
  ~InDetServMatFactoryFS();
  
  // Creation of geometry:
  virtual void create(GeoPhysVol *world);
  // manager
  virtual const InDetDD::InDetServMatManager* getDetectorManager () const;

 private:  
  
  // Illegal operations:
  const InDetServMatFactoryFS & operator=(const InDetServMatFactoryFS &right);
  InDetServMatFactoryFS(const InDetServMatFactoryFS &right);

  // private data
  StoreGateSvc                   *m_detStore;
  ServiceHandle<IRDBAccessSvc>    m_rdbAccess;
  InDetDD::InDetServMatManager   *m_manager;
};

#endif
