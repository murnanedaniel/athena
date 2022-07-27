/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef INDETSERVMATGEOMODEL_SCT_SERVMATFACTORYDC2_H
#define INDETSERVMATGEOMODEL_SCT_SERVMATFACTORYDC2_H

#include "AthenaBaseComps/AthMessaging.h"
#include "GaudiKernel/ServiceHandle.h"

class StoreGateSvc;
class GeoPhysVol;
class GeoShape;
class GeoMaterial;
class StoredMaterialManager;
class IRDBAccessSvc;

#include <string>

class SCT_ServMatFactoryDC2 : public AthMessaging {

 public:
  
  // Constructor:
  SCT_ServMatFactoryDC2(StoreGateSvc  *pDetStore,
			ServiceHandle<IRDBAccessSvc>& pRDBAccess);  
  // Destructor:
  virtual ~SCT_ServMatFactoryDC2() = default;
  
  // Creation of geometry:
  void create(GeoPhysVol *mother);

 private:  
  const GeoShape* createShape(int volType, 
			      double rmin1, 
			      double rmax1, 
			      double halflength,
			      double rmin2,
			      double rmax2);
  
  const GeoMaterial* createMaterial(const std::string & name,
				    int volType, 
				    double fractionRL,
				    double rmin1, 
				    double rmax1, 
				    double halflength,
				    double rmin2,
				    double rmax2);

 // Illegal operations:
  const SCT_ServMatFactoryDC2 & operator=(const SCT_ServMatFactoryDC2 &right);
  SCT_ServMatFactoryDC2(const SCT_ServMatFactoryDC2 &right);

  // private data
  StoreGateSvc                   *m_detStore;
  ServiceHandle<IRDBAccessSvc>    m_rdbAccess;
  StoredMaterialManager          *m_materialManager;
};

#endif //  INDETSERVMATGEOMODEL_SCT_SERVMATFACTORYDC2_H
