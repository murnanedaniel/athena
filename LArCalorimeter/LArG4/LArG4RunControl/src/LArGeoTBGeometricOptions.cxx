/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#include "LArG4RunControl/LArGeoTBGeometricOptions.h"

#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/Bootstrap.h"
#include "StoreGate/StoreGateSvc.h"

void LArGeoTBGeometricOptions::saveMe()
{
  IService* pSvc;
  ISvcLocator* svcLocator = Gaudi::svcLocator();
  StatusCode result = svcLocator->service("DetectorStore",pSvc);

  if(result.isSuccess())
  {
    StoreGateSvc* m_detStore = dynamic_cast<StoreGateSvc*>(pSvc);
    result=m_detStore->record(this,"LArGeoTBGeometricOptions");
    if(!result.isSuccess())
      std::cout << "Can not record LArGeoTBGeometricOptions" << std::endl;
  }
}

void LArGeoTBGeometricOptions::printMe()
{
  std::cout << " *** *** This is the object of type LArGeoTBGeometricOptions *** *** \n";
  std::cout << " ** CryoEtaPosition = " << m_CryoEtaPosition << "\n *** *** \n";
  std::cout << " ** CryoPhiPosition = " << m_CryoPhiPosition << "\n *** *** \n";
}
