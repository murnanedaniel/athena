/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

/**
 * @file StoreDump.cxx
 * @brief Implementation of class StoreDump
 */
  
#include <string>
#include "StoreDump.h"
#include "AthenaPoolUtilities/AthenaAttributeList.h"
    
StoreDump::StoreDump(const std::string& name, 
			ISvcLocator* pSvcLocator) :
              AthAlgorithm(name, pSvcLocator)
{
   declareProperty("StoreName", m_storename="DetectorStore");
   declareProperty("Key", m_key="*");
}

StoreDump::~StoreDump()
{}

StatusCode StoreDump::initialize() { return StatusCode::SUCCESS; }

StatusCode StoreDump::execute() 
{
   StoreGateSvc* detStore = 0;
   ATH_CHECK( service(m_storename,detStore) );
   ATH_MSG_INFO( detStore->dump()  );

   AthenaAttributeList* test = 0;
   ATH_CHECK( detStore->retrieve(test,m_key) );
   ATH_MSG_INFO( "AttributeList has size " << test->size()  );

   return StatusCode::SUCCESS; 
}

StatusCode StoreDump::finalize() { return StatusCode::SUCCESS; }
