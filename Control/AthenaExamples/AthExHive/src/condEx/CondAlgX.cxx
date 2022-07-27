/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "CondAlgX.h"
#include "AthExHive/IASCIICondDbSvc.h"

#include "StoreGate/WriteCondHandle.h"

#include "GaudiKernel/ServiceHandle.h"

#include "GaudiKernel/EventIDBase.h"
#include "GaudiKernel/EventIDRange.h"

#include <thread>
#include <chrono>
#include <memory>

CondAlgX::CondAlgX( const std::string& name, 
            ISvcLocator* pSvcLocator ) : 
  ::AthAlgorithm( name, pSvcLocator ),
  m_cds("ASCIICondDbSvc",name)
{
}

CondAlgX::~CondAlgX() {}

StatusCode CondAlgX::initialize() {
  ATH_MSG_DEBUG("initialize " << name());

  ATH_CHECK( m_evt.initialize() );

  ATH_CHECK( m_cds.retrieve() );

  m_wchk.setDbKey(m_dbKey);
  ATH_CHECK( m_wchk.initialize() );

  return StatusCode::SUCCESS;
}

StatusCode CondAlgX::execute() {
  ATH_MSG_DEBUG("execute " << name());

  SG::ReadHandle<xAOD::EventInfo> evt( m_evt );
  if (!evt.isValid()) {
    ATH_MSG_ERROR ("Could not retrieve EventInfo");
    return StatusCode::FAILURE;
  }

  ATH_MSG_DEBUG("   EventInfo:  r: " << evt->runNumber()
                << " e: " << evt->eventNumber() );


  SG::WriteCondHandle<CondDataObj> wch(m_wchk);

  EventIDBase now(getContext().eventID());

  if (evt->eventNumber() == 10) {
      std::this_thread::sleep_for(std::chrono::milliseconds( 500 ));
  }

  // do we have a valid m_wch for current time?
  if ( wch.isValid(now) ) {
    ATH_MSG_DEBUG("CondHandle is already valid for " << now
		  << ". In theory this should not be called, but may happen"
		  << " if multiple concurrent events are being processed out of order.");

  } else {

    ATH_MSG_DEBUG("  CondHandle " << wch.key() 
                  << " not valid. Getting new info for dbKey \"" 
                  << wch.dbKey() << "\" from CondDb");

    EventIDRange r;
    IASCIICondDbSvc::dbData_t val;
    if (m_cds->getRange(wch.dbKey(), getContext(), r, val).isFailure()) {
      ATH_MSG_ERROR("  could not find dbKey \"" << wch.dbKey() 
		    << "\" in CondSvc registry");
      return StatusCode::FAILURE;
    }

    CondDataObj* cdo = new CondDataObj( val );
    if (wch.record(r, cdo).isFailure()) {
      ATH_MSG_ERROR("could not record CondDataObj " << wch.key() 
		    << " = " << *cdo
                    << " with EventRange " << r);
      return StatusCode::FAILURE;
    }
    ATH_MSG_INFO("recorded new CDO " << wch.key() << " = " << *cdo 
		 << " with range " << r);
  }



  return StatusCode::SUCCESS;

}

