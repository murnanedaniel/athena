/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#include "AthenaKernel/errorcheck.h"

// AFP_Raw2Digi includes
#include "AFP_Raw2Digi/AFP_Raw2Digi.h"
#include "xAODForward/AFPSiHit.h"
#include "xAODForward/AFPSiHitContainer.h"
#include "xAODForward/AFPSiHitAuxContainer.h"

AFP_Raw2Digi::AFP_Raw2Digi(const std::string &name, ISvcLocator *pSvcLocator)
  : AthReentrantAlgorithm(name, pSvcLocator)
{
}

AFP_Raw2Digi::~AFP_Raw2Digi() {}

StatusCode AFP_Raw2Digi::initialize() {
  ATH_MSG_INFO("Initializing " << name() << "...");

  CHECK( m_DigiTool.retrieve() );

  return StatusCode::SUCCESS;
}

StatusCode AFP_Raw2Digi::finalize() {
  ATH_MSG_INFO("Finalizing " << name() << "...");

  return StatusCode::SUCCESS;
}

StatusCode AFP_Raw2Digi::execute(const EventContext &ctx) const{

  ATH_MSG_DEBUG("Executing " << name() << "...");

  CHECK (m_DigiTool->recoAll(ctx) );

  return StatusCode::SUCCESS;
}
