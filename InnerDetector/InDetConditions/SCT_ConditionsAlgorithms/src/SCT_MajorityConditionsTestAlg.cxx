/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
*/

/**
 * @file SCT_MajorityConditionsTestAlg.cxx
 *
 * @brief 
 * 
 *
 * @author gwilliam@mail.cern.ch
 **/

#include "SCT_MajorityConditionsTestAlg.h"

//Athena includes
#include "Identifier/Identifier.h"

SCT_MajorityConditionsTestAlg::SCT_MajorityConditionsTestAlg(const std::string& name, ISvcLocator* pSvcLocator ) : 
  AthAlgorithm( name, pSvcLocator ),
  m_majoritySvc("SCT_MajorityConditionsSvc", name)
{
  declareProperty("MajoritySvc",      m_majoritySvc);
}

//Initialize
StatusCode SCT_MajorityConditionsTestAlg::initialize(){
  ATH_MSG_INFO("Calling initialize");
  
  // Retrieve link masking service
  if (m_majoritySvc.retrieve().isFailure()) {
    ATH_MSG_ERROR("Could not retrieve the link masking service");
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

//Execute
StatusCode SCT_MajorityConditionsTestAlg::execute(){
  ATH_MSG_INFO("Calling execute");

  ATH_MSG_INFO("Detector is " << (m_majoritySvc->isGood()   ? "GOOD" : "BAD"));
  ATH_MSG_INFO("ECC is      " << (m_majoritySvc->isGood(-2) ? "GOOD" : "BAD"));
  ATH_MSG_INFO("Barrel is   " << (m_majoritySvc->isGood(0)  ? "GOOD" : "BAD"));
  ATH_MSG_INFO("ECA is      " << (m_majoritySvc->isGood(2)  ? "GOOD" : "BAD"));

  return StatusCode::SUCCESS;
}


//Finalize
StatusCode SCT_MajorityConditionsTestAlg::finalize(){
  ATH_MSG_INFO("Calling finalize");
  return StatusCode::SUCCESS;
}
