/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "GaudiKernel/ServiceHandle.h"
#include "StoreGate/StoreGateSvc.h"

#include "TrigConfIO/JsonFileLoader.h"
#include "TrigConfIO/TrigDBL1BunchGroupSetLoader.h"
#include "TrigConfIO/TrigDBMenuLoader.h"
#include "TrigConfInterfaces/IJobOptionsSvc.h"

#include "LVL1ConfigSvc.h"
#include "TrigConfMD5.h"

#include <memory>

TrigConf::LVL1ConfigSvc::LVL1ConfigSvc(const std::string& name, ISvcLocator* pSvcLocator) :
  AthService(name, pSvcLocator)
{}

StatusCode TrigConf::LVL1ConfigSvc::loadRun3StyleMenu()
{
  auto l1menu = std::make_unique<TrigConf::L1Menu>();

  if (m_inputType == "DB") {
    // load l1menu
    TrigConf::TrigDBMenuLoader dbmenuloader(m_dbConnection);
    dbmenuloader.setLevel(TrigConf::MSGTC::WARNING);
    ATH_CHECK( dbmenuloader.loadL1Menu(m_smk, *l1menu) );

  }
  else if (m_inputType == "FILE") {
    // json file menu loader
    TrigConf::JsonFileLoader fileLoader;
    fileLoader.setLevel(TrigConf::MSGTC::WARNING);

    ATH_CHECK( fileLoader.loadFile(m_l1FileName, *l1menu) );

    uint32_t smk =  m_smk.value();
    if (!m_hltFileName.empty() && smk == 0u) {
      auto hltmenu = std::make_unique<TrigConf::HLTMenu>();
      const bool status = fileLoader.loadFile(m_hltFileName, *hltmenu);
      if (status) {
        smk = TrigConf::truncatedHash(*l1menu, *hltmenu);
      } else {
        ATH_MSG_DEBUG("No HLT menu created, cannot compute a MC-SMK in this job");
      }
    }
    ATH_MSG_INFO("Setting file-loaded L1 Menu SMK to:" << smk);
    l1menu->setSMK(smk); // allow assigning a specified or hashed SMK when running from FILE


  }
  else {
    ATH_MSG_ERROR("Unknown input type '" << m_inputType
                  << "'. Allowed values: " << m_inputType.documentation());
    return StatusCode::FAILURE;
  }

  ServiceHandle<StoreGateSvc> detStore("StoreGateSvc/DetectorStore", name());
  ATH_CHECK(detStore.retrieve());
  if (detStore->record(std::move(l1menu), "L1TriggerMenu").isSuccess()) {
    ATH_MSG_INFO("Recorded L1 menu as 'L1TriggerMenu' in detector store");
  }

  return StatusCode::SUCCESS;
}

StatusCode TrigConf::LVL1ConfigSvc::initialize()
{
  // Handle to JobOptionsSvc to retrieve configuration keys
  if (auto joSvc = serviceLocator()->service<TrigConf::IJobOptionsSvc>("JobOptionsSvc")) {
    if (joSvc->superMasterKey() > 0) {
      m_inputType = "DB";
      m_smk = joSvc->superMasterKey();
      m_dbConnection = joSvc->server();
    }
  }
  else {
    ATH_MSG_DEBUG("Did not locate TrigConf::JobOptionsSvc, not running athenaHLT");
  }

  ATH_MSG_INFO("Loading L1 trigger menu from:");
  ATH_MSG_INFO(m_inputType);
  if (m_inputType == "FILE") {
    ATH_MSG_INFO(m_l1FileName);
  }
  else if (m_inputType == "DB") {
    ATH_MSG_INFO(m_dbConnection);
    ATH_MSG_INFO(m_smk);
  }

  ATH_CHECK(loadRun3StyleMenu());

  return StatusCode::SUCCESS;
}
