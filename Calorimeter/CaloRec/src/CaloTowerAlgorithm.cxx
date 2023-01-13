/*
  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
*/


#include "AthenaKernel/errorcheck.h"

#include "CaloEvent/CaloTowerSeg.h"
#include "CaloEvent/CaloTowerContainer.h"
#include "CaloTowerAlgorithm.h"

#include <string>
#include <vector>
#include <iomanip>

CaloTowerAlgorithm::CaloTowerAlgorithm(const std::string& name, 
				       ISvcLocator* pSvcLocator) 
  : AthReentrantAlgorithm(name,pSvcLocator)
  , m_nEtaTowers(50)
  , m_nPhiTowers(64)
  , m_minEta(-2.5)
  , m_maxEta(2.5)
  , m_genericLink(true) 
  , m_ptools( this )
  , m_towerContainerKey("")
{
  // tool names
  //declareProperty("TowerBuilderTools",m_toolNames);
  declareProperty("TowerBuilderTools",m_ptools);
  // output data
  declareProperty("TowerContainerName",m_towerContainerKey);
  // tower grid
  declareProperty("NumberOfEtaTowers",m_nEtaTowers);
  declareProperty("NumberOfPhiTowers",m_nPhiTowers);
  declareProperty("EtaMin",m_minEta);
  declareProperty("EtaMax",m_maxEta);
  // linkable
  declareProperty("GenericLinked",m_genericLink);
}

CaloTowerAlgorithm::~CaloTowerAlgorithm()
= default;

////////////////
// Initialize //
////////////////

StatusCode CaloTowerAlgorithm::initialize()
{
  ATH_CHECK(m_towerContainerKey.initialize());
  ////////////////////
  // Allocate Tools //
  ////////////////////

  // check tool names
  if (m_ptools.empty()) {
    ATH_MSG_ERROR(" no tools given for this algorithm.");
    return StatusCode::FAILURE;
  }

  // find tools


  CaloTowerSeg theTowerSeg(m_nEtaTowers,m_nPhiTowers,m_minEta,m_maxEta);


  unsigned int toolCtr = 0;
  ATH_MSG_INFO(" ");
  ATH_MSG_INFO("List of tools in execution sequence:");
  ATH_MSG_INFO("------------------------------------");

  ATH_CHECK(m_ptools.retrieve());

  for (ToolHandle<ICaloTowerBuilderToolBase>& tool : m_ptools) {
    toolCtr++;

    ATH_MSG_INFO(std::setw(2) << toolCtr << ".) " << tool->type()
        << "::name() = \042" << tool->name() << "\042");

    ATH_MSG_INFO("------------------------------------");
    ATH_MSG_INFO(" ");

    ATH_MSG_DEBUG(" set correct tower seg for this tool "
        << tool->name());

    tool->setTowerSeg(theTowerSeg);

//    if (tool->initializeTool().isFailure()) {
//      ATH_MSG_WARNING(" Tool failed to initialize");
//    }

  } //close iteration over tools
  return StatusCode::SUCCESS;
}

/////////////
// Execute //
/////////////

StatusCode CaloTowerAlgorithm::execute (const EventContext& ctx) const
{

  //////////////////////////
  // Re-allocate Services //
  //////////////////////////

  //timing
  IChronoStatSvc* theTicker = chronoSvc();

  /////////////////////
  // Tool Processing //
  /////////////////////

  CaloTowerSeg theTowerSeg(m_nEtaTowers,m_nPhiTowers,m_minEta,m_maxEta);

  SG::WriteHandle<CaloTowerContainer> theTowers(m_towerContainerKey, ctx);
  ATH_CHECK( theTowers.record(std::make_unique<CaloTowerContainer>(theTowerSeg)) );
  

  ToolHandleArray<ICaloTowerBuilderToolBase>::const_iterator firstITool  = m_ptools.begin();
  ToolHandleArray<ICaloTowerBuilderToolBase>::const_iterator lastITool   = m_ptools.end();
  StatusCode processStatus = StatusCode::SUCCESS;
  //
  // loop stops only when Failure indicated by one of the tools
  //

  ATH_MSG_DEBUG("In execute() ");
  
  while (!processStatus.isFailure() && firstITool != lastITool) {

    if (theTicker != nullptr) {
      theTicker->chronoStart((*firstITool)->name());
    }

    processStatus = (*firstITool)->execute(ctx, theTowers.ptr());

    if (theTicker != nullptr) {
      theTicker->chronoStop((*firstITool)->name());
    }
    if (!processStatus.isFailure()) {
      ATH_MSG_DEBUG((*firstITool)->name()
          << ": CaloTowerContainer::size() = " << theTowers->size());

      ++firstITool;
    } else {
      // some problem - but do not skip event loop!
      ATH_MSG_ERROR("problems while or after processing tool \042"
          << (*firstITool)->name()
          << "\042 - cross-check CaloTowerContainer::size() = "
          << theTowers->size());

      ++firstITool;
    }
  }

  return StatusCode::SUCCESS;
}

//////////////
// Finalize //
//////////////

StatusCode CaloTowerAlgorithm::finalize()
{
  return StatusCode::SUCCESS;
}
