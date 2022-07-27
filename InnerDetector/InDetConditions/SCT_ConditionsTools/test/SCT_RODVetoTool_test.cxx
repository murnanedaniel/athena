/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

/**
 * @file SCT_RODVetoTool_test.cxx
 * Unit tests to test the methods of 
 * SCT_RodVetoTool
 *
 * @author fschenck@cern.ch
 *
 * */

// This code was updated on 2018-02-22 by following
// Simulation/ISF/ISF_HepMC/ISF_HepMC_Tools/test/GenParticleGenericFilter_test.cxx

#include "CxxUtils/checker_macros.h"
ATLAS_NO_CHECK_FILE_THREAD_SAFETY;

// Google Test
#include "gtest/gtest.h"

// Framework
#include "TestTools/initGaudi.h"

//Gaudi includes
#include "GaudiKernel/IAppMgrUI.h"
#include "GaudiKernel/SmartIF.h"
#include "GaudiKernel/IEvtSelector.h"
#include "GaudiKernel/IToolSvc.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/IProperty.h"
#include "GaudiKernel/EventContext.h"

// Tested AthAlgTool
#include "../src/SCT_RODVetoTool.h"

// ATLAS C++
#include "AthenaKernel/ExtendedEventContext.h"
#include "SCT_ConditionsTools/ISCT_ConditionsTool.h"
#include "SCT_ConditionsData/IdentifierSet.h"
#include "InDetIdentifier/SCT_ID.h"
#include "StoreGate/StoreGateSvc.h"
#include "AthenaKernel/DummyRCUSvc.h"

namespace SCT_test {
   EventIDBase timestamp (int t)
   {
      return EventIDBase (EventIDBase::UNDEFNUM,  // run
                          EventIDBase::UNDEFEVT,  // event
                          t);
   }


// Gaudi Test fixture that provides a clean Gaudi environment for
// each individual test case
class GaudiFixture {

protected:
  GaudiFixture() {
    SetUpGaudi();
  }
  ~GaudiFixture() {
    TearDownGaudi();
  }

  void SetUpGaudi() {
    m_appMgr = Gaudi::createApplicationMgr();
    ASSERT_TRUE( m_appMgr!=nullptr );

    m_svcLoc = m_appMgr;
    ASSERT_TRUE( m_svcLoc.isValid() );

    m_svcMgr = m_appMgr;
    ASSERT_TRUE( m_svcMgr.isValid() );

    m_sg = nullptr;
    ASSERT_TRUE( m_svcLoc->service ("StoreGateSvc", m_sg).isSuccess() );

    m_evtSel = m_svcLoc->service("EventSelector");
    ASSERT_TRUE( m_evtSel.isValid() );

    m_propMgr = m_appMgr;
    ASSERT_TRUE( m_propMgr.isValid() );
    ASSERT_TRUE( m_propMgr->setProperty("EvtSel", "EventSelector").isSuccess() );
    ASSERT_TRUE( m_propMgr->setProperty("JobOptionsType", "NONE").isSuccess() );

    m_toolSvc = m_svcLoc->service("ToolSvc");
    ASSERT_TRUE( m_toolSvc.isValid() );

    ASSERT_TRUE( m_appMgr->configure().isSuccess() );
    ASSERT_TRUE( m_appMgr->initialize().isSuccess() );

    // Easy way to set up services
    const ServiceLocatorHelper helper{*m_svcLoc, "HELPER"};

    // Set up a new StoreGateSvc instance named DetectorStore, as SCT_RODVetoTool needs
    IService* i_svcD{helper.service("StoreGateSvc/DetectorStore", true /*quiet*/ , true /*createIf*/)};
    StoreGateSvc* detStore{dynamic_cast<StoreGateSvc*>(i_svcD)};
    if (detStore) {
      if (not detStore->contains<SCT_ID>("SCT_ID")) {
        // Needed for SCT_RODVetoTool
        const SCT_ID* pHelper{new SCT_ID()};
        ASSERT_TRUE( detStore->record(pHelper, "SCT_ID").isSuccess() );
      }
    }
    StoreGateSvc* conditionStore = nullptr;
    ASSERT_TRUE( m_svcLoc->service ("StoreGateSvc/ConditionStore", conditionStore).isSuccess() );

  }

  void TearDownGaudi() {
    ASSERT_TRUE( m_appMgr->finalize().isSuccess() );
    ASSERT_TRUE( m_appMgr->terminate().isSuccess() );
    Gaudi::setInstance(static_cast<IAppMgrUI*>(nullptr));
  }

  // protected member variables for Core Gaudi components
  IAppMgrUI* m_appMgr{nullptr};
  SmartIF<ISvcLocator> m_svcLoc;
  SmartIF<ISvcManager> m_svcMgr;
  SmartIF<IToolSvc> m_toolSvc;
  SmartIF<IEvtSelector> m_evtSel;
  SmartIF<IProperty> m_propMgr;
  StoreGateSvc* m_sg{nullptr};
  Athena_test::DummyRCUSvc m_rcu;
};

class SCT_RODVetoTool_test: public ::testing::Test, public GaudiFixture {

protected:
  virtual void SetUp() override {
    IAlgTool* tool{nullptr};
    ASSERT_TRUE( m_toolSvc->retrieveTool("SCT_RODVetoTool/SCT_RODVetoTool", tool) );
    m_tool = dynamic_cast<SCT_RODVetoTool*>(tool);
  }

  virtual void TearDown() override {
    for (size_t refCount{m_tool->refCount()}; refCount>0; refCount--) {
      ASSERT_TRUE( m_toolSvc->releaseTool(m_tool) );
    }
  }   

  // the tested AthAlgTool
  SCT_RODVetoTool* m_tool{nullptr};
};

TEST_F(SCT_RODVetoTool_test, Initialization) {
  ASSERT_TRUE( m_tool->initialize().isSuccess() );
}
  
TEST_F(SCT_RODVetoTool_test, Finalization) {
  ASSERT_TRUE( m_tool->finalize().isSuccess() );
}

TEST_F(SCT_RODVetoTool_test, queryInterface) {
  //It wants a pointer to a pointer, I give it a pointer to a pointer
  void* b{nullptr};
  void** ppvInterface{&b};
  const InterfaceID rrid{"ISCT_ConditionsTool", 1, 0};
  ASSERT_TRUE( m_tool->queryInterface(rrid, ppvInterface).isSuccess() );
}

TEST_F(SCT_RODVetoTool_test, canReportAbout) {
  ASSERT_TRUE( m_tool->canReportAbout(InDetConditions::DEFAULT) );
}


/*
  This segfaults the tests.
  This test should probably be done as well,
  but to do it you have to create the cabling
  service and properly initialize it with
  rods etc.

TEST_F(SCT_RODVetoTool_test, isGood_Hash) {
  m_tool->initialize();
  IdentifierHash hashId{1234};
  ASSERT_FALSE( m_tool->isGood(hashId) );
}
*/

TEST_F(SCT_RODVetoTool_test, isGood_Id) {
  ASSERT_TRUE( m_tool->initialize() );
  const Identifier elementId{0};
  EventContext ctx{0, 0};
  ctx.setExtension( Atlas::ExtendedEventContext( m_sg, 0 ) );
  EventIDBase eid (1, 0, 0, 0, 20);
  ctx.setEventID (eid);

  StoreGateSvc* conditionStore = nullptr;
  ASSERT_TRUE( m_svcLoc->service ("StoreGateSvc/ConditionStore", conditionStore).isSuccess() );
  std::unique_ptr<IdentifierSet> dummyData{std::make_unique<IdentifierSet>()};

  DataObjID id2 ("BadSCTModuleIds_RODVeto");
  auto cc2 = std::make_unique<CondCont<IdentifierSet> > (m_rcu, id2);
  const EventIDRange range2 (timestamp (0), timestamp (80));
  assert( cc2->insert (range2, std::move(dummyData), ctx).isSuccess() );
  assert( conditionStore->record (std::move (cc2), "BadSCTModuleIds_RODVeto") );

  ASSERT_TRUE( m_tool->isGood(elementId, ctx, InDetConditions::DEFAULT) );
}

} // namespace SCT_test

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
