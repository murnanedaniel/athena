//        Copyright Iowa State University 2017.
//                  Author: Nils Krumnack
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

// Please feel free to contact me (nils.erik.krumnack@cern.ch) for bug
// reports, feature suggestions, praise and complaints.


//
// includes
//

#include <AnaAlgorithm/AnaAlgorithm.h>
#include <AnaAlgorithm/AnaAlgorithmConfig.h>

#include <EventLoopTest/UnitTestAlg2.h>
#include <AsgTools/ToolHandle.h>
#include <AsgTools/MessageCheck.h>
#include <AsgTools/UnitTest.h>
#include <cmath>
#include <gtest/gtest.h>

#ifdef ROOTCORE
#include <xAODRootAccess/Init.h>
#endif

#include <RootCoreUtils/Assert.h>

using namespace EL;

//
// unit test
//

TEST (AnaAlgorithmTest, create_basic)
{
  AnaAlgorithmConfig config;
  config.setName ("name");
  config.setType ("EL::AnaAlgorithm");
  std::vector<std::shared_ptr<void> > cleanup;
  std::unique_ptr<AnaAlgorithm> alg;
  ASSERT_SUCCESS (config.makeAlgorithm (alg, cleanup));
  ASSERT_NE (nullptr, alg.get());
  ASSERT_EQ ("name", alg->name());
}

TEST (AnaAlgorithmTest, newAlg)
{
  std::unique_ptr<UnitTestAlg2> alg (new UnitTestAlg2 ("name", nullptr));
}

TEST (AnaAlgorithmTest, create)
{
  AnaAlgorithmConfig config;
  config.setName ("name");
  config.setType ("EL::UnitTestAlg2");
  std::vector<std::shared_ptr<void> > cleanup;
  std::unique_ptr<AnaAlgorithm> alg;
  ASSERT_SUCCESS (config.makeAlgorithm (alg, cleanup));
  ASSERT_NE (nullptr, alg.get());
  ASSERT_EQ ("name", alg->name());
}

TEST (AnaAlgorithmTest, setProperty_string)
{
  AnaAlgorithmConfig config;
  config.setName ("name");
  config.setType ("EL::UnitTestAlg2");
  ASSERT_SUCCESS (config.setProperty ("string_property", "42"));
  std::unique_ptr<AnaAlgorithm> alg;
  std::vector<std::shared_ptr<void> > cleanup;
  ASSERT_SUCCESS (config.makeAlgorithm (alg, cleanup));
  UnitTestAlg2 *myalg = dynamic_cast<UnitTestAlg2*>(alg.get());
  ASSERT_NE (nullptr, myalg);
  ASSERT_EQ ("42", myalg->m_string_property);
}

TEST (AnaAlgorithmTest, setProperty)
{
  AnaAlgorithmConfig config;
  config.setName ("name");
  config.setType ("EL::UnitTestAlg2");
  ASSERT_SUCCESS (config.setProperty ("property", 42));
  std::vector<std::shared_ptr<void> > cleanup;
  std::unique_ptr<AnaAlgorithm> alg;
  ASSERT_SUCCESS (config.makeAlgorithm (alg, cleanup));
  UnitTestAlg2 *myalg = dynamic_cast<UnitTestAlg2*>(alg.get());
  ASSERT_NE (nullptr, myalg);
  ASSERT_EQ (42, myalg->m_property);
}

TEST (AnaAlgorithmTest, setSubTool)
{
  AnaAlgorithmConfig config;
  config.setName ("name");
  config.setType ("EL::UnitTestAlg2");
  ASSERT_SUCCESS (config.createPrivateTool ("toolHandle", "EL::UnitTestTool"));
  ASSERT_SUCCESS (config.setProperty ("toolHandle.propertyInt", 17));
  std::vector<std::shared_ptr<void> > cleanup;
  std::unique_ptr<AnaAlgorithm> alg;
  ASSERT_SUCCESS (config.makeAlgorithm (alg, cleanup));
  UnitTestAlg2 *myalg = dynamic_cast<UnitTestAlg2*>(alg.get());
  ASSERT_NE (nullptr, myalg);
  ASSERT_NE (nullptr, &*myalg->m_toolHandle);
  ASSERT_EQ (17, myalg->m_toolHandle->getPropertyInt());
  ASSERT_EQ (nullptr, myalg->m_toolHandle->getSubtool());
}

TEST (AnaAlgorithmTest, setSubSubTool)
{
  AnaAlgorithmConfig config;
  config.setName ("name");
  config.setType ("EL::UnitTestAlg2");
  ASSERT_SUCCESS (config.createPrivateTool ("toolHandle", "EL::UnitTestTool"));
  ASSERT_SUCCESS (config.createPrivateTool ("toolHandle.subtool", "EL::UnitTestTool"));
  ASSERT_SUCCESS (config.setProperty ("toolHandle.subtool.propertyInt", 17));
  std::vector<std::shared_ptr<void> > cleanup;
  std::unique_ptr<AnaAlgorithm> alg;
  ASSERT_SUCCESS (config.makeAlgorithm (alg, cleanup));
  UnitTestAlg2 *myalg = dynamic_cast<UnitTestAlg2*>(alg.get());
  ASSERT_NE (nullptr, myalg);
  ASSERT_NE (nullptr, &*myalg->m_toolHandle);
  ASSERT_EQ (0, myalg->m_toolHandle->getPropertyInt());
  ASSERT_NE (nullptr, myalg->m_toolHandle->getSubtool());
  ASSERT_EQ (17, myalg->m_toolHandle->getSubtool()->getPropertyInt());
}

ATLAS_GOOGLE_TEST_MAIN
