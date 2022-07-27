/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/
/*
 */
/**
 * @file IOVDbSvc/test/IOVDbFolder_test.cxx
 * @author Shaun Roe
 * @date Feb, 2019
 * @brief Some tests for IOVDbFolder interface in the Boost framework
 */
 
 
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE TEST_IOVDBSVC

#include <boost/test/unit_test.hpp>
//
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/IMessageSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IAddressCreator.h"
#include "GaudiKernel/IOpaqueAddress.h"
//
#include "../src/FolderTypes.h"
#include "../src/IOVDbParser.h"
#include "../src/IOVDbConn.h"
#include "../src/IOVDbFolder.h"
//
#include "GaudiKernelFixtureBase.h"
#include "TestFolderFixture.h"
//
#include "CxxUtils/checker_macros.h"
ATLAS_NO_CHECK_FILE_THREAD_SAFETY;

struct GaudiKernelFixture:public GaudiKernelFixtureBase{
  GaudiKernelFixture():GaudiKernelFixtureBase(){
    //nop, everything in base.
  }
};

struct TestFolderFixture:public TestFolderFixtureBase{
  TestFolderFixture():TestFolderFixtureBase("IOVDbFolderTest.db"){
    //nop, everything in base.
  }
};


struct IOVDbConnFixture{
  ServiceHandle<IMessageSvc> msgSvc;
  const std::string connectionString;
  MsgStream log;
  IOVDbConn connection;
  IOVDbConnFixture():msgSvc("msgSvc","test"),
   connectionString("sqlite://;schema=IOVDbFolderTest.db;dbname=OFLP200"),
   log(msgSvc.get(), "IOVDbFolder_test"),
   connection(connectionString, true, log){
  }
};

struct IOVDbParserFixture{
  ServiceHandle<IMessageSvc> msgSvc;
  const std::string descriptionString;
  MsgStream log;
  IOVDbParser parser;
  IOVDbParserFixture():msgSvc("msgSvc","test"),
   descriptionString{"/key1<timeStamp>run-lumi</timeStamp><addrHeader><address_header service_type=\"71\" clid=\"1238547719\" /></addrHeader><typeName>CondAttrListCollection</typeName>"},
   log(msgSvc.get(), "IOVDbFolder_test"),
   parser(descriptionString, log){
  }
};

BOOST_FIXTURE_TEST_SUITE(IOVDbFolderTest , GaudiKernelFixture)
  //pre-requisites
  IOVDbConnFixture connectionFixture;
  IOVDbParserFixture parserFixture;
  //can use the parser msgstream?
  //need IClassIDSvc
  ServiceHandle<IClassIDSvc> clidSvc("ClassIDSvc","test");
  //tests construction
  BOOST_AUTO_TEST_CASE(IOVDbFolderConstruction){
    //note:default construction is not explicitly deleted; perhaps it should be.
    BOOST_CHECK_NO_THROW(IOVDbFolder f(&(connectionFixture.connection), parserFixture.parser, parserFixture.log, clidSvc.get(), nullptr, false));
  }
  BOOST_FIXTURE_TEST_SUITE(IOVDbFolderMethods, TestFolderFixture)
    BOOST_AUTO_TEST_CASE(PublicMethods){
      //preload tests
      IOVDbConn connection("sqlite://;schema=IOVDbFolderTest.db;dbname=OFLP200", true, parserFixture.log);
      IOVDbFolder iovDbFolder(&connection, parserFixture.parser, parserFixture.log, clidSvc.get(), nullptr, false, true);
      BOOST_TEST_CHECKPOINT("After instantiation, but before any loading method call");
      BOOST_TEST(iovDbFolder.folderName() == "/key1");
      BOOST_TEST(iovDbFolder.key() == "/key1");
      BOOST_TEST(iovDbFolder.conn() != nullptr);
      BOOST_TEST(iovDbFolder.multiVersion() == false);
      const bool isEpochTimestamp(true);
      BOOST_TEST(iovDbFolder.timeStamp() != isEpochTimestamp);//before looking
      BOOST_TEST(iovDbFolder.tagOverride() == false);
      BOOST_TEST(iovDbFolder.retrieved() == false);
      BOOST_TEST(iovDbFolder.noOverride() == false);
      BOOST_TEST(iovDbFolder.folderType() == IOVDbNamespace::AttrList);
      BOOST_TEST(iovDbFolder.readMeta() == false);
      BOOST_TEST(iovDbFolder.writeMeta() == false);
      BOOST_TEST(iovDbFolder.fromMetaDataOnly() == false);
      BOOST_TEST(iovDbFolder.extensible() == false);
      BOOST_TEST(iovDbFolder.dropped() == false);
      BOOST_TEST(iovDbFolder.joTag().empty());
      BOOST_TEST(iovDbFolder.resolvedTag().empty());
      BOOST_TEST(iovDbFolder.eventStore() == "StoreGateSvc");
      BOOST_TEST(iovDbFolder.clid() == 0);
      BOOST_TEST(iovDbFolder.bytesRead() == 0);
      BOOST_TEST(iovDbFolder.readTime() == 0.0);
      BOOST_TEST(iovDbFolder.cacheValid(cool::ValidityKey(600)) == false );
      BOOST_TEST_CHECKPOINT("After preLoadFolder method call");
      const std::string tag("");
      auto addr=iovDbFolder.preLoadFolder(tagInfoMgr.get(),0,600);
      BOOST_TEST(iovDbFolder.loadCache(cool::ValidityKey(50),600, tag, true) == true);
      BOOST_TEST_MESSAGE("After loadCache method call...");
      BOOST_TEST(iovDbFolder.timeStamp() == isEpochTimestamp);//after looking
      BOOST_TEST(iovDbFolder.retrieved() == false);//only changed by getAddress method
      BOOST_TEST(iovDbFolder.bytesRead() == 8);
      BOOST_TEST(iovDbFolder.readTime() > 0.0);
      BOOST_TEST(iovDbFolder.clid() == 1238547719);
      auto theValidityKey=iovDbFolder.iovTime(IOVTime(cool::ValidityKey(600)));
      BOOST_TEST(theValidityKey == cool::ValidityKey(600));
      IOVRange zeroRange{IOVTime(0),IOVTime(0)};
      IOVRange returnRange{zeroRange};
      BOOST_TEST(iovDbFolder.currentRange() ==  zeroRange);//why?
      BOOST_TEST(iovDbFolder.cacheValid(cool::ValidityKey(600)) );
      std::unique_ptr<IOpaqueAddress> returnedAddress;
      ServiceHandle<IAddressCreator> persistencySvc("EventPersistencySvc", "test");
      bool poolRequested{};
      BOOST_TEST( iovDbFolder.getAddress(cool::ValidityKey(600), &(*persistencySvc),0,returnedAddress,returnRange, poolRequested));
      BOOST_TEST_CHECKPOINT("After getAddress method call");
      BOOST_TEST(iovDbFolder.retrieved() == true);
      IOVRange actualRange{IOVTime(100),IOVTime(cool::ValidityKeyMax)};
      BOOST_TEST(returnRange == actualRange);
      BOOST_CHECK_NO_THROW(iovDbFolder.resetCache());
    }
  BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
