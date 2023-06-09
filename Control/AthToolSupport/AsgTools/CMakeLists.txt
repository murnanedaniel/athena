# $Id: CMakeLists.txt 745006 2016-05-05 10:40:32Z will $
################################################################################
# Package: AsgTools
################################################################################


# Declare the package name:
atlas_subdir( AsgTools )

# Dependencies are taken based on what environment we are in:
if( XAOD_STANDALONE )
   set( deps Control/xAODRootAccessInterfaces Control/xAODRootAccess )
   set( libs xAODRootAccessInterfaces xAODRootAccess )
else()
   set( deps Control/AthenaBaseComps Control/SGTools Database/IOVDbDataModel
      GaudiKernel )
   set( libs AthenaBaseComps SGTools IOVDbDataModel GaudiKernel )
endif()

# Declare the package's dependencies:
atlas_depends_on_subdirs( PUBLIC ${deps} )

# External dependencies:
find_package( Boost COMPONENTS regex )
find_package( ROOT COMPONENTS Core )
find_package( GTest )

# Decide which sources to use:
if( XAOD_STANDALONE )
   set( sources Root/*.cxx )
else()
   set( sources Root/MsgLevel.cxx Root/AsgTool.cxx Root/AsgMetadataTool.cxx
      Root/AsgMessaging.cxx Root/ToolStore.cxx Root/AnaToolHandle.cxx
      Root/MessageCheck.cxx Root/UnitTest.cxx src/*.cxx )
endif()

# Component(s) in the package:
atlas_add_library( AsgTools
   AsgTools/*.h AsgTools/*.icc ${sources}
   PUBLIC_HEADERS AsgTools
   INCLUDE_DIRS ${GTEST_INCLUDE_DIRS}
   PRIVATE_INCLUDE_DIRS ${Boost_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS}
   LINK_LIBRARIES ${libs}
   PRIVATE_LINK_LIBRARIES ${Boost_LIBRARIES} ${ROOT_LIBRARIES} )

atlas_add_dictionary( AsgToolsDict
   AsgTools/AsgToolsDict.h AsgTools/selection.xml
   LINK_LIBRARIES AsgTools )

# Test(s) in the package:
atlas_add_test( gt_UnitTest_test
   SOURCES test/gt_UnitTest_test.cxx
   INCLUDE_DIRS ${GTEST_INCLUDE_DIRS}
   LINK_LIBRARIES ${GTEST_LIBRARIES} AsgTools
   EXTRA_PATTERNS "\([0-9]+ ms\)|\([0-9]+ ms total\)" )

atlas_add_test( ut_asgtools_statuscode
   SOURCES test/ut_asgtools_statuscode.cxx
   LINK_LIBRARIES AsgTools )

atlas_add_test( ut_asgtools_toolstore
   SOURCES test/ut_asgtools_toolstore.cxx
   LINK_LIBRARIES AsgTools )

if( XAOD_STANDALONE )

   atlas_add_test( gt_asgtools_toolhandle_test
      SOURCES test/gt_asgtools_toolhandle_test.cxx
      INCLUDE_DIRS ${GTEST_INCLUDE_DIRS}
      LINK_LIBRARIES ${GTEST_LIBRARIES} AsgTools
      EXTRA_PATTERNS "\([0-9]+ ms\)|\([0-9]+ ms total\)" )

   atlas_add_test( ut_MessageCheck
      SOURCES test/ut_MessageCheck.cxx
      LINK_LIBRARIES AsgTools )

endif()
