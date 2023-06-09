################################################################################
# Package: AthenaROOTAccess
################################################################################

# Declare the package name:
atlas_subdir( AthenaROOTAccess )

# Declare the package's dependencies:
atlas_depends_on_subdirs( PUBLIC
                          Control/AthContainersInterfaces
                          Control/AthenaKernel
                          Control/CxxUtils
                          Control/SGTools
                          GaudiKernel
                          PRIVATE
                          AtlasTest/TestTools
                          Control/AthAllocators
                          Control/AthContainers
                          Control/AthLinks
                          Control/RootUtils
                          Database/AthenaPOOL/AthenaPoolCnvSvc
                          Database/PersistentDataModel
                          Database/PersistentDataModelTPCnv
                          Event/xAOD/xAODHIEvent
                          Event/xAOD/xAODMuon
                          Event/xAOD/xAODTracking )

# External dependencies:
find_package( ROOT COMPONENTS PyROOT Core Tree MathCore Hist RIO pthread )

# Component(s) in the package:
atlas_add_library( AthenaROOTAccess
                   src/*.cxx
                   PUBLIC_HEADERS AthenaROOTAccess
                   PRIVATE_INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                   LINK_LIBRARIES AthenaKernel CxxUtils SGTools GaudiKernel RootUtilsPyROOT AthenaPoolCnvSvcLib
                   PRIVATE_LINK_LIBRARIES ${ROOT_LIBRARIES} TestTools AthAllocators AthContainers AthLinks RootUtils PersistentDataModel PersistentDataModelTPCnv xAODHIEvent xAODMuon xAODTracking )

atlas_add_dictionary( AthenaROOTAccessDict
                      AthenaROOTAccess/AthenaROOTAccessDict.h
                      AthenaROOTAccess/selection.xml
                      INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                      LINK_LIBRARIES ${ROOT_LIBRARIES} AthenaKernel CxxUtils SGTools GaudiKernel TestTools AthAllocators AthContainers AthLinks RootUtils RootUtilsPyROOT AthenaPoolCnvSvcLib PersistentDataModel PersistentDataModelTPCnv xAODHIEvent xAODMuon xAODTracking AthenaROOTAccess )

atlas_add_dictionary( AthenaROOTAccessTestDict
                      AthenaROOTAccess/AthenaROOTAccessTestDict.h
                      AthenaROOTAccess/test_selection.xml
                      INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                      LINK_LIBRARIES ${ROOT_LIBRARIES} AthenaKernel CxxUtils SGTools GaudiKernel TestTools AthAllocators AthContainers AthLinks RootUtils RootUtilsPyROOT AthenaPoolCnvSvcLib PersistentDataModel PersistentDataModelTPCnv xAODHIEvent xAODMuon xAODTracking AthenaROOTAccess )

atlas_add_test( DVCollectionProxy_test
                SOURCES
                test/DVCollectionProxy_test.cxx
                INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                LINK_LIBRARIES ${ROOT_LIBRARIES} AthenaKernel CxxUtils SGTools GaudiKernel TestTools AthAllocators AthContainers AthLinks RootUtils RootUtilsPyROOT AthenaPoolCnvSvcLib PersistentDataModel PersistentDataModelTPCnv xAODHIEvent xAODMuon xAODTracking AthenaROOTAccess
                EXTRA_PATTERNS "is already in" )

atlas_add_test( DataBucketVoid_test
                SOURCES
                test/DataBucketVoid_test.cxx
                INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                LINK_LIBRARIES ${ROOT_LIBRARIES} AthenaKernel CxxUtils SGTools GaudiKernel TestTools AthAllocators AthContainers AthLinks RootUtils RootUtilsPyROOT AthenaPoolCnvSvcLib PersistentDataModel PersistentDataModelTPCnv xAODHIEvent xAODMuon xAODTracking AthenaROOTAccess
                EXTRA_PATTERNS "is already in" )

# Install files from the package:
atlas_install_python_modules( python/*.py )
atlas_install_joboptions( share/*.py )

