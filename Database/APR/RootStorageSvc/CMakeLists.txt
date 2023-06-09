################################################################################
# Package: RootStorageSvc
################################################################################

# Declare the package name:
atlas_subdir( RootStorageSvc )

# Declare the package's dependencies:
atlas_depends_on_subdirs( PRIVATE
                          AtlasTest/TestTools
                          Control/AthContainers
                          Control/AthContainersInterfaces
                          Control/CxxUtils
                          Control/RootUtils
                          Database/APR/POOLCore
                          Database/APR/StorageSvc
                          Database/PersistentDataModel
                          GaudiKernel )

# External dependencies:
find_package( ROOT COMPONENTS Core RIO TreePlayer Tree MathCore Hist pthread )

# tag ROOTCintexLibs was not recognized in automatic conversion in cmt2cmake

# Component(s) in the package:
atlas_add_library( RootStorageSvc
                   src/*.cpp
                   NO_PUBLIC_HEADERS
                   INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                   LINK_LIBRARIES ${ROOT_LIBRARIES} RootUtilsPyROOT
                   PRIVATE_LINK_LIBRARIES TestTools AthContainers CxxUtils RootUtils POOLCore StorageSvc PersistentDataModel GaudiKernel )

atlas_add_test( RootAuxVectorFactory_test
                SOURCES
                test/RootAuxVectorFactory_test.cxx
                INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                LINK_LIBRARIES ${ROOT_LIBRARIES} TestTools AthContainers CxxUtils RootUtils RootUtilsPyROOT POOLCore StorageSvc PersistentDataModel GaudiKernel RootStorageSvc
                EXTRA_PATTERNS "is already in" )

# Component list generation:
atlas_generate_componentslist( RootStorageSvc )
