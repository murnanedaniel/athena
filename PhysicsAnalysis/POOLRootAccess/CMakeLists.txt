################################################################################
# Package: POOLRootAccess
################################################################################

# Declare the package name:
atlas_subdir( POOLRootAccess )

# Declare the package's dependencies:
atlas_depends_on_subdirs( PUBLIC
                          Control/AthenaKernel
                          Control/StoreGate
                          GaudiKernel
			  Control/StoreGateBindings
                          PRIVATE
                          Control/xAODRootAccess
			  Control/AthAnalysisBaseComps
                          Event/xAOD/xAODBase
                          Event/xAOD/xAODEventInfo )

# External dependencies:
find_package( ROOT COMPONENTS Core Tree MathCore Hist RIO pthread )

# Component(s) in the package:
atlas_add_library( POOLRootAccess
                   src/*.cxx
                   PUBLIC_HEADERS POOLRootAccess
                   PRIVATE_INCLUDE_DIRS ${ROOT_INCLUDE_DIRS} 
		   LINK_LIBRARIES ${ROOT_LIBRARIES} xAODRootAccess AthAnalysisBaseCompsLib StoreGateLib 
                   PRIVATE_LINK_LIBRARIES ${ROOT_LIBRARIES} GaudiKernel SGtests StoreGateBindings  )

atlas_add_dictionary( POOLRootAccessDict
                      POOLRootAccess/POOLRootAccessDict.h
                      POOLRootAccess/selection.xml
                      INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                      LINK_LIBRARIES ${ROOT_LIBRARIES} POOLRootAccess )

atlas_add_executable( ut_basicRead_test
                      test/ut_basicRead_test.cxx
                      INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                      LINK_LIBRARIES ${ROOT_LIBRARIES} POOLRootAccess xAODEventInfo )

atlas_add_executable( ut_basicxAODRead_test
                      test/ut_basicxAODRead_test.cxx
                      INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                      LINK_LIBRARIES ${ROOT_LIBRARIES} POOLRootAccess xAODEventInfo )

# Install files from the package:
atlas_install_joboptions( share/*.opts )

