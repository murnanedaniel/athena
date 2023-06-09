################################################################################
# Package: PileupReweighting
################################################################################

# Declare the package name:
atlas_subdir( PileupReweighting )

# Declare the package's dependencies:
atlas_depends_on_subdirs( PUBLIC
                          Control/AthToolSupport/AsgTools
                          DataQuality/GoodRunsLists
                          Event/xAOD/xAODEventInfo
                          GaudiKernel
                          PhysicsAnalysis/AnalysisCommon/PATInterfaces
                          PRIVATE
                          AtlasTest/TestTools
                          Control/AthenaBaseComps
                          Tools/PathResolver )

# External dependencies:
find_package( ROOT COMPONENTS MathCore RIO Cint Core Tree Hist pthread MathMore Minuit Minuit2 Matrix Physics HistPainter Rint )

# tag ROOTBasicLibs was not recognized in automatic conversion in cmt2cmake

# tag ROOTCintexLibs was not recognized in automatic conversion in cmt2cmake

# Component(s) in the package:
atlas_add_root_dictionary( PileupReweightingLib
                           PileupReweightingLibDictSource
                           ROOT_HEADERS PileupReweighting/TPileupReweighting.h Root/LinkDef.h
                           EXTERNAL_PACKAGES ROOT )

atlas_add_library( PileupReweightingLib
                   Root/*.cxx
                   ${PileupReweightingLibDictSource}
                   PUBLIC_HEADERS PileupReweighting
                   PRIVATE_INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                   LINK_LIBRARIES AsgTools xAODEventInfo GaudiKernel PATInterfaces GoodRunsListsLib
                   PRIVATE_LINK_LIBRARIES ${ROOT_LIBRARIES} TestTools AthenaBaseComps PathResolver )

atlas_add_component( PileupReweighting
                     src/components/*.cxx
                     INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                     LINK_LIBRARIES ${ROOT_LIBRARIES} AsgTools GoodRunsListsLib xAODEventInfo GaudiKernel PATInterfaces TestTools AthenaBaseComps PathResolver PileupReweightingLib )

atlas_add_dictionary( PileupReweightingDict
                      PileupReweighting/PileupReweightingDict.h
                      PileupReweighting/selection.xml
                      INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                      LINK_LIBRARIES ${ROOT_LIBRARIES} AsgTools GoodRunsListsLib xAODEventInfo GaudiKernel PATInterfaces TestTools AthenaBaseComps PathResolver PileupReweightingLib )

atlas_add_executable( testPRW
                      src/testPRW.C
                      INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                      LINK_LIBRARIES ${ROOT_LIBRARIES} AsgTools GoodRunsListsLib xAODEventInfo GaudiKernel PATInterfaces TestTools AthenaBaseComps PathResolver PileupReweightingLib )

atlas_add_executable( testPRWTool
                      src/testPRWTool.C
                      INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                      LINK_LIBRARIES ${ROOT_LIBRARIES} AsgTools GoodRunsListsLib xAODEventInfo GaudiKernel PATInterfaces TestTools AthenaBaseComps PathResolver PileupReweightingLib )

atlas_add_executable( checkPRWConfigFile
                      src/checkPRWConfigFile.C
                      INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                      LINK_LIBRARIES ${ROOT_LIBRARIES} AsgTools GoodRunsListsLib xAODEventInfo GaudiKernel PATInterfaces TestTools AthenaBaseComps PathResolver PileupReweightingLib )

atlas_add_test( ut_PRWExample_test
                SOURCES
                test/ut_PRWExample_test.cxx
                INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                LINK_LIBRARIES ${ROOT_LIBRARIES} AsgTools GoodRunsListsLib xAODEventInfo GaudiKernel PATInterfaces TestTools AthenaBaseComps PathResolver PileupReweightingLib )

# Install files from the package:
atlas_install_python_modules( python/*.py )
atlas_install_joboptions( share/*.py )
atlas_install_data( share/*.root )

