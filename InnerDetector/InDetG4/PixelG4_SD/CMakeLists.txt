################################################################################
# Package: PixelG4_SD
################################################################################

# Declare the package name:
atlas_subdir( PixelG4_SD )

# Declare the package's dependencies:
atlas_depends_on_subdirs( PRIVATE
                          Control/CxxUtils
                          Control/StoreGate
                          GaudiKernel
                          InnerDetector/InDetSimEvent
                          Simulation/G4Atlas/G4AtlasTools
                          Simulation/G4Sim/MCTruth
                          Tools/AtlasDoxygen )

# External dependencies:
find_package( CLHEP )
find_package( Geant4 )
find_package( XercesC )

# Component(s) in the package:
atlas_add_component( PixelG4_SD
                     src/*.cxx
                     src/components/*.cxx
                     INCLUDE_DIRS ${GEANT4_INCLUDE_DIRS} ${XERCESC_INCLUDE_DIRS} ${CLHEP_INCLUDE_DIRS}
                     LINK_LIBRARIES ${GEANT4_LIBRARIES} ${XERCESC_LIBRARIES} ${CLHEP_LIBRARIES} CxxUtils StoreGateLib SGtests GaudiKernel InDetSimEvent G4AtlasToolsLib MCTruth )

# Install files from the package:
atlas_install_python_modules( python/*.py )

