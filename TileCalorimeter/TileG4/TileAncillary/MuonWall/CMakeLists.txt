################################################################################
# Package: MuonWall
################################################################################

# Declare the package name:
atlas_subdir( MuonWall )

# Declare the package's dependencies:
atlas_depends_on_subdirs( PUBLIC
                          Tools/AtlasDoxygen
                          PRIVATE
                          Calorimeter/CaloIdentifier
                          Control/CxxUtils
                          Control/StoreGate
                          DetectorDescription/GeoModel/GeoModelInterfaces
                          DetectorDescription/GeoModel/GeoModelKernel
                          GaudiKernel
                          Simulation/G4Atlas/G4AtlasTools
                          Simulation/G4Utilities/Geo2G4
                          TileCalorimeter/TileSimEvent )

# External dependencies:
find_package( CLHEP )
find_package( Geant4 )
find_package( XercesC )

# Component(s) in the package:
atlas_add_component( MuonWall
                     src/*.cxx
                     src/components/*.cxx
                     INCLUDE_DIRS ${GEANT4_INCLUDE_DIRS} ${XERCESC_INCLUDE_DIRS} ${CLHEP_INCLUDE_DIRS}
                     LINK_LIBRARIES ${GEANT4_LIBRARIES} ${XERCESC_LIBRARIES} ${CLHEP_LIBRARIES} CaloIdentifier CxxUtils StoreGateLib SGtests GeoModelKernel GaudiKernel G4AtlasToolsLib Geo2G4Lib TileSimEvent )

# Install files from the package:
atlas_install_python_modules( python/*.py )
atlas_install_runtime( share/muonwall.dtd share/muonwall.xml )

