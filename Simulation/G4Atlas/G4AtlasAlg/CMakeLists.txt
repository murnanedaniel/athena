################################################################################
# Package: G4AtlasAlg
################################################################################

# Declare the package name:
atlas_subdir( G4AtlasAlg )

# Declare the package's dependencies:
atlas_depends_on_subdirs( PUBLIC
                          GaudiKernel
                          PRIVATE
                          Control/AthenaBaseComps
                          Control/AthenaKernel
                          Control/CxxUtils
                          Control/SGTools
                          Control/StoreGate
                          DetectorDescription/GeoModel/GeoModelInterfaces
                          Event/EventInfo
                          Generators/GeneratorObjects
                          Simulation/G4Atlas/G4AtlasInterfaces
                          Simulation/G4Atlas/G4AtlasTools
                          Simulation/G4Sim/FADS/FadsKinematics
                          Simulation/G4Sim/MCTruth
                          Simulation/G4Sim/MCTruthBase
                          Simulation/G4Sim/SimHelpers )

# External dependencies:
find_package( CLHEP )
find_package( Geant4 )
find_package( HepMC )
find_package( XercesC )

# Component(s) in the package:
atlas_add_component( G4AtlasAlg
                     src/*.cxx
                     src/components/*.cxx
                     INCLUDE_DIRS ${GEANT4_INCLUDE_DIRS} ${XERCESC_INCLUDE_DIRS} ${CLHEP_INCLUDE_DIRS} ${HEPMC_INCLUDE_DIRS}
                     LINK_LIBRARIES ${GEANT4_LIBRARIES} ${XERCESC_LIBRARIES} ${CLHEP_LIBRARIES} ${HEPMC_LIBRARIES} AthenaBaseComps AthenaKernel GaudiKernel G4AtlasInterfaces G4AtlasToolsLib FadsKinematics CxxUtils SGTools StoreGateLib SGtests EventInfo GeneratorObjects MCTruth MCTruthBaseLib SimHelpers )

# Install files from the package:
atlas_install_python_modules( python/*.py )

