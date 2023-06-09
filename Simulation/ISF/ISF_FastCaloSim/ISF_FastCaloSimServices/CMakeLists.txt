################################################################################
# Package: ISF_FastCaloSimServices
################################################################################

# Declare the package name:
atlas_subdir( ISF_FastCaloSimServices )

# Declare the package's dependencies:
atlas_depends_on_subdirs( PUBLIC
                          GaudiKernel
                          PRIVATE
                          Control/AthenaBaseComps
                          Simulation/Barcode/BarcodeInterfaces
                          Simulation/ISF/ISF_Core/ISF_Interfaces
                          Tracking/TrkEvent/TrkEventPrimitives
                          Tracking/TrkExtrapolation/TrkExInterfaces
                          Calorimeter/CaloInterface
                          Calorimeter/CaloEvent
                          Control/StoreGate
                          Event/NavFourMom
                          Generators/GeneratorObjects
                          Simulation/FastShower/FastCaloSim
                          Simulation/ISF/ISF_Core/ISF_Event
                          Simulation/ISF/ISF_FastCaloSim/ISF_FastCaloSimEvent
                          Simulation/ISF/ISF_FastCaloSim/ISF_FastCaloSimInterfaces )

# External dependencies:
find_package( CLHEP )
find_package( HepMC )

# Component(s) in the package:
atlas_add_component( ISF_FastCaloSimServices
                     src/*.cxx
                     src/components/*.cxx
                     INCLUDE_DIRS ${CLHEP_INCLUDE_DIRS} ${HEPMC_INCLUDE_DIRS}
                     LINK_LIBRARIES ${CLHEP_LIBRARIES} ${HEPMC_LIBRARIES} AthenaBaseComps GaudiKernel ISF_Interfaces TrkEventPrimitives TrkExInterfaces CaloEvent StoreGateLib SGtests NavFourMom GeneratorObjects FastCaloSimLib ISF_Event ISF_FastCaloSimEvent ISF_FastCaloSimInterfaces )

# Install files from the package:
atlas_install_headers( ISF_FastCaloSimServices )
atlas_install_python_modules( python/*.py )

