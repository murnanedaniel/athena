################################################################################
# Package: MCTruthSimAlgs
################################################################################

# Declare the package name:
atlas_subdir( MCTruthSimAlgs )

# Declare the package's dependencies:
atlas_depends_on_subdirs( PRIVATE
                          Calorimeter/CaloSimEvent
                          Control/AthenaBaseComps
                          Control/AthenaKernel
                          Control/CxxUtils
                          Control/PileUpTools
                          Control/StoreGate
                          Event/xAOD/xAODJet
                          GaudiKernel
                          Generators/GeneratorObjects
                          Generators/TruthUtils
                          MuonSpectrometer/MuonSimEvent
                          Reconstruction/RecEvent
                          Simulation/G4Sim/TrackRecord )

# Component(s) in the package:
atlas_add_component( MCTruthSimAlgs
                     src/*.cxx
                     src/components/MCTruthSimAlgs_entries.cxx
                     src/components/MCTruthSimAlgs_load.cxx
                     LINK_LIBRARIES CaloSimEvent AthenaBaseComps AthenaKernel CxxUtils PileUpToolsLib StoreGateLib SGtests xAODJet GaudiKernel GeneratorObjects TruthUtils MuonSimEvent RecEvent )

# Install files from the package:
atlas_install_python_modules( python/*.py )

