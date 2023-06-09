################################################################################
# Package: LArG4Validation
################################################################################

# Declare the package name:
atlas_subdir( LArG4Validation )

# Declare the package's dependencies:
atlas_depends_on_subdirs( PUBLIC
                          GaudiKernel
                          PRIVATE
                          Calorimeter/CaloDetDescr
                          Calorimeter/CaloIdentifier
                          Control/AthenaBaseComps
                          Control/StoreGate
                          DetectorDescription/GeoModel/GeoAdaptors
                          DetectorDescription/GeoModel/GeoModelKernel
                          Event/EventInfo
                          Generators/GeneratorObjects
                          LArCalorimeter/LArSimEvent
                          MagneticField/MagFieldInterfaces
                          Reconstruction/egamma/egammaEvent
                          TestPolicy )

# External dependencies:
find_package( CLHEP )
find_package( HepPDT )
find_package( ROOT COMPONENTS Core Tree MathCore Hist RIO pthread )

# Component(s) in the package:
atlas_add_component( LArG4Validation
                     src/*.cxx
                     src/components/*.cxx
                     INCLUDE_DIRS ${HEPPDT_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS} ${CLHEP_INCLUDE_DIRS}
                     LINK_LIBRARIES ${HEPPDT_LIBRARIES} ${ROOT_LIBRARIES} ${CLHEP_LIBRARIES} GaudiKernel CaloDetDescrLib CaloIdentifier AthenaBaseComps StoreGateLib SGtests GeoAdaptors GeoModelKernel EventInfo GeneratorObjects LArSimEvent MagFieldInterfaces egammaEvent )

# Install files from the package:
atlas_install_headers( LArG4Validation )
atlas_install_python_modules( python/*.py )
atlas_install_joboptions( share/*.py )
atlas_install_runtime( test/*.xml )
atlas_install_scripts( share/AODElectronContainerReader.py share/LArG4ValidationGenerate.py share/LArG4AODNtuplePlotter.py share/LArG4ValidationPlotter.py )

