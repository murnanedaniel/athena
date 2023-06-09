################################################################################
# Package: Geo2G4
################################################################################

# Declare the package name:
atlas_subdir( Geo2G4 )

# Declare the package's dependencies:
atlas_depends_on_subdirs( PUBLIC
                          Control/AthenaKernel
                          DetectorDescription/GeoModel/GeoModelKernel
                          DetectorDescription/GeoModel/GeoModelUtilities
                          GaudiKernel
                          PRIVATE
                          Control/AthenaBaseComps
                          Control/SGTools
                          Control/StoreGate
                          DetectorDescription/GeoModel/GeoModelInterfaces
                          DetectorDescription/GeoModel/GeoSpecialShapes
                          Simulation/G4Sim/SimHelpers )

# External dependencies:
find_package( Boost COMPONENTS filesystem thread system )
find_package( CLHEP )
find_package( Geant4 )
find_package( ROOT COMPONENTS MathCore RIO Core Tree Hist pthread )
find_package( XercesC )

# Component(s) in the package:
atlas_add_library( Geo2G4Lib
                   src/*.cxx
                   PUBLIC_HEADERS Geo2G4
                   INCLUDE_DIRS ${XERCESC_INCLUDE_DIRS} ${GEANT4_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS}
                   PRIVATE_INCLUDE_DIRS ${Boost_INCLUDE_DIRS} ${CLHEP_INCLUDE_DIRS}
                   PRIVATE_DEFINITIONS ${CLHEP_DEFINITIONS}
                   LINK_LIBRARIES ${XERCESC_LIBRARIES} ${GEANT4_LIBRARIES} ${ROOT_LIBRARIES} AthenaKernel GeoModelKernel GeoModelUtilities GaudiKernel StoreGateLib SGtests
                   PRIVATE_LINK_LIBRARIES ${Boost_LIBRARIES} ${CLHEP_LIBRARIES} AthenaBaseComps SGTools GeoSpecialShapes SimHelpers )

atlas_add_component( Geo2G4
                     src/components/*.cxx
                     INCLUDE_DIRS ${Boost_INCLUDE_DIRS} ${XERCESC_INCLUDE_DIRS} ${GEANT4_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS} ${CLHEP_INCLUDE_DIRS}
                     LINK_LIBRARIES ${Boost_LIBRARIES} ${XERCESC_LIBRARIES} ${GEANT4_LIBRARIES} ${ROOT_LIBRARIES} ${CLHEP_LIBRARIES} AthenaKernel GeoModelKernel GeoModelUtilities GaudiKernel AthenaBaseComps SGTools StoreGateLib SGtests GeoSpecialShapes SimHelpers Geo2G4Lib )

atlas_add_dictionary( LArWheelSolidCheckerDict
                      src/LArWheelSolidDDProxy.h
                      src/lcg_dict/selection.xml
                      INCLUDE_DIRS ${Boost_INCLUDE_DIRS} ${XERCESC_INCLUDE_DIRS} ${GEANT4_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS} ${CLHEP_INCLUDE_DIRS}
                      LINK_LIBRARIES ${Boost_LIBRARIES} ${XERCESC_LIBRARIES} ${GEANT4_LIBRARIES} ${ROOT_LIBRARIES} ${CLHEP_LIBRARIES} AthenaKernel GeoModelKernel GeoModelUtilities GaudiKernel AthenaBaseComps SGTools StoreGateLib SGtests GeoSpecialShapes SimHelpers Geo2G4Lib )

