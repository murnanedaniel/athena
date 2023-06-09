################################################################################
# Package: LArG4EC
################################################################################

# Declare the package name:
atlas_subdir( LArG4EC )

# Declare the package's dependencies:
atlas_depends_on_subdirs( PUBLIC
                          Calorimeter/CaloG4Sim
                          DetectorDescription/GeoModel/GeoSpecialShapes
                          LArCalorimeter/LArG4/LArG4Code
                          LArCalorimeter/LArGeoModel/LArGeoCode
                          PRIVATE
                          Control/StoreGate
                          Database/AthenaPOOL/RDBAccessSvc
                          DetectorDescription/GeoModel/GeoModelInterfaces
                          DetectorDescription/GeoModel/GeoModelUtilities
                          GaudiKernel
                          LArCalorimeter/LArG4/LArG4RunControl
                          LArCalorimeter/LArGeoModel/LArHV
                          Tools/PathResolver )

# External dependencies:
find_package( Boost COMPONENTS filesystem thread system )
find_package( CLHEP )
find_package( CORAL COMPONENTS CoralBase CoralKernel RelationalAccess )
find_package( Geant4 )
find_package( XercesC )

# Component(s) in the package:
atlas_add_library( LArG4EC
                   src/*.cc
                   PUBLIC_HEADERS LArG4EC
                   INCLUDE_DIRS ${XERCESC_INCLUDE_DIRS} ${GEANT4_INCLUDE_DIRS}
                   PRIVATE_INCLUDE_DIRS ${Boost_INCLUDE_DIRS} ${CORAL_INCLUDE_DIRS} ${CLHEP_INCLUDE_DIRS}
                   PRIVATE_DEFINITIONS ${CLHEP_DEFINITIONS}
                   LINK_LIBRARIES ${XERCESC_LIBRARIES} ${GEANT4_LIBRARIES} GeoSpecialShapes LArG4Code LArGeoCode CaloG4SimLib StoreGateLib SGtests
                   PRIVATE_LINK_LIBRARIES ${Boost_LIBRARIES} ${CORAL_LIBRARIES} ${CLHEP_LIBRARIES} GeoModelUtilities GaudiKernel LArG4RunControl LArHV PathResolver )

