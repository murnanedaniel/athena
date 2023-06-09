################################################################################
# Package: TrigExPartialEB
################################################################################

# Declare the package name:
atlas_subdir( TrigExPartialEB )

# Declare the package's dependencies:
atlas_depends_on_subdirs( PUBLIC
                          Control/AthenaBaseComps
                          Event/EventInfo
                          GaudiKernel
                          Trigger/TrigDataAccess/TrigROBDataProviderSvc
                          PRIVATE
                          Control/AthenaKernel
                          Control/StoreGate
                          Event/ByteStreamCnvSvcBase
                          Trigger/TrigDataAccess/TrigDataAccessMonitoring
                          Trigger/TrigEvent/TrigSteeringEvent
                          Trigger/TrigT1/TrigT1Result )

# External dependencies:
find_package( Boost COMPONENTS filesystem thread system )
find_package( tdaq-common )

# Component(s) in the package:
atlas_add_component( TrigExPartialEB
                     src/*.cxx
                     src/components/*.cxx
                     INCLUDE_DIRS ${Boost_INCLUDE_DIRS} ${TDAQ-COMMON_INCLUDE_DIRS}
                     LINK_LIBRARIES ${Boost_LIBRARIES} ${TDAQ-COMMON_LIBRARIES} AthenaBaseComps EventInfo GaudiKernel TrigROBDataProviderSvcLib AthenaKernel StoreGateLib SGtests ByteStreamCnvSvcBaseLib TrigDataAccessMonitoringLib TrigSteeringEvent TrigT1Result )

# Install files from the package:
atlas_install_headers( TrigExPartialEB )
atlas_install_joboptions( share/*.py )

