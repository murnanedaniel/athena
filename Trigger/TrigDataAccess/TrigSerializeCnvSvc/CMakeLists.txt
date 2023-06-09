################################################################################
# Package: TrigSerializeCnvSvc
################################################################################

# Declare the package name:
atlas_subdir( TrigSerializeCnvSvc )

# Declare the package's dependencies:
atlas_depends_on_subdirs( PUBLIC
                          Control/AthContainers
                          Control/AthenaBaseComps
                          Control/SGTools
                          Control/StoreGate
                          GaudiKernel
                          PRIVATE
                          Control/AthenaKernel
                          Control/DataModelRoot
                          Database/APR/CollectionBase
                          Database/APR/CollectionUtilities
                          Database/APR/POOLCore
                          Database/APR/PersistencySvc
                          Database/APR/StorageSvc
                          Database/PersistentDataModel
                          Trigger/TrigDataAccess/TrigSerializeResult
                          Trigger/TrigDataAccess/TrigSerializeTP )

# External dependencies:
find_package( ROOT COMPONENTS Core Tree MathCore Hist RIO pthread )

# Component(s) in the package:
atlas_add_library( TrigSerializeCnvSvcLib
                   src/*.cxx
                   PUBLIC_HEADERS TrigSerializeCnvSvc
                   PRIVATE_INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                   LINK_LIBRARIES AthContainers AthenaBaseComps SGTools GaudiKernel StoreGateLib SGtests TrigSerializeResultLib TrigSerializeTPLib
                   PRIVATE_LINK_LIBRARIES ${ROOT_LIBRARIES} AthenaKernel DataModelRoot CollectionBase CollectionUtilities POOLCore PersistencySvc StorageSvc PersistentDataModel )

atlas_add_component( TrigSerializeCnvSvc
                     src/components/*.cxx
                     INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                     LINK_LIBRARIES ${ROOT_LIBRARIES} AthContainers AthenaBaseComps SGTools StoreGateLib SGtests GaudiKernel AthenaKernel DataModelRoot CollectionBase CollectionUtilities POOLCore PersistencySvc StorageSvc PersistentDataModel TrigSerializeResultLib TrigSerializeTPLib TrigSerializeCnvSvcLib )

