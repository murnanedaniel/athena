# $Id: CMakeLists.txt 753795 2016-06-08 23:37:56Z ssnyder $
################################################################################
# Package: AthenaServices
################################################################################

# Declare the package name:
atlas_subdir( AthenaServices )

# Declare the package's dependencies:
atlas_depends_on_subdirs(
   PRIVATE
   AtlasTest/TestTools
   Control/AthContainersInterfaces
   Control/AthenaBaseComps
   Control/AthenaKernel
   Control/CxxUtils
   Control/DataModel
   Control/DataModelRoot
   Control/Navigation
   Control/PerformanceMonitoring/PerfMonEvent
   Control/PerformanceMonitoring/PerfMonKernel
   Control/RngComps
   Control/SGTools
   Control/StoreGate
   Database/PersistentDataModel
   Event/EventInfo
   Event/xAOD/xAODEventInfo
   Event/EventInfoUtils
   GaudiKernel )

# External dependencies:
find_package( Boost COMPONENTS thread )
find_package( CLHEP )
find_package( PythonLibs )
find_package( yampl )

# Component(s) in the package:
atlas_add_component( AthenaServices src/*.cxx src/components/*.cxx
   INCLUDE_DIRS ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS}
   ${YAMPL_INCLUDE_DIRS} ${CLHEP_INCLUDE_DIRS}
   LINK_LIBRARIES ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ${YAMPL_LIBRARIES}
   ${CLHEP_LIBRARIES} z TestTools AthenaBaseComps AthenaKernel CxxUtils
   DataModel DataModelRoot Navigation PerfMonEvent PerfMonKernel SGTools
   StoreGateLib SGtests PersistentDataModel EventInfo xAODEventInfo EventInfoUtils GaudiKernel )

# Test library checking the ability to build T/P converters:
atlas_add_tpcnv_library( AthenaServicesTest src/test/*.cxx
   NO_PUBLIC_HEADERS
   PRIVATE_LINK_LIBRARIES AthenaKernel )

# The test(s) of the package:
atlas_add_test( AthenaOutputStream_test
   SOURCES test/AthenaOutputStream_test.cxx src/AthenaOutputStream.cxx
   src/OutputStreamSequencerSvc.cxx src/MetaDataSvc.cxx
   INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
   LINK_LIBRARIES TestTools AthenaKernel SGTools StoreGateLib GaudiKernel
   AthenaBaseComps PersistentDataModel
   EXTRA_PATTERNS "^AthenaRootStrea... +(INFO|DEBUG)|^AthenaSealSvc +DEBUG|^SGAudSvc +(INFO|DEBUG)|of type DataHistory|DEBUG Recorded object|object modifiable when retrieved|^ToolSvc +DEBUG Service base class initialized|^IncidentSvc *DEBUG Adding .* listener|^Py:Athena *INFO using release|DecisionSvc +DEBUG|: duplicate entry|^IoComponentMgr +(INFO|DEBUG)|DBReplicaSvc|^HistogramPersis.*DEBUG|^JobOptionsSvc +INFO|^//GP:=|^ItemListSvc +(INFO|DEBUG)|Info File PoolFileCatalog.xml does not exist|DataModelCompatSvc::initialize|^ProxyProviderSvc +DEBUG|^DataModelCompatSvc +DEBUG|^AthenaOutputStreamVERBOSE|^AthenaOutputStream +DEBUG|Service base class initialized successfully|^TimelineSvc +DEBUG|DEBUG Property update for OutputLevel"
   ENVIRONMENT "JOBOPTSEARCHPATH=${CMAKE_CURRENT_SOURCE_DIR}/share" )

atlas_add_test( FPEControlSvc_test
   SOURCES test/FPEControlSvc_test.cxx
   LINK_LIBRARIES TestTools AthenaKernel GaudiKernel
   EXTRA_PATTERNS "^JobOptionsSvc +INFO|^//GP:="
   ENVIRONMENT "JOBOPTSEARCHPATH=${CMAKE_CURRENT_SOURCE_DIR}/share" )

atlas_add_test( AthenaEventLoopMgr_test
   SOURCES test/AthenaEventLoopMgr_test.cxx
   LINK_LIBRARIES TestTools AthenaKernel GaudiKernel EventInfo AthenaBaseComps
   EXTRA_PATTERNS "^JobOptionsSvc +INFO|^//GP:="
   ENVIRONMENT "JOBOPTSEARCHPATH=${CMAKE_CURRENT_SOURCE_DIR}/share" )

atlas_add_test( TestStopRun
   SCRIPT test/TestStopRun.sh
   EXTRA_PATTERNS "SGAudSvc +INFO Finalizing|SGAudSvc +INFO Initializing|^Py:Athena +INFO executing ROOT6Setup|No duplicates have been found|duplicate entry.*ignored|^Py:ConfigurableDb WARNING|Read module info for|^ToolSvc.* INFO( finalize:)? Removing all tools|^CoreDumpSvc *INFO Handling signals|types added|including file|^(Mon|Tue|Wed|Thu|Fri|Sat|Sun) (Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)|local .* libGaudiKernelDict.so|^Number of My|^Py:Athena +INFO using release|^StoreGateSvc +INFO Start|^Py:AthenaDsoDb +INFO could not install alias|Bindings.py:660newobj = object.__new__|Updating ROOT::Reflex::PluginService::SetDebug|^ApplicationMgr +INFO|^StoreGateSvc +INFO"
   PROPERTIES TIMEOUT 300 )

atlas_add_test( AthTPCnvSvc
   SCRIPT test/AthTPCnvSvc.sh
   PROPERTIES TIMEOUT 300
   EXTRA_PATTERNS "SGAudSvc +INFO Finalizing|SGAudSvc +INFO Initializing|^Py:Athena +INFO executing ROOT6Setup|No duplicates have been found|duplicate entry.*ignored|^Py:ConfigurableDb WARNING|Read module info for|^ToolSvc.* INFO( finalize:)? Removing all tools|^CoreDumpSvc *INFO Handling signals|types added|including file|^(Mon|Tue|Wed|Thu|Fri|Sat|Sun) (Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)|local .* libGaudiKernelDict.so|^Number of My|^Py:Athena +INFO using release|^StoreGateSvc +INFO Start|^Py:AthenaDsoDb +INFO could not install alias|Bindings.py:660newobj = object.__new__|Updating ROOT::Reflex::PluginService::SetDebug|^ApplicationMgr +INFO|^StoreGateSvc +INFO" )

# Install files from the package:
atlas_install_python_modules( python/*.py )
atlas_install_joboptions( share/AthTPCnvSvc_test.py
   share/SimplePOSIXTimeKeeperOptions.py
   share/SimplePOSIXTimeKeeperSvc_test.py
   share/MixingEventSelector_test.py
   share/MultiplePassBootstrap.py
   share/MultiplePass_test.py
   share/MultiplePassWithAlgFilter_test.py
   share/OverrideEventNumber.py
   share/TestStopRun.py
   share/TestSeedRunEvent.py
   share/AthDictLoaderSvc_test.py
   share/ReadAthenaPoolSeek_jobOptions.py
   share/AthenaJobOptionsSvc_jobOptions.py
   share/test_tpcnvdb.py )

# Need to make sure that python file installation happens before genconf
# in order for our custom __init__.py to be installed correctly.
add_dependencies( AthenaServicesConfigurables AthenaServicesPythonInstall )
