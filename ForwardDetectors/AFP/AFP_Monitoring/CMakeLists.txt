# $Id: CMakeLists.txt 766963 2016-08-08 20:07:09Z ggach $
################################################################################
# Package: AFP_Monitoring
################################################################################

# Declare the package name:
atlas_subdir( AFP_Monitoring )

# Declare the package's dependencies:
atlas_depends_on_subdirs( PUBLIC
   Control/AthenaMonitoring
   Event/xAOD/xAODForward
   ForwardDetectors/AFP/AFP_RawEv
   PRIVATE
   External/AtlasROOT
   Tools/LWHists
   Event/xAOD/xAODEventInfo
   )

# Component(s) in the package:
atlas_add_library( AFP_Monitoring
   AFP_Monitoring/*.h src/*.cxx
   PUBLIC_HEADERS AFP_Monitoring
   LINK_LIBRARIES ${ROOT_LIBRARIES} AthenaMonitoringLib LWHists xAODForward xAODEventInfo AFP_RawEv)
