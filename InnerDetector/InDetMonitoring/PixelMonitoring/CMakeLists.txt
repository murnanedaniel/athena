################################################################################
# Package: PixelMonitoring
################################################################################

# Declare the package name:
atlas_subdir( PixelMonitoring )

# Declare the package's dependencies:
atlas_depends_on_subdirs( PUBLIC
                          Control/AthenaMonitoring
                          DetectorDescription/AtlasDetDescr
                          GaudiKernel
                          InnerDetector/InDetDetDescr/PixelGeoModel
                          InnerDetector/InDetRawEvent/InDetRawData
                          InnerDetector/InDetRecEvent/InDetPrepRawData
                          LumiBlock/LumiBlockComps
                          Tracking/TrkEvent/TrkTrack
                          PRIVATE
                          Database/AthenaPOOL/AthenaPoolUtilities
                          DetectorDescription/GeoPrimitives
                          Event/EventInfo
                          Event/EventPrimitives
                          InnerDetector/InDetConditions/InDetConditionsSummaryService
                          InnerDetector/InDetConditions/PixelConditionsServices
                          InnerDetector/InDetDetDescr/InDetIdentifier
                          InnerDetector/InDetDetDescr/InDetReadoutGeometry
                          InnerDetector/InDetDetDescr/PixelCabling
                          InnerDetector/InDetRecEvent/InDetRIO_OnTrack
                          Tools/LWHists
                          Tracking/TrkEvent/TrkParameters
                          Tracking/TrkEvent/TrkSpacePoint
                          Tracking/TrkEvent/TrkTrackSummary
                          Tracking/TrkTools/TrkToolInterfaces )

# External dependencies:
find_package( Eigen )
find_package( ROOT COMPONENTS Core Tree MathCore Hist RIO pthread )

# Component(s) in the package:
atlas_add_component( PixelMonitoring
                     src/*.cxx
                     src/componets/*.cxx
                     INCLUDE_DIRS ${ROOT_INCLUDE_DIRS} ${EIGEN_INCLUDE_DIRS}
                     LINK_LIBRARIES ${ROOT_LIBRARIES} ${EIGEN_LIBRARIES} AthenaMonitoringLib AtlasDetDescr GaudiKernel InDetRawData InDetPrepRawData LumiBlockCompsLib TrkTrack AthenaPoolUtilities GeoPrimitives EventInfo EventPrimitives InDetIdentifier InDetReadoutGeometry InDetRIO_OnTrack LWHists TrkParameters TrkSpacePoint TrkTrackSummary TrkToolInterfaces )

# Install files from the package:
atlas_install_headers( PixelMonitoring )
atlas_install_runtime( share/*.csv )

