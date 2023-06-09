################################################################################
# Package: InDetAmbiTrackSelectionTool
################################################################################

# Declare the package name:
atlas_subdir( InDetAmbiTrackSelectionTool )

# Declare the package's dependencies:
atlas_depends_on_subdirs( PUBLIC
                          Control/AthenaBaseComps
                          GaudiKernel
                          Tracking/TrkEvent/TrkRIO_OnTrack
                          Tracking/TrkEvent/TrkTrack
                          Tracking/TrkTools/TrkToolInterfaces
                          Tracking/TrkValidation/TrkValInterfaces
                          PRIVATE
                          Control/DataModel
                          Control/StoreGate
                          InnerDetector/InDetDetDescr/InDetIdentifier
                          InnerDetector/InDetDetDescr/PixelGeoModel
                          InnerDetector/InDetRecEvent/InDetPrepRawData
                          InnerDetector/InDetRecTools/InDetRecToolInterfaces
                          Tracking/TrkDetDescr/TrkDetElementBase
                          Tracking/TrkDetDescr/TrkSurfaces
                          Tracking/TrkEvent/TrkCaloClusterROI
                          Tracking/TrkEvent/TrkMeasurementBase
                          Tracking/TrkEvent/TrkPrepRawData
                          Tracking/TrkEvent/TrkPseudoMeasurementOnTrack
                          Tracking/TrkEvent/TrkTrackSummary 
                          Tracking/TrkEvent/TrkParameters )

# External dependencies:
find_package( ROOT COMPONENTS Core Tree MathCore Hist RIO pthread )

# Component(s) in the package:
atlas_add_component( InDetAmbiTrackSelectionTool
                     src/*.cxx
                     src/components/*.cxx
                     INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                     LINK_LIBRARIES ${ROOT_LIBRARIES} AthenaBaseComps GaudiKernel TrkRIO_OnTrack TrkTrack TrkToolInterfaces TrkValInterfaces DataModel StoreGateLib SGtests InDetIdentifier InDetPrepRawData InDetRecToolInterfaces TrkDetElementBase TrkSurfaces TrkCaloClusterROI TrkMeasurementBase TrkPrepRawData TrkPseudoMeasurementOnTrack TrkTrackSummary TrkParameters )

# Install files from the package:
atlas_install_headers( InDetAmbiTrackSelectionTool )

