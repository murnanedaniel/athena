################################################################################
# Package: FTK_RecTools
################################################################################

# Declare the package name:
atlas_subdir( FTK_RecTools )

# Declare the package's dependencies:
atlas_depends_on_subdirs( PUBLIC
                          AtlasPolicy
                          Control/AthenaBaseComps
                          Event/xAOD/xAODTracking
                          GaudiKernel
                          Tools/PyJobTransforms
                          DetectorDescription/GeoPrimitives
                          Tracking/TrkEvent/TrkTrack
                          Tracking/TrkEvent/VxVertex
                          Tracking/TrkVertexFitter/TrkVxEdmCnv
                          Trigger/TrigFTK/FTK_DataProviderInterfaces
                          Trigger/TrigFTK/FTK_RecToolInterfaces
                          Trigger/TrigFTK/TrigFTK_RawData
                          InnerDetector/InDetDetDescr/SCT_ModuleDistortions
                          InnerDetector/InDetRecEvent/InDetPrepRawData
                          InnerDetector/InDetRecEvent/InDetRIO_OnTrack
                          Tracking/TrkEvent/TrkParameters
                          Tracking/TrkTools/TrkToolInterfaces
                          Tracking/TrkTools/TrkAmbiguityProcessor
                          PRIVATE
                          Event/EventPrimitives
                          InnerDetector/InDetConditions/PixelConditionsServices
                          InnerDetector/InDetConditions/PixelConditionsTools
                          InnerDetector/InDetDetDescr/InDetIdentifier
                          InnerDetector/InDetDetDescr/InDetReadoutGeometry
                          InnerDetector/InDetDetDescr/PixelGeoModel
                          InnerDetector/InDetRecTools/SiClusterizationTool
                          Tracking/TrkDetDescr/TrkSurfaces )


# External dependencies:
find_package( ROOT COMPONENTS Core Tree MathCore Hist RIO pthread )

# Component(s) in the package:
atlas_add_library( FTK_RecToolsLib
                   src/*.cxx
                   PUBLIC_HEADERS FTK_RecTools
                   INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                   LINK_LIBRARIES ${ROOT_LIBRARIES} AthenaBaseComps GeoPrimitives xAODTracking TrkTrack VxVertex FTK_DataProviderInterfaces TrigFTK_RawData TrkVxEdmCnvLib  InDetPrepRawData InDetRIO_OnTrack TrkParameters TrkToolInterfaces TrkAmbiguityProcessorLib StoreGateLib SGtests EventPrimitives InDetIdentifier InDetReadoutGeometry SiClusterizationToolLib TrkSurfaces)

atlas_add_component( FTK_RecTools
                     src/*.cxx src/components/*.cxx
                     INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                     LINK_LIBRARIES ${ROOT_LIBRARIES} AthenaBaseComps GeoPrimitives xAODTracking TrkTrack VxVertex TrkVxEdmCnvLib FTK_DataProviderInterfaces TrigFTK_RawData TrkParameters FTK_RecToolsLib InDetPrepRawData InDetRIO_OnTrack TrkParameters TrkToolInterfaces TrkAmbiguityProcessorLib StoreGateLib SGtests EventPrimitives InDetIdentifier InDetReadoutGeometry SiClusterizationToolLib TrkSurfaces)

