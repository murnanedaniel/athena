################################################################################
# Package: InDetV0Finder
################################################################################

# Declare the package name:
atlas_subdir( InDetV0Finder )

# Declare the package's dependencies:
atlas_depends_on_subdirs( PUBLIC
                          Control/AthenaBaseComps
                          Control/DataModel
                          GaudiKernel
                          InnerDetector/InDetRecTools/InDetConversionFinderTools
                          Reconstruction/RecoTools/ITrackToVertex
                          PRIVATE
                          Event/xAOD/xAODTracking
                          InnerDetector/InDetConditions/InDetBeamSpotService
                          Tracking/TrkExtrapolation/TrkExInterfaces
                          Tracking/TrkTools/TrkToolInterfaces
                          Tracking/TrkVertexFitter/TrkV0Fitter
                          Tracking/TrkVertexFitter/TrkVKalVrtFitter
                          Tracking/TrkVertexFitter/TrkVertexAnalysisUtils
                          Tracking/TrkVertexFitter/TrkVertexFitterInterfaces )

# External dependencies:
find_package( CLHEP )
find_package( HepPDT )
find_package( ROOT COMPONENTS Core Tree MathCore Hist RIO pthread )

# Component(s) in the package:
atlas_add_component( InDetV0Finder
                     src/*.cxx
                     src/components/*.cxx
                     INCLUDE_DIRS ${HEPPDT_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS} ${CLHEP_INCLUDE_DIRS}
                     LINK_LIBRARIES ${HEPPDT_LIBRARIES} ${ROOT_LIBRARIES} ${CLHEP_LIBRARIES} AthenaBaseComps DataModel GaudiKernel InDetConversionFinderToolsLib ITrackToVertex xAODTracking TrkExInterfaces TrkToolInterfaces TrkV0FitterLib TrkVKalVrtFitterLib TrkVertexAnalysisUtilsLib TrkVertexFitterInterfaces )

# Install files from the package:
atlas_install_headers( InDetV0Finder )

