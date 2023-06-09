################################################################################
# Package: eflowRec
################################################################################

# Declare the package name:
atlas_subdir( eflowRec )

# Declare the package's dependencies:
atlas_depends_on_subdirs( PUBLIC
                          Calorimeter/CaloEvent
                          Calorimeter/CaloIdentifier
                          Calorimeter/CaloRec
                          Control/AthContainers
                          Control/AthLinks
                          Control/AthenaBaseComps
                          Control/CLIDSvc
                          Control/CxxUtils
                          Control/SGTools
                          DetectorDescription/GeoPrimitives
                          Event/xAOD/xAODCaloEvent
                          Event/xAOD/xAODEgamma
                          Event/xAOD/xAODMuon
                          Event/xAOD/xAODPFlow
                          Event/xAOD/xAODTracking
                          GaudiKernel
                          InnerDetector/InDetRecTools/InDetTrackSelectionTool
                          InnerDetector/InDetRecTools/TrackVertexAssociationTool
                          Reconstruction/Particle
                          Reconstruction/RecoTools/RecoToolInterfaces
                          Tracking/TrkEvent/TrkParameters
                          PRIVATE
                          Calorimeter/CaloDetDescr
                          Calorimeter/CaloUtils
                          Control/StoreGate
                          Event/FourMomUtils
                          Tools/PathResolver
                          Tracking/TrkEvent/TrkCaloExtension
                          Tracking/TrkEvent/TrkParametersIdentificationHelpers )

# External dependencies:
find_package( Eigen )
find_package( ROOT COMPONENTS Core Tree MathCore Hist RIO pthread )

# Component(s) in the package:
atlas_add_component( eflowRec
                     src/*.cxx
                     src/components/*.cxx
                     INCLUDE_DIRS ${ROOT_INCLUDE_DIRS} ${EIGEN_INCLUDE_DIRS}
                     LINK_LIBRARIES ${ROOT_LIBRARIES} ${EIGEN_LIBRARIES} CaloEvent CaloIdentifier CaloRecLib AthContainers AthLinks AthenaBaseComps CxxUtils SGTools GeoPrimitives xAODCaloEvent xAODEgamma xAODMuon xAODPFlow xAODTracking GaudiKernel TrackVertexAssociationToolLib Particle RecoToolInterfaces TrkParameters CaloDetDescrLib CaloUtilsLib StoreGateLib SGtests FourMomUtils PathResolver TrkCaloExtension TrkParametersIdentificationHelpers )

# Install files from the package:
atlas_install_headers( eflowRec )
atlas_install_python_modules( python/*.py )
atlas_install_joboptions( share/*.py )
atlas_install_runtime( share/rmsvsdr.dat )

