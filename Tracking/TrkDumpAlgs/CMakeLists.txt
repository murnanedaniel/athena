################################################################################
# Package: TrkDumpAlgs
################################################################################

# Declare the package name:
atlas_subdir( TrkDumpAlgs )

# Declare the package's dependencies:
atlas_depends_on_subdirs( PUBLIC
                          Control/AthenaBaseComps
                          Control/AthenaKernel
                          DetectorDescription/Identifier
                          InnerDetector/InDetDetDescr/InDetIdentifier
                          InnerDetector/InDetDetDescr/PixelReadoutGeometry
                          InnerDetector/InDetDetDescr/SCT_ReadoutGeometry
                          InnerDetector/InDetDetDescr/InDetReadoutGeometry
                          InnerDetector/InDetRawEvent/InDetRawData
			  InnerDetector/InDetRecEvent/InDetRIO_OnTrack
			  InnerDetector/InDetSimEvent
			  Tracking/TrkEvent/TrkTrack
                          PRIVATE
                          InnerDetector/InDetRawEvent/InDetSimData
                          InnerDetector/InDetRecEvent/InDetPrepRawData
                          Tracking/TrkEvent/TrkSpacePoint
                          Tracking/TrkEvent/TrkTrack
			  Tracking/TrkEvent/TrkRIO_OnTrack
			  Tracking/TrkEvent/TrkTruthData
                          Generators/GeneratorObjects
                          Event/xAOD/xAODTruth
                          )

# External dependencies:
find_package( Eigen )
find_package( ROOT COMPONENTS Core Tree MathCore Hist RIO pthread )
find_package( HepMC )
find_package( HepPDT )

# tag use_new_tracking_edm was not recognized in automatic conversion in cmt2cmake

# Component(s) in the package:
atlas_add_component( TrkDumpAlgs
                     src/*.cxx
                     src/components/*.cxx
                     INCLUDE_DIRS ${ROOT_INCLUDE_DIRS} ${EIGEN_INCLUDE_DIRS} ${HEPMC_INCLUDE_DIRS} ${HEPPDT_INCLUDE_DIRS}
                     LINK_LIBRARIES ${HEPPDT_LIBRARIES} ${ROOT_LIBRARIES} ${EIGEN_LIBRARIES} ${HEPMC_LIBRARIES} AthenaBaseComps Identifier InDetReadoutGeometry InDetIdentifier InDetRIO_OnTrack PixelReadoutGeometry SCT_ReadoutGeometry InDetRawData InDetSimData InDetPrepRawData TrkSpacePoint TrkTrack TrkRIO_OnTrack InDetSimEvent)

# Install files from the package:
atlas_install_headers( TrkDumpAlgs )
atlas_install_joboptions( share/*.py )
