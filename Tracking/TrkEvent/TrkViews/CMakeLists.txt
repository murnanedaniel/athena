################################################################################
# Package: TrkViews
################################################################################

# Declare the package name:
atlas_subdir( TrkViews )

# Declare the package's dependencies:
atlas_depends_on_subdirs( PUBLIC
                          Event/xAOD/xAODCore
                          Tracking/TrkEvent/TrkTrack
                          PRIVATE
                          Control/AthContainers )

# External dependencies:
find_package( ROOT COMPONENTS Core Tree MathCore Hist RIO pthread )

# Component(s) in the package:
atlas_add_dictionary( TrkViewsDict
                      TrkViews/TrkViewsDict.h
                      TrkViews/selection.xml
                      INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                      LINK_LIBRARIES ${ROOT_LIBRARIES} xAODCore TrkTrack AthContainers )

# Install files from the package:
atlas_install_headers( TrkViews )

