################################################################################
# Package: Lhef_i
################################################################################

# Declare the package name:
atlas_subdir( Lhef_i )

# Declare the package's dependencies:
atlas_depends_on_subdirs( PRIVATE
                          Generators/GeneratorFortranCommon )

# External dependencies:
find_package( Pythia6 )

# tag allowUndefined was not recognized in automatic conversion in cmt2cmake

# Component(s) in the package:
atlas_add_library( Lhef_i
                   src/*.F
                   NO_PUBLIC_HEADERS
                   PRIVATE_INCLUDE_DIRS ${PYTHIA6_INCLUDE_DIRS}
                   LINK_LIBRARIES GeneratorFortranCommonLib
                   PRIVATE_LINK_LIBRARIES ${PYTHIA6_LIBRARIES} )

# Install files from the package:
atlas_install_joboptions( share/*.py )
atlas_install_runtime( share/lhefmadgraph.unweighted_events )

