################################################################################
# Package: EventIndexProducer
################################################################################

# Declare the package name:
atlas_subdir( EventIndexProducer )

# Declare the package's dependencies:
atlas_depends_on_subdirs( PUBLIC
                          Tools/PyJobTransforms
                          Tools/PyJobTransformsCore )

# Install files from the package:
atlas_install_headers( EventIndexProducer )
atlas_install_python_modules( python/*.py )
atlas_install_joboptions( share/*.py )
atlas_install_runtime( scripts/*.py )

