################################################################################
# Package: DQDefects
################################################################################

# Declare the package name:
atlas_subdir( DQDefects )

# External dependencies:
find_package( Oracle )
find_package( pytools )

# tag NEEDS_PYCOOL was not recognized in automatic conversion in cmt2cmake

# tag NEEDS_COOL_FACTORY was not recognized in automatic conversion in cmt2cmake

# Install files from the package:
atlas_install_python_modules( python/*.py python/tests python/data )
atlas_install_scripts( share/*.py share/*.sh )

