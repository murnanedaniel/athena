#!/bin/bash
#
# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
#
# Script building all the externals necessary for AthDerivation.
#

# Set up the variables necessary for the script doing the heavy lifting.
ATLAS_PROJECT_DIR=$(cd $(dirname ${BASH_SOURCE[0]}) && pwd)
ATLAS_EXT_PROJECT_NAME="AthDerivationExternals"
ATLAS_BUILDTYPE="RelWithDebInfo"
ATLAS_EXTRA_CMAKE_ARGS=(-DLCG_VERSION_NUMBER=101
                        -DLCG_VERSION_POSTFIX="_ATLAS_26"
                        -DATLAS_GAUDI_TAG="v36r6.002"
                        -DATLAS_ACTS_TAG="v19.2.0"
                        -DATLAS_GEOMODEL_TAG="4.2.8.1")
ATLAS_EXTRA_MAKE_ARGS=()

# Let "the common script" do all the heavy lifting.
source "${ATLAS_PROJECT_DIR}/../../Build/AtlasBuildScripts/build_project_externals.sh"
