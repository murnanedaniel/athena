# Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
#
# This module is used to set up the environment for Sherpa
# 
#

# Set the environment variable(s):
find_package( Sherpa )

if( SHERPA_FOUND )
  set( SHERPAENVIRONMENT_ENVIRONMENT 
        FORCESET SHERPAVER ${SHERPA_VERSION} )
endif()

# Silently declare the module found:
set( SHERPAENVIRONMENT_FOUND TRUE )


