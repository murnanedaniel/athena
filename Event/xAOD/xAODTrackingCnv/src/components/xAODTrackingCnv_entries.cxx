// $Id: xAODCreatorAlgs_entries.cxx 298282 2013-12-03 13:05:16Z emoyse $

// Gaudi/Athena include(s):
#include "GaudiKernel/DeclareFactoryEntries.h"

// Local include(s):
#include "../TrackParticleCnvAlg.h"
#include "../VertexCnvAlg.h"

DECLARE_NAMESPACE_ALGORITHM_FACTORY( xAODMaker, TrackParticleCnvAlg )
DECLARE_NAMESPACE_ALGORITHM_FACTORY( xAODMaker, VertexCnvAlg )

DECLARE_FACTORY_ENTRIES( xAODCreatorAlgs ) {
  DECLARE_NAMESPACE_ALGORITHM( xAODMaker, TrackParticleCnvAlg )
    DECLARE_NAMESPACE_ALGORITHM( xAODMaker, VertexCnvAlg )
}
