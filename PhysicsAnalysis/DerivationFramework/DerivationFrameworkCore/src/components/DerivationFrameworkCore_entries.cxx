#include "DerivationFrameworkCore/DerivationKernel.h"
#include "DerivationFrameworkCore/CommonAugmentation.h"
#include "GaudiKernel/DeclareFactoryEntries.h"
#include "DerivationFrameworkCore/TriggerMatchingAugmentation.h"

DECLARE_NAMESPACE_ALGORITHM_FACTORY( DerivationFramework , DerivationKernel )
DECLARE_NAMESPACE_ALGORITHM_FACTORY( DerivationFramework , CommonAugmentation )

#include "../GoodRunsListFilterAlgorithm.h"
DECLARE_NAMESPACE_ALGORITHM_FACTORY( DerivationFramework, GoodRunsListFilterAlgorithm )

DECLARE_FACTORY_ENTRIES(DerivationFrameworkCore) {
 DECLARE_NAMESPACE_ALGORITHM( DerivationFramework, GoodRunsListFilterAlgorithm )
 DECLARE_NAMESPACE_ALGORITHM( DerivationFramework , DerivationKernel )
 DECLARE_NAMESPACE_ALGORITHM( DerivationFramework , CommonAugmentation )
}

DECLARE_TOOL_FACTORY( TriggerMatchingAugmentation  )
DECLARE_FACTORY_ENTRIES( DerivationFrameworkCore ) {
    DECLARE_TOOL( TriggerMatchingAugmentation )
}
