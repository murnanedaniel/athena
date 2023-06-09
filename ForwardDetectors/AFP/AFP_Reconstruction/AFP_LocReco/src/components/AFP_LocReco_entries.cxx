#include "AFP_LocReco/AFP_TDLocReco.h"
#include "AFP_LocReco/AFP_SIDLocReco.h"
#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_ALGORITHM_FACTORY(AFP_TDLocReco)
DECLARE_ALGORITHM_FACTORY(AFP_SIDLocReco)

DECLARE_FACTORY_ENTRIES(AFP_LocReco) {
	DECLARE_ALGORITHM  (AFP_TDLocReco)
	DECLARE_ALGORITHM  (AFP_SIDLocReco)
}
