#include "TrigLongLivedParticlesHypo/TrigL2HVJetHypoAllCuts.h"
#include "TrigLongLivedParticlesHypo/TrigL2HVJetHypo.h"
#include "TrigLongLivedParticlesHypo/TrigL2HVJetHypoTrk.h"
#include "TrigLongLivedParticlesHypo/MuonClusterHypo.h"
#include "TrigLongLivedParticlesHypo/TrigLoFRemovalHypo.h"
#include "TrigLongLivedParticlesHypo/TrigCaloRatioHypo.h"
#include "GaudiKernel/DeclareFactoryEntries.h"


DECLARE_ALGORITHM_FACTORY( TrigL2HVJetHypoAllCuts )
DECLARE_ALGORITHM_FACTORY( TrigL2HVJetHypo )
DECLARE_ALGORITHM_FACTORY( TrigL2HVJetHypoTrk )
DECLARE_ALGORITHM_FACTORY( MuonClusterHypo )
DECLARE_ALGORITHM_FACTORY( TrigLoFRemovalHypo )
DECLARE_ALGORITHM_FACTORY( TrigCaloRatioHypo )

DECLARE_FACTORY_ENTRIES( TrigLongLivedParticlesHypo ) {
  DECLARE_ALGORITHM( TrigL2HVJetHypoAllCuts )
  DECLARE_ALGORITHM( TrigL2HVJetHypo )
  DECLARE_ALGORITHM( TrigL2HVJetHypoTrk )
  DECLARE_ALGORITHM( MuonClusterHypo )
  DECLARE_ALGORITHM( TrigLoFRemovalHypo )
  DECLARE_ALGORITHM( TrigCaloRatioHypo )
}

