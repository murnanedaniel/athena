// $Id$
       
// Gaudi/Athena include(s):
#include "GaudiKernel/DeclareFactoryEntries.h"       
// Local include(s):

#include "../CMMCPHitsCnvTool.h"
#include "../CMMCPHitsCnvAlg.h"
#include "../CMMEtSumsCnvTool.h"
#include "../CMMEtSumsCnvAlg.h"
#include "../CMMJetHitsCnvTool.h"
#include "../CMMJetHitsCnvAlg.h" 
#include "../CMMRoICnvTool.h"
#include "../CMMRoICnvAlg.h" 

#include "../CPMHitsCnvTool.h"
#include "../CPMHitsCnvAlg.h"
#include "../CPMTowerCnvTool.h"
#include "../CPMTowerCnvAlg.h"
#include "../CPMRoICnvTool.h"
#include "../CPMRoICnvAlg.h"

#include "../JEMHitsCnvTool.h"
#include "../JEMHitsCnvAlg.h"
#include "../JEMEtSumsCnvTool.h"
#include "../JEMEtSumsCnvAlg.h"
#include "../JEMRoICnvTool.h"
#include "../JEMRoICnvAlg.h"

#include "../JetElementCnvTool.h"
#include "../JetElementCnvAlg.h"
#include "../RODHeaderCnvTool.h"
#include "../RODHeaderCnvAlg.h"
#include "../TriggerTowerCnvTool.h"
#include "../TriggerTowerCnvAlg.h"

DECLARE_NAMESPACE_TOOL_FACTORY( xAODMaker, CMMCPHitsCnvTool )
DECLARE_NAMESPACE_ALGORITHM_FACTORY( xAODMaker, CMMCPHitsCnvAlg )
DECLARE_NAMESPACE_TOOL_FACTORY( xAODMaker, CMMEtSumsCnvTool )
DECLARE_NAMESPACE_ALGORITHM_FACTORY( xAODMaker, CMMEtSumsCnvAlg )
DECLARE_NAMESPACE_TOOL_FACTORY( xAODMaker, CMMJetHitsCnvTool )
DECLARE_NAMESPACE_ALGORITHM_FACTORY( xAODMaker, CMMJetHitsCnvAlg )
DECLARE_NAMESPACE_TOOL_FACTORY( xAODMaker, CMMRoICnvTool )
DECLARE_NAMESPACE_ALGORITHM_FACTORY( xAODMaker, CMMRoICnvAlg )

DECLARE_NAMESPACE_TOOL_FACTORY( xAODMaker, CPMHitsCnvTool )
DECLARE_NAMESPACE_ALGORITHM_FACTORY( xAODMaker, CPMHitsCnvAlg )
DECLARE_NAMESPACE_TOOL_FACTORY( xAODMaker, CPMTowerCnvTool )
DECLARE_NAMESPACE_ALGORITHM_FACTORY( xAODMaker, CPMTowerCnvAlg )
DECLARE_NAMESPACE_TOOL_FACTORY( xAODMaker, CPMRoICnvTool )
DECLARE_NAMESPACE_ALGORITHM_FACTORY( xAODMaker, CPMRoICnvAlg )

DECLARE_NAMESPACE_TOOL_FACTORY( xAODMaker, JEMHitsCnvTool )
DECLARE_NAMESPACE_ALGORITHM_FACTORY( xAODMaker, JEMHitsCnvAlg )
DECLARE_NAMESPACE_TOOL_FACTORY( xAODMaker, JEMEtSumsCnvTool )
DECLARE_NAMESPACE_ALGORITHM_FACTORY( xAODMaker, JEMEtSumsCnvAlg )
DECLARE_NAMESPACE_TOOL_FACTORY( xAODMaker, JEMRoICnvTool )
DECLARE_NAMESPACE_ALGORITHM_FACTORY( xAODMaker, JEMRoICnvAlg )

DECLARE_NAMESPACE_TOOL_FACTORY( xAODMaker, JetElementCnvTool )
DECLARE_NAMESPACE_ALGORITHM_FACTORY( xAODMaker, JetElementCnvAlg )
DECLARE_NAMESPACE_TOOL_FACTORY( xAODMaker, RODHeaderCnvTool )
DECLARE_NAMESPACE_ALGORITHM_FACTORY( xAODMaker, RODHeaderCnvAlg )
DECLARE_NAMESPACE_TOOL_FACTORY( xAODMaker, TriggerTowerCnvTool )
DECLARE_NAMESPACE_ALGORITHM_FACTORY( xAODMaker, TriggerTowerCnvAlg )

DECLARE_FACTORY_ENTRIES( xAODTrigL1CaloCnv ) {      

  DECLARE_NAMESPACE_TOOL( xAODMaker, CMMCPHitsCnvTool )
  DECLARE_NAMESPACE_ALGORITHM( xAODMaker, CMMCPHitsCnvAlg )
  DECLARE_NAMESPACE_TOOL( xAODMaker, CMMEtSumsCnvTool )
  DECLARE_NAMESPACE_ALGORITHM( xAODMaker, CMMEtSumsCnvAlg )
  DECLARE_NAMESPACE_TOOL( xAODMaker, CMMJetHitsCnvTool )
  DECLARE_NAMESPACE_ALGORITHM( xAODMaker, CMMJetHitsCnvAlg ) 
  DECLARE_NAMESPACE_TOOL( xAODMaker, CMMRoICnvTool )
  DECLARE_NAMESPACE_ALGORITHM( xAODMaker, CMMRoICnvAlg )
  
  DECLARE_NAMESPACE_TOOL( xAODMaker, CPMHitsCnvTool )
  DECLARE_NAMESPACE_ALGORITHM( xAODMaker, CPMHitsCnvAlg )
  DECLARE_NAMESPACE_TOOL( xAODMaker, CPMTowerCnvTool )
  DECLARE_NAMESPACE_ALGORITHM( xAODMaker, CPMTowerCnvAlg )
  DECLARE_NAMESPACE_TOOL( xAODMaker, CPMRoICnvTool )
  DECLARE_NAMESPACE_ALGORITHM( xAODMaker, CPMRoICnvAlg )  
  
  DECLARE_NAMESPACE_TOOL( xAODMaker, JEMHitsCnvTool )
  DECLARE_NAMESPACE_ALGORITHM( xAODMaker, JEMHitsCnvAlg )
  DECLARE_NAMESPACE_TOOL( xAODMaker, JEMEtSumsCnvTool )
  DECLARE_NAMESPACE_ALGORITHM( xAODMaker, JEMEtSumsCnvAlg )
  DECLARE_NAMESPACE_TOOL( xAODMaker, JEMRoICnvTool )
  DECLARE_NAMESPACE_ALGORITHM( xAODMaker, JEMRoICnvAlg )   
  
  DECLARE_NAMESPACE_TOOL( xAODMaker, JetElementCnvTool )
  DECLARE_NAMESPACE_ALGORITHM( xAODMaker, JetElementCnvAlg )
  DECLARE_NAMESPACE_TOOL( xAODMaker, RODHeaderCnvTool )
  DECLARE_NAMESPACE_ALGORITHM( xAODMaker, RODHeaderCnvAlg )  
  DECLARE_NAMESPACE_TOOL( xAODMaker, TriggerTowerCnvTool )
  DECLARE_NAMESPACE_ALGORITHM( xAODMaker, TriggerTowerCnvAlg )   
}
