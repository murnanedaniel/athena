#include "../CPMon.h"
#include "../CPSimMon.h"
#include "../JEPJEMMon.h"
#include "../JEPCMXMon.h"
#include "../JEPSimMon.h"
#include "../PPrMon.h"
#include "../PPrStabilityMon.h"
#include "../PPrSpareMon.h"
#include "../PPMSimBSMon.h"
#include "../RODMon.h"
#include "../OverviewMon.h"

//Run 1
#include "../CMMMon.h"
#include "../CPMSimBSMon.h"
#include "../JEMMon.h"
//#include "../JEPSimBSMon.h"
#include "../TrigT1CaloBSMon.h"
#include "../TrigT1CaloCpmMonTool.h"
#include "../TrigT1CaloGlobalMonTool.h"
#include "../EmEfficienciesMonTool.h"
#include "../JetEfficienciesMonTool.h"
#include "../RODMonV1.h"

#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_NAMESPACE_TOOL_FACTORY(LVL1, OverviewMon)
DECLARE_NAMESPACE_TOOL_FACTORY(LVL1, CPMon)
DECLARE_NAMESPACE_TOOL_FACTORY(LVL1, CPSimMon)
DECLARE_NAMESPACE_TOOL_FACTORY(LVL1, JEPJEMMon)
DECLARE_NAMESPACE_TOOL_FACTORY(LVL1, JEPCMXMon)
DECLARE_NAMESPACE_TOOL_FACTORY(LVL1, JEPSimMon)

DECLARE_NAMESPACE_TOOL_FACTORY(LVL1, PPrMon)
DECLARE_NAMESPACE_TOOL_FACTORY(LVL1, PPrStabilityMon)
DECLARE_NAMESPACE_TOOL_FACTORY(LVL1, PPrSpareMon)
DECLARE_NAMESPACE_TOOL_FACTORY(LVL1, PPMSimBSMon)
DECLARE_NAMESPACE_TOOL_FACTORY(LVL1, RODMon)

//Run 1
DECLARE_NAMESPACE_TOOL_FACTORY(LVL1, CMMMon)
DECLARE_NAMESPACE_TOOL_FACTORY(LVL1, CPMSimBSMon)
DECLARE_NAMESPACE_TOOL_FACTORY(LVL1, JEMMon)
//DECLARE_NAMESPACE_TOOL_FACTORY(LVL1, JEPSimBSMon)
DECLARE_NAMESPACE_TOOL_FACTORY(LVL1, TrigT1CaloBSMon)
DECLARE_NAMESPACE_TOOL_FACTORY(LVL1, TrigT1CaloCpmMonTool)
DECLARE_NAMESPACE_TOOL_FACTORY(LVL1, TrigT1CaloGlobalMonTool)
DECLARE_NAMESPACE_TOOL_FACTORY(LVL1, EmEfficienciesMonTool)
DECLARE_NAMESPACE_TOOL_FACTORY(LVL1, JetEfficienciesMonTool)
DECLARE_NAMESPACE_TOOL_FACTORY(LVL1, RODMonV1)


DECLARE_FACTORY_ENTRIES(TrigT1CaloMonitoring) {
  DECLARE_NAMESPACE_ALGTOOL(LVL1, OverviewMon)
  DECLARE_NAMESPACE_ALGTOOL(LVL1, CPMon)
  DECLARE_NAMESPACE_ALGTOOL(LVL1, CPSimMon)
  DECLARE_NAMESPACE_ALGTOOL(LVL1, JEPJEMMon)
  DECLARE_NAMESPACE_ALGTOOL(LVL1, JEPCMXMon)
  DECLARE_NAMESPACE_ALGTOOL(LVL1, JEPSimMon)
  DECLARE_NAMESPACE_ALGTOOL(LVL1, PPrMon)
  DECLARE_NAMESPACE_ALGTOOL(LVL1, PPrStabilityMon)  
  DECLARE_NAMESPACE_ALGTOOL(LVL1, PPrSpareMon)
  DECLARE_NAMESPACE_ALGTOOL(LVL1, PPMSimBSMon)
  DECLARE_NAMESPACE_ALGTOOL(LVL1, RODMon)
  
  //Run 1
  DECLARE_NAMESPACE_ALGTOOL(LVL1, CMMMon)
  DECLARE_NAMESPACE_ALGTOOL(LVL1, CPMSimBSMon)
  DECLARE_NAMESPACE_ALGTOOL(LVL1, JEMMon)
  //DECLARE_NAMESPACE_ALGTOOL(LVL1, JEPSimBSMon)
  DECLARE_NAMESPACE_ALGTOOL(LVL1, TrigT1CaloBSMon)
  DECLARE_NAMESPACE_ALGTOOL(LVL1, TrigT1CaloCpmMonTool)
  DECLARE_NAMESPACE_ALGTOOL(LVL1, TrigT1CaloGlobalMonTool)
  DECLARE_NAMESPACE_ALGTOOL(LVL1, EmEfficienciesMonTool)
  DECLARE_NAMESPACE_ALGTOOL(LVL1, JetEfficienciesMonTool)
  DECLARE_NAMESPACE_ALGTOOL(LVL1, RODMonV1)
}

