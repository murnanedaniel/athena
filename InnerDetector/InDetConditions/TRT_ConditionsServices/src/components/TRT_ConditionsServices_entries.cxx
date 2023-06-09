#include "GaudiKernel/DeclareFactoryEntries.h"
#include "src/TRT_ConditionsSummarySvc.h"
#include "src/TRT_ConditionsTestSvc.h"
#include "src/TRT_AlignDbSvc.h"
#include "src/TRT_CalDbSvc.h"
#include "src/TRT_StrawAlignDbSvc.h"
#include "src/TRT_DCS_ConditionsSvc.h"
#include "src/TRT_HWMappingSvc.h"
#include "src/TRT_StrawNeighbourSvc.h"
#include "src/TRT_StrawStatusSummarySvc.h"
#include "src/TRT_ByteStream_ConditionsSvc.h"
#include "src/TRT_DAQ_ConditionsSvc.h"
#include "src/TRT_ActiveFractionSvc.h"

DECLARE_SERVICE_FACTORY( TRT_ConditionsSummarySvc )
DECLARE_SERVICE_FACTORY( TRT_ConditionsTestSvc )
DECLARE_SERVICE_FACTORY( TRT_AlignDbSvc )
DECLARE_SERVICE_FACTORY( TRT_CalDbSvc )
DECLARE_SERVICE_FACTORY( TRT_StrawAlignDbSvc )
DECLARE_SERVICE_FACTORY( TRT_DCS_ConditionsSvc )
DECLARE_SERVICE_FACTORY( TRT_HWMappingSvc )
DECLARE_SERVICE_FACTORY( TRT_StrawNeighbourSvc )
DECLARE_SERVICE_FACTORY( TRT_StrawStatusSummarySvc )
DECLARE_SERVICE_FACTORY( TRT_ByteStream_ConditionsSvc )
DECLARE_SERVICE_FACTORY( TRT_DAQ_ConditionsSvc )
DECLARE_SERVICE_FACTORY( TRT_ActiveFractionSvc )

DECLARE_FACTORY_ENTRIES( TRT_ConditionsServices ) {
  DECLARE_SERVICE( TRT_ConditionsSummarySvc );
  DECLARE_SERVICE( TRT_ConditionsTestSvc );
  DECLARE_SERVICE( TRT_AlignDbSvc );
  DECLARE_SERVICE( TRT_CalDbSvc );
  DECLARE_SERVICE( TRT_StrawAlignDbSvc );
  DECLARE_SERVICE( TRT_DCS_ConditionsSvc );
  DECLARE_SERVICE( TRT_HWMappingSvc );
  DECLARE_SERVICE( TRT_StrawNeighbourSvc );
  DECLARE_SERVICE( TRT_StrawStatusSummarySvc );
  DECLARE_SERVICE( TRT_ByteStream_ConditionsSvc );
  DECLARE_SERVICE( TRT_DAQ_ConditionsSvc );
  DECLARE_SERVICE( TRT_ActiveFractionSvc );
}
