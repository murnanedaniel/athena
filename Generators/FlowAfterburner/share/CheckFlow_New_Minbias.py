##########################################################
#
# 28/02/2006 Andrzej Olszewski
# jobOptions to run CheckFlow
# on a file with flow modified generated events
# modified by Cano Ay aycano@cern.ch 14/05/2008 
#
##########################################################
#import AthenaCommon.AtlasUnixGeneratorJob
import AthenaCommon.AtlasUnixStandardJob

import AthenaPoolCnvSvc.ReadAthenaPool
svcMgr.EventSelector.InputCollections= ["hijing.flow.10-15.pool.root"]
#svcMgr.EventSelector.InputCollections= ["/tmp/soumya/soumya.pool.root"]

from AthenaCommon.GlobalFlags import globalflags
from PartPropSvc.PartPropSvcConf import PartPropSvc
svcMgr += PartPropSvc()

#include( "AthenaCommon/AthenaCommonFlags.py" )


#------------------------------------------------------------
# Load main libraries
#------------------------------------------------------------
from AthenaCommon.AlgSequence import AlgSequence
job=AlgSequence()
from FlowAfterburner.FlowAfterburnerConf import CheckFlow_New_Minbias
job += CheckFlow_New_Minbias()

theApp.Dlls   += [ "GeneratorObjectsAthenaPoolPoolCnv" ]
#------------------------------------------------------------
# Number of events and OutputLevel
#------------------------------------------------------------
theApp.EvtMax                = -1
svcMgr.MessageSvc.OutputLevel               = 3

#------------------------------------------------------------
# Persistency services ntuples etc...
#------------------------------------------------------------
# THistService for native root in Athena
from GaudiSvc.GaudiSvcConf import THistSvc
svcMgr += THistSvc()
svcMgr.THistSvc.Output = ["FlowOutPut DATAFILE='CheckFlow.hist.root' OPT='RECREATE'"]

#theApp.Dlls += [ "GaudiSvc" ]
#THistSvc = Algorithm( "THistSvc" )
#THistSvc.Output = ["FlowOutPut DATAFILE='CheckFlow.hist.root' OPT='RECREATE'"]
        
#CheckFlow = Algorithm( "CheckFlow" )
CheckFlow_New_Minbias = job.CheckFlow_New_Minbias
CheckFlow_New_Minbias.McEventKey = "GEN_EVENT"
CheckFlow_New_Minbias.RapidityCutMax = 5.5
#==============================================================
#
# End of job options file
#
###############################################################
