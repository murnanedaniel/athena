
from AthenaCommon.Logging import logging 
log = logging.getLogger("TrigInDetValidation")

log.info( "preinclude: TIDAvtx_preinclude.py" ) 

from TrigInDetConfig.ConfigSettings import getInDetTrigConfig

getInDetTrigConfig("jet")._addSingleTrackVertices = True
getInDetTrigConfig("jet")._minNSiHits_vtx = 8

log.info( "ID Trigger addSingleVertices: "+str(getInDetTrigConfig("jet").addSingleTrackVertices) )
log.info( "ID Trigger minNSiHits:        "+str(getInDetTrigConfig("jet").minNSiHits_vtx) )


