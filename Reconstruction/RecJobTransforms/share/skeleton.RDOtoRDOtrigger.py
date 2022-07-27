# Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
#
# Transform skeleton for RDO -> RDO_TRIG (running trigger and adding its output to the RDO file).
# This is only a wrapper interfacing transform arguments into the main job options file
# TriggerJobOpts/runHLT_standalone.py
#

from AthenaCommon.Include import include
from AthenaCommon.Logging import logging
from AthenaCommon.AthenaCommonFlags import athenaCommonFlags
from AthenaConfiguration.AllConfigFlags import ConfigFlags
from AthenaConfiguration.ComponentAccumulator import CAtoGlobalWrapper
from PerfMonComps.PerfMonCompsConfig import PerfMonMTSvcCfg
from PerfMonComps.PerfMonConfigHelpers import setPerfmonFlagsFromRunArgs

log = logging.getLogger('skeleton.RDOtoRDOtrigger')

##################################################
# Ensure AthenaMT is used as running Trigger requires the MT scheduler
##################################################
from AthenaCommon.ConcurrencyFlags import jobproperties as cfjp
if cfjp.ConcurrencyFlags.NumThreads() == 0 and ConfigFlags.Concurrency.NumThreads == 0:
    raise RuntimeError('RDOtoRDOTrigger requires AthenaMT, but serial Athena was used. Please use threads=1 or more')

##################################################
# Options read by runHLT_standalone, can be overwritten in runArgs/preExec/preInclude
##################################################
doWriteRDOTrigger = True
doWriteBS         = False

##################################################
# Parse runArgs
##################################################
if 'runArgs' not in globals():
    raise RuntimeError('runArgs not defined')

ra = globals()['runArgs']
def getFromRunArgs(propName, failIfMissing=True):
    if not hasattr(ra, propName):
        if failIfMissing:
            raise RuntimeError(propName + ' not defined in runArgs')
        else:
            return None
    return getattr(ra, propName)

# Input/Output
athenaCommonFlags.FilesInput = getFromRunArgs('inputRDOFile')
ConfigFlags.Input.Files = getFromRunArgs('inputRDOFile')
ConfigFlags.Output.RDOFileName = getFromRunArgs('outputRDO_TRIGFile')

# Max/skip events
maxEvents = getFromRunArgs('maxEvents', False) or -1
skipEvents = getFromRunArgs('skipEvents', False) or 0
athenaCommonFlags.EvtMax = maxEvents
athenaCommonFlags.SkipEvents = skipEvents

#conditions/geometry setup for runHLT_standalone
setGlobalTag = getFromRunArgs("conditionsTag", False)
setDetDescr = getFromRunArgs("geometryVersion", False)

# PerfMon setup (ATR-25439)
setPerfmonFlagsFromRunArgs(ConfigFlags, ra)

##################################################
# Parse preExec / preInclude
##################################################
preExec = getFromRunArgs('preExec', False)
if preExec:
    log.info('Executing transform preExec')
    for cmd in preExec:
        log.info(cmd)
        exec(cmd)

preInclude = getFromRunArgs('preInclude', False)
if preInclude:
    log.info('Executing transform preInclude')
    for fragment in preInclude:
        include(fragment)

##################################################
# Include the main job options
##################################################
include("TriggerJobOpts/runHLT_standalone.py")

##################################################
# Include PerfMon configuration (ATR-25439)
##################################################
# Translate old concurrency flags to new for PerfMon
# - do that only here as runHLT_standalone must be independent of these flags
flagsForPerfMon = ConfigFlags.clone()
flagsForPerfMon.Concurrency.NumProcs = cfjp.ConcurrencyFlags.NumProcs()
flagsForPerfMon.Concurrency.NumThreads = cfjp.ConcurrencyFlags.NumThreads()
flagsForPerfMon.Concurrency.NumConcurrentEvents = cfjp.ConcurrencyFlags.NumConcurrentEvents()
flagsForPerfMon.lock()
CAtoGlobalWrapper(PerfMonMTSvcCfg, flagsForPerfMon)

##################################################
# Parse postExec / postInclude
##################################################
postExec = getFromRunArgs('postExec', failIfMissing=False)
if postExec:
    log.info('Executing transform postExec')
    for cmd in postExec:
        log.info(cmd)
        exec(cmd)

postInclude = getFromRunArgs('postInclude', failIfMissing=False)
if postInclude:
    log.info('Executing transform postInclude')
    for fragment in postInclude:
        include(fragment)
