# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

## @file TriggerUnixStandardSetup.py
## @brief py-module to configure the Athena AppMgr for trigger
## @author Werner Wiedenmann <Werner.Wiedenmann@cern.ch>
###############################################################

class _Conf:
    """Some configuration flags for this module with defaults"""
    useOnlineTHistSvc = True    # set in athenaHLT.py

def setupCommonServices(flags):
    from AthenaCommon import CfgMgr
    from AthenaCommon.Logging import logging
    from AthenaCommon.Constants import INFO
    from AthenaCommon.AppMgr import ServiceMgr as svcMgr, theApp
    from AthenaCommon.ConcurrencyFlags import jobproperties as jps

    # Setup messaging for Python and C++
    from AthenaCommon.Logging import log
    log.setFormat("%(asctime)s  Py:%(name)-31s %(levelname)7s %(message)s")

    # Create our own logger
    log = logging.getLogger( 'TriggerUnixStandardSetup::setupCommonServices:' )

    # Do the default Atlas job configuration first
    import AthenaCommon.AtlasUnixStandardJob   # noqa: F401

    # Now do HLT/thread specific configuration (see e.g. AtlasThreadedJob.py)
    from StoreGate.StoreGateConf import SG__HiveMgrSvc
    svcMgr += SG__HiveMgrSvc("EventDataSvc",
                             NSlots = jps.ConcurrencyFlags.NumConcurrentEvents())

    import StoreGate.StoreGateConf as StoreGateConf
    svcMgr += StoreGateConf.StoreGateSvc("ConditionStore")

    # Configure the CoreDumpSvc
    if not hasattr(svcMgr, "CoreDumpSvc"):
        from AthenaServices.Configurables import CoreDumpSvc
        svcMgr += CoreDumpSvc()

    # ThreadPoolService thread local initialization
    from GaudiHive.GaudiHiveConf import ThreadPoolSvc
    svcMgr += ThreadPoolSvc("ThreadPoolSvc")
    svcMgr.ThreadPoolSvc.ThreadInitTools = ["ThreadInitTool"]

    from GaudiHive.GaudiHiveConf import AlgResourcePool
    svcMgr += AlgResourcePool( OutputLevel = INFO,
                               TopAlg=["AthSequencer/AthMasterSeq"] )

    from AthenaCommon.AlgSequence import AlgSequence
    from SGComps.SGCompsConf import SGInputLoader
    topSequence = AlgSequence()
    topSequence += SGInputLoader(FailIfNoProxy = True)

    from AthenaCommon.AlgScheduler import AlgScheduler
    AlgScheduler.ShowDataDependencies(False)
    AlgScheduler.ShowControlFlow(False)

    # Setup SGCommitAuditor to sweep new DataObjects at end of Alg execute
    theApp.AuditAlgorithms = True
    from SGComps.SGCompsConf import SGCommitAuditor
    svcMgr.AuditorSvc += SGCommitAuditor()

    # setup ROOT6
    from PyUtils.Helpers import ROOT6Setup
    ROOT6Setup(batch=True)

    # Setup online THistSvc unless specifically configured otherwise
    #    setup the THistSvc early and force the creation of the THistSvc 
    #    so that it can be used by infrastructure services to book histograms  
    #    (to avoid problems e.g. with histograms in ROBDataProviderSvc)
    if _Conf.useOnlineTHistSvc:
        if hasattr(svcMgr, 'THistSvc'):
            log.fatal("The offline histogramming THistSvc is already in place.")
            raise RuntimeError("Cannot setup online histogramming TrigMonTHistSvc")
        log.debug("Using online histogramming service (TrigMonTHistSvc)")
        from TrigServices.TrigServicesConf import TrigMonTHistSvc
        svcMgr += TrigMonTHistSvc("THistSvc")
    else:
        log.debug("Using offline histogramming service (THistSvc)")
        from GaudiSvc.GaudiSvcConf import THistSvc
        svcMgr += THistSvc()

    # Online event loop manager
    from AthenaConfiguration.ComponentAccumulator import CAtoGlobalWrapper
    from TrigServices.TrigServicesConfig import TrigServicesCfg
    CAtoGlobalWrapper(TrigServicesCfg, flags)
    svcMgr.HltEventLoopMgr.WhiteboardSvc = "EventDataSvc"
    svcMgr.HltEventLoopMgr.SchedulerSvc = AlgScheduler.getScheduler().getName()

    # StoreGateSvc
    svcMgr.StoreGateSvc.ActivateHistory = False
    
    # Initialization of DetDescrCnvSvc
    svcMgr += CfgMgr.DetDescrCnvSvc(
        # specify primary Identifier dictionary to be used
        IdDictName = "IdDictParser/ATLAS_IDS.xml")

    theApp.CreateSvc += [ svcMgr.DetDescrCnvSvc.getFullName() ]
    svcMgr.EventPersistencySvc.CnvServices += [ "DetDescrCnvSvc" ]

    # Configuration of Interval of Validity Service
    svcMgr += CfgMgr.IOVSvc()
    
    # Explicitly set a few OutputLevels (needed because some services are created in
    # different order when running with the PSC)
    svcMgr.IncidentSvc.OutputLevel = theApp.OutputLevel
    svcMgr.ProxyProviderSvc.OutputLevel = theApp.OutputLevel
    svcMgr.StoreGateSvc.OutputLevel = theApp.OutputLevel
    
    return


def setupCommonServicesEnd():
    from AthenaCommon.AppMgr import ServiceMgr as svcMgr, athCondSeq
    from AthenaCommon.Logging import logging
    from AthenaCommon.AlgSequence import AlgSequence

    log = logging.getLogger( 'TriggerUnixStandardSetup::setupCommonServicesEnd:' )
    topSequence = AlgSequence()

    # --- create the ByteStreamCnvSvc after the Detector Description otherwise
    # --- the initialization of converters fails
    #from AthenaCommon.AppMgr import theApp
    #theApp.CreateSvc += [ svcMgr.ByteStreamCnvSvc.getFullName() ]

    # Make sure no THistSvc output/input stream is defined for online running
    if _Conf.useOnlineTHistSvc:
        svcMgr.THistSvc.Output = []
        if len(svcMgr.THistSvc.Input)>0:
            log.error('THistSvc.Input = %s. Input not allowed for online running. Disabling input.', svcMgr.THistSvc.Input)
            svcMgr.THistSvc.Input = []

    # For offline running make sure at least the EXPERT stream is defined
    else:
        if 1 not in [ o.count('EXPERT') for o in svcMgr.THistSvc.Output ]:
            svcMgr.THistSvc.Output += ["EXPERT DATAFILE='expert-monitoring.root' OPT='RECREATE'"]

    # Basic operational monitoring
    from TrigOnlineMonitor.TrigOnlineMonitorConfig import TrigOpMonitor
    topSequence += TrigOpMonitor()

    # Set default properties for some important services after all user job options
    log.info('Configure core services for online running')

    from TrigServices.TrigServicesConfig import setupMessageSvc
    setupMessageSvc()

    from TrigServices.TrigServicesConfig import enableCOOLFolderUpdates
    enableCOOLFolderUpdates()

    svcMgr.CoreDumpSvc.CoreDumpStream = "stdout"
    svcMgr.CoreDumpSvc.CallOldHandler = False  # avoid calling e.g. ROOT signal handler
    svcMgr.CoreDumpSvc.FastStackTrace = True   # first produce a fast stacktrace
    svcMgr.CoreDumpSvc.StackTrace = True       # then produce full stacktrace using gdb
    svcMgr.CoreDumpSvc.DumpCoreFile = True     # also produce core file (if allowed by ulimit -c)
    svcMgr.CoreDumpSvc.FatalHandler = 0        # no extra fatal handler
    svcMgr.CoreDumpSvc.TimeOut = 120000000000   # timeout for stack trace generation changed to 120s (ATR-17112,ATR-25404)

    svcMgr.IOVSvc.updateInterval = "RUN"
    svcMgr.IOVSvc.preLoadData = True
    svcMgr.IOVSvc.preLoadExtensibleFolders = False  # ATR-19392
    svcMgr.IOVSvc.forceResetAtBeginRun = False

    if hasattr(svcMgr,'IOVDbSvc'):
        svcMgr.IOVDbSvc.CacheAlign = 0  # VERY IMPORTANT to get unique queries for folder udpates (see Savannah #81092)
        svcMgr.IOVDbSvc.CacheRun = 0
        svcMgr.IOVDbSvc.CacheTime = 0

    if hasattr(athCondSeq, 'AtlasFieldMapCondAlg'):
        athCondSeq.AtlasFieldMapCondAlg.LoadMapOnStart = True

    return
