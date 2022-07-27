from AthenaCommon.Logging import logging
BSFilterLog = logging.getLogger('BSFilter')
BSFilterLog.info( '****************** STARTING BSFilter *****************' )

BSFilterLog.info( '**** Transformation run arguments' )
BSFilterLog.info( str(runArgs) )

if hasattr( runArgs, 'trfSubstepName'):
    if runArgs.trfSubstepName=="BSFilter" and hasattr(runArgs, 'InputLbnMapFile') and hasattr( runArgs, "triggerBit"):
        delattr(runArgs, 'InputLbnMapFile')
        BSFilterLog.info( '**** Removed InputLbnMapFile argument: transformation run arguments are now...' )
        BSFilterLog.info( str(runArgs) )

#---------------------------
theApp.EvtMax = runArgs.maxEvents
svcMgr.MessageSvc.infoLimit=100000
svcMgr.MessageSvc.defaultLimit=100000

from AthenaCommon.GlobalFlags  import globalflags
globalflags.InputFormat = "bytestream"
globalflags.DataSource = "data"
from AthenaCommon.AthenaCommonFlags  import athenaCommonFlags
if hasattr( runArgs, 'inputZeroBiasBSFile'):
    athenaCommonFlags.BSRDOInput=runArgs.inputZeroBiasBSFile
else:
    athenaCommonFlags.BSRDOInput=runArgs.inputBS_SKIMFile

#---------------------------
## Run performance monitoring (memory logging)
from PerfMonComps.PerfMonFlags import jobproperties as perfmonjp
perfmonjp.PerfMonFlags.doMonitoring = True # to enable monitoring
perfmonjp.PerfMonFlags.doFastMon = True    # to only enable a lightweight monitoring
#perfmonjp.PerfMonFlags.doSemiDetailedMonitoring = True

#---------------------------
BSFilterLog.info( '**** ByteStreamInputSvc configuration' )

include( "ByteStreamCnvSvc/BSEventStorageEventSelector_jobOptions.py" )
ByteStreamInputSvc = svcMgr.ByteStreamInputSvc
# svcMgr.EventSelector.Input = open(runArgs.InputFileMapFile).readline().rstrip().split(',')
if hasattr( runArgs, 'inputZeroBiasBSFile'):
    svcMgr.EventSelector.Input=runArgs.inputZeroBiasBSFile
else:
    svcMgr.EventSelector.Input=runArgs.inputBS_SKIMFile
print (ByteStreamInputSvc)

# ---------------------------
# Service to write out BS events
BSFilterLog.info( '**** ByteStreamOutputSvc configuration' )
for n in range(0,runArgs.noutputs):
    from ByteStreamCnvSvc.ByteStreamCnvSvcConf import ByteStreamEventStorageOutputSvc
    if hasattr( runArgs, 'outputBS_SKIMFile'):
        bsOutputSvc=ByteStreamEventStorageOutputSvc("BSESOutputSvc",OutputDirectory='./',SimpleFileName=runArgs.outputBS_SKIMFile)
    else:
        if n==0: myn=""
        else: myn=str(n)
        bsOutputSvc=ByteStreamEventStorageOutputSvc("BSESOutputSvc"+myn,OutputDirectory='./',SimpleFileName=getattr(runArgs,"outputBS_TRIGSKIM"+myn+"File"))
    svcMgr += bsOutputSvc
    print (bsOutputSvc)

# ---------------------------
BSFilterLog.info( '**** ByteStreamFilter configuration' )

from AthenaCommon.AlgSequence import AlgSequence
topSequence = AlgSequence()

# main alg
from OverlayByteStreamUtils.OverlayByteStreamUtilsConf import BSFilter
filAlg = BSFilter("BSFilter")

# convert flags to ConfigFlags
from AthenaConfiguration.OldFlags2NewFlags import getNewConfigFlags
ConfigFlags = getNewConfigFlags()

if hasattr( runArgs, "triggerBit"):
    filAlg.L1Trigger = runArgs.triggerBit
    # We don't need any Calo/Muon RoI information. The following hack
    # disables the decoding of it in RoIBResultToxAOD:
    ConfigFlags.Detector.EnableMuon = False
    ConfigFlags.Detector.EnableCalo = False
    ConfigFlags.lock()
    from TriggerJobOpts.TriggerRecoGetter import TriggerRecoGetter
    triggerGetter = TriggerRecoGetter()
else:
    filAlg.L1Trigger = ""

if hasattr(runArgs, 'outputTXT_EVENTIDFile'):
    filAlg.EventIdFile=runArgs.outputTXT_EVENTIDFile # The name of the file to write to for EventIdModifierSvc lines
else:
    filAlg.EventIdFile=""

if hasattr( runArgs, "inputFilterFile") and not hasattr( runArgs, 'InputLbnMapFile'): #don't filter if you're skimming using the map file
    filAlg.filterfile=runArgs.inputFilterFile # The thing made from the TAG files via "root -l -b -q HITAGprinter_run.C"
else:
    filAlg.filterfile = ""

topSequence += filAlg

# ---------------------------
BSFilterLog.info( '**** ByteStreamCopyTool configuration' )

if hasattr( runArgs, 'InputLbnMapFile'):
    from OverlayByteStreamUtils.OverlayByteStreamUtilsConf import ByteStreamMultipleOutputStreamCopyTool
    bsCopyTool = ByteStreamMultipleOutputStreamCopyTool("MultipleOutputStreamBSCopyTool")
    bsCopyTool.lbn_map_file = runArgs.InputLbnMapFile
    if hasattr( runArgs, "inputFilterFile"): bsCopyTool.trigfile = runArgs.inputFilterFile
    else: bsCopyTool.trigfile = ""
    bsCopyTool.NoutputSvc = runArgs.noutputs
    for n in range(0,runArgs.noutputs):
        if n==0: myn=""
        else: myn=str(n)
        setattr(bsCopyTool,"ByteStreamOutputSvc"+str(n),getattr(svcMgr,"BSESOutputSvc"+myn))
else:
    from ByteStreamCnvSvc.ByteStreamCnvSvcConf import ByteStreamOutputStreamCopyTool
    bsCopyTool = ByteStreamOutputStreamCopyTool("OutputStreamBSCopyTool")
    bsCopyTool.ByteStreamOutputSvc=bsOutputSvc

svcMgr.ToolSvc += bsCopyTool

bsCopyTool.ByteStreamInputSvc=svcMgr.ByteStreamInputSvc

# ---------------------------
BSFilterLog.info( '**** AthenaOutputStream configuration' )

# create AthenaOutputStream for BS Copy 
from AthenaServices.AthenaServicesConf import AthenaOutputStream
OutputStreamBSCopy = AthenaOutputStream( "OutputStreamBSCopy", WritingTool=bsCopyTool )
topSequence += OutputStreamBSCopy

OutputStreamBSCopy.AcceptAlgs =["BSFilter"] 

print (topSequence)

# ---------------------------
# Post-include/exec
if hasattr(runArgs,"postInclude"):
    for fragment in runArgs.postInclude:
        include(fragment)
if hasattr(runArgs, "postExec") and runArgs.postExec != 'NONE':
    for cmd in runArgs.postExec:
        exec(cmd)
