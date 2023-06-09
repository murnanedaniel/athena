from AthenaCommon.Logging import logging
overlaylog = logging.getLogger('overlay')
overlaylog.info( '****************** STARTING OVERLAY *****************' )

overlaylog.info( '**** Transformation run arguments' )
overlaylog.info( str(runArgs) )

from PerfMonComps.PerfMonFlags import jobproperties as jp
jp.PerfMonFlags.doMonitoring = True # to enable monitoring
jp.PerfMonFlags.doFastMon = True    # to only enable a lightweight monitoring

if hasattr(runArgs, "preExec") and runArgs.preExec != 'NONE':
    for cmd in runArgs.preExec:
        exec(cmd)

import AthenaCommon.AtlasUnixStandardJob

if hasattr(runArgs, 'preInclude'):
for cf in runArgs.preInclude:
    include(cf)

#==============================================================
# Job definition parameters:
#==============================================================
from AthenaCommon.GlobalFlags  import globalflags

globalflags.isOverlay.set_Value_and_Lock(True)

if hasattr(runArgs,"maxEvents"): from AthenaCommon.AthenaCommonFlags  import athenaCommonFlags
athenaCommonFlags.EvtMax = runArgs.maxEvents
if hasattr(runArgs,"skipEvents"): athenaCommonFlags.SkipEvents= runArgs.skipEvents

if hasattr(runArgs,"inputHITSFile"):
    athenaCommonFlags.PoolHitsInput.set_Value_and_Lock( runArgs.inputHITSFile )
    SignalInputCollections = runArgs.inputHITSFile
else:
    raise RuntimeError ("No input HITS file defined")

if hasattr(runArgs,"outputRDOFile"): 
    athenaCommonFlags.PoolRDOOutput.set_Value_and_Lock( runArgs.outputRDOFile )
    OverlayCollection = runArgs.outputRDOFile

if hasattr(runArgs,"geometryVersion"): globalFlags.DetDescrVersion=runArgs.geometryVersion

from Digitization.DigitizationFlags import digitizationFlags
if hasattr(runArgs,"digiSeedOffset1"): digitizationFlags.rndmSeedOffset1=int(runArgs.digiSeedOffset1)
if hasattr(runArgs,"digiSeedOffset2"): digitizationFlags.rndmSeedOffset2=int(runArgs.digiSeedOffset2)
if hasattr(runArgs,"samplingFractionDbTag"): digitizationFlags.physicsList=runArgs.samplingFractionDbTag
if hasattr(runArgs,"digiRndmSvc"): digitizationFlags.rndmSvc=runArgs.digiRndmSvc

from AthenaCommon.DetFlags import DetFlags
if not hasattr(runArgs, 'outputRDO_SGNLFile') or runArgs.outputRDO_SGNLFile=="NONE":
   overlayFlags.doSignal=False
   SignalCollection = "NONE"
else:
    overlayFlags.doSignal=True
    SignalCollection = runArgs.outputRDO_SGNLFile

DetFlags.overlay.BCM_setOff()
DetFlags.overlay.Lucid_setOff()

readBS = False
if hasattr(runArgs, 'ReadByteStream'):
    readBS = runArgs.ReadByteStream
isRealData = False

if readBS:
   globalflags.InputFormat.set_Value_and_Lock('bytestream')
   DataInputCollections=runArgs.inputPileUpBSFile
else:
   DataInputCollections=runArgs.inputRDO_BKGFile
   athenaCommonFlags.PoolRDOInput=runArgs.inputRDO_BKGFile

if hasattr(runArgs, 'conditionsTag') and runArgs.conditionsTag!='NONE' and runArgs.conditionsTag!='':
   globalflags.ConditionsTag=runArgs.conditionsTag
   if len(globalflags.ConditionsTag())!=0:
      from IOVDbSvc.CondDB import conddb
      conddb.setGlobalTag(globalflags.ConditionsTag())

DetFlags.Print()

# LVL1 Trigger Menu
if hasattr(runArgs, "triggerConfig") and runArgs.triggerConfig!="NONE":
    # LVL1 Trigger Menu
    # PJB 9/2/2009 Setup the new triggerConfig flags here
    from TriggerJobOpts.TriggerFlags import TriggerFlags
    triggerArg = runArgs.triggerConfig
    #if not prefixed with LVL1: add it here
    Args = triggerArg.split(":")
    if Args[0] != "LVL1":
        TriggerFlags.triggerConfig ="LVL1:"+triggerArg
    else:
        TriggerFlags.triggerConfig =triggerArg
    overlaylog.info( 'triggerConfig argument is: %s ', TriggerFlags.triggerConfig.get_Value() )
    from TriggerJobOpts.TriggerConfigGetter import TriggerConfigGetter
    cfg = TriggerConfigGetter("HIT2RDO")

from AthenaCommon.DetFlags import DetFlags
DetFlags.ID_setOn()
DetFlags.Muon_setOn()
DetFlags.LAr_setOn()
DetFlags.Tile_setOn()
if hasattr(runArgs, "triggerConfig") and runArgs.triggerConfig=="NONE":
  DetFlags.LVL1_setOff()
else:
  DetFlags.LVL1_setOn()

DetFlags.BCM_setOn()
DetFlags.Lucid_on()

DetFlags.simulateLVL1.Lucid_setOff()

print "================ DetFlags ================ "
DetFlags.Print()

# Geometry, controlled by DetFlags
#GlobalFlags.DetGeo.set_atlas()
#GlobalFlags.DataSource.set_geant4()
globalflags.DataSource.set_Value_and_Lock('geant4')

#from AtlasGeoModel import SetGeometryVersion
#from AtlasGeoModel import GeoModelInit

include ( "RecExCond/AllDet_detDescr.py" )

#from AtlasGeoModel import SetGeometryVersion
#from AtlasGeoModel import GeoModelInit
#from AtlasGeoModel import SetupRecoGeometry

#include( "BFieldAth/BFieldAth_jobOptions.py" )

#--------
from AthenaCommon.AppMgr import theApp
theApp.EventLoop = "PileUpEventLoopMgr"
if hasattr(runArgs,"maxEvents"): theApp.EvtMax = runArgs.maxEvents

include ( "EventOverlayJobTransforms/ConfiguredOverlay_jobOptions.py" )

if DetFlags.overlay.Truth_on():
   include ( "EventOverlayJobTransforms/TruthOverlay_jobOptions.py" )

if DetFlags.overlay.BCM_on() or DetFlags.overlay.Lucid_on():
   include ( "EventOverlayJobTransforms/BeamOverlay_jobOptions.py" )

if DetFlags.overlay.pixel_on() or DetFlags.overlay.SCT_on() or DetFlags.overlay.TRT_on():
   include ( "EventOverlayJobTransforms/InnerDetectorOverlay_jobOptions.py" )

if DetFlags.overlay.LAr_on() or DetFlags.overlay.Tile_on():
   include ( "EventOverlayJobTransforms/CaloOverlay_jobOptions.py" )

if DetFlags.overlay.CSC_on() or DetFlags.overlay.MDT_on() or DetFlags.overlay.RPC_on() or DetFlags.overlay.TGC_on():
   include ( "EventOverlayJobTransforms/MuonOverlay_jobOptions.py" )

if DetFlags.overlay.LVL1_on():
   include ( "EventOverlayJobTransforms/Level1Overlay_jobOptions.py" )

# save the overlay output first
include ( "EventOverlayJobTransforms/OverlayOutputItemList_jobOptions.py" )

# now save the signal information in the same job
if overlayFlags.doSignal==True:
   include ( "EventOverlayJobTransforms/SignalOutputItemList_jobOptions.py" )

# For random number initialization
from AthenaCommon.ConfigurableDb import getConfigurable
ServiceMgr += getConfigurable(digitizationFlags.rndmSvc.get_Value())()
digitizationFlags.rndmSeedList.addtoService()
digitizationFlags.rndmSeedList.printSeeds()

# To not overwrite the BCID
#from AthenaCommon.AppMgr import ServiceMgr
#ServiceMgr.PileUpEventLoopMgr.IsEventOverlayJob=True

#================================================================
print "overlay_trf: final outStream = ", outStream

#================================================================
# Logging
#--------------------------------------------------------------
#ServiceMgr.MessageSvc.OutputLevel = INFO
ServiceMgr.MessageSvc.OutputLevel = INFO
ServiceMgr.MessageSvc.Format = "% F%45W%S%7W%R%T %0W%M"

print "overlay_trf: at the end. job=\n", job
print "\n\noverlay_trf: at the end. ServiceMgr=\n", ServiceMgr

# Post-include
if hasattr(runArgs,"postInclude"):
    for fragment in runArgs.postInclude:
        include(fragment)

if hasattr(runArgs, "postExec") and runArgs.postExec != 'NONE':
    for cmd in runArgs.postExec:
        exec(cmd)
