include.block ( "EventOverlayJobTransforms/CaloOverlay_jobOptions.py" )

from Digitization.DigitizationFlags import jobproperties
from AthenaCommon.DetFlags import DetFlags
from OverlayCommonAlgs.OverlayFlags import overlayFlags


if DetFlags.overlay.LAr_on() or DetFlags.overlay.Tile_on():

   jobproperties.Digitization.doCaloNoise=False

   if overlayFlags.isDataOverlay():
      include("LArConditionsCommon/LArIdMap_comm_jobOptions.py")
      include("LArConditionsCommon/LArConditionsCommon_comm_jobOptions.py")
   else:
      include("LArConditionsCommon/LArIdMap_MC_jobOptions.py")
      include("LArConditionsCommon/LArConditionsCommon_MC_jobOptions.py")

   include( "LArIdCnv/LArIdCnv_joboptions.py" )
   include( "CaloDetMgrDetDescrCnv/CaloDetMgrDetDescrCnv_joboptions.py" )

if DetFlags.overlay.LAr_on():

    # FIXME this is for doing Overlay with MC RDO + MC hits
    #   real data RDO will require some changes in the setup for proper db acces
    #include( "LArDetDescr/LArDetDescr_joboptions.py" )
    #include( "LArAthenaPool/LArAthenaPool_joboptions.py" )
    # We also need the conditions svc for MC constants:
    if overlayFlags.isDataOverlay():
       theApp.Dlls += [ "LArByteStream"]
       LArDigitKey = "FREE"
       ServiceMgr.ByteStreamAddressProviderSvc.TypeNames += ["LArDigitContainer/"+LArDigitKey]
       ServiceMgr.ByteStreamAddressProviderSvc.TypeNames += ["LArFebHeaderContainer/LArFebHeader"]
       ServiceMgr.ByteStreamAddressProviderSvc.TypeNames+=["LArDigitContainer/LArDigitContainer_MC"]

    from LArROD.LArRawChannelGetter import LArRawChannelGetter
    LArRawChannelGetter()

    from LArROD.LArDigits import DefaultLArDigitThinner
    LArDigitThinner = DefaultLArDigitThinner('LArDigitThinner') # automatically added to topSequence
    if overlayFlags.isDataOverlay():
       job.LArDigitThinner.InputContainerName = overlayFlags.dataStore()+"+FREE"
       job.LArDigitThinner.RawChannelContainerName = "LArRawChannels_FromDigits"
       #job.digitmaker1.LArPileUpTool.OutputLevel=DEBUG
       #MessageSvc.debugLimit = 100000
       #job.digitmaker1.LArPileUpTool.useLArFloat=False
       job.digitmaker1.LArPileUpTool.PedestalKey = "LArPedestal"
       job.LArRawChannelBuilder.DataLocation = "LArDigitContainer_MC"

#----------------------------------------------------------------
if DetFlags.overlay.Tile_on():

    include( "TileIdCnv/TileIdCnv_jobOptions.py" )
    include( "TileConditions/TileConditions_jobOptions.py" )

    include( "TileSimAlgs/TileDigitization_jobOptions.py" )
    include( "TileL2Algs/TileL2Algs_jobOptions.py" )

    job.TileHitVecToCnt.DigitizationTool.RndmEvtOverlay = True
    theTileDigitsMaker.RndmEvtOverlay = True
    if overlayFlags.isDataOverlay():
       theApp.Dlls += [ "TileByteStream"]
       ServiceMgr.ByteStreamAddressProviderSvc.TypeNames += [ "TileBeamElemContainer/TileBeamElemCnt"]
       ServiceMgr.ByteStreamAddressProviderSvc.TypeNames += [ "TileRawChannelContainer/TileRawChannelCnt"]
       ServiceMgr.ByteStreamAddressProviderSvc.TypeNames += [ "TileDigitsContainer/TileDigitsCnt"]
       ServiceMgr.ByteStreamAddressProviderSvc.TypeNames += [ "TileL2Container/TileL2Cnt"]
       ServiceMgr.ByteStreamAddressProviderSvc.TypeNames += [ "TileLaserObject/TileLaserObj"]

#--------------------
