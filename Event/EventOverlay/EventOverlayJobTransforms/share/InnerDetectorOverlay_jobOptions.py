include.block ( "EventOverlayJobTransforms/InnerDetectorOverlay_jobOptions.py" )

from Digitization.DigitizationFlags import digitizationFlags
from AthenaCommon.DetFlags import DetFlags
from AthenaCommon import CfgGetter
from OverlayCommonAlgs.OverlayFlags import overlayFlags


if DetFlags.overlay.pixel_on() or DetFlags.overlay.SCT_on() or DetFlags.overlay.TRT_on():

    digitizationFlags.doInDetNoise=False # FIXME THIS SHOULD BE SET EARLIER IN THE CONFIGURATION

    #if overlayFlags.isDataOverlay():
    #   include( "InDetCosmicRecExample/InDetCosmicFlags_jobOptions.py" )

    if DetFlags.overlay.pixel_on():
        job += CfgGetter.getAlgorithm("PixelOverlayDigitization")
        if not overlayFlags.doTrackOverlay():
            job += CfgGetter.getAlgorithm("PixelOverlay")
            if DetFlags.overlay.Truth_on():
                job += CfgGetter.getAlgorithm("PixelSDOOverlay")
        else:
            job.PixelOverlayDigitization.DigitizationTool.RDOCollName="PixelRDOs"
            job.PixelOverlayDigitization.DigitizationTool.SDOCollName="PixelSDO_Map"

        if overlayFlags.isDataOverlay():
            if overlayFlags.isOverlayMT():
                job.InDetPixelRawDataProvider.RDOKey = overlayFlags.bkgPrefix()+"PixelRDOs"
            else:
                job.InDetPixelRawDataProvider.RDOKey = overlayFlags.dataStore()+"+PixelRDOs"

            from RegionSelector.RegSelToolConfig import makeRegSelTool_Pixel
            job.InDetPixelRawDataProvider.RegSelTool = makeRegSelTool_Pixel()

            #ServiceMgr.ByteStreamAddressProviderSvc.TypeNames += [ "PixelRDO_Container/PixelRDOs" ]
            #ServiceMgr.ByteStreamAddressProviderSvc.TypeNames += [ "Trk::PixelClusterContainer/PixelOnlineClusters" ]

    if DetFlags.overlay.SCT_on():

        # Setup the ReadCalibChip folders and Svc
        if overlayFlags.isDataOverlay():
            #conddb.blockFolder("/SCT/DAQ/Calibration/ChipGain")
            #conddb.blockFolder("/SCT/DAQ/Calibration/ChipNoise")
            #conddb.addFolder("SCT_OFL","/SCT/DAQ/Calibration/ChipGain",forceMC=True)
            #conddb.addFolder("SCT_OFL","/SCT/DAQ/Calibration/ChipNoise",forceMC=True)
            conddb.addFolder("SCT_OFL","/SCT/DAQ/Calibration/ChipGain<tag>SctDaqCalibrationChipGain-Apr10-01</tag>",forceMC=True, className="CondAttrListCollection")
            conddb.addFolder("SCT_OFL","/SCT/DAQ/Calibration/ChipNoise<tag>SctDaqCalibrationChipNoise-Apr10-01</tag>",forceMC=True, className="CondAttrListCollection")

            #if not conddb.folderRequested('/SCT/DAQ/Calibration/ChipGain'):
            #   conddb.addFolderSplitOnline("SCT","/SCT/DAQ/Calibration/ChipGain","/SCT/DAQ/Calibration/ChipGain",forceMC=True)
            #if not conddb.folderRequested('/SCT/DAQ/Calibration/ChipNoise'):
            #   conddb.addFolderSplitOnline("SCT","/SCT/DAQ/Calibration/ChipNoise","/SCT/DAQ/Calibration/ChipNoise",forceMC=True)

        job += CfgGetter.getAlgorithm("SCT_OverlayDigitization")
        if not overlayFlags.doTrackOverlay():
            job += CfgGetter.getAlgorithm("SCTOverlay")
            if DetFlags.overlay.Truth_on():
                job += CfgGetter.getAlgorithm("SCTSDOOverlay")
        else:
            job.SCT_OverlayDigitization.DigitizationTool.OutputObjectName="SCT_RDOs"
            job.SCT_OverlayDigitization.DigitizationTool.OutputSDOName="SCT_SDO_Map"            

        if overlayFlags.isDataOverlay():
            include("InDetRecExample/InDetRecConditionsAccess.py")

            if overlayFlags.isOverlayMT():
                job.InDetSCTRawDataProvider.RDOKey = overlayFlags.bkgPrefix()+"SCT_RDOs"
                job.InDetSCTRawDataProvider.LVL1IDKey = overlayFlags.bkgPrefix()+"SCT_LVL1ID"
                job.InDetSCTRawDataProvider.BCIDKey = overlayFlags.bkgPrefix()+"SCT_BCID"
            else:
                job.InDetSCTRawDataProvider.RDOKey = overlayFlags.dataStore()+"+SCT_RDOs"
                job.InDetSCTRawDataProvider.LVL1IDKey = overlayFlags.dataStore()+"+SCT_LVL1ID"
                job.InDetSCTRawDataProvider.BCIDKey = overlayFlags.dataStore()+"+SCT_BCID"
            #ServiceMgr.ByteStreamAddressProviderSvc.TypeNames += [ "SCT_RDO_Container/SCT_RDOs" ]
            #ServiceMgr.ByteStreamAddressProviderSvc.TypeNames += [ "Trk::SCT_ClusterContainer/SCT_OnlineClusters" ]

    if DetFlags.overlay.TRT_on():
        if overlayFlags.isDataOverlay():
            conddb.blockFolder("/TRT/Cond/DigVers")
            #conddb.addFolderWithTag("TRT_OFL","/TRT/Cond/DigVers","TRTCondDigVers-Collisions-01",force=True,forceMC=True)
            conddb.addFolder("TRT_OFL","/TRT/Cond/DigVers",forceMC=True,
                             className = 'AthenaAttributeList')

        from TRT_ElectronPidTools.TRT_ElectronPidToolsConf import InDet__TRT_LocalOccupancy
        TRT_LocalOccupancy = InDet__TRT_LocalOccupancy(name="TRT_LocalOccupancy", isTrigger= False,
                                                       TRT_DriftCircleCollection="")

        job += CfgGetter.getAlgorithm("TRT_OverlayDigitization")
        if not overlayFlags.doTrackOverlay():
            job += CfgGetter.getAlgorithm("TRTOverlay")
            job.TRTOverlay.TRT_LocalOccupancyTool  = TRT_LocalOccupancy
            if DetFlags.overlay.Truth_on():
                job += CfgGetter.getAlgorithm("TRTSDOOverlay")
        else:
            job.TRT_OverlayDigitization.DigitizationTool.OutputObjectName="TRT_RDOs"
            job.TRT_OverlayDigitization.DigitizationTool.OutputSDOName="TRT_SDO_Map"
            job.TRT_OverlayDigitization.DigitizationTool.Override_isOverlay=0


        from InDetRecExample.InDetJobProperties import InDetFlags
        include("InDetRecExample/InDetRecConditionsAccess.py")
        
        if overlayFlags.isDataOverlay():
            if overlayFlags.isOverlayMT():
                job.InDetTRTRawDataProvider.RDOKey = overlayFlags.bkgPrefix()+"TRT_RDOs"
            else:
                job.InDetTRTRawDataProvider.RDOKey = overlayFlags.dataStore()+"+TRT_RDOs"
            #ServiceMgr.ByteStreamAddressProviderSvc.TypeNames += [ "TRT_RDO_Container/TRT_RDOs" ]

            from RegionSelector.RegSelToolConfig import makeRegSelTool_TRT
            InDetTRTRawDataProvider.RegSelTool = makeRegSelTool_TRT()

            #from IOVDbSvc.CondDB import conddb
            #conddb.addFolder("TRT","/TRT/Calib/T0","<tag>TrtCalibt0-UPD2-FDR2-01</tag>")
            #conddb.addFolder("TRT","/TRT/Calib/RT","<tag>TrtCalibRt-UPD2-FDR2-01</tag>")
            #conddb.addFolder("TRT","/TRT/Calib/T0","<tag>TrtCalibRt-HLT-UPD1-01</tag>")
            #conddb.addFolder("TRT","/TRT/Calib/RT","<tag>TrtCalibT0-HLT-UPD1-01</tag>")
            conddb.addFolder("TRT_ONL","/TRT/Onl/ROD/Compress",className='CondAttrListCollection')
