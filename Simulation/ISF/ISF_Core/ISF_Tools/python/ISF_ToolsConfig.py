# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration

"""
Tools configurations for ISF
KG Tan, 17/06/2012
"""

from AthenaCommon.CfgGetter import getPrivateTool,getPrivateToolClone,getPublicTool,getPublicToolClone,\
        getService,getServiceClone,getAlgorithm,getAlgorithmClone

from AthenaCommon.Constants import *  # FATAL,ERROR etc.
from AthenaCommon.SystemOfUnits import *
from AthenaCommon.DetFlags import DetFlags

from ISF_Config.ISF_jobProperties import ISF_Flags # IMPORTANT: Flags must be set before tools are retrieved


def getParticleHelper(name="ISF_ParticleHelper", **kwargs):
    kwargs.setdefault("BarcodeSvc"             , ISF_Flags.BarcodeService())
    from ISF_Tools.ISF_ToolsConf import ISF__ParticleHelper
    return ISF__ParticleHelper(name, **kwargs)

def getMemoryMonitor(name="ISF_MemoryMonitor", **kwargs):
    from ISF_Tools.ISF_ToolsConf import ISF__MemoryMonitoringTool
    return ISF__MemoryMonitoringTool(name, **kwargs)

def getMC12EntryLayerFilter(name="ISF_MC12EntryLayerFilter", **kwargs):
    kwargs.setdefault('AllowOnlyDefinedBarcodes'          , True    )
    kwargs.setdefault('AllowOnlyLegacyPrimaries'          , False   )
    kwargs.setdefault('LegacyParticleGenerationIncrement' , 1000000 )
    kwargs.setdefault('LegacyFirstSecondaryBarcode'       , 200001  )
    from ISF_Tools.ISF_ToolsConf import ISF__GenericBarcodeFilter
    return ISF__GenericBarcodeFilter(name, **kwargs)

def getMC12LLPEntryLayerFilter(name="ISF_MC12LLPEntryLayerFilter", **kwargs):
    return getMC12EntryLayerFilter(name, **kwargs)

def getMC12PlusEntryLayerFilter(name="ISF_MC12PlusEntryLayerFilter", **kwargs):
    return getMC12EntryLayerFilter(name, **kwargs)

def getMC15EntryLayerFilter(name="ISF_MC15EntryLayerFilter", **kwargs):
    kwargs.setdefault('MinEkinCharged'      , 100.*MeV  )
    kwargs.setdefault('MinEkinNeutral'      , -1.       )
    from ISF_Tools.ISF_ToolsConf import ISF__EntryLayerFilter
    return ISF__EntryLayerFilter(name, **kwargs  )

def getMC15aEntryLayerFilter(name="ISF_MC15aEntryLayerFilter", **kwargs):
    return getMC15EntryLayerFilter(name, **kwargs)

def getMC15aPlusEntryLayerFilter(name="ISF_MC15aPlusEntryLayerFilter", **kwargs):
    return getMC15EntryLayerFilter(name, **kwargs)

def getValidationEntryLayerFilter(name="ISF_ValidationEntryLayerFilter", **kwargs):
    return getMC12EntryLayerFilter(name, **kwargs)

def getCosmicEventFilterTool(name="ISF_CosmicEventFilter", **kwargs):
    from G4AtlasApps.SimFlags import simFlags
    volumeNames = ["CaloEntryLayer"]
    useAndFilter = True
    if simFlags.CosmicFilterVolumeName=="Muon":
        volumeNames=["MuonExitLayer"]
    if simFlags.CosmicFilterVolumeName=="Calo":
        volumeNames=["MuonEntryLayer"]
    if simFlags.CosmicFilterVolumeName=="TRT_Barrel":
        volumeNames=["TRTBarrelEntryLayer"]
    if simFlags.CosmicFilterVolumeName=="TRT_EC":
        volumeNames=["TRTECAEntryLayer"]
    if simFlags.CosmicFilterVolumeName=="SCT_Barrel":
        volumeNames=["SCTBarrelEntryLayer"]
    if simFlags.CosmicFilterVolumeName=="Pixel":
        volumeNames=["PixelEntryLayer"]
    print 'G4CosmicFilter: Filter volume is %s' % volumeNames

    if simFlags.CosmicFilterVolumeName=="TRT_EC":
        print 'G4CosmicFilter: adding TRTECBEntryLayer'
        useAndFilter = True
        volumeNames += ["TRTECBEntryLayer"]

    if simFlags.CosmicFilterVolumeName2.statusOn:
        if simFlags.CosmicFilterVolumeName2=="Muon":
            useAndFilter = True
            volumeNames += ["MuonExitLayer"]
        if simFlags.CosmicFilterVolumeName=="Calo":
            useAndFilter = True
            volumeNames += ["MuonEntryLayer"]
        if simFlags.CosmicFilterVolumeName2=="InnerDetector":
            useAndFilter = True
            volumeNames += ["CaloEntryLayer"]
        if simFlags.CosmicFilterVolumeName2=="TRT_Barrel":
            useAndFilter = True
            volumeNames += ["TRTBarrelEntryLayer"]
        if simFlags.CosmicFilterVolumeName2=="TRT_EC":
            useAndFilter = False
            volumeNames += ["TRTECAEntryLayer","TRTECBEntryLayer"]
        if simFlags.CosmicFilterVolumeName2=="SCT_Barrel":
            useAndFilter = True
            volumeNames += ["SCTBarrelEntryLayer"]
        if simFlags.CosmicFilterVolumeName2=="Pixel":
            useAndFilter = True
            volumeNames += ["PixelEntryLayer"]
    ##
    kwargs.setdefault('UseAndFilter' , useAndFilter )
    kwargs.setdefault('VolumeNames'  , volumeNames  )
    if simFlags.CosmicFilterID.statusOn:
         kwargs.setdefault('PDG_ID', simFlags.CosmicFilterID.get_Value() )
    if simFlags.CosmicFilterPTmin.statusOn:
        kwargs.setdefault('ptMin'  , simFlags.CosmicFilterPTmin.get_Value() )
    if simFlags.CosmicFilterPTmax.statusOn:
        kwargs.setdefault('ptMax'  , simFlags.CosmicFilterPTmax.get_Value() )
    from ISF_Tools.ISF_ToolsConf import ISF__CosmicEventFilterTool
    return ISF__CosmicEventFilterTool(name, **kwargs)

def getInToOutSubDetOrderingTool(name="ISF_InToOutSubDetOrderingTool", **kwargs):
    # higher ordered particles will be simulated first
    kwargs.setdefault('OrderID'          , 100000000   )
    kwargs.setdefault('OrderBeamPipe'    , 1000000     )
    kwargs.setdefault('OrderCalo'        , 10000       )
    kwargs.setdefault('OrderMS'          , 100         )
    kwargs.setdefault('OrderCavern'      , 1           )
    from ISF_Tools.ISF_ToolsConf import ISF__GenericParticleOrderingTool
    return ISF__GenericParticleOrderingTool(name, **kwargs)

def getParticleOrderingTool(name="ISF_ParticleOrderingTool", **kwargs):
    kwargs.setdefault('OrderID'          , 1           )
    kwargs.setdefault('OrderBeamPipe'    , 1           )
    kwargs.setdefault('OrderCalo'        , 1           )
    kwargs.setdefault('OrderMS'          , 1           )
    kwargs.setdefault('OrderCavern'      , 1           )
    from ISF_Tools.ISF_ToolsConf import ISF__GenericParticleOrderingTool
    return ISF__GenericParticleOrderingTool(name, **kwargs)

