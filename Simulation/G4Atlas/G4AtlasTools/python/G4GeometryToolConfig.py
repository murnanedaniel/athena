# Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
from __future__ import print_function

from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory

from AthenaCommon import Logging

#the physics region tools
from G4AtlasTools.G4PhysicsRegionConfigNew import SX1PhysicsRegionToolCfg, BedrockPhysicsRegionToolCfg, CavernShaftsConcretePhysicsRegionToolCfg, PixelPhysicsRegionToolCfg, SCTPhysicsRegionToolCfg, TRTPhysicsRegionToolCfg, TRT_ArPhysicsRegionToolCfg,ITkPixelPhysicsRegionToolCfg,ITkStripPhysicsRegionToolCfg,BeampipeFwdCutPhysicsRegionToolCfg, FWDBeamLinePhysicsRegionToolCfg, EMBPhysicsRegionToolCfg, EMECPhysicsRegionToolCfg, HECPhysicsRegionToolCfg, FCALPhysicsRegionToolCfg, FCAL2ParaPhysicsRegionToolCfg, EMECParaPhysicsRegionToolCfg, FCALParaPhysicsRegionToolCfg, PreSampLArPhysicsRegionToolCfg, DeadMaterialPhysicsRegionToolCfg
from G4AtlasTools.G4PhysicsRegionConfigNew import DriftWallPhysicsRegionToolCfg, DriftWall1PhysicsRegionToolCfg, DriftWall2PhysicsRegionToolCfg

#the field config tools
from G4AtlasTools.G4FieldConfigNew import ATLASFieldManagerToolCfg, TightMuonsATLASFieldManagerToolCfg, BeamPipeFieldManagerToolCfg, InDetFieldManagerToolCfg, MuonsOnlyInCaloFieldManagerToolCfg, MuonFieldManagerToolCfg, Q1FwdFieldManagerToolCfg, Q2FwdFieldManagerToolCfg, Q3FwdFieldManagerToolCfg, D1FwdFieldManagerToolCfg, D2FwdFieldManagerToolCfg, Q4FwdFieldManagerToolCfg, Q5FwdFieldManagerToolCfg, Q6FwdFieldManagerToolCfg, Q7FwdFieldManagerToolCfg, Q1HKickFwdFieldManagerToolCfg, Q1VKickFwdFieldManagerToolCfg, Q2HKickFwdFieldManagerToolCfg, Q2VKickFwdFieldManagerToolCfg, Q3HKickFwdFieldManagerToolCfg, Q3VKickFwdFieldManagerToolCfg, Q4VKickAFwdFieldManagerToolCfg, Q4HKickFwdFieldManagerToolCfg, Q4VKickBFwdFieldManagerToolCfg, Q5HKickFwdFieldManagerToolCfg,  Q6VKickFwdFieldManagerToolCfg, FwdRegionFieldManagerToolCfg

from G4AtlasTools.G4AtlasToolsConfigNew import SensitiveDetectorMasterToolCfg

GeoDetectorTool=CompFactory.GeoDetectorTool
from BeamPipeGeoModel.BeamPipeGMConfig import BeamPipeGeometryCfg
from AtlasGeoModel.InDetGMConfig import InDetGeometryCfg, InDetServiceMaterialCfg
from AtlasGeoModel.ITkGMConfig import ITkGeometryCfg
from HGTD_GeoModel.HGTD_GeoModelConfig import HGTD_GeometryCfg
from LArGeoAlgsNV.LArGMConfig import LArGMCfg
from TileGeoModel.TileGMConfig import TileGMCfg
from MuonConfig.MuonGeometryConfig import MuonGeoModelCfg
from AtlasGeoModel.ForDetGeoModelConfig import ForDetGeometryCfg

CylindricalEnvelope, PolyconicalEnvelope, MaterialDescriptionTool,SmartlessnessTool,G4AtlasDetectorConstructionTool=CompFactory.getComps("CylindricalEnvelope","PolyconicalEnvelope","MaterialDescriptionTool","SmartlessnessTool","G4AtlasDetectorConstructionTool",)

from AthenaCommon.SystemOfUnits import mm, cm, m

#ToDo - finish migrating this (dnoel)
#to still migrate: getCavernWorld, getCavernInfraGeoDetectorTool
#from ForwardRegionProperties.ForwardRegionPropertiesToolConfig import ForwardRegionPropertiesCfg


#put it here to avoid circular import?
G4GeometryNotifierSvc=CompFactory.G4GeometryNotifierSvc
def G4GeometryNotifierSvcCfg(ConfigFlags, name="G4GeometryNotifierSvc", **kwargs):
    kwargs.setdefault("ActivateLVNotifier", True)
    kwargs.setdefault("ActivatePVNotifier", False)
    return G4GeometryNotifierSvc(name, **kwargs)


def BeamPipeGeoDetectorToolCfg(ConfigFlags, name='BeamPipe', **kwargs):
    #set up geometry
    result=BeamPipeGeometryCfg(ConfigFlags)
    kwargs.setdefault("DetectorName", "BeamPipe")
    #add the GeometryNotifierSvc
    result.addService(G4GeometryNotifierSvcCfg(ConfigFlags))
    kwargs.setdefault("GeometryNotifierSvc", result.getService("G4GeometryNotifierSvc"))
    result.setPrivateTools(GeoDetectorTool(name, **kwargs))
    return result


def PixelGeoDetectorToolCfg(ConfigFlags, name='Pixel', **kwargs):
    #set up geometry
    result=InDetGeometryCfg(ConfigFlags)
    kwargs.setdefault("DetectorName", "Pixel")
    #add the GeometryNotifierSvc
    result.addService(G4GeometryNotifierSvcCfg(ConfigFlags))
    kwargs.setdefault("GeometryNotifierSvc", result.getService("G4GeometryNotifierSvc"))
    result.setPrivateTools(GeoDetectorTool(name, **kwargs))
    return result


def SCTGeoDetectorToolCfg(ConfigFlags, name='SCT', **kwargs):
    #set up geometry
    result=InDetGeometryCfg(ConfigFlags)
    kwargs.setdefault("DetectorName", "SCT")
    #add the GeometryNotifierSvc
    result.addService(G4GeometryNotifierSvcCfg(ConfigFlags))
    kwargs.setdefault("GeometryNotifierSvc", result.getService("G4GeometryNotifierSvc"))
    result.setPrivateTools(GeoDetectorTool(name, **kwargs))
    return result

def ITkPixelGeoDetectorToolCfg(ConfigFlags, name='ITkPixel', **kwargs):
    #set up geometry
    result=ITkGeometryCfg(ConfigFlags)
    kwargs.setdefault("DetectorName", "ITkPixel")
    #add the GeometryNotifierSvc
    result.addService(G4GeometryNotifierSvcCfg(ConfigFlags))
    kwargs.setdefault("GeometryNotifierSvc", result.getService("G4GeometryNotifierSvc"))
    result.setPrivateTools(GeoDetectorTool(name, **kwargs))
    return result


def ITkStripGeoDetectorToolCfg(ConfigFlags, name='ITkStrip', **kwargs):
    #set up geometry
    result=ITkGeometryCfg(ConfigFlags)
    kwargs.setdefault("DetectorName", "ITkStrip")
    #add the GeometryNotifierSvc
    result.addService(G4GeometryNotifierSvcCfg(ConfigFlags))
    kwargs.setdefault("GeometryNotifierSvc", result.getService("G4GeometryNotifierSvc"))
    result.setPrivateTools(GeoDetectorTool(name, **kwargs))
    return result


def HGTDGeoDetectorToolCfg(ConfigFlags, name='HGTD', **kwargs):
    #set up geometry
    result=HGTD_GeometryCfg(ConfigFlags)
    kwargs.setdefault("DetectorName", "HGTD")
    #add the GeometryNotifierSvc
    result.addService(G4GeometryNotifierSvcCfg(ConfigFlags))
    kwargs.setdefault("GeometryNotifierSvc", result.getService("G4GeometryNotifierSvc"))
    result.setPrivateTools(GeoDetectorTool(name, **kwargs))
    return result


def TRTGeoDetectorToolCfg(ConfigFlags, name='TRT', **kwargs):
    #set up geometry
    result=InDetGeometryCfg(ConfigFlags)
    kwargs.setdefault("DetectorName", "TRT")
    #add the GeometryNotifierSvc
    result.addService(G4GeometryNotifierSvcCfg(ConfigFlags))
    kwargs.setdefault("GeometryNotifierSvc", result.getService("G4GeometryNotifierSvc"))
    result.setPrivateTools(GeoDetectorTool(name, **kwargs))
    return result


def IDetServicesMatGeoDetectorToolCfg(ConfigFlags, name='IDetServicesMat', **kwargs):
    #set up geometry
    result=InDetServiceMaterialCfg(ConfigFlags)
    kwargs.setdefault("DetectorName", "InDetServMat")
    #add the GeometryNotifierSvc
    result.addService(G4GeometryNotifierSvcCfg(ConfigFlags))
    kwargs.setdefault("GeometryNotifierSvc", result.getService("G4GeometryNotifierSvc"))
    result.setPrivateTools(GeoDetectorTool(name, **kwargs))
    return result


def LArMgrGeoDetectorToolCfg(ConfigFlags, name='LArMgr', **kwargs):
    #set up geometry
    result=LArGMCfg(ConfigFlags)
    kwargs.setdefault("DetectorName", "LArMgr")
    #add the GeometryNotifierSvc
    result.addService(G4GeometryNotifierSvcCfg(ConfigFlags))
    kwargs.setdefault("GeometryNotifierSvc", result.getService("G4GeometryNotifierSvc"))
    result.setPrivateTools(GeoDetectorTool(name, **kwargs))
    return result


def TileGeoDetectorToolCfg(ConfigFlags, name='Tile', **kwargs):
    #set up geometry
    result=TileGMCfg(ConfigFlags)
    kwargs.setdefault("DetectorName", "Tile")
    #add the GeometryNotifierSvc
    result.addService(G4GeometryNotifierSvcCfg(ConfigFlags))
    kwargs.setdefault("GeometryNotifierSvc", result.getService("G4GeometryNotifierSvc"))
    result.setPrivateTools(GeoDetectorTool(name, **kwargs))
    return result


def LucidGeoDetectorToolCfg(ConfigFlags, name='Lucid', **kwargs):
    #set up geometry
    result=ForDetGeometryCfg(ConfigFlags)
    kwargs.setdefault("DetectorName", "LUCID")
    #add the GeometryNotifierSvc
    result.addService(G4GeometryNotifierSvcCfg(ConfigFlags))
    kwargs.setdefault("GeometryNotifierSvc", result.getService("G4GeometryNotifierSvc"))
    result.setPrivateTools(GeoDetectorTool(name, **kwargs))
    return result


def ALFAGeoDetectorToolCfg(ConfigFlags, name='ALFA', **kwargs):
    #set up geometry
    result=ForDetGeometryCfg(ConfigFlags)
    kwargs.setdefault("DetectorName", "ALFA")
    #add the GeometryNotifierSvc
    result.addService(G4GeometryNotifierSvcCfg(ConfigFlags))
    kwargs.setdefault("GeometryNotifierSvc", result.getService("G4GeometryNotifierSvc"))
    result.setPrivateTools(GeoDetectorTool(name, **kwargs))
    return result


def ZDCGeoDetectorToolCfg(ConfigFlags, name='ZDC', **kwargs):
    #set up geometry
    result=ForDetGeometryCfg(ConfigFlags)
    kwargs.setdefault("DetectorName", "ZDC")
    #add the GeometryNotifierSvc
    result.addService(G4GeometryNotifierSvcCfg(ConfigFlags))
    kwargs.setdefault("GeometryNotifierSvc", result.getService("G4GeometryNotifierSvc"))
    result.setPrivateTools(GeoDetectorTool(name, **kwargs))
    return result


def AFPGeoDetectorToolCfg(ConfigFlags, name='AFP', **kwargs):
    #set up geometry
    result=ForDetGeometryCfg(ConfigFlags)
    kwargs.setdefault("DetectorName", "AFP")
    kwargs.setdefault("GeoDetectorName", "AFP_GeoModel")
    #add the GeometryNotifierSvc
    result.addService(G4GeometryNotifierSvcCfg(ConfigFlags))
    kwargs.setdefault("GeometryNotifierSvc", result.getService("G4GeometryNotifierSvc"))
    result.setPrivateTools(GeoDetectorTool(name, **kwargs))
    return result


def FwdRegionGeoDetectorToolCfg(ConfigFlags, name='FwdRegion', **kwargs):
    #set up geometry
    result=ForDetGeometryCfg(ConfigFlags)
    kwargs.setdefault("DetectorName", "FwdRegion")
    kwargs.setdefault("GeoDetectorName", "ForwardRegionGeoModel")
    #add the GeometryNotifierSvc
    result.addService(G4GeometryNotifierSvcCfg(ConfigFlags))
    kwargs.setdefault("GeometryNotifierSvc", result.getService("G4GeometryNotifierSvc"))
    result.setPrivateTools(GeoDetectorTool(name, **kwargs))
    return result


def MuonGeoDetectorToolCfg(ConfigFlags, name='Muon', **kwargs):
    #set up geometry
    result=MuonGeoModelCfg(ConfigFlags)
    kwargs.setdefault("DetectorName", "Muon")
    #add the GeometryNotifierSvc
    result.addService(G4GeometryNotifierSvcCfg(ConfigFlags))
    kwargs.setdefault("GeometryNotifierSvc", result.getService("G4GeometryNotifierSvc"))
    result.setPrivateTools(GeoDetectorTool(name, **kwargs))
    return result


#todo - set this up
def getCavernInfraGeoDetectorTool(ConfigFlags, name='CavernInfra', **kwargs):
    result = ComponentAccumulator() #needs geometry setting up!
    kwargs.setdefault("DetectorName", "CavernInfra")
    #add the GeometryNotifierSvc
    result.addService(G4GeometryNotifierSvcCfg(ConfigFlags))
    kwargs.setdefault("GeometryNotifierSvc", result.getService("G4GeometryNotifierSvc"))
    result.setPrivateTools(GeoDetectorTool(name, **kwargs))
    return result


def ITKEnvelopeCfg(ConfigFlags, name="ITK", **kwargs):
    result = ComponentAccumulator()

    kwargs.setdefault("DetectorName", "ITK")
    kwargs.setdefault("InnerRadius", 32.15*mm)
    kwargs.setdefault("OuterRadius", 1.148*m)
    if ConfigFlags.Detector.GeometryHGTD:
        # ITk should include the HGTD (3420 mm < |z| < 3545 mm) when turned on
        kwargs.setdefault("dZ", 354.5*cm)
    else:
        kwargs.setdefault("dZ", 347.5*cm)

    SubDetectorList=[]
    if ConfigFlags.Detector.GeometryITkPixel:
        toolITkPixel = result.popToolsAndMerge(ITkPixelGeoDetectorToolCfg(ConfigFlags))
        SubDetectorList += [toolITkPixel]
    if ConfigFlags.Detector.GeometryITkStrip:
        toolITkStrip = result.popToolsAndMerge(ITkStripGeoDetectorToolCfg(ConfigFlags))
        SubDetectorList += [toolITkStrip]

    kwargs.setdefault("SubDetectors", SubDetectorList)
    result.setPrivateTools(CylindricalEnvelope(name, **kwargs))
    return result

def IDETEnvelopeCfg(ConfigFlags, name="IDET", **kwargs):
    result = ComponentAccumulator()
    isRUN2 = ConfigFlags.GeoModel.Run in ["RUN2", "RUN3"]
    kwargs.setdefault("DetectorName", "IDET")
    innerRadius = 37.*mm # RUN1 default
    if isRUN2:
        innerRadius = 28.9*mm #29.15*mm
    kwargs.setdefault("InnerRadius", innerRadius)
    kwargs.setdefault("OuterRadius", 1.148*m)
    kwargs.setdefault("dZ", 347.5*cm)

    SubDetectorList=[]
    if ConfigFlags.Detector.GeometryPixel:
        toolPixel = result.popToolsAndMerge(PixelGeoDetectorToolCfg(ConfigFlags))
        SubDetectorList += [toolPixel]
    if ConfigFlags.Detector.GeometrySCT:
        toolSCT = result.popToolsAndMerge(SCTGeoDetectorToolCfg(ConfigFlags))
        SubDetectorList += [toolSCT]
    if ConfigFlags.Detector.GeometryTRT:
        toolTRT = result.popToolsAndMerge(TRTGeoDetectorToolCfg(ConfigFlags))
        SubDetectorList += [toolTRT]

    toolIDetServices = result.popToolsAndMerge(IDetServicesMatGeoDetectorToolCfg(ConfigFlags))
    SubDetectorList += [toolIDetServices]
    kwargs.setdefault("SubDetectors", SubDetectorList)
    result.setPrivateTools(CylindricalEnvelope(name, **kwargs))
    return result


def CALOEnvelopeCfg(ConfigFlags, name="CALO", **kwargs):
    result = ComponentAccumulator()

    kwargs.setdefault("DetectorName", "CALO")
    kwargs.setdefault("NSurfaces", 18)
    kwargs.setdefault("InnerRadii", [41.,41.,41.,41.,41.,41.,120.,120.,1148.,1148.,120.,120.,41.,41.,41.,41.,41.,41.]) #FIXME Units?
    kwargs.setdefault("OuterRadii", [415.,415.,3795.,3795.,4251.,4251.,4251.,4251.,4251.,4251.,4251.,4251.,4251.,4251.,3795.,3795.,415.,415.]) #FIXME Units?
    if ConfigFlags.Detector.GeometryHGTD:
        # Make room for HGTD (3420 mm < |z| < 3545 mm) when turned on
        kwargs.setdefault("ZSurfaces", [-6781.,-6747.,-6747.,-6530.,-6530.,-4587.,-4587.,-3545.,-3545.,3545.,3545.,4587.,4587.,6530.,6530.,6747.,6747.,6781.]) #FIXME Units?
    else:
        kwargs.setdefault("ZSurfaces", [-6781.,-6747.,-6747.,-6530.,-6530.,-4587.,-4587.,-3475.,-3475.,3475.,3475.,4587.,4587.,6530.,6530.,6747.,6747.,6781.]) #FIXME Units?
    SubDetectorList=[]
    if ConfigFlags.Detector.GeometryLAr:
        toolLArMgr = result.popToolsAndMerge(LArMgrGeoDetectorToolCfg(ConfigFlags))
        SubDetectorList += [ toolLArMgr ]
    if ConfigFlags.Detector.GeometryTile:
        toolTile = result.popToolsAndMerge(TileGeoDetectorToolCfg(ConfigFlags))
        SubDetectorList += [ toolTile ]
    kwargs.setdefault("SubDetectors", SubDetectorList)
    result.setPrivateTools(PolyconicalEnvelope(name, **kwargs))
    return result


def ForwardRegionEnvelopeCfg(ConfigFlags, name='ForwardRegion', **kwargs):
    result = ComponentAccumulator()

    kwargs.setdefault("DetectorName", "ForDetEnvelope")
    SubDetectorList=[]

    if ConfigFlags.Detector.GeometryFwdRegion: # I.e. fully simulate the FwdRegion rather than using BeamTransport to get to Forward Detectors
        toolFwdRegion = result.popToolsAndMerge(FwdRegionGeoDetectorToolCfg(ConfigFlags))
        SubDetectorList += [ toolFwdRegion ]

        #TODO - migrate this over (WIP at the moment) (dnoel)
        #toolFwdRegionProperties = ForwardRegionPropertiesCfg(ConfigFlags)
        #result.addPublicTool(toolFwdRegionProperties) #add this as a service later?
    if ConfigFlags.Detector.GeometryZDC:
        toolZDC = result.popToolsAndMerge(ZDCGeoDetectorToolCfg(ConfigFlags))
        SubDetectorList += [ toolZDC ]
    if ConfigFlags.Detector.GeometryALFA:
        toolALFA = result.popToolsAndMerge(ALFAGeoDetectorToolCfg(ConfigFlags))
        SubDetectorList += [ toolALFA ]
    if ConfigFlags.Detector.GeometryAFP:
        toolAFP = result.popToolsAndMerge(AFPGeoDetectorToolCfg(ConfigFlags))
        SubDetectorList += [ toolAFP ]
    kwargs.setdefault("SubDetectors", SubDetectorList)
    ##FIXME Should this really be a GeoDetectorTool???
    result.setPrivateTools(GeoDetectorTool(name, **kwargs))
    return result


def MUONEnvelopeCfg(ConfigFlags, name="MUONQ02", **kwargs): #FIXME rename to MUON when safe (IS IT SAFE?))
    result = ComponentAccumulator()

    kwargs.setdefault("DetectorName", "MUONQ02") #FIXME rename to MUON when safe
    kwargs.setdefault("NSurfaces", 34)
    kwargs.setdefault("InnerRadii", [1050.,1050.,1050.,1050.,436.7,436.7,279.,279.,70.,70.,420.,420.,3800.,3800.,4255.,4255.,4255.,4255.,4255.,4255.,3800.,3800.,420.,420.,70.,70.,279.,279.,436.7,436.7,1050.,1050.,1050.,1050.]) #FIXME Units?
    kwargs.setdefault("OuterRadii", [1500.,1500.,2750.,2750.,12650.,12650.,13400.,13400.,14200.,14200.,14200.,14200.,14200.,14200.,14200.,14200.,13000.,13000.,14200.,14200.,14200.,14200.,14200.,14200.,14200.,14200.,13400.,13400.,12650.,12650.,2750.,2750.,1500.,1500.]) #FIXME Units?
    kwargs.setdefault("ZSurfaces", [-26046.,-23001.,-23001.,-22030.,-22030.,-18650.,-18650.,-12900.,-12900.,-6783.,-6783.,-6748.,-6748.,-6550.,-6550.,-4000.,-4000.,4000.,4000.,6550.,6550.,6748.,6748.,6783.,6783.,12900.,12900.,18650.,18650.,22030.,22030.,23001.,23001.,26046.]) #FIXME Units?
    SubDetectorList=[]
    if ConfigFlags.Detector.GeometryMuon:
        toolMuon = result.popToolsAndMerge(MuonGeoDetectorToolCfg(ConfigFlags))
        SubDetectorList += [ toolMuon ]

    kwargs.setdefault("SubDetectors", SubDetectorList)
    result.setPrivateTools(PolyconicalEnvelope(name, **kwargs))
    return result


def CosmicShortCutCfg(ConfigFlags, name="CosmicShortCut", **kwargs):
    kwargs.setdefault("DetectorName", "TTR_BARREL")
    kwargs.setdefault("NSurfaces", 14)
    kwargs.setdefault("InnerRadii", [70.,70.,12500.,12500.,12500.,12500.,13000.,13000.,12500.,12500.,12500.,12500.,70.,70.]) #FIXME Units?
    kwargs.setdefault("OuterRadii", [12501.,12501.,12501.,12501.,13001.,13001.,13001.,13001.,13001.,13001.,12501.,12501.,12501.,12501.]) #FIXME Units?
    kwargs.setdefault("ZSurfaces", [-22031.,-22030.,-22030.,-12901.,-12901.,-12900.,-12900., 12900.,12900.,12901.,12901.,22030.,22030.,22031.]) #FIXME Units?
    SubDetectorList=[]
    kwargs.setdefault("SubDetectors", SubDetectorList)
    return PolyconicalEnvelope(name, **kwargs)


def generateSubDetectorList(ConfigFlags):
    result = ComponentAccumulator()
    SubDetectorList=[]

    if ConfigFlags.Beam.Type == 'cosmics' or ConfigFlags.Sim.CavernBG != 'Signal':
        if ConfigFlags.Beam.Type == 'cosmics' and ConfigFlags.Sim.ReadTR:
            SubDetectorList += [ CosmicShortCutCfg(ConfigFlags) ]

    if ConfigFlags.Detector.GeometryMuon:
        accMuon = MUONEnvelopeCfg(ConfigFlags)
        toolMuon = accMuon.popPrivateTools()
        SubDetectorList += [ toolMuon ] #FIXME rename to MUON when safe
    if ConfigFlags.Detector.GeometryID:
        toolIDET = result.popToolsAndMerge(IDETEnvelopeCfg(ConfigFlags))
        SubDetectorList += [ toolIDET ]
    if ConfigFlags.Detector.GeometryITk:
        toolITK = result.popToolsAndMerge(ITKEnvelopeCfg(ConfigFlags))
        SubDetectorList += [ toolITK ]
    if ConfigFlags.Detector.GeometryHGTD:
        toolHGTD = result.popToolsAndMerge(HGTDGeoDetectorToolCfg(ConfigFlags))
        SubDetectorList += [ toolHGTD ]
    if ConfigFlags.Detector.GeometryCalo:
        toolCALO = result.popToolsAndMerge(CALOEnvelopeCfg(ConfigFlags))
        SubDetectorList += [ toolCALO ]
    if ConfigFlags.Detector.GeometryMuon:
        result.merge(accMuon) #add the acc later to match the old style config
    if ConfigFlags.Detector.GeometryBpipe:
        toolBpipe = result.popToolsAndMerge(BeamPipeGeoDetectorToolCfg(ConfigFlags))
        SubDetectorList += [ toolBpipe ]
    if ConfigFlags.Detector.GeometryLucid:
        toolLucid = result.popToolsAndMerge(LucidGeoDetectorToolCfg(ConfigFlags))
        SubDetectorList += [ toolLucid ]
    if ConfigFlags.Detector.GeometryForward:
        toolForward = result.popToolsAndMerge(ForwardRegionEnvelopeCfg(ConfigFlags))
        SubDetectorList += [ toolForward ]

    #if DetFlags.Muon_on(): #HACK
    #    SubDetectorList += ['MUONQ02'] #FIXME rename to MUON when safe #HACK
    #SubDetectorList += generateFwdSubDetectorList() #FIXME Fwd Detectors not supported yet.
    result.setPrivateTools(SubDetectorList)
    return result


def ATLASEnvelopeCfg(ConfigFlags, name="Atlas", **kwargs):
    result = ComponentAccumulator()

    kwargs.setdefault("DetectorName", "Atlas")
    kwargs.setdefault("NSurfaces", 18)
    ## InnerRadii
    innerRadii = [0.0] * 18
    kwargs.setdefault("InnerRadii", innerRadii)

    ## Shrink the global ATLAS envelope to the activated detectors,
    ## except when running on special setups.

    ## OuterRadii
    AtlasForwardOuterR = 2751.
    AtlasOuterR1 = 14201.
    AtlasOuterR2 = 14201.
    if ConfigFlags.Beam.Type != 'cosmics' and not ConfigFlags.Detector.GeometryMuon and not \
       (ConfigFlags.Sim.CavernBG != 'Signal'):
        AtlasOuterR1 = 4251.
        AtlasOuterR2 = 4251.
        if not ConfigFlags.Detector.GeometryCalo:
            AtlasOuterR1 = 1150.
            AtlasOuterR2 = 1150.

    outerRadii = [0.0] * 18
    for i in (0, 1, 16, 17):
        outerRadii[i] = 1501.
    for i in (2, 3, 14, 15):
        outerRadii[i] = AtlasForwardOuterR
    for i in (4, 5, 12, 13):
        outerRadii[i] = AtlasOuterR2
    for i in range(6, 12):
        outerRadii[i] = AtlasOuterR1

    ## World R range
    if ConfigFlags.Sim.WorldRRange:
        routValue = ConfigFlags.Sim.WorldRRange
        if ConfigFlags.Sim.WorldRRange > max(AtlasOuterR1, AtlasOuterR2):
            for i in range(4, 14):
                outerRadii[i] = routValue
        else:
            raise RuntimeError('getATLASEnvelope: ERROR ConfigFlags.Sim.WorldRRange must be > %f. Current value %f' % (max(AtlasOuterR1, AtlasOuterR2), routValue) )
    kwargs.setdefault("OuterRadii", outerRadii)

    ## ZSurfaces
    zSurfaces = [-26046., -23001., -23001., -22031., -22031., -12899., -12899., -6741., -6741.,  6741.,  6741.,  12899., 12899., 22031., 22031., 23001., 23001., 26046.] # FIXME units mm??

    if ConfigFlags.Detector.GeometryForward:
        zSurfaces[0]  = -400000.
        zSurfaces[17] =  400000.

    #leave a check in for WorldRrange and WorldZrange?
    if ConfigFlags.Sim.WorldZRange:
        print (ConfigFlags.Sim.WorldZRange)
        if ConfigFlags.Sim.WorldZRange < 26046.:
              raise RuntimeError('getATLASEnvelope: ERROR ConfigFlags.Sim.WorldZRange must be > 26046. Current value: %f' % ConfigFlags.Sim.WorldZRange)
        zSurfaces[17] =  ConfigFlags.Sim.WorldZRange + 100.
        zSurfaces[16] =  ConfigFlags.Sim.WorldZRange + 50.
        zSurfaces[15] =  ConfigFlags.Sim.WorldZRange + 50.
        zSurfaces[14] =  ConfigFlags.Sim.WorldZRange
        zSurfaces[13] =  ConfigFlags.Sim.WorldZRange
        zSurfaces[0] =  -ConfigFlags.Sim.WorldZRange - 100.
        zSurfaces[1] =  -ConfigFlags.Sim.WorldZRange - 50.
        zSurfaces[2] =  -ConfigFlags.Sim.WorldZRange - 50.
        zSurfaces[3] =  -ConfigFlags.Sim.WorldZRange
        zSurfaces[4] =  -ConfigFlags.Sim.WorldZRange

    kwargs.setdefault("ZSurfaces", zSurfaces)
    SubDetectorList = result.popToolsAndMerge(generateSubDetectorList(ConfigFlags))
    kwargs.setdefault("SubDetectors", SubDetectorList)
    result.setPrivateTools(PolyconicalEnvelope(name, **kwargs))
    return result


def MaterialDescriptionToolCfg(ConfigFlags, name="MaterialDescriptionTool", **kwargs):
    ## kwargs.setdefault("SomeProperty", aValue)
    result = ComponentAccumulator()
    result.setPrivateTools(MaterialDescriptionTool(name, **kwargs))
    return result


def SmartlessnessToolCfg(ConfigFlags, name="SmartlessnessTool", **kwargs):
    ## kwargs.setdefault("SomeProperty", aValue)
    result = ComponentAccumulator()
    result.setPrivateTools(SmartlessnessTool(name, **kwargs))
    return result


def getATLAS_RegionCreatorList(ConfigFlags):
    regionCreatorList = []

    isRUN2 = (ConfigFlags.GeoModel.Run in ["RUN2", "RUN3"]) or (ConfigFlags.GeoModel.Run=="UNDEFINED" and ConfigFlags.GeoModel.IBLLayout not in ["noIBL", "UNDEFINED"])

    if ConfigFlags.Beam.Type == 'cosmics' or ConfigFlags.Sim.CavernBG != 'Signal':
        regionCreatorList += [SX1PhysicsRegionToolCfg(ConfigFlags), BedrockPhysicsRegionToolCfg(ConfigFlags), CavernShaftsConcretePhysicsRegionToolCfg(ConfigFlags)]
        #regionCreatorList += ['CavernShaftsAirPhysicsRegionTool'] # Not used currently
    if ConfigFlags.Detector.GeometryID:
        if ConfigFlags.Detector.GeometryPixel:
            regionCreatorList += [PixelPhysicsRegionToolCfg(ConfigFlags)]
        if ConfigFlags.Detector.GeometrySCT:
            regionCreatorList += [SCTPhysicsRegionToolCfg(ConfigFlags)]
        if ConfigFlags.Detector.GeometryTRT:
            regionCreatorList += [TRTPhysicsRegionToolCfg(ConfigFlags)]
            if isRUN2:
                regionCreatorList += [TRT_ArPhysicsRegionToolCfg(ConfigFlags)] #'TRT_KrPhysicsRegionTool'
        # FIXME dislike the ordering here, but try to maintain the same ordering as in the old configuration.
        if ConfigFlags.Detector.GeometryBpipe:
            if ConfigFlags.Sim.BeamPipeSimMode != "Normal":
                regionCreatorList += [BeampipeFwdCutPhysicsRegionToolCfg(ConfigFlags)]
            #if simFlags.ForwardDetectors.statusOn and simFlags.ForwardDetectors() == 2:
            if False:
                regionCreatorList += [FWDBeamLinePhysicsRegionToolCfg(ConfigFlags)]
    if ConfigFlags.Detector.GeometryITk:
        if ConfigFlags.Detector.GeometryITkPixel:
            regionCreatorList += [ITkPixelPhysicsRegionToolCfg(ConfigFlags)] #TODO: add dedicated config
        if ConfigFlags.Detector.GeometryITkStrip:
            regionCreatorList += [ITkStripPhysicsRegionToolCfg(ConfigFlags)] #TODO: And here...
        # FIXME dislike the ordering here, but try to maintain the same ordering as in the old configuration.
        if ConfigFlags.Detector.GeometryBpipe:
            if ConfigFlags.Sim.BeamPipeSimMode != "Normal":
                regionCreatorList += [BeampipeFwdCutPhysicsRegionToolCfg(ConfigFlags)]
            #if simFlags.ForwardDetectors.statusOn and simFlags.ForwardDetectors() == 2:
            if False:
                regionCreatorList += [FWDBeamLinePhysicsRegionToolCfg(ConfigFlags)]
    if ConfigFlags.Detector.GeometryCalo:
        if ConfigFlags.Detector.GeometryLAr:
            # Shower parameterization overrides the calibration hit flag
            if ConfigFlags.Sim.LArParameterization > 0 \
               and ConfigFlags.Sim.CalibrationRun in ['LAr','LAr+Tile','DeadLAr']:
                Logging.log.info('You requested both calibration hits and frozen showers / parameterization in the LAr.')
                Logging.log.info('  Such a configuration is not allowed, and would give junk calibration hits where the showers are modified.')
                Logging.log.info('  Please try again with a different value of either ConfigFlags.Sim.LArParameterization (' + str(ConfigFlags.Sim.LArParameterization) + ') or ConfigFlags.Sim.CalibrationRun ('+str(ConfigFlags.Sim.CalibrationRun)+')')
                raise RuntimeError('Configuration not allowed')
            regionCreatorList += [EMBPhysicsRegionToolCfg(ConfigFlags),
                                  EMECPhysicsRegionToolCfg(ConfigFlags),
                                  HECPhysicsRegionToolCfg(ConfigFlags),
                                  FCALPhysicsRegionToolCfg(ConfigFlags)]
            if ConfigFlags.Sim.LArParameterization > 0:
                # FIXME 'EMBPhysicsRegionTool' used for parametrization also - do we need a second instance??
                regionCreatorList += [EMECParaPhysicsRegionToolCfg(ConfigFlags),
                                      FCALParaPhysicsRegionToolCfg(ConfigFlags),
                                      FCAL2ParaPhysicsRegionToolCfg(ConfigFlags)]
                if ConfigFlags.Sim.LArParameterization > 1:
                    pass
                    #todo - add the line below
                    regionCreatorList += [PreSampLArPhysicsRegionToolCfg(ConfigFlags), DeadMaterialPhysicsRegionToolCfg(ConfigFlags)]
    ## FIXME _initPR never called for FwdRegion??
    #if simFlags.ForwardDetectors.statusOn:
    #    if DetFlags.geometry.FwdRegion_on():
    #        regionCreatorList += ['FwdRegionPhysicsRegionTool']
    if ConfigFlags.Detector.GeometryMuon:
        #todo - add the line below
        regionCreatorList += [DriftWallPhysicsRegionToolCfg(ConfigFlags), DriftWall1PhysicsRegionToolCfg(ConfigFlags), DriftWall2PhysicsRegionToolCfg(ConfigFlags)]
        #if ConfigFlags.Sim.CavernBG != 'Read' and not (simFlags.RecordFlux.statusOn and simFlags.RecordFlux()):
            #pass
            #todo - add the line below
            #regionCreatorList += [MuonSystemFastPhysicsRegionTool(ConfigFlags)]
    return regionCreatorList


def getTB_RegionCreatorList(ConfigFlags):
    regionCreatorList = []
    #from G4AtlasApps.SimFlags import simFlags

    #TODO - migrate below>>
    #if (ConfigFlags.GeoModel.AtlasVersion=="tb_LArH6_2003"):
    #    if (ConfigFlags.Detector.GeometryLAr):
    #        regionCreatorList += [FCALPhysicsRegionTool(ConfigFlags)]
    #elif (ConfigFlags.GeoModel.AtlasVersion=="tb_LArH6_2002"):
    #    if (ConfigFlags.Detector.GeometryLAr):
    #        regionCreatorList += [HECPhysicsRegionTool(ConfigFlags)]
    #elif (ConfigFlags.GeoModel.AtlasVersion=="tb_LArH6EC_2002"):
    #    if (ConfigFlags.Detector.GeometryLAr):
    #        regionCreatorList += [EMECPhysicsRegionTool(ConfigFlags)]
    #elif (ConfigFlags.GeoModel.AtlasVersion=="tb_LArH6_2004"):
    #    if (simFlags.LArTB_H6Hec.get_Value()):
    #        regionCreatorList += [HECPhysicsRegionTool(ConfigFlags)]
    #    if (simFlags.LArTB_H6Emec.get_Value()):
    #        regionCreatorList += [EMECPhysicsRegionTool(ConfigFlags)]
    #    if (simFlags.LArTB_H6Fcal.get_Value()):
    #        regionCreatorList += [FCALPhysicsRegionTool(ConfigFlags)]
    #<<migrate above
    return regionCreatorList


#########################################################################
def ATLAS_FieldMgrListCfg(ConfigFlags):
    result = ComponentAccumulator()
    fieldMgrList = []

    if ConfigFlags.Sim.TightMuonStepping:
        acc   = TightMuonsATLASFieldManagerToolCfg(ConfigFlags)
        tool  = result.popToolsAndMerge(acc)
        fieldMgrList += [tool]
    else:
        acc   = ATLASFieldManagerToolCfg(ConfigFlags)
        tool  = result.popToolsAndMerge(acc)
        fieldMgrList += [tool]

    if ConfigFlags.Detector.GeometryBpipe:
        acc = BeamPipeFieldManagerToolCfg(ConfigFlags)
        tool  = result.popToolsAndMerge(acc)
        fieldMgrList += [tool]
    if ConfigFlags.Detector.GeometryID:
        acc = InDetFieldManagerToolCfg(ConfigFlags)
        tool  = result.popToolsAndMerge(acc)
        fieldMgrList += [tool]
    #if ConfigFlags.Detector.GeometryCalo and simFlags.MuonFieldOnlyInCalo.statusOn and simFlags.MuonFieldOnlyInCalo():
    if False:
        acc = MuonsOnlyInCaloFieldManagerToolCfg(ConfigFlags)
        tool  = result.popToolsAndMerge(acc)
        fieldMgrList += [tool]
    if ConfigFlags.Detector.GeometryMuon:
        acc = MuonFieldManagerToolCfg(ConfigFlags)
        tool  = result.popToolsAndMerge(acc)
        fieldMgrList += [tool]

    #sort these forward ones later
    if ConfigFlags.Detector.GeometryForward: #needed?
        if ConfigFlags.Detector.GeometryFwdRegion: #or forward?
          accQ1FwdRegionFieldManager = Q1FwdFieldManagerToolCfg(ConfigFlags)
          accQ2FwdRegionFieldManager = Q2FwdFieldManagerToolCfg(ConfigFlags)
          accQ3FwdRegionFieldManager = Q3FwdFieldManagerToolCfg(ConfigFlags)
          accD1FwdRegionFieldManager = D1FwdFieldManagerToolCfg(ConfigFlags)
          accD2FwdRegionFieldManager = D2FwdFieldManagerToolCfg(ConfigFlags)
          accQ4FwdRegionFieldManager = Q4FwdFieldManagerToolCfg(ConfigFlags)
          accQ5FwdRegionFieldManager = Q5FwdFieldManagerToolCfg(ConfigFlags)
          accQ6FwdRegionFieldManager = Q6FwdFieldManagerToolCfg(ConfigFlags)
          accQ7FwdRegionFieldManager = Q7FwdFieldManagerToolCfg(ConfigFlags)
          accQ1HKickFwdRegionFieldManager = Q1HKickFwdFieldManagerToolCfg(ConfigFlags)
          accQ1VKickFwdRegionFieldManager = Q1VKickFwdFieldManagerToolCfg(ConfigFlags)
          accQ2HKickFwdRegionFieldManager = Q2HKickFwdFieldManagerToolCfg(ConfigFlags)
          accQ2VKickFwdRegionFieldManager = Q2VKickFwdFieldManagerToolCfg(ConfigFlags)
          accQ3HKickFwdRegionFieldManager = Q3HKickFwdFieldManagerToolCfg(ConfigFlags)
          accQ3VKickFwdRegionFieldManager = Q3VKickFwdFieldManagerToolCfg(ConfigFlags)
          accQ4VKickAFwdRegionFieldManager = Q4VKickAFwdFieldManagerToolCfg(ConfigFlags)
          accQ4HKickFwdRegionFieldManager = Q4HKickFwdFieldManagerToolCfg(ConfigFlags)
          accQ4VKickBFwdRegionFieldManager = Q4VKickBFwdFieldManagerToolCfg(ConfigFlags)
          accQ5HKickFwdRegionFieldManager = Q5HKickFwdFieldManagerToolCfg(ConfigFlags)
          accQ6VKickFwdRegionFieldManager = Q6VKickFwdFieldManagerToolCfg(ConfigFlags)
          accFwdRegionFieldManager = FwdRegionFieldManagerToolCfg(ConfigFlags)

          toolQ1FwdRegionFieldManager = result.popToolsAndMerge(accQ1FwdRegionFieldManager)
          toolQ2FwdFieldManager = result.popToolsAndMerge(accQ2FwdRegionFieldManager)
          toolQ3FwdFieldManager = result.popToolsAndMerge(accQ3FwdRegionFieldManager)
          toolD1FwdFieldManager = result.popToolsAndMerge(accD1FwdRegionFieldManager)
          toolD2FwdFieldManager = result.popToolsAndMerge(accD2FwdRegionFieldManager)
          toolQ4FwdFieldManager = result.popToolsAndMerge(accQ4FwdRegionFieldManager)
          toolQ5FwdFieldManager = result.popToolsAndMerge(accQ5FwdRegionFieldManager)
          toolQ6FwdFieldManager = result.popToolsAndMerge(accQ6FwdRegionFieldManager)
          toolQ7FwdFieldManager = result.popToolsAndMerge(accQ7FwdRegionFieldManager)
          toolQ1HKickFwdFieldManager = result.popToolsAndMerge(accQ1HKickFwdRegionFieldManager)
          toolQ1VKickFwdFieldManager = result.popToolsAndMerge(accQ1VKickFwdRegionFieldManager)
          toolQ2HKickFwdFieldManager = result.popToolsAndMerge(accQ2HKickFwdRegionFieldManager)
          toolQ2VKickFwdFieldManager = result.popToolsAndMerge(accQ2VKickFwdRegionFieldManager)
          toolQ3HKickFwdFieldManager = result.popToolsAndMerge(accQ3HKickFwdRegionFieldManager)
          toolQ3VKickFwdFieldManager = result.popToolsAndMerge(accQ3VKickFwdRegionFieldManager)
          toolQ4VKickAFwdFieldManager = result.popToolsAndMerge(accQ4VKickAFwdRegionFieldManager)
          toolQ4HKickFwdFieldManager = result.popToolsAndMerge(accQ4HKickFwdRegionFieldManager)
          toolQ4VKickBFwdFieldManager = result.popToolsAndMerge(accQ4VKickBFwdRegionFieldManager)
          toolQ5HKickFwdFieldManager = result.popToolsAndMerge(accQ5HKickFwdRegionFieldManager)
          toolQ6VKickFwdFieldManager = result.popToolsAndMerge(accQ6VKickFwdRegionFieldManager)
          toolFwdRegionFieldManager = result.popToolsAndMerge(accFwdRegionFieldManager)

          fieldMgrList+=[toolQ1FwdRegionFieldManager,
                         toolQ2FwdFieldManager,
                         toolQ3FwdFieldManager,
                         toolD1FwdFieldManager,
                         toolD2FwdFieldManager,
                         toolQ4FwdFieldManager,
                         toolQ5FwdFieldManager,
                         toolQ6FwdFieldManager,
                         toolQ7FwdFieldManager,
                         toolQ1HKickFwdFieldManager,
                         toolQ1VKickFwdFieldManager,
                         toolQ2HKickFwdFieldManager,
                         toolQ2VKickFwdFieldManager,
                         toolQ3HKickFwdFieldManager,
                         toolQ3VKickFwdFieldManager,
                         toolQ4VKickAFwdFieldManager,
                         toolQ4HKickFwdFieldManager,
                         toolQ4VKickBFwdFieldManager,
                         toolQ5HKickFwdFieldManager,
                         toolQ6VKickFwdFieldManager,
                         toolFwdRegionFieldManager]

    result.setPrivateTools(fieldMgrList)
    return result


def getTB_FieldMgrList(ConfigFlags):
    fieldMgrList = []
    return fieldMgrList


def getGeometryConfigurationTools(ConfigFlags):
    geoConfigToolList = []
    # The methods for these tools should be defined in the
    # package containing each tool, so G4AtlasTools in this case
    result =ComponentAccumulator()
    geoConfigToolList += [result.popToolsAndMerge(MaterialDescriptionToolCfg(ConfigFlags))]
    geoConfigToolList += [result.popToolsAndMerge(SmartlessnessToolCfg(ConfigFlags))]
    return result, geoConfigToolList


def G4AtlasDetectorConstructionToolCfg(ConfigFlags, name="G4AtlasDetectorConstructionTool", **kwargs):
    result = ComponentAccumulator()

    ## For now just have the same geometry configurations tools loaded for ATLAS and TestBeam
    geoConfAcc, listOfGeoConfTools = getGeometryConfigurationTools(ConfigFlags)
    result.merge(geoConfAcc)
    kwargs.setdefault("GeometryConfigurationTools", listOfGeoConfTools)

    if "SenDetMasterTool" not in kwargs:
        tool = result.popToolsAndMerge(SensitiveDetectorMasterToolCfg(ConfigFlags))
        result.addPublicTool(tool)
        kwargs.setdefault("SenDetMasterTool", result.getPublicTool(tool.name))

    #if hasattr(simFlags,"Eta"): #FIXME ugly hack
    if False:
        kwargs.setdefault("World", 'TileTB_World') # NEED TO ADD THIS
        kwargs.setdefault("RegionCreators", getTB_RegionCreatorList(ConfigFlags))
        kwargs.setdefault("FieldManagers", getTB_FieldMgrList(ConfigFlags))
    #elif hasattr(simFlags,"LArTB_H1TableYPos"): #FIXME ugly hack
    elif False:
        kwargs.setdefault("World", 'LArTB_World')
        kwargs.setdefault("RegionCreators", getTB_RegionCreatorList(ConfigFlags))
        kwargs.setdefault("FieldManagers", getTB_FieldMgrList(ConfigFlags))
    else:
        #if ConfigFlags.Beam.Type == 'cosmics' or ConfigFlags.Sim.CavernBG != 'Signal':
        if False:
            kwargs.setdefault("World", 'Cavern')
        else:
            toolGeo = result.popToolsAndMerge(ATLASEnvelopeCfg(ConfigFlags))
            kwargs.setdefault("World", toolGeo)
        kwargs.setdefault("RegionCreators", getATLAS_RegionCreatorList(ConfigFlags))
        #if hasattr(simFlags, 'MagneticField') and simFlags.MagneticField.statusOn:
        if True:
            acc = ATLAS_FieldMgrListCfg(ConfigFlags)
            fieldMgrList = result.popToolsAndMerge(acc)
            kwargs.setdefault("FieldManagers", fieldMgrList)
    result.setPrivateTools(G4AtlasDetectorConstructionTool(name, **kwargs))
    return result
