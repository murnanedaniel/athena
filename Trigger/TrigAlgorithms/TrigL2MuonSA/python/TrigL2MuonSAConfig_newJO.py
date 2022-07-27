#
#  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
#
#  This file configs the L2MuonSA reco alg in the newJO way, 
#      but now is located here temporarily until newJO migrations are done in all trigger signatures.
#  This should be moved at somewhere in offline.

from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory

### Output Name ###
muFastInfo = "MuonL2SAInfo"


# Get Rpc data decoder for MuFast data preparator 
def RpcDataPreparatorCfg( flags, roisKey ):

    acc = ComponentAccumulator()

    # Set Rpc data preparator for MuFast data preparator
    TrigL2MuonSA__RpcDataPreparator=CompFactory.getComp("TrigL2MuonSA::RpcDataPreparator")
    from TrigT1MuonRecRoiTool.TrigT1MuonRecRoiToolConfig import getRun3RPCRecRoiTool
    RpcDataPreparator = TrigL2MuonSA__RpcDataPreparator()
    RpcDataPreparator.TrigT1RPCRecRoiTool = getRun3RPCRecRoiTool(name="RPCRecRoiTool",useRun3Config=flags.Trigger.enableL1MuonPhase1)
    from RegionSelector.RegSelToolConfig import regSelTool_RPC_Cfg
    RpcDataPreparator.RegSel_RPC = acc.popToolsAndMerge( regSelTool_RPC_Cfg( flags ) )

    return acc, RpcDataPreparator

# Get Tgc data decoder for MuFast data preparator 
def TgcDataPreparatorCfg( flags, roisKey ):

    acc = ComponentAccumulator()

    # Set Tgc data preparator for MuFast data preparator
    TrigL2MuonSA__TgcDataPreparator=CompFactory.getComp("TrigL2MuonSA::TgcDataPreparator")
    TgcDataPreparator = TrigL2MuonSA__TgcDataPreparator()
 
    return acc, TgcDataPreparator

# Get Mdt data decoder for MuFast data preparator 
def MdtDataPreparatorCfg( flags, roisKey ):

    acc = ComponentAccumulator()

    # Set Mdt data preparator for MuFast data preparator
    TrigL2MuonSA__MdtDataPreparator=CompFactory.getComp("TrigL2MuonSA::MdtDataPreparator")
    MdtDataPreparator = TrigL2MuonSA__MdtDataPreparator()
    from RegionSelector.RegSelToolConfig import regSelTool_MDT_Cfg
    MdtDataPreparator.RegSel_MDT = acc.popToolsAndMerge( regSelTool_MDT_Cfg( flags ) )
 
    return acc, MdtDataPreparator

# Get Csc data decoder for MuFast data preparator 
def CscDataPreparatorCfg( flags, roisKey ):

    acc = ComponentAccumulator()

    # Set Csc data preparator for MuFast data preparator
    TrigL2MuonSA__CscDataPreparator=CompFactory.getComp("TrigL2MuonSA::CscDataPreparator")
    CscDataPreparator = TrigL2MuonSA__CscDataPreparator()


    return acc, CscDataPreparator

def StgcDataPreparatorCfg( flags, roisKey ):

    acc = ComponentAccumulator()

    # Set Stgc data preparator for MuFast data preparator
    TrigL2MuonSA__StgcDataPreparator=CompFactory.getComp("TrigL2MuonSA::StgcDataPreparator")
    StgcDataPreparator = TrigL2MuonSA__StgcDataPreparator()
    from RegionSelector.RegSelToolConfig import regSelTool_STGC_Cfg
    StgcDataPreparator.RegSel_STGC = acc.popToolsAndMerge( regSelTool_STGC_Cfg( flags ) )

    return acc, StgcDataPreparator

def MmDataPreparatorCfg( flags, roisKey ):

    acc = ComponentAccumulator()

    # Set Mm data preparator for MuFast data preparator
    TrigL2MuonSA__MmDataPreparator=CompFactory.getComp("TrigL2MuonSA::MmDataPreparator")
    MmDataPreparator = TrigL2MuonSA__MmDataPreparator()
    from RegionSelector.RegSelToolConfig import regSelTool_MM_Cfg
    MmDataPreparator.RegSel_MM = acc.popToolsAndMerge( regSelTool_MM_Cfg( flags ) )

    return acc, MmDataPreparator

def RpcRoadDefinerCfg( flags, roisKey ):

    acc = ComponentAccumulator()

    # Set RPC road definer for MuFast data preparator
    TrigL2MuonSA__RpcRoadDefiner=CompFactory.getComp("TrigL2MuonSA::RpcRoadDefiner")
    RpcRoadDefiner = TrigL2MuonSA__RpcRoadDefiner( )
    from RegionSelector.RegSelToolConfig import regSelTool_MDT_Cfg
    RpcRoadDefiner.RegionSelectionTool = acc.popToolsAndMerge( regSelTool_MDT_Cfg( flags ) )

    return acc, RpcRoadDefiner

def TgcRoadDefinerCfg( flags, roisKey ):

    acc = ComponentAccumulator()

    # Set TGC road definer for MuFast data preparator
    TrigL2MuonSA__TgcRoadDefiner=CompFactory.getComp("TrigL2MuonSA::TgcRoadDefiner")
    TgcRoadDefiner = TrigL2MuonSA__TgcRoadDefiner( )
    from RegionSelector.RegSelToolConfig import regSelTool_MDT_Cfg
    TgcRoadDefiner.RegionSelectionTool = acc.popToolsAndMerge( regSelTool_MDT_Cfg( flags ) )

    return acc, TgcRoadDefiner

def ClusterRoadDefinerCfg(flags):

    acc = ComponentAccumulator()
    TrigL2MuonSA__ClusterRoadDefiner = CompFactory.getComp("TrigL2MuonSA::ClusterRoadDefiner")
    from RegionSelector.RegSelToolConfig import regSelTool_MDT_Cfg
    ClusterRoadDefiner = TrigL2MuonSA__ClusterRoadDefiner()
    ClusterRoadDefiner.RegionSelectionTool = acc.popToolsAndMerge(regSelTool_MDT_Cfg(flags))
    return acc, ClusterRoadDefiner

# Based on TrigL2MuonSAConfig at TrigL2MuonSA/TrigL2MuonSAConfig.py
def muFastSteeringCfg( flags, roisKey, setup="" ):
    from MuonConfig.MuonCalibrationConfig import MdtCalibrationToolCfg

    acc = ComponentAccumulator()

    # Get RPC decoder
    rpcAcc, RpcDataPreparator = RpcDataPreparatorCfg( flags, roisKey )
    acc.merge( rpcAcc )

    # Get TGC decoder
    tgcAcc, TgcDataPreparator = TgcDataPreparatorCfg( flags, roisKey )
    acc.merge( tgcAcc )

    # Get MDT decoder
    mdtAcc, MdtDataPreparator = MdtDataPreparatorCfg( flags, roisKey )
    acc.merge( mdtAcc )

    # Get CSC decoder
    if flags.Detector.GeometryCSC:
        cscAcc, CscDataPreparator = CscDataPreparatorCfg( flags, roisKey )
        acc.merge( cscAcc )
    else:
        CscDataPreparator = ""

    # Get sTGC decoder
    if flags.Detector.GeometrysTGC:
        stgcAcc, StgcDataPreparator = StgcDataPreparatorCfg( flags, roisKey )
        acc.merge( stgcAcc )
    else:
        StgcDataPreparator = ""

    # Get MM decoder
    if flags.Detector.GeometryMM:
        mmAcc, MmDataPreparator = MmDataPreparatorCfg( flags, roisKey )
        acc.merge( mmAcc )
    else:
        MmDataPreparator = ""

    # Get RPC road definer
    rpcRDAcc, RpcRoadDefiner = RpcRoadDefinerCfg( flags, roisKey )
    acc.merge( rpcRDAcc )

    # Get TGC road definer
    tgcRDAcc, TgcRoadDefiner = TgcRoadDefinerCfg( flags, roisKey )
    acc.merge( tgcRDAcc )

    #Get Cluster Road Definer
    clusRDAcc, ClusterRoadDefiner = ClusterRoadDefinerCfg(flags)
    acc.merge(clusRDAcc)

    # Set MuFast data preparator
    TrigL2MuonSA__MuFastDataPreparator=CompFactory.getComp("TrigL2MuonSA::MuFastDataPreparator")

    from TrigT1MuonRecRoiTool.TrigT1MuonRecRoiToolConfig import getRun3RPCRecRoiTool
    MuFastDataPreparator = TrigL2MuonSA__MuFastDataPreparator( CSCDataPreparator = CscDataPreparator,
                                                               MDTDataPreparator = MdtDataPreparator,
                                                               RPCDataPreparator = RpcDataPreparator,
                                                               TGCDataPreparator = TgcDataPreparator,
                                                               STGCDataPreparator = StgcDataPreparator,
                                                               MMDataPreparator = MmDataPreparator,
                                                               RpcRoadDefiner = RpcRoadDefiner,
                                                               TgcRoadDefiner = TgcRoadDefiner,
                                                               ClusterRoadDefiner = ClusterRoadDefiner,
                                                               TrigT1RPCRecRoiTool = getRun3RPCRecRoiTool(name="RPCRecRoiTool",useRun3Config=flags.Trigger.enableL1MuonPhase1) )

    # Setup the station fitter
    TrigL2MuonSA__MuFastStationFitter,TrigL2MuonSA__PtFromAlphaBeta=CompFactory.getComps("TrigL2MuonSA::MuFastStationFitter","TrigL2MuonSA::PtFromAlphaBeta")
    PtFromAlphaBeta = TrigL2MuonSA__PtFromAlphaBeta()
    PtFromAlphaBeta.useCscPt = True
    PtFromAlphaBeta.AvoidMisalignedCSCs = False

    MuFastStationFitter = TrigL2MuonSA__MuFastStationFitter( PtFromAlphaBeta = PtFromAlphaBeta )
    TrigL2MuonSA__MuFastPatternFinder,TrigL2MuonSA__MuFastTrackFitter,TrigL2MuonSA__MuFastTrackExtrapolator,TrigL2MuonSA__MuCalStreamerTool,TrigL2MuonSA__CscSegmentMaker=CompFactory.getComps("TrigL2MuonSA::MuFastPatternFinder","TrigL2MuonSA::MuFastTrackFitter","TrigL2MuonSA::MuFastTrackExtrapolator","TrigL2MuonSA::MuCalStreamerTool","TrigL2MuonSA::CscSegmentMaker")
    MuFastPatternFinder     = TrigL2MuonSA__MuFastPatternFinder(CalibrationTool=acc.popToolsAndMerge( MdtCalibrationToolCfg(flags)))
    MuFastTrackFitter       = TrigL2MuonSA__MuFastTrackFitter()
    MuFastTrackExtrapolator = TrigL2MuonSA__MuFastTrackExtrapolator()
    MuCalStreamerTool       = TrigL2MuonSA__MuCalStreamerTool()
    CscSegmentMaker         = TrigL2MuonSA__CscSegmentMaker()

    if not flags.Detector.GeometrysTGC and not flags.Detector.GeometryMM:
        MuFastStationFitter.NswStationFitter=""

    # Set Reco alg of muFast step
    #from TrigL2MuonSA.TrigL2MuonSAMonitoring import TrigL2MuonSAMonitoring
    from TrkConfig.AtlasExtrapolatorConfig import AtlasExtrapolatorCfg
    MuFastSteering=CompFactory.MuFastSteering
    muFastAlg = MuFastSteering( name                   = "MuFastSteering_Muon"+setup,
                                DataPreparator         = MuFastDataPreparator,
                                StationFitter          = MuFastStationFitter,
                                PatternFinder          = MuFastPatternFinder,
                                TrackFitter            = MuFastTrackFitter,
                                TrackExtrapolator      = MuFastTrackExtrapolator,
                                FtfRoadDefiner         = 
                                    CompFactory.TrigL2MuonSA.FtfRoadDefiner( IOExtrapolator=acc.popToolsAndMerge(AtlasExtrapolatorCfg(flags))),
                                CalibrationStreamer    = MuCalStreamerTool, 
                                CscSegmentMaker        = CscSegmentMaker,
                                R_WIDTH_TGC_FAILED     = 200,
                                R_WIDTH_RPC_FAILED     = 400,
                                DoCalibrationStream    = False,
                                USE_ROIBASEDACCESS_CSC = True,
                                RpcErrToDebugStream    = True,
                                topoRoad=True,
                                dEtasurrRoI = 0.14,
                                dPhisurrRoI = 0.14,
                                MonTool                = None,
                                UseRun3Config          = flags.Trigger.enableL1MuonPhase1 )

    # Default backextrapolator is for MC Misaligned Detector
    # Based on MuonBackExtrapolatorForMisalignedDet at TrigMuonBackExtrapolator/TrigMuonBackExtrapolatorConfig.py
    TrigMuonBackExtrapolator=CompFactory.TrigMuonBackExtrapolator
    muFastAlg.BackExtrapolator = TrigMuonBackExtrapolator( name        = "MisalignedBackExtrapolator",
                                                           Aligned     = False,
                                                           DataSet     = False )

    # set the flag whether to use NSW or not
    muFastAlg.USE_STGC = True
    muFastAlg.USE_MM = True

    muFastAlg.UseEndcapInnerFromBarrel = True

    if setup == '900GeV':
        muFastAlg.WinPt = 4.0
        muFastAlg.Scale_Road_BarrelInner  = 3
        muFastAlg.Scale_Road_BarrelMiddle = 3
        muFastAlg.Scale_Road_BarrelOuter  = 3
    else:
        muFastAlg.WinPt = 6.0
        muFastAlg.Scale_Road_BarrelInner  = 1
        muFastAlg.Scale_Road_BarrelMiddle = 1
        muFastAlg.Scale_Road_BarrelOuter  = 1

    if setup == 'MuonCalib':
        muFastAlg.DoCalibrationStream = True
        muFastAlg.MuonCalDataScouting = False
        muFastAlg.MuonCalBufferSize   = 1024*1024

    if setup == 'MuonCalibDataScouting':
        muFastAlg.DoCalibrationStream = True
        muFastAlg.MuonCalDataScouting = True
        muFastAlg.MuonCalBufferSize   = 1024*1024

    return acc, muFastAlg

def PtBarrelLUTSvcCfg( flags ):

    acc = ComponentAccumulator()
    ptBarrelLUTSvc = CompFactory.getComp("TrigL2MuonSA::PtBarrelLUTSvc")(name = 'PtBarrelLUTSvc')
    ptBarrelLUTSvc.LUTfile = "pt_barrel.lut"
    ptBarrelLUTSvc.SP_LUTfile = "pt_barrelSP_new.lut"
    
    acc.addService( ptBarrelLUTSvc )

    return acc, ptBarrelLUTSvc

def PtBarrelLUTSvcCfg_MC( flags ):

    acc = ComponentAccumulator()
    ptBarrelLUTSvc_MC = CompFactory.getComp("TrigL2MuonSA::PtBarrelLUTSvc")(name = 'PtBarrelLUTSvc_MC')
    ptBarrelLUTSvc_MC.LUTfile = "pt_barrel.mc10.lut"
    acc.addService( ptBarrelLUTSvc_MC )

    return acc, ptBarrelLUTSvc_MC

def PtEndcapLUTSvcCfg( flags ):

    acc = ComponentAccumulator()
    ptEndcapLUTSvc = CompFactory.getComp("TrigL2MuonSA::PtEndcapLUTSvc")(name = 'PtEndcapLUTSvc')
    ptEndcapLUTSvc.FileName = "pt_endcap.lut"
    ptEndcapLUTSvc.EMeanLUT = "pt_comb_mean.lut"
    ptEndcapLUTSvc.ESigmaLUT = "pt_comb_sigma.lut"
    if flags.Detector.GeometrysTGC or flags.Detector.GeometryMM:
        ptEndcapLUTSvc.UseRun3LUT = True
    else:
        ptEndcapLUTSvc.UseRun3LUT = False
    acc.addService( ptEndcapLUTSvc )

    return acc, ptEndcapLUTSvc

def PtEndcapLUTSvcCfg_MC( flags ):

    acc = ComponentAccumulator()
    ptEndcapLUTSvc_MC = CompFactory.getComp("TrigL2MuonSA::PtEndcapLUTSvc")(name = 'PtEndcapLUTSvc_MC')
    ptEndcapLUTSvc_MC.FileName = "pt_endcap.mc10.lut"
    ptEndcapLUTSvc_MC.EMeanLUT = "pt_comb_mean.lut"
    ptEndcapLUTSvc_MC.ESigmaLUT = "pt_comb_sigma.lut"
    if flags.Detector.GeometrysTGC or flags.Detector.GeometryMM:
        ptEndcapLUTSvc_MC.UseRun3LUT = True
    else:
        ptEndcapLUTSvc_MC.UseRun3LUT = False
    acc.addService( ptEndcapLUTSvc_MC )

    return acc, ptEndcapLUTSvc_MC


def AlignmentBarrelLUTSvcCfg( flags ):

    acc = ComponentAccumulator()
    alignmentBarrelLUTSvc = CompFactory.getComp("TrigL2MuonSA::AlignmentBarrelLUTSvc")(name = 'AlignmentBarrelLUTSvc')
    alignmentBarrelLUTSvc.LUTfile = "dZ_barrel.lut"
    acc.addService( alignmentBarrelLUTSvc )

    return acc, alignmentBarrelLUTSvc

# In the future, above functions should be moved to TrigL2MuonSA package(?)

def l2MuFastAlgCfg( flags, roisKey="" ):

    acc = ComponentAccumulator()

    if not roisKey:
        from HLTSeeding.HLTSeedingConfig import mapThresholdToL1RoICollection
        roisKey = mapThresholdToL1RoICollection("MU")

    # Get Reco alg of muFast step
    muFastAcc, muFastFex = muFastSteeringCfg( flags, roisKey )  
    muFastFex.MuRoIs = roisKey
    muFastFex.Run2RecMuonRoI = "HLT_RecMURoIs"
    muFastFex.RecMuonRoI = "LVL1MuonRoIs"
    muFastFex.MuonL2SAInfo = muFastInfo
    muFastFex.forID = "forID"
    muFastFex.forMS = "forMS"
    acc.merge( muFastAcc )

    # Get services of the Reco alg
    acc.merge( PtBarrelLUTSvcCfg(flags)[0] )   
    acc.merge( PtBarrelLUTSvcCfg_MC(flags)[0] )   
    acc.merge( PtEndcapLUTSvcCfg(flags)[0] )   
    acc.merge( PtEndcapLUTSvcCfg_MC(flags)[0] )   
    acc.merge( AlignmentBarrelLUTSvcCfg(flags)[0] )
    acc.addEventAlgo(muFastFex)

    return acc



def l2MuFastHypoCfg( flags, name="UNSPECIFIED", muFastInfo="UNSPECIFIED" ):

    TrigMufastHypoAlg=CompFactory.TrigMufastHypoAlg
    muFastHypo = TrigMufastHypoAlg( name )
    muFastHypo.MuonL2SAInfoFromMuFastAlg = muFastInfo 

    return muFastHypo
 
