# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory

#############################################################################################
#TRTSegmentFinder
#############################################################################################

def InDetTRT_RoadMakerCfg(flags, name='InDetTRT_RoadMaker', **kwargs):
    from MagFieldServices.MagFieldServicesConfig import MagneticFieldSvcCfg
    acc = MagneticFieldSvcCfg(flags)
    acc.merge(TRT_DetElementsRoadCondAlgCfg(flags)) # To produce the input TRT_DetElementsRoadData_xk CondHandle
    
    from TrkConfig.TrkExRungeKuttaPropagatorConfig import RungeKuttaPropagatorCfg
    InDetPatternPropagator = acc.popToolsAndMerge(RungeKuttaPropagatorCfg(flags, name="InDetPatternPropagator"))
    kwargs.setdefault("RoadWidth", 20.)
    kwargs.setdefault("PropagatorTool", InDetPatternPropagator)
    acc.setPrivateTools(CompFactory.InDet.TRT_DetElementsRoadMaker_xk(name, **kwargs))
    return acc

def TRT_DetElementsRoadCondAlgCfg(flags, name="InDet__TRT_DetElementsRoadCondAlg_xk", **kwargs):
    acc = ComponentAccumulator()
    acc.addCondAlgo(CompFactory.InDet.TRT_DetElementsRoadCondAlg_xk(name, **kwargs))
    return acc

#############################################################################################
# BackTracking
#############################################################################################

def InDetAmbiScoringToolBaseCfg(flags, name='InDetAmbiScoringTool', **kwargs) :
    acc = ComponentAccumulator()

    if "Extrapolator" not in kwargs:
        from TrkConfig.AtlasExtrapolatorConfig import InDetExtrapolatorCfg
        kwargs.setdefault("Extrapolator", acc.popToolsAndMerge(InDetExtrapolatorCfg(flags)))

    if "SummaryTool" not in kwargs:
        from TrkConfig.TrkTrackSummaryToolConfig import InDetTrackSummaryToolCfg
        kwargs.setdefault("SummaryTool", acc.popToolsAndMerge(InDetTrackSummaryToolCfg(flags)))

    if 'DriftCircleCutTool' not in kwargs:
        from InDetConfig.InDetTrackSelectorToolConfig import InDetTRTDriftCircleCutToolCfg
        kwargs.setdefault("DriftCircleCutTool", acc.popToolsAndMerge(InDetTRTDriftCircleCutToolCfg(flags)))

    have_calo_rois = flags.InDet.Tracking.doBremRecovery and flags.InDet.Tracking.doCaloSeededBrem and flags.Detector.EnableCalo
    if have_calo_rois:
        from InDetConfig.InDetCaloClusterROISelectorConfig import CaloClusterROIPhiRZContainerMakerCfg
        acc.merge(CaloClusterROIPhiRZContainerMakerCfg(flags))
        kwargs.setdefault("EMROIPhiRZContainer", "InDetCaloClusterROIPhiRZ5GeV")

    kwargs.setdefault("useAmbigFcn", True )
    kwargs.setdefault("useTRT_AmbigFcn", False )
    kwargs.setdefault("maxZImp", flags.InDet.Tracking.ActivePass.maxZImpact )
    kwargs.setdefault("maxEta", flags.InDet.Tracking.ActivePass.maxEta )
    kwargs.setdefault("usePixel", flags.InDet.Tracking.ActivePass.usePixel )
    kwargs.setdefault("useSCT", flags.InDet.Tracking.ActivePass.useSCT )
    kwargs.setdefault("doEmCaloSeed", have_calo_rois )
    acc.setPrivateTools(CompFactory.InDet.InDetAmbiScoringTool(name, **kwargs))
    return acc

def InDetCosmicsScoringToolBaseCfg(flags, name='InDetCosmicsScoringTool', **kwargs) :
    acc = ComponentAccumulator()

    from TrkConfig.TrkTrackSummaryToolConfig import InDetTrackSummaryToolCfg
    InDetTrackSummaryTool = acc.popToolsAndMerge(InDetTrackSummaryToolCfg(flags))
    acc.addPublicTool(InDetTrackSummaryTool)

    kwargs.setdefault("nWeightedClustersMin", flags.InDet.Tracking.ActivePass.nWeightedClustersMin )
    kwargs.setdefault("minTRTHits", 0 )
    kwargs.setdefault("SummaryTool", InDetTrackSummaryTool )

    acc.setPrivateTools(CompFactory.InDet.InDetCosmicScoringTool(name, **kwargs ))
    return acc

def InDetCosmicExtenScoringToolCfg(flags, name='InDetCosmicExtenScoringTool',**kwargs):
    kwargs.setdefault("nWeightedClustersMin", 0)
    kwargs.setdefault("minTRTHits", flags.InDet.Tracking.ActivePass.minTRTonTrk )
    return InDetCosmicsScoringToolBaseCfg(flags, name, **kwargs)

def InDetCosmicScoringTool_TRTCfg(flags, name='InDetCosmicScoringTool_TRT',**kwargs):
    acc = ComponentAccumulator()
    from TrkConfig.TrkTrackSummaryToolConfig import InDetTrackSummaryToolNoHoleSearchCfg
    InDetTrackSummaryToolNoHoleSearch = acc.popToolsAndMerge(InDetTrackSummaryToolNoHoleSearchCfg(flags))
    acc.addPublicTool(InDetTrackSummaryToolNoHoleSearch)

    kwargs.setdefault("minTRTHits", flags.InDet.Tracking.ActivePass.minSecondaryTRTonTrk)
    kwargs.setdefault("SummaryTool", InDetTrackSummaryToolNoHoleSearch)

    acc.setPrivateTools(acc.popToolsAndMerge(InDetCosmicExtenScoringToolCfg(flags, name, **kwargs)))
    return acc

def InDetTRT_SeededScoringToolCfg(flags, name='InDetTRT_SeededScoringTool', **kwargs) :
    acc = ComponentAccumulator()

    kwargs.setdefault("useAmbigFcn", False)
    kwargs.setdefault("useTRT_AmbigFcn", True)
    kwargs.setdefault("minTRTonTrk", flags.InDet.Tracking.ActivePass.minSecondaryTRTonTrk)
    kwargs.setdefault("minTRTPrecisionFraction", flags.InDet.Tracking.ActivePass.minSecondaryTRTPrecFrac)
    kwargs.setdefault("minPt", flags.InDet.Tracking.ActivePass.minSecondaryPt)
    kwargs.setdefault("maxRPhiImp", flags.InDet.Tracking.ActivePass.maxSecondaryImpact)
    kwargs.setdefault("minSiClusters", flags.InDet.Tracking.ActivePass.minSecondaryClusters)
    kwargs.setdefault("maxSiHoles", flags.InDet.Tracking.ActivePass.maxSecondaryHoles)
    kwargs.setdefault("maxPixelHoles", flags.InDet.Tracking.ActivePass.maxSecondaryPixelHoles)
    kwargs.setdefault("maxSCTHoles", flags.InDet.Tracking.ActivePass.maxSecondarySCTHoles)
    kwargs.setdefault("maxDoubleHoles", flags.InDet.Tracking.ActivePass.maxSecondaryDoubleHoles)

    acc.setPrivateTools(acc.popToolsAndMerge(InDetAmbiScoringToolBaseCfg(flags, name=name, **kwargs)))
    return acc

#########################################################################################################
#TRTExtension
#########################################################################################################

def InDetAmbiScoringToolCfg(flags, name='InDetAmbiScoringTool', **kwargs) :
    acc = ComponentAccumulator()
    kwargs.setdefault("useAmbigFcn", True )
    kwargs.setdefault("useTRT_AmbigFcn", False )
    kwargs.setdefault("minTRTonTrk", 0 )
    kwargs.setdefault("minTRTPrecisionFraction", 0 )
    kwargs.setdefault("minPt", flags.InDet.Tracking.ActivePass.minPT )
    kwargs.setdefault("maxRPhiImp", flags.InDet.Tracking.ActivePass.maxPrimaryImpact )
    kwargs.setdefault("minSiClusters", flags.InDet.Tracking.ActivePass.minClusters )
    kwargs.setdefault("minPixel", flags.InDet.Tracking.ActivePass.minPixel )
    kwargs.setdefault("maxSiHoles", flags.InDet.Tracking.ActivePass.maxHoles )
    kwargs.setdefault("maxPixelHoles", flags.InDet.Tracking.ActivePass.maxPixelHoles )
    kwargs.setdefault("maxSCTHoles", flags.InDet.Tracking.ActivePass.maxSctHoles )
    kwargs.setdefault("maxDoubleHoles", flags.InDet.Tracking.ActivePass.maxDoubleHoles )
    acc.setPrivateTools(acc.popToolsAndMerge(InDetAmbiScoringToolBaseCfg(flags, name + flags.InDet.Tracking.ActivePass.extension, **kwargs)))
    return acc

def InDetAmbiScoringToolSiCfg(flags, name='InDetAmbiScoringToolSi', **kwargs) :
    kwargs.setdefault('DriftCircleCutTool', None)
    return InDetAmbiScoringToolCfg(flags, name, **kwargs)

def InDetExtenScoringToolCfg(flags, name='InDetExtenScoringTool', **kwargs) :
    acc = ComponentAccumulator()
    kwargs.setdefault("minTRTonTrk", flags.InDet.Tracking.ActivePass.minTRTonTrk)
    kwargs.setdefault("minTRTPrecisionFraction", flags.InDet.Tracking.ActivePass.minTRTPrecFrac)
    acc.setPrivateTools(acc.popToolsAndMerge(InDetAmbiScoringToolCfg(flags, name = name,  **kwargs)))
    return acc

#############################################################################################
#TrackingSiPatternConfig
#############################################################################################

def PRDtoTrackMapToolCfg(name='PRDtoTrackMapTool',**kwargs) :
    acc = ComponentAccumulator()
    acc.setPrivateTools(CompFactory.Trk.PRDtoTrackMapTool( name, **kwargs))
    return acc

def InDetNNScoringToolBaseCfg(flags, name='InDetNNScoringTool', **kwargs) :
    acc = ComponentAccumulator()

    have_calo_rois = flags.InDet.Tracking.doBremRecovery and flags.InDet.Tracking.doCaloSeededBrem and flags.Detector.EnableCalo
    if have_calo_rois :
        from InDetConfig.InDetCaloClusterROISelectorConfig import CaloClusterROIPhiRZContainerMakerCfg
        acc.merge(CaloClusterROIPhiRZContainerMakerCfg(flags))
        kwargs.setdefault("EMROIPhiRZContainer", "InDetCaloClusterROIPhiRZ5GeV")

    if "Extrapolator" not in kwargs:
        from TrkConfig.AtlasExtrapolatorConfig import InDetExtrapolatorCfg
        kwargs.setdefault("Extrapolator", acc.popToolsAndMerge(InDetExtrapolatorCfg(flags)))

    if "SummaryTool" not in kwargs:
        from TrkConfig.TrkTrackSummaryToolConfig import InDetTrackSummaryToolCfg
        kwargs.setdefault("SummaryTool", acc.popToolsAndMerge(InDetTrackSummaryToolCfg(flags)))

    if 'DriftCircleCutTool' not in kwargs:
        from InDetConfig.InDetTrackSelectorToolConfig import InDetTRTDriftCircleCutToolCfg
        kwargs.setdefault("DriftCircleCutTool", acc.popToolsAndMerge(InDetTRTDriftCircleCutToolCfg(flags)))

    kwargs.setdefault("nnCutConfig", "dev/TrackingCP/LRTAmbiNetwork/20200727_225401/nn-config.json" )
    kwargs.setdefault("nnCutThreshold", flags.InDet.Tracking.nnCutLargeD0Threshold )

    kwargs.setdefault("useAmbigFcn", True )
    kwargs.setdefault("useTRT_AmbigFcn", False )
    kwargs.setdefault("maxZImp", flags.InDet.Tracking.ActivePass.maxZImpact )
    kwargs.setdefault("maxEta", flags.InDet.Tracking.ActivePass.maxEta )
    kwargs.setdefault("usePixel", flags.InDet.Tracking.ActivePass.usePixel )
    kwargs.setdefault("useSCT", flags.InDet.Tracking.ActivePass.useSCT )
    kwargs.setdefault("doEmCaloSeed", have_calo_rois )

    acc.setPrivateTools(CompFactory.InDet.InDetNNScoringTool(name, **kwargs ))
    return acc

def InDetNNScoringToolCfg(flags, name='InDetNNScoringTool', **kwargs) :
    kwargs.setdefault("useAmbigFcn", True )
    kwargs.setdefault("useTRT_AmbigFcn", False )
    kwargs.setdefault("minTRTonTrk", 0 )
    kwargs.setdefault("minTRTPrecisionFraction", 0 )
    kwargs.setdefault("minPt", flags.InDet.Tracking.ActivePass.minPT )
    kwargs.setdefault("maxRPhiImp", flags.InDet.Tracking.ActivePass.maxPrimaryImpact )
    kwargs.setdefault("minSiClusters", flags.InDet.Tracking.ActivePass.minClusters )
    kwargs.setdefault("minPixel", flags.InDet.Tracking.ActivePass.minPixel )
    kwargs.setdefault("maxSiHoles", flags.InDet.Tracking.ActivePass.maxHoles )
    kwargs.setdefault("maxPixelHoles", flags.InDet.Tracking.ActivePass.maxPixelHoles )
    kwargs.setdefault("maxSCTHoles", flags.InDet.Tracking.ActivePass.maxSctHoles )
    kwargs.setdefault("maxDoubleHoles", flags.InDet.Tracking.ActivePass.maxDoubleHoles)

    return InDetNNScoringToolBaseCfg(flags, name=name+flags.InDet.Tracking.ActivePass.extension, **kwargs )

def InDetNNScoringToolSiCfg(flags, name='InDetNNScoringToolSi', **kwargs) :
    kwargs.setdefault('DriftCircleCutTool', None)
    return InDetNNScoringToolCfg(flags,name, **kwargs)

def InDetCosmicsScoringToolCfg(flags, name='InDetCosmicsScoringTool', **kwargs) :
    return InDetCosmicsScoringToolBaseCfg(flags,
                                          name=name+flags.InDet.Tracking.ActivePass.extension)

def FullLinearizedTrackFactoryCfg(flags, name='TrackToVertexIPEstimator', **kwargs):
    from TrkConfig.AtlasExtrapolatorConfig import InDetExtrapolatorCfg

    acc = InDetExtrapolatorCfg(flags)
    kwargs.setdefault("Extrapolator", acc.getPrimary()) #private
    # FIXME Configure AtlasFieldCacheCondObj
    acc.setPrivateTools(CompFactory.Trk.FullLinearizedTrackFactory( **kwargs))
    return acc

def TrackToVertexIPEstimatorCfg(flags, name='TrackToVertexIPEstimator', **kwargs):
    from TrkConfig.AtlasExtrapolatorConfig import InDetExtrapolatorCfg

    # All sub-tools are currently public
    acc = FullLinearizedTrackFactoryCfg(flags, 'FullLinearizedTrkFactory')
    linearizedTrackFactory = acc.popPrivateTools()
    acc.addPublicTool(linearizedTrackFactory)
    kwargs.setdefault("LinearizedTrackFactory", linearizedTrackFactory )
    if 'Extrapolator' not in kwargs:
        extrapolator = acc.getPrimaryAndMerge( InDetExtrapolatorCfg(flags) )
        kwargs.setdefault("Extrapolator", extrapolator)
    # Think it's okay not to configure VertexUpdator
    acc.setPrivateTools(CompFactory.Trk.TrackToVertexIPEstimator( **kwargs))
    return acc
