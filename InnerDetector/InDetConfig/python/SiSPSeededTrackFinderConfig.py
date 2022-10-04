# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
# Configuration of SiSPSeededTrackFinder package
from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory

def SiSPSeededTrackFinderCfg(flags, name="InDetSiSpTrackFinder", **kwargs) :
    acc = ComponentAccumulator()

    if "TrackTool" not in kwargs:
        from InDetConfig.SiTrackMakerConfig import SiTrackMaker_xkCfg
        kwargs.setdefault("TrackTool", acc.popToolsAndMerge(
            SiTrackMaker_xkCfg(flags)))

    if "PropagatorTool" not in kwargs:
        from TrkConfig.TrkExRungeKuttaPropagatorConfig import InDetPropagatorCfg
        InDetPropagator = acc.popToolsAndMerge(InDetPropagatorCfg(flags))
        acc.addPublicTool(InDetPropagator)
        kwargs.setdefault("PropagatorTool", InDetPropagator)

    if "TrackSummaryTool" not in kwargs:
        from TrkConfig.TrkTrackSummaryToolConfig import InDetTrackSummaryToolNoHoleSearchCfg
        kwargs.setdefault("TrackSummaryTool", acc.popToolsAndMerge(
            InDetTrackSummaryToolNoHoleSearchCfg(flags)))

    if "SeedsTool" not in kwargs:
        from InDetConfig.SiSpacePointsSeedToolConfig import SiSpacePointsSeedMakerCfg
        kwargs.setdefault("SeedsTool", acc.popToolsAndMerge(
            SiSpacePointsSeedMakerCfg(flags)))

    kwargs.setdefault("useMBTSTimeDiff", flags.Reco.EnableHI) # Heavy-ion config

    if flags.InDet.Tracking.ActivePass.usePrdAssociationTool:
        # not all classes have that property !!!
        kwargs.setdefault("PRDtoTrackMap", 'InDetPRDtoTrackMap'+ flags.InDet.Tracking.ActivePass.extension)

    if flags.InDet.Tracking.ActivePass.extension == "Forward":
        kwargs.setdefault("useZvertexTool", flags.Reco.EnableHI) # For heavy-ion
        kwargs.setdefault("useZBoundFinding", False)
    else:
        kwargs.setdefault("useZvertexTool", flags.Reco.EnableHI) # For heavy-ion
        kwargs.setdefault("useZBoundFinding", flags.InDet.Tracking.ActivePass.doZBoundary)
    
    #
    # --- Z-coordinates primary vertices finder (only for collisions)
    #
    if kwargs["useZvertexTool"] and "ZvertexTool" not in kwargs:
        from InDetConfig.SiZvertexToolConfig import SiZvertexMaker_xkCfg
        kwargs.setdefault("ZvertexTool", acc.popToolsAndMerge(
            SiZvertexMaker_xkCfg(flags)))

    if flags.Reco.EnableHI:
        kwargs.setdefault("FreeClustersCut",2) #Heavy Ion optimization from Igor

    acc.addEventAlgo(CompFactory.InDet.SiSPSeededTrackFinder(name+flags.InDet.Tracking.ActivePass.extension, **kwargs))
    return acc

def ITkSiSPSeededTrackFinderCfg(flags, name="ITkSiSpTrackFinder", **kwargs) :
    acc = ComponentAccumulator()

    if "TrackTool" not in kwargs:
        from InDetConfig.SiTrackMakerConfig import ITkSiTrackMaker_xkCfg
        kwargs.setdefault("TrackTool", acc.popToolsAndMerge(
            ITkSiTrackMaker_xkCfg(flags)))

    if "PropagatorTool" not in kwargs:
        from TrkConfig.TrkExRungeKuttaPropagatorConfig import ITkPropagatorCfg
        ITkPropagator = acc.popToolsAndMerge(ITkPropagatorCfg(flags))
        acc.addPublicTool(ITkPropagator)
        kwargs.setdefault("PropagatorTool", ITkPropagator)

    if "TrackSummaryTool" not in kwargs:
        from TrkConfig.TrkTrackSummaryToolConfig import ITkTrackSummaryToolNoHoleSearchCfg
        kwargs.setdefault("TrackSummaryTool", acc.popToolsAndMerge(
            ITkTrackSummaryToolNoHoleSearchCfg(flags)))

    if "SeedsTool" not in kwargs:
        ITkSiSpacePointsSeedMaker = None
        if flags.ITk.Tracking.ActivePass.extension != "ConversionFinding" and flags.Acts.TrackFinding.useSiSpacePointSeedMaker:
            from ActsTrkSeedingTool.ActsTrkSeedingToolConfig import ActsTrkSiSpacePointsSeedMakerCfg
            ITkSiSpacePointsSeedMaker = acc.popToolsAndMerge(ActsTrkSiSpacePointsSeedMakerCfg(flags))
        else:
            from InDetConfig.SiSpacePointsSeedToolConfig import ITkSiSpacePointsSeedMakerCfg
            ITkSiSpacePointsSeedMaker = acc.popToolsAndMerge(ITkSiSpacePointsSeedMakerCfg(flags))

        kwargs.setdefault("SeedsTool", ITkSiSpacePointsSeedMaker)

    if flags.ITk.Tracking.ActivePass.usePrdAssociationTool:
        # not all classes have that property !!!
        kwargs.setdefault("PRDtoTrackMap", 'ITkPRDtoTrackMap'+ flags.ITk.Tracking.ActivePass.extension)

    kwargs.setdefault("useZvertexTool", False)
    kwargs.setdefault("useZBoundFinding", flags.ITk.Tracking.ActivePass.doZBoundary)
    kwargs.setdefault("ITKGeometry", True)
    kwargs.setdefault("SpacePointsSCTName", "ITkStripSpacePoints")
    kwargs.setdefault("SpacePointsPixelName", "ITkPixelSpacePoints")

    if flags.ITk.Tracking.doFastTracking :
        kwargs.setdefault("doFastTracking", True)

        if 'InDetEtaDependentCutsSvc' not in kwargs :
            from InDetConfig.InDetEtaDependentCutsConfig import ITkEtaDependentCutsSvcCfg
            acc.merge(ITkEtaDependentCutsSvcCfg(flags))
            kwargs.setdefault("InDetEtaDependentCutsSvc", acc.getService("ITkEtaDependentCutsSvc"+flags.ITk.Tracking.ActivePass.extension))

    acc.addEventAlgo(CompFactory.InDet.SiSPSeededTrackFinder(name+flags.ITk.Tracking.ActivePass.extension, **kwargs))
    return acc

def ITkSiSPSeededTrackFinderROIConvCfg(flags, name="ITkSiSpTrackFinderROIConv", **kwargs) :
    from InDetConfig.InDetCaloClusterROISelectorConfig import ITkCaloClusterROIPhiRZContainerMakerCfg
    acc = ITkCaloClusterROIPhiRZContainerMakerCfg(flags)

    if "RegSelTool_Strip" not in kwargs:
        from RegionSelector.RegSelToolConfig import regSelTool_ITkStrip_Cfg
        kwargs.setdefault("RegSelTool_Strip", acc.popToolsAndMerge(
            regSelTool_ITkStrip_Cfg(flags)))

    kwargs.setdefault("useITkConvSeeded", True)
    kwargs.setdefault("EMROIPhiRZContainer", "ITkCaloClusterROIPhiRZ15GeVUnordered")

    acc.merge(ITkSiSPSeededTrackFinderCfg(flags, name, **kwargs))
    return acc
