# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory

def ActsTrkPixelSpacePointToolCfg(flags, name = "ActsTrkPixelSpacePointTool", **kwargs):
    acc = ComponentAccumulator()
    acc.setPrivateTools(CompFactory.ActsTrk.PixelSpacePointFormationTool(name, **kwargs))
    return acc

def ActsTrkStripSpacePointToolCfg(flags, name = "ActsTrkStripSpacePointTool", **kwargs):
    acc = ComponentAccumulator()

    from SiLorentzAngleTool.ITkStripLorentzAngleConfig import ITkStripLorentzAngleToolCfg
    kwargs.setdefault("LorentzAngleTool", acc.popToolsAndMerge(ITkStripLorentzAngleToolCfg(flags)) )
    kwargs.setdefault("AllClusters", False)

    acc.setPrivateTools(CompFactory.ActsTrk.StripSpacePointFormationTool(name, **kwargs))
    return acc

def ActsTrkPixelSpacePointFormationCfg(ConfigFlags,
                                       name = "ActsTrkPixelSpacePointFormation",
                                       **kwargs):

    from PixelGeoModelXml.ITkPixelGeoModelConfig import ITkPixelReadoutGeometryCfg
    acc = ITkPixelReadoutGeometryCfg(ConfigFlags)

    ActsTrkPixelSpacePointTool = acc.popToolsAndMerge(ActsTrkPixelSpacePointToolCfg(ConfigFlags))
    kwargs.setdefault("SpacePointFormationTool", ActsTrkPixelSpacePointTool)

    kwargs.setdefault("PixelClusters", "ITkPixelClusters")
    kwargs.setdefault("PixelDetectorElements", "ITkPixelDetectorElementCollection")

    kwargs.setdefault("PixelSpacePoints", "ITkPixelSpacePoints")
    kwargs.setdefault("PixelSpacePointData", "ITkPixelSpacePointData")

    if ConfigFlags.Acts.doMonitoring:
        from ActsTrkAnalysis.ActsTrkLiveMonitoringConfig import ActsTrkPixelSpacePointFormatioLiveMonitoringToolCfg
        kwargs.setdefault("MonTool", acc.popToolsAndMerge(ActsTrkPixelSpacePointFormatioLiveMonitoringToolCfg(ConfigFlags)))

    acc.addEventAlgo(CompFactory.ActsTrk.PixelSpacePointFormationAlg(name, **kwargs))
    return acc

def ActsTrkStripSpacePointFormationCfg(ConfigFlags,
                                       name = "ActsTrkStripSpacePointFormation",
                                       **kwargs):

    from StripGeoModelXml.ITkStripGeoModelConfig import ITkStripReadoutGeometryCfg
    acc = ITkStripReadoutGeometryCfg(ConfigFlags)

    ActsTrkStripSpacePointTool = acc.popToolsAndMerge(ActsTrkStripSpacePointToolCfg(ConfigFlags))
    kwargs.setdefault("SpacePointFormationTool", ActsTrkStripSpacePointTool)

    kwargs.setdefault("StripClusters", "ITkStripClusters")
    kwargs.setdefault("StripDetectorElements", "ITkStripDetectorElementCollection")
    kwargs.setdefault("StripElementPropertiesTable", "ITkStripElementPropertiesTable")

    kwargs.setdefault("StripSpacePoints", "ITkStripSpacePoints")
    kwargs.setdefault("StripSpacePointData", "ITkStripSpacePointData")
    kwargs.setdefault("StripOverlapSpacePoints", "ITkStripOverlapSpacePoints")
    kwargs.setdefault("StripOverlapSpacePointData", "ITkStripOverlapSpacePointData")
    kwargs.setdefault("ProcessOverlapForStrip", True)

    if ConfigFlags.Acts.doMonitoring:
        from ActsTrkAnalysis.ActsTrkLiveMonitoringConfig import ActsTrkStripSpacePointFormatioLiveMonitoringToolCfg
        kwargs.setdefault("MonTool", acc.popToolsAndMerge(ActsTrkStripSpacePointFormatioLiveMonitoringToolCfg(ConfigFlags)))

    acc.addEventAlgo(CompFactory.ActsTrk.StripSpacePointFormationAlg(name, **kwargs))
    return acc

def ActsTrkSpacePointFormationCfg(flags):
    acc = ComponentAccumulator()
    if flags.Detector.EnableITkPixel:
        acc.merge(ActsTrkPixelSpacePointFormationCfg(flags))
    if flags.Detector.EnableITkStrip:
        # Need to schedule this here in case the Athena space point formation is not schedule
        # This is because as of now requires at least ITkSiElementPropertiesTableCondAlgCfg
        # This may be because the current strip space point formation algorithm is not using Acts
        # May be not necessary once the Acts-based strip space point maker is ready
        from StripGeoModelXml.ITkStripGeoModelConfig import ITkStripReadoutGeometryCfg
        acc.merge(ITkStripReadoutGeometryCfg(flags))
        
        from BeamSpotConditions.BeamSpotConditionsConfig import BeamSpotCondAlgCfg
        acc.merge(BeamSpotCondAlgCfg(flags))
        
        from InDetConfig.SiSpacePointFormationConfig import ITkSiElementPropertiesTableCondAlgCfg
        acc.merge(ITkSiElementPropertiesTableCondAlgCfg(flags))
        
        acc.merge(ActsTrkStripSpacePointFormationCfg(flags))

    if flags.Acts.doAnalysis:
        from ActsTrkAnalysis.ActsTrkAnalysisConfig import ActsTrkSpacePointAnalysisCfg
        acc.merge(ActsTrkSpacePointAnalysisCfg(flags))
        
    return acc
