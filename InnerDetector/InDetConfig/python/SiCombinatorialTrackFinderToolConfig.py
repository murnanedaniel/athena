# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
# Configuration of SiCombinatorialTrackFinderTool_xk package

from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory

def SiDetElementBoundaryLinksCondAlg_xk_Pixel_Cfg(flags, name = "InDetSiDetElementBoundaryLinksPixelCondAlg", **kwargs):
    from PixelGeoModel.PixelGeoModelConfig import PixelReadoutGeometryCfg
    acc = PixelReadoutGeometryCfg(flags) # To produce PixelDetectorElementCollection

    kwargs.setdefault("ReadKey", "PixelDetectorElementCollection")
    kwargs.setdefault("WriteKey", "PixelDetElementBoundaryLinks_xk")

    SiDetElementBoundaryLinksCondAlg = CompFactory.InDet.SiDetElementBoundaryLinksCondAlg_xk(name, **kwargs)
    acc.addEventAlgo(SiDetElementBoundaryLinksCondAlg)
    return acc

def SiDetElementBoundaryLinksCondAlg_xk_SCT_Cfg(flags, name = "InDetSiDetElementBoundaryLinksSCTCondAlg", **kwargs):
    from SCT_GeoModel.SCT_GeoModelConfig import SCT_ReadoutGeometryCfg
    acc = SCT_ReadoutGeometryCfg(flags) # To produce SCT_DetectorElementCollection

    kwargs.setdefault("ReadKey", "SCT_DetectorElementCollection")
    kwargs.setdefault("WriteKey", "SCT_DetElementBoundaryLinks_xk")

    SiDetElementBoundaryLinksCondAlg = CompFactory.InDet.SiDetElementBoundaryLinksCondAlg_xk(name, **kwargs)
    acc.addEventAlgo(SiDetElementBoundaryLinksCondAlg)
    return acc

def SiDetElementBoundaryLinksCondAlg_xk_ITkPixel_Cfg(flags, name = "ITkSiDetElementBoundaryLinksPixelCondAlg", **kwargs):
    from PixelGeoModelXml.ITkPixelGeoModelConfig import ITkPixelReadoutGeometryCfg
    acc = ITkPixelReadoutGeometryCfg(flags) # To produce ITkPixelDetectorElementCollection

    kwargs.setdefault("ReadKey", "ITkPixelDetectorElementCollection")
    kwargs.setdefault("WriteKey", "ITkPixelDetElementBoundaryLinks_xk")
    kwargs.setdefault("ITkGeometry", True)

    SiDetElementBoundaryLinksCondAlg = CompFactory.InDet.SiDetElementBoundaryLinksCondAlg_xk(name, **kwargs)
    acc.addEventAlgo(SiDetElementBoundaryLinksCondAlg)
    return acc

def SiDetElementBoundaryLinksCondAlg_xk_ITkStrip_Cfg(flags, name = "ITkSiDetElementBoundaryLinksStripCondAlg", **kwargs):
    from StripGeoModelXml.ITkStripGeoModelConfig import ITkStripReadoutGeometryCfg
    acc = ITkStripReadoutGeometryCfg(flags) # To produce ITkStripDetectorElementCollection

    kwargs.setdefault("ReadKey", "ITkStripDetectorElementCollection")
    kwargs.setdefault("WriteKey", "ITkStripDetElementBoundaryLinks_xk")
    kwargs.setdefault("ITkGeometry", True)

    SiDetElementBoundaryLinksCondAlg = CompFactory.InDet.SiDetElementBoundaryLinksCondAlg_xk(name, **kwargs)
    acc.addEventAlgo(SiDetElementBoundaryLinksCondAlg)
    return acc

def SiCombinatorialTrackFinder_xkCfg(flags, name="InDetSiComTrackFinder", **kwargs) :
    acc = ComponentAccumulator()

    # For SiDetElementBoundaryLinks_xk ReadCondHandle
    if flags.InDet.Tracking.ActivePass.usePixel:
        acc.merge(SiDetElementBoundaryLinksCondAlg_xk_Pixel_Cfg(flags))

    if flags.InDet.Tracking.ActivePass.useSCT:
        acc.merge(SiDetElementBoundaryLinksCondAlg_xk_SCT_Cfg(flags))

    #
    # --- Local track finding using sdCaloSeededSSSpace point seed
    #
    if flags.InDet.Tracking.doDBMstandalone:
        from InDetConfig.TrackingCommonConfig import InDetRotCreatorDBMCfg
        RotCreator = acc.popToolsAndMerge(InDetRotCreatorDBMCfg(flags))
        kwargs.setdefault("useSCT", False)
        kwargs.setdefault("MagneticFieldMode", "NoField")
        kwargs.setdefault("TrackQualityCut", 9.3)
    else:
        from InDetConfig.TrackingCommonConfig import InDetRotCreatorDigitalCfg
        RotCreator = acc.popToolsAndMerge(InDetRotCreatorDigitalCfg(flags))
        kwargs.setdefault("useSCT", flags.Detector.EnableSCT)

    acc.addPublicTool(RotCreator)
    kwargs.setdefault("RIOonTrackTool", RotCreator)

    from TrkConfig.TrkExRungeKuttaPropagatorConfig import RungeKuttaPropagatorCfg
    InDetPatternPropagator = acc.popToolsAndMerge(RungeKuttaPropagatorCfg(flags, name="InDetPatternPropagator"))
    acc.addPublicTool(InDetPatternPropagator)
    kwargs.setdefault("PropagatorTool", InDetPatternPropagator)

    from InDetConfig.TrackingCommonConfig import InDetPatternUpdatorCfg
    kwargs.setdefault("UpdatorTool", acc.getPrimaryAndMerge(InDetPatternUpdatorCfg()))

    from InDetConfig.InDetBoundaryCheckToolConfig import InDetBoundaryCheckToolCfg
    kwargs.setdefault("BoundaryCheckTool", acc.popToolsAndMerge(InDetBoundaryCheckToolCfg(flags)))
    
    kwargs.setdefault("usePixel", flags.Detector.EnablePixel)
    kwargs.setdefault("PixelClusterContainer", "PixelClusters")
    kwargs.setdefault("SCT_ClusterContainer", "SCT_Clusters")

    if flags.Detector.EnablePixel:
        from PixelConditionsTools.PixelConditionsSummaryConfig import PixelConditionsSummaryCfg
        kwargs.setdefault("PixelSummaryTool", acc.popToolsAndMerge(PixelConditionsSummaryCfg(flags)))
    else:
        kwargs.setdefault("PixelSummaryTool", "")

    if flags.Detector.EnableSCT:
        from SCT_ConditionsTools.SCT_ConditionsToolsConfig import SCT_ConditionsSummaryToolCfg
        kwargs.setdefault("SctSummaryTool", acc.popToolsAndMerge(SCT_ConditionsSummaryToolCfg(flags)))
    else:
        kwargs.setdefault("SctSummaryTool", "")

    track_finder = CompFactory.InDet.SiCombinatorialTrackFinder_xk(name = name+flags.InDet.Tracking.ActivePass.extension, **kwargs)
    acc.setPrivateTools(track_finder)
    return acc

def SiCombinatorialTrackFinder_xk_Trig_Cfg( flags, name="InDetTrigSiComTrackFinder", **kwargs ):
  """
  based  on: InnerDetector/InDetExample/InDetTrigRecExample/python/InDetTrigConfigRecLoadTools.py
  """
  acc = ComponentAccumulator()

  # For SiDetElementBoundaryLinks_xk ReadCondHandle
  if flags.InDet.Tracking.ActivePass.usePixel:
      acc.merge(SiDetElementBoundaryLinksCondAlg_xk_Pixel_Cfg(flags))

  if flags.InDet.Tracking.ActivePass.useSCT:
      acc.merge(SiDetElementBoundaryLinksCondAlg_xk_SCT_Cfg(flags))

  from TrkConfig.TrkExRungeKuttaPropagatorConfig import RungeKuttaPropagatorCfg
  propagatorTool = acc.popToolsAndMerge( RungeKuttaPropagatorCfg( flags, name="InDetTrigPatternPropagator" ) )
  acc.addPublicTool(propagatorTool)

  from TrigInDetConfig.TrigInDetConfig import KalmanxkUpdatorCfg, RIO_OnTrackCreatorCfg
  patternUpdatorTool = acc.getPrimaryAndMerge( KalmanxkUpdatorCfg( flags ) )
  rioOnTrackTool = acc.getPrimaryAndMerge( RIO_OnTrackCreatorCfg( flags ) )

  from PixelConditionsTools.PixelConditionsSummaryConfig import PixelConditionsSummaryCfg
  pixelCondSummaryTool = acc.popToolsAndMerge( PixelConditionsSummaryCfg(flags) )

  from SCT_ConditionsTools.SCT_ConditionsToolsConfig import SCT_ConditionsSummaryToolCfg
  sctCondSummaryTool = acc.popToolsAndMerge( SCT_ConditionsSummaryToolCfg( flags, withFlaggedCondTool=False, withTdaqTool=False ) )

  kwargs.setdefault("PropagatorTool", propagatorTool)
  kwargs.setdefault("UpdatorTool", patternUpdatorTool)
  kwargs.setdefault("RIOonTrackTool", rioOnTrackTool)
  kwargs.setdefault("usePixel", flags.Detector.EnablePixel)
  kwargs.setdefault("useSCT", flags.Detector.EnableSCT)
  kwargs.setdefault("PixelClusterContainer", 'PixelTrigClusters')
  kwargs.setdefault("SCT_ClusterContainer", 'SCT_TrigClusters')
  kwargs.setdefault("PixelSummaryTool", pixelCondSummaryTool)
  kwargs.setdefault("SctSummaryTool", sctCondSummaryTool)

  SiCombinatorialTrackFinder = CompFactory.InDet.SiCombinatorialTrackFinder_xk(name, **kwargs)
  acc.setPrivateTools( SiCombinatorialTrackFinder )
  return acc

def ITkSiCombinatorialTrackFinder_xkCfg(flags, name="ITkSiComTrackFinder", **kwargs) :
    acc = ComponentAccumulator()

    # For SiDetElementBoundaryLinks_xk ReadCondHandle
    if flags.ITk.Tracking.ActivePass.useITkPixel:
        acc.merge(SiDetElementBoundaryLinksCondAlg_xk_ITkPixel_Cfg(flags))

    if flags.ITk.Tracking.ActivePass.useITkStrip:
        acc.merge(SiDetElementBoundaryLinksCondAlg_xk_ITkStrip_Cfg(flags))

    #
    # --- Local track finding using sdCaloSeededSSSpace point seed
    #
    from InDetConfig.ITkTrackingCommonConfig import ITkRotCreatorDigitalCfg
    ITkRotCreatorDigital = acc.getPrimaryAndMerge(ITkRotCreatorDigitalCfg(flags))
    kwargs.setdefault("RIOonTrackTool", ITkRotCreatorDigital)

    from TrkConfig.TrkExRungeKuttaPropagatorConfig import RungeKuttaPropagatorCfg
    ITkPatternPropagator = acc.popToolsAndMerge(RungeKuttaPropagatorCfg(flags, name="ITkPatternPropagator"))
    acc.addPublicTool(ITkPatternPropagator)
    kwargs.setdefault("PropagatorTool", ITkPatternPropagator)

    from InDetConfig.ITkRecToolConfig import ITkPatternUpdatorCfg
    ITkPatternUpdator = acc.popToolsAndMerge(ITkPatternUpdatorCfg(flags))
    from InDetConfig.InDetBoundaryCheckToolConfig import ITkBoundaryCheckToolCfg
    ITkBoundaryCheckTool = acc.popToolsAndMerge(ITkBoundaryCheckToolCfg(flags))

    kwargs.setdefault("UpdatorTool", ITkPatternUpdator)
    kwargs.setdefault("BoundaryCheckTool", ITkBoundaryCheckTool)
    kwargs.setdefault("usePixel", flags.Detector.EnableITkPixel)
    kwargs.setdefault("useSCT", flags.Detector.EnableITkStrip)
    kwargs.setdefault("PixelClusterContainer", 'ITkPixelClusters')
    kwargs.setdefault("SCT_ClusterContainer", 'ITkStripClusters')
    kwargs.setdefault("PixelDetElementBoundaryLinks_xk", "ITkPixelDetElementBoundaryLinks_xk")
    kwargs.setdefault("SCT_DetElementBoundaryLinks_xk", "ITkStripDetElementBoundaryLinks_xk")
    kwargs.setdefault("SCTDetEleCollKey","ITkStripDetectorElementCollection")
    kwargs.setdefault("ITkGeometry", True)
    kwargs.setdefault("doFastTracking", flags.ITk.Tracking.doFastTracking)

    if flags.Detector.EnableITkPixel:
        from PixelConditionsTools.ITkPixelConditionsSummaryConfig import ITkPixelConditionsSummaryCfg
        kwargs.setdefault("PixelSummaryTool", acc.popToolsAndMerge(ITkPixelConditionsSummaryCfg(flags)))
    else:
        kwargs.setdefault("PixelSummaryTool", None)

    if flags.Detector.EnableITkStrip:
        from SCT_ConditionsTools.ITkStripConditionsToolsConfig import ITkStripConditionsSummaryToolCfg
        kwargs.setdefault("SctSummaryTool", acc.popToolsAndMerge(ITkStripConditionsSummaryToolCfg(flags)))
    else:
        kwargs.setdefault("SctSummaryTool", None)

    ITkSiComTrackFinder = CompFactory.InDet.SiCombinatorialTrackFinder_xk(name = name+flags.ITk.Tracking.ActivePass.extension, **kwargs)
    acc.setPrivateTools(ITkSiComTrackFinder)
    return acc
