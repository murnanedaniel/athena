# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
# Configuration of pixel tools of SiClusterOnTrackTool package
from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory
from AthenaConfiguration.Enums import BeamType

#############################################
### InDet PixelClusterOnTrackTool offline ###
#############################################

def InDetPixelClusterOnTrackToolBaseCfg(flags, name="PixelClusterOnTrackTool", **kwargs):
    from PixelConditionsAlgorithms.PixelConditionsConfig import PixelDistortionAlgCfg, PixelOfflineCalibCondAlgCfg
    acc = PixelOfflineCalibCondAlgCfg(flags)    # To produce PixelOfflineCalibData
    if not flags.InDet.Tracking.doDBMstandalone:
        acc.merge(PixelDistortionAlgCfg(flags)) # To produce PixelDistortionData

    from TrkConfig.TrkRIO_OnTrackCreatorConfig import RIO_OnTrackErrorScalingCondAlgCfg
    acc.merge(RIO_OnTrackErrorScalingCondAlgCfg(flags)) # To produce RIO_OnTrackErrorScaling

    if 'LorentzAngleTool' not in kwargs:
        from SiLorentzAngleTool.PixelLorentzAngleConfig import PixelLorentzAngleToolCfg
        kwargs.setdefault("LorentzAngleTool", acc.popToolsAndMerge(
            PixelLorentzAngleToolCfg(flags)))

    if flags.Beam.Type is BeamType.Cosmics or flags.InDet.Tracking.doDBMstandalone:
        kwargs.setdefault("ErrorStrategy", 0)
        kwargs.setdefault("PositionStrategy", 0)

    kwargs.setdefault("DisableDistortions", flags.InDet.Tracking.doDBMstandalone)
    kwargs.setdefault("applyNNcorrection", flags.InDet.Tracking.doPixelClusterSplitting and flags.InDet.Tracking.pixelClusterSplittingType == "NeuralNet")
    kwargs.setdefault("NNIBLcorrection", flags.InDet.Tracking.doPixelClusterSplitting and flags.InDet.Tracking.pixelClusterSplittingType == "NeuralNet")
    split_cluster_map_extension = flags.InDet.Tracking.ActivePass.extension if flags.InDet.Tracking.ActivePass.useTIDE_Ambi else ""
    kwargs.setdefault("SplitClusterAmbiguityMap", f"SplitClusterAmbiguityMap{split_cluster_map_extension}")
    kwargs.setdefault("RunningTIDE_Ambi", flags.InDet.Tracking.doTIDE_Ambi)

    acc.setPrivateTools(CompFactory.InDet.PixelClusterOnTrackTool(name, **kwargs))
    return acc

def InDetPixelClusterOnTrackToolDigitalCfg(flags, name="InDetPixelClusterOnTrackToolDigital", **kwargs):
    kwargs.setdefault("SplitClusterAmbiguityMap", "")

    if flags.InDet.Tracking.doDigitalROTCreation:
        kwargs.setdefault("applyNNcorrection", False)
        kwargs.setdefault("NNIBLcorrection", False)
        kwargs.setdefault("ErrorStrategy", 2)
        kwargs.setdefault("PositionStrategy", 1)

    return InDetPixelClusterOnTrackToolBaseCfg(flags, name, **kwargs)

def InDetPixelClusterOnTrackToolNNSplittingCfg(flags, name="InDetPixelClusterOnTrackToolNNSplitting", **kwargs):
    acc = ComponentAccumulator()

    if flags.InDet.Tracking.doPixelClusterSplitting \
       and flags.InDet.Tracking.pixelClusterSplittingType == "NeuralNet" \
       and "NnClusterizationFactory" not in kwargs:
        from InDetConfig.SiClusterizationToolConfig import NnClusterizationFactoryCfg
        kwargs.setdefault("NnClusterizationFactory", acc.popToolsAndMerge(
            NnClusterizationFactoryCfg(flags)))

    acc.setPrivateTools(acc.popToolsAndMerge(
        InDetPixelClusterOnTrackToolBaseCfg(flags, name, **kwargs)))
    return acc

def InDetPixelClusterOnTrackToolCfg(flags, name="InDetPixelClusterOnTrackTool", **kwargs):
    if flags.InDet.Tracking.doDigitalROTCreation:
        return InDetPixelClusterOnTrackToolDigitalCfg(flags, name, **kwargs)
    else:
        return InDetPixelClusterOnTrackToolNNSplittingCfg(flags, name, **kwargs)


def InDetBroadPixelClusterOnTrackToolCfg(flags, name='InDetBroadPixelClusterOnTrackTool', **kwargs):
    kwargs.setdefault("ErrorStrategy", 0)
    return InDetPixelClusterOnTrackToolCfg(flags, name, **kwargs)

def InDetPixelClusterOnTrackToolDBMCfg(flags, name='InDetPixelClusterOnTrackToolDBM', **kwargs):
    kwargs.setdefault("DisableDistortions", True )
    kwargs.setdefault("applyNNcorrection", False )
    kwargs.setdefault("NNIBLcorrection", False )
    kwargs.setdefault("RunningTIDE_Ambi", False )
    kwargs.setdefault("ErrorStrategy", 0 )
    kwargs.setdefault("PositionStrategy", 0 )
    return InDetPixelClusterOnTrackToolBaseCfg(flags, name, **kwargs)

#############################################
### InDet PixelClusterOnTrackTool trigger ###
#############################################

def TrigPixelClusterOnTrackToolBaseCfg(flags, name="InDetTrigPixelClusterOnTrackTool", **kwargs):
    from PixelConditionsAlgorithms.PixelConditionsConfig import PixelDistortionAlgCfg, PixelOfflineCalibCondAlgCfg
    acc = PixelOfflineCalibCondAlgCfg(flags) # To produce PixelOfflineCalibData
    acc.merge(PixelDistortionAlgCfg(flags))  # To produce PixelDistortionData

    from TrkConfig.TrkRIO_OnTrackCreatorConfig import RIO_OnTrackErrorScalingCondAlgCfg
    acc.merge(RIO_OnTrackErrorScalingCondAlgCfg(flags)) # To produce RIO_OnTrackErrorScaling

    if 'LorentzAngleTool' not in kwargs:
        from SiLorentzAngleTool.PixelLorentzAngleConfig import PixelLorentzAngleToolCfg
        kwargs.setdefault("LorentzAngleTool", acc.popToolsAndMerge(
            PixelLorentzAngleToolCfg(flags)) )

    if 'NnClusterizationFactory' not in kwargs:
        from InDetConfig.SiClusterizationToolConfig import TrigNnClusterizationFactoryCfg
        kwargs.setdefault("NnClusterizationFactory", acc.popToolsAndMerge(
            TrigNnClusterizationFactoryCfg(flags)))

    kwargs.setdefault("ErrorStrategy", 2)
    kwargs.setdefault("SplitClusterAmbiguityMap", "TrigPixelClusterAmbiguitiesMap")

    acc.setPrivateTools(CompFactory.InDet.PixelClusterOnTrackTool(name, **kwargs))
    return acc

###########################################
### ITk PixelClusterOnTrackTool offline ###
###########################################

def ITkPixelClusterOnTrackToolBaseCfg(flags, name="ITkPixelClusterOnTrackTool", **kwargs):
    from PixelConditionsAlgorithms.ITkPixelConditionsConfig import ITkPixelOfflineCalibCondAlgCfg
    acc = ITkPixelOfflineCalibCondAlgCfg(flags) # To produce PixelOfflineCalibData

    if 'LorentzAngleTool' not in kwargs:
        from SiLorentzAngleTool.ITkPixelLorentzAngleConfig import ITkPixelLorentzAngleToolCfg
        kwargs.setdefault("LorentzAngleTool", acc.popToolsAndMerge(
            ITkPixelLorentzAngleToolCfg(flags)))

    if flags.Beam.Type is BeamType.Cosmics:
        kwargs.setdefault("ErrorStrategy", 0)
        kwargs.setdefault("PositionStrategy", 0)

    kwargs.setdefault("applyNNcorrection", False )
    split_cluster_map_extension = flags.ITk.Tracking.ActivePass.extension
    kwargs.setdefault("SplitClusterAmbiguityMap", f"SplitClusterAmbiguityMap{split_cluster_map_extension}")
    kwargs.setdefault("RunningTIDE_Ambi", True )

    kwargs.setdefault("PixelErrorScalingKey", "")

    acc.setPrivateTools(CompFactory.ITk.PixelClusterOnTrackTool(name, **kwargs))
    return acc

def ITkPixelClusterOnTrackToolTruthSplittingCfg(flags, name='ITkPixelClusterOnTrackToolTruthSplitting', **kwargs):
    acc = ComponentAccumulator()

    if flags.ITk.Tracking.doPixelClusterSplitting \
       and flags.ITk.Tracking.pixelClusterSplittingType == "Truth":
        if 'NnClusterizationFactory' not in kwargs:
            from InDetConfig.SiClusterizationToolConfig import ITkTruthClusterizationFactoryCfg
            kwargs.setdefault("NnClusterizationFactory", acc.popToolsAndMerge(
                ITkTruthClusterizationFactoryCfg(flags))) #Truth-based for ITk for now

    acc.setPrivateTools(acc.popToolsAndMerge(
        ITkPixelClusterOnTrackToolBaseCfg(flags, name, **kwargs)))
    return acc

def ITkPixelClusterOnTrackToolCfg(flags, name='ITkPixelClusterOnTrackTool', **kwargs):
    if flags.ITk.Tracking.doDigitalClustering:
        kwargs.setdefault("PositionStrategy", 0)
        kwargs.setdefault("ErrorStrategy", 1)

    return ITkPixelClusterOnTrackToolTruthSplittingCfg(flags, name, **kwargs)

def ITkBroadPixelClusterOnTrackToolCfg(flags, name='ITkBroadPixelClusterOnTrackTool', **kwargs):
    kwargs.setdefault("ErrorStrategy", 0)
    return ITkPixelClusterOnTrackToolCfg(flags, name, **kwargs)
