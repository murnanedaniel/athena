"""Define methods to construct configured BCM Digitization tools and algs

Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
"""
from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory
from AthenaConfiguration.Enums import ProductionStep
from RngComps.RandomServices import AthRNGSvcCfg
from PixelGeoModel.PixelGeoModelConfig import PixelReadoutGeometryCfg
from OutputStreamAthenaPool.OutputStreamConfig import OutputStreamCfg
from Digitization.TruthDigitizationOutputConfig import TruthDigitizationOutputCfg
from Digitization.PileUpToolsConfig import PileUpToolsCfg
from Digitization.PileUpMergeSvcConfigNew import PileUpMergeSvcCfg, PileUpXingFolderCfg


# The earliest and last bunch crossing times for which interactions will be sent
# to the BCM Digitization code.
def BCM_FirstXing():
    return -25


def BCM_LastXing():
    return 0


def BCM_RangeCfg(flags, name="BCM_Range", **kwargs):
    """Return a BCM configured PileUpXingFolder tool"""
    kwargs.setdefault("FirstXing", BCM_FirstXing())
    kwargs.setdefault("LastXing",  BCM_LastXing())
    # Default 0 no dataproxy reset
    kwargs.setdefault("CacheRefreshFrequency", 1.0)
    kwargs.setdefault("ItemList", ["SiHitCollection#BCMHits"])
    return PileUpXingFolderCfg(flags, name, **kwargs)


def BCM_DigitizationToolCommonCfg(flags, name="BCM_DigitizationTool", **kwargs):
    """Return a ComponentAccumulator with configured BCM_DigitizationTool"""
    acc = ComponentAccumulator()
    kwargs.setdefault("RndmSvc",
                      acc.getPrimaryAndMerge(AthRNGSvcCfg(flags)).name)
    kwargs.setdefault("HitCollName", "BCMHits")
    if flags.Common.ProductionStep == ProductionStep.PileUpPresampling:
        kwargs.setdefault("OutputRDOKey", flags.Overlay.BkgPrefix + "BCM_RDOs")
        kwargs.setdefault("OutputSDOKey", flags.Overlay.BkgPrefix + "BCM_SDO_Map")
    elif flags.Common.isOverlay:
        kwargs.setdefault("OnlyUseContainerName", False)
        kwargs.setdefault("OutputRDOKey", flags.Overlay.SigPrefix + "BCM_RDOs")
        kwargs.setdefault("OutputSDOKey", flags.Overlay.SigPrefix + "BCM_SDO_Map")
    else:
        kwargs.setdefault("OutputRDOKey", "BCM_RDOs")
        kwargs.setdefault("OutputSDOKey", "BCM_SDO_Map")
    if flags.Digitization.DoInnerDetectorNoise:
        kwargs.setdefault("ModNoise", [90.82] * 8)
    else:
        kwargs.setdefault("ModNoise", [0.0] * 8)
    kwargs.setdefault("ModSignal", [450.0] * 8)
    kwargs.setdefault("NinoThr", [330.0] * 8)
    # BCM with diamond
    kwargs.setdefault("MIPDeposit", 0.33)
    # Alternative BCM with graphite
    # kwargs.setdefault("MIPDeposit", 0.25)
    kwargs.setdefault("EffDistanceParam", 4.0)
    kwargs.setdefault("EffSharpnessParam", 0.11)
    kwargs.setdefault("TimeDelay", 9.0)
    if flags.Digitization.DoXingByXingPileUp:
        kwargs.setdefault("FirstXing", BCM_FirstXing())
        kwargs.setdefault("LastXing",  BCM_LastXing())
    kwargs.setdefault("RndmSvc", acc.getPrimaryAndMerge(AthRNGSvcCfg(flags)).name)

    BCM_DigitizationTool = CompFactory.BCM_DigitizationTool
    acc.setPrivateTools(BCM_DigitizationTool(name, **kwargs))
    return acc


def BCM_DigitizationToolCfg(flags, name="BCM_DigitizationTool", **kwargs):
    """Return ComponentAccumulator with BCM_DigitizationTool for non-overlay"""
    acc = ComponentAccumulator()
    if flags.Digitization.PileUp:
        intervals = []
        if not flags.Digitization.DoXingByXingPileUp:
            intervals += [acc.popToolsAndMerge(BCM_RangeCfg(flags))]
        kwargs.setdefault("MergeSvc", acc.getPrimaryAndMerge(PileUpMergeSvcCfg(flags, Intervals=intervals)).name)
    else:
        kwargs.setdefault("MergeSvc", '')
    kwargs.setdefault("OnlyUseContainerName", flags.Digitization.PileUp)
    tool = acc.popToolsAndMerge(BCM_DigitizationToolCommonCfg(flags, name, **kwargs))
    acc.setPrivateTools(tool)
    return acc


def BCM_OverlayDigitizationToolCfg(flags, name="BCM_OverlayDigitizationTool", **kwargs):
    """Return ComponentAccumulator with BCM_DigitizationTool for Overlay"""
    acc = ComponentAccumulator()
    kwargs.setdefault("MergeSvc", '')
    kwargs.setdefault("OnlyUseContainerName", False)
    tool = acc.popToolsAndMerge(BCM_DigitizationToolCommonCfg(flags, name, **kwargs))
    acc.setPrivateTools(tool)
    return acc


def BCM_OutputCfg(flags):
    """Return ComponentAccumulator with Output for BCM. Not standalone."""
    acc = ComponentAccumulator()
    if flags.Output.doWriteRDO:
        ItemList = ["BCM_RDO_Container#*"]
        if flags.Digitization.TruthOutput:
            ItemList += ["InDetSimDataCollection#*"]
            acc.merge(TruthDigitizationOutputCfg(flags))
        acc.merge(OutputStreamCfg(flags, "RDO", ItemList))
    return acc


def BCM_DigitizationBasicCfg(flags, **kwargs):
    """Return ComponentAccumulator for BCM digitization"""
    acc = PixelReadoutGeometryCfg(flags)
    if "PileUpTools" not in kwargs:
        PileUpTools = acc.popToolsAndMerge(BCM_DigitizationToolCfg(flags))
        kwargs["PileUpTools"] = PileUpTools
    acc.merge(PileUpToolsCfg(flags, **kwargs))
    return acc


def BCM_OverlayDigitizationBasicCfg(flags, name="BCM_OverlayDigitization", **kwargs):
    """Return ComponentAccumulator with BCM Overlay digitization"""
    acc = PixelReadoutGeometryCfg(flags)
    if flags.Common.ProductionStep != ProductionStep.FastChain:
        from SGComps.SGInputLoaderConfig import SGInputLoaderCfg
        acc.merge(SGInputLoaderCfg(flags, ["SiHitCollection#BCMHits"]))

    if "DigitizationTool" not in kwargs:
        tool = acc.popToolsAndMerge(BCM_OverlayDigitizationToolCfg(flags))
        kwargs["DigitizationTool"] = tool

    if flags.Concurrency.NumThreads > 0:
        kwargs.setdefault("Cardinality", flags.Concurrency.NumThreads)

    # Set common overlay extra inputs
    kwargs.setdefault("ExtraInputs", flags.Overlay.ExtraInputs)

    BCM_Digitization = CompFactory.BCM_Digitization
    acc.addEventAlgo(BCM_Digitization(name, **kwargs))
    return acc


# with output defaults
def BCM_DigitizationCfg(flags, **kwargs):
    """Return ComponentAccumulator for BCM digitization and Output"""
    acc = BCM_DigitizationBasicCfg(flags, **kwargs)
    acc.merge(BCM_OutputCfg(flags))
    return acc


def BCM_OverlayDigitizationCfg(flags, **kwargs):
    """Return ComponentAccumulator with BCM Overlay digitization and Output"""
    acc = BCM_OverlayDigitizationBasicCfg(flags, **kwargs)
    acc.merge(BCM_OutputCfg(flags))
    return acc
