"""Define methods to construct configured BCM overlay algorithms

Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
"""

from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory


def BCMRawDataProviderAlgCfg(flags, name="BCMRawDataProvider", **kwargs):
    """Return a ComponentAccumulator for BCM raw data provider"""
    # Temporary until available in the central location
    acc = ComponentAccumulator()

    kwargs.setdefault("RDOKey", f"{flags.Overlay.BkgPrefix}BCM_RDOs")

    acc.addEventAlgo(CompFactory.BCM_RawDataProvider(name, **kwargs))

    return acc


def BCMDataOverlayExtraCfg(flags, **kwargs):
    """Return a ComponentAccumulator with BCM data overlay specifics"""
    acc = ComponentAccumulator()

    # We need to convert BS to RDO for data overlay
    acc.merge(BCMRawDataProviderAlgCfg(flags))

    return acc


def BCMOverlayAlgCfg(flags, name="BCMOverlay", **kwargs):
    """Return a ComponentAccumulator for BCMOverlay algorithm"""
    acc = ComponentAccumulator()

    kwargs.setdefault("BkgInputKey", f"{flags.Overlay.BkgPrefix}BCM_RDOs")
    kwargs.setdefault("SignalInputKey", f"{flags.Overlay.SigPrefix}BCM_RDOs")
    kwargs.setdefault("OutputKey", "BCM_RDOs")

    if not flags.Overlay.DataOverlay:
        from SGComps.SGInputLoaderConfig import SGInputLoaderCfg
        acc.merge(SGInputLoaderCfg(flags, [f'BCM_RDO_Container#{kwargs["BkgInputKey"]}']))

    kwargs.setdefault("isDataOverlay", flags.Overlay.DataOverlay)

    # Do BCM overlay
    acc.addEventAlgo(CompFactory.BCMOverlay(name, **kwargs))

    # Setup output
    if flags.Output.doWriteRDO:
        from OutputStreamAthenaPool.OutputStreamConfig import OutputStreamCfg
        acc.merge(OutputStreamCfg(flags, "RDO", ItemList=[
            "BCM_RDO_Container#BCM_RDOs"
        ]))

    if flags.Output.doWriteRDO_SGNL:
        from OutputStreamAthenaPool.OutputStreamConfig import OutputStreamCfg
        acc.merge(OutputStreamCfg(flags, "RDO_SGNL", ItemList=[
            f"BCM_RDO_Container#{flags.Overlay.SigPrefix}BCM_RDOs"
        ]))

    return acc


def BCMTruthOverlayCfg(flags, name="BCMSDOOverlay", **kwargs):
    """Return a ComponentAccumulator for the BCM SDO overlay algorithm"""
    acc = ComponentAccumulator()

    # We do not need background BCM SDOs
    kwargs.setdefault("BkgInputKey", "")

    kwargs.setdefault("SignalInputKey", f"{flags.Overlay.SigPrefix}BCM_SDO_Map")
    kwargs.setdefault("OutputKey", "BCM_SDO_Map")

    # Do BCM truth overlay
    acc.addEventAlgo(CompFactory.InDetSDOOverlay(name, **kwargs))

    # Setup output
    if flags.Output.doWriteRDO:
        from OutputStreamAthenaPool.OutputStreamConfig import OutputStreamCfg
        acc.merge(OutputStreamCfg(flags, "RDO", ItemList=[
            "InDetSimDataCollection#BCM_SDO_Map"
        ]))
    
    if flags.Output.doWriteRDO_SGNL:
        from OutputStreamAthenaPool.OutputStreamConfig import OutputStreamCfg
        acc.merge(OutputStreamCfg(flags, "RDO_SGNL", ItemList=[
            f"InDetSimDataCollection#{flags.Overlay.SigPrefix}BCM_SDO_Map"
        ]))

    return acc


def BCMOverlayCfg(flags):
    """Configure and return a ComponentAccumulator for BCM overlay"""
    acc = ComponentAccumulator()

    # Add data overlay specifics
    if flags.Overlay.DataOverlay:
        acc.merge(BCMDataOverlayExtraCfg(flags))

    # Add BCM overlay digitization algorithm
    from BCM_Digitization.BCM_DigitizationConfigNew import BCM_OverlayDigitizationBasicCfg
    acc.merge(BCM_OverlayDigitizationBasicCfg(flags))
    # Add BCM overlay algorithm
    acc.merge(BCMOverlayAlgCfg(flags))
    # Add BCM truth overlay
    acc.merge(BCMTruthOverlayCfg(flags))

    return acc
