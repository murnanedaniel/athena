"""Define methods to construct configured EventInfo conversion algorithms

Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
"""

from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory
from AthenaConfiguration.Enums import LHCPeriod


def EventInfoCnvAlgCfg(flags, name="EventInfoCnvAlg",
                       inputKey="McEventInfo",
                       outputKey="EventInfo",
                       disableBeamSpot=False,
                       **kwargs):
    """Return a ComponentAccumulator for EventInfoCnvAlg algorithm"""

    acc = ComponentAccumulator()

    kwargs.setdefault("AODKey", inputKey)
    kwargs.setdefault("xAODKey", outputKey)

    # TODO: luminosity

    if not disableBeamSpot:
        from BeamSpotConditions.BeamSpotConditionsConfig import BeamSpotCondAlgCfg
        acc.merge(BeamSpotCondAlgCfg(flags))

    xAODMaker__EventInfoCnvAlg = CompFactory.xAODMaker.EventInfoCnvAlg
    alg = xAODMaker__EventInfoCnvAlg(name, **kwargs)
    acc.addEventAlgo(alg)

    return acc


def EventInfoOverlayAlgCfg(flags, name="EventInfoOverlay", **kwargs):
    """Return a ComponentAccumulator for EventInfoOverlay algorithm"""
    acc = ComponentAccumulator()

    # Add beam spot conditions
    from BeamSpotConditions.BeamSpotConditionsConfig import BeamSpotCondAlgCfg
    acc.merge(BeamSpotCondAlgCfg(flags))

    kwargs.setdefault("BkgInputKey", flags.Overlay.BkgPrefix + "EventInfo")
    kwargs.setdefault("SignalInputKey", flags.Overlay.SigPrefix + "EventInfo")
    kwargs.setdefault("OutputKey", "EventInfo")

    kwargs.setdefault("DataOverlay", flags.Overlay.DataOverlay)
    kwargs.setdefault("ValidateBeamSpot", not flags.Overlay.DataOverlay and flags.GeoModel.Run is LHCPeriod.Run3)

    if flags.Input.MCChannelNumber > 0:
        kwargs.setdefault("MCChannelNumber", flags.Input.MCChannelNumber)

    # Do the xAOD::EventInfo overlay
    xAODMaker__EventInfoOverlay = CompFactory.xAODMaker.EventInfoOverlay
    alg = xAODMaker__EventInfoOverlay(name, **kwargs)
    acc.addEventAlgo(alg)

    # Re-map signal address
    from SGComps.AddressRemappingConfig import AddressRemappingCfg
    acc.merge(AddressRemappingCfg([
        "xAOD::EventInfo#EventInfo->" + flags.Overlay.SigPrefix + "EventInfo",
        "xAOD::EventAuxInfo#EventInfoAux.->" + flags.Overlay.SigPrefix + "EventInfoAux.",
    ]))

    return acc


def EventInfoOverlayOutputCfg(flags, **kwargs):
    """Return event info overlay output configuration"""
    acc = ComponentAccumulator()

    from OutputStreamAthenaPool.OutputStreamConfig import OutputStreamCfg
    if flags.Output.doWriteRDO:
        acc.merge(OutputStreamCfg(flags, "RDO"))

    # Add signal output
    if flags.Output.doWriteRDO_SGNL:
        acc.merge(OutputStreamCfg(flags, "RDO_SGNL",
                                  ItemList=[
                                      "xAOD::EventInfo#Sig_EventInfo",
                                      "xAOD::EventAuxInfo#Sig_EventInfoAux."
                                  ]))

    return acc


def EventInfoOverlayCfg(flags, **kwargs):
    """Return a ComponentAccumulator for the full EventInfoOverlay algorithm accumulator"""

    acc = ComponentAccumulator()

    # Check if running on legacy HITS
    if "EventInfo" not in flags.Input.Collections and "EventInfo" not in flags.Input.SecondaryCollections:
        acc.merge(EventInfoCnvAlgCfg(flags,
                                     inputKey=flags.Overlay.SigPrefix+"McEventInfo",
                                     outputKey=flags.Overlay.SigPrefix+"EventInfo",
                                     **kwargs))
        # Re-map signal address
        from SGComps.AddressRemappingConfig import AddressRemappingCfg
        acc.merge(AddressRemappingCfg([
            "EventInfo#McEventInfo->" + flags.Overlay.SigPrefix + "McEventInfo",
        ]))

    acc.merge(EventInfoOverlayAlgCfg(flags, **kwargs))
    acc.merge(EventInfoOverlayOutputCfg(flags, **kwargs))
    return acc





def EventInfoUpdateFromContextAlgCfg(flags, name="EventInfoUpdateFromContextAlg", **kwargs):
    """Return a ComponentAccumulator for EventInfoUpdateFromContext algorithm"""
    acc = ComponentAccumulator()

    # Add beam spot conditions
    from BeamSpotConditions.BeamSpotConditionsConfig import BeamSpotCondAlgCfg
    acc.merge(BeamSpotCondAlgCfg(flags))

    kwargs.setdefault("SignalInputKey", "Input_EventInfo")
    kwargs.setdefault("OutputKey", "EventInfo")

    if flags.Input.MCChannelNumber > 0:
        kwargs.setdefault("MCChannelNumber", flags.Input.MCChannelNumber)

    # Do the xAOD::EventInfo overlay
    acc.addEventAlgo(CompFactory.xAODMaker.EventInfoUpdateFromContextAlg(name, **kwargs))

    # Re-map signal address
    from SGComps.AddressRemappingConfig import AddressRemappingCfg
    acc.merge(AddressRemappingCfg([
        "xAOD::EventInfo#EventInfo->" + "Input_EventInfo",
        "xAOD::EventAuxInfo#EventInfoAux.->" + "Input_EventInfoAux.",
    ]))

    return acc


def EventInfoUpdateFromContextOutputCfg(flags, **kwargs):
    """Return event info overlay output configuration"""
    acc = ComponentAccumulator()

    from OutputStreamAthenaPool.OutputStreamConfig import OutputStreamCfg
    if flags.Output.doWriteRDO:
        acc.merge(OutputStreamCfg(flags, "RDO"))

    return acc


def EventInfoUpdateFromContextCfg(flags, **kwargs):
    """Return a ComponentAccumulator for the full EventInfoOverlay algorithm accumulator"""

    acc = ComponentAccumulator()

    # Check if running on legacy HITS
    if "EventInfo" not in flags.Input.Collections and "EventInfo" not in flags.Input.SecondaryCollections:
        acc.merge(EventInfoCnvAlgCfg(flags,
                                     inputKey="McEventInfo",
                                     outputKey="Input_EventInfo",
                                     **kwargs))

    acc.merge(EventInfoUpdateFromContextAlgCfg(flags, **kwargs))
    acc.merge(EventInfoUpdateFromContextOutputCfg(flags, **kwargs))
    return acc
