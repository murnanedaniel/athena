"""ComponentAccumulator configuration for PileUpMergeSvc

Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
"""
from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory


def PileUpMergeSvcCfg(flags, name="PileUpMergeSvc", Intervals=[], **kwargs):
    """Return ComponentAccumulator with PileUpMergeSvc

    If doing XingByXingPileUp, Intervals should contian PileUpXingFolder tools.
    Otherwise it should be empty, and we enforce that here.
    """
    acc = ComponentAccumulator()

    if not flags.Digitization.DoXingByXingPileUp:
        # handle input type variety
        if not isinstance(Intervals, list):
            Intervals = [Intervals]
        kwargs["Intervals"] = Intervals

    kwargs.setdefault("EventInfoKeyName", "EventInfo") # FIXME Make default?
    acc.addService(CompFactory.PileUpMergeSvc(name, **kwargs), primary = True)
    return acc


def PileUpXingFolderCfg(flags, name="PileUpXingFolder" , **kwargs):
    """Return a configured PileUpXingFolder tool"""
    acc = ComponentAccumulator()
    acc.setPrivateTools(CompFactory.PileUpXingFolder(name, **kwargs))
    return acc
