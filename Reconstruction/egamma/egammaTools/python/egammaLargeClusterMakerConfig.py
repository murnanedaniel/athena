# Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration

__doc__ = """ Configure egammaLargeClusterMaker,
               which chooses cells to store in the AOD"""
__author__ = "Jovan Mitrevski"

from AthenaCommon.Logging import logging
from AthenaConfiguration.ComponentFactory import CompFactory
from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
egammaLargeClusterMaker = CompFactory.egammaLargeClusterMaker


def egammaLargeClusterMakerCfg(flags, name="egammaLCMakerTool",  **kwargs):

    acc = ComponentAccumulator()

    kwargs.setdefault("CellsName", flags.Egamma.Keys.Input.CaloCells)
    kwargs.setdefault("InputClusterCollection",
                      flags.Egamma.Keys.Output.CaloClusters)

    acc.setPrivateTools(egammaLargeClusterMaker(name, **kwargs))

    return acc


if __name__ == "__main__":

    from AthenaConfiguration.AllConfigFlags import ConfigFlags
    from AthenaConfiguration.ComponentAccumulator import printProperties
    from AthenaCommon.Configurable import Configurable
    from AthenaConfiguration.TestDefaults import defaultTestFiles
    Configurable.configurableRun3Behavior = True

    ConfigFlags.Input.Files = defaultTestFiles.RDO
    ConfigFlags.fillFromArgs()
    ConfigFlags.lock()
    ConfigFlags.dump()

    cfg = ComponentAccumulator()
    mlog = logging.getLogger("egammaLargeClusterMakerConfigTest")
    mlog.info("Configuring  egammaLargeClusterMaker: ")
    printProperties(mlog, cfg.popToolsAndMerge(
        egammaLargeClusterMakerCfg(ConfigFlags)),
        nestLevel=1,
        printDefaults=True)

    f = open("egammalargeclustermaker.pkl", "wb")
    cfg.store(f)
    f.close()
