# Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration

from AthenaCommon.Logging import logging
from AthenaConfiguration.ComponentFactory import CompFactory
from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
egammaTrackSlimmer = CompFactory.egammaTrackSlimmer


def egammaTrackSlimmerCfg(
        flags,
        name='egammaTrackSlimmer',
        **kwargs):

    mlog = logging.getLogger(name)
    mlog.info('Starting configuration')

    acc = ComponentAccumulator()

    kwargs.setdefault("StreamName", 'StreamAOD')
    kwargs.setdefault("InputElectronContainerName",
                      flags.Egamma.Keys.Output.Electrons)
    kwargs.setdefault("InputPhotonContainerName",
                      flags.Egamma.Keys.Output.Photons)
    kwargs.setdefault("TrackParticleContainerName",
                      flags.Egamma.Keys.Output.GSFTrackParticles)
    kwargs.setdefault("VertexContainerName",
                      flags.Egamma.Keys.Output.ConversionVertices)

    egtrkslimmerAlg = egammaTrackSlimmer(name, **kwargs)

    acc.addEventAlgo(egtrkslimmerAlg)
    return acc


if __name__ == "__main__":
    from AthenaCommon.Configurable import Configurable
    Configurable.configurableRun3Behavior = True
    from AthenaConfiguration.AllConfigFlags import ConfigFlags as flags
    from AthenaConfiguration.TestDefaults import defaultTestFiles
    from AthenaConfiguration.ComponentAccumulator import printProperties
    from AthenaConfiguration.MainServicesConfig import MainServicesCfg
    flags.Input.Files = defaultTestFiles.ESD

    acc = MainServicesCfg(flags)
    acc.merge(egammaTrackSlimmerCfg(flags))
    mlog = logging.getLogger("egammaTrackSlimmerConfigTest")
    mlog.info("Configuring  egammaTrackSlimmer: ")
    printProperties(mlog,
                    acc.getEventAlgo("egammaTrackSlimmer"),
                    nestLevel=1,
                    printDefaults=True)
    with open("egammatrackslimmer.pkl", "wb") as f:
        acc.store(f)
