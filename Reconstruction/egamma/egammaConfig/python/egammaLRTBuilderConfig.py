# Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration

__doc__ = """
          Instantiate the EGamma LRT reconstruction.
          Note that
          egammaTopoClusterCopier is scheduled in TrackRecoConfig
          """

from AthenaCommon.Logging import logging
from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from egammaTrackTools.egammaTrackToolsConfig import (
    EMExtrapolationToolsLRTCommonCacheCfg, EMExtrapolationToolsLRTCacheCfg)


def EGammaLRTReconstructionCfg(flags, name="EGammaLRTReconstruction"):

    mlog = logging.getLogger(name)
    mlog.info('Starting EGamma LRT reconstruction configuration')

    acc = ComponentAccumulator()

    # if large radius tracking not enabled add nothing
    if not flags.InDet.doR3LargeD0 or not flags.Egamma.enabled:
        if not flags.InDet.doR3LargeD0 and flags.Egamma.enabled:
            mlog.info('Large radius tracking not enabled. Do nothing')
        return acc

    # Add e/gamma tracking algorithms
    if flags.Egamma.doGSF:

        from egammaAlgs.egammaSelectedTrackCopyConfig import (
            egammaSelectedTrackCopyCfg)
        emextLRTCommonCache = acc.popToolsAndMerge(
            EMExtrapolationToolsLRTCommonCacheCfg(flags))
        acc.merge(egammaSelectedTrackCopyCfg(
            flags,
            name="egammaSelectedLRTTrackCopy",
            TrackParticleContainerName="InDetLargeD0TrackParticles",
            OutputTrkPartContainerName="LRTegammaSelectedTrackParticles",
            ExtrapolationToolCommonCache=emextLRTCommonCache)
        )

        from egammaAlgs.EMBremCollectionBuilderConfig import (
            EMBremCollectionBuilderCfg)
        acc.merge(EMBremCollectionBuilderCfg(
            flags,
            name='EMLRTBremCollectionBuilder',
            TrackParticleContainerName='InDetLargeD0TrackParticles',
            SelectedTrackParticleContainerName='LRTegammaSelectedTrackParticles',
            OutputTrkPartContainerName='LRT'+flags.Egamma.Keys.Output.GSFTrackParticles,
            OutputTrackContainerName='LRT'+flags.Egamma.Keys.Output.GSFTracks)
        )

        from egammaAlgs.EMGSFCaloExtensionBuilderConfig import (
            EMGSFCaloExtensionBuilderCfg)
        acc.merge(EMGSFCaloExtensionBuilderCfg(
            flags,
            name='EMGSFLRTCaloExtensionBuilder',
            GSFPerigeeCache='LRTGSFPerigeeCaloExtension',
            GSFLastCache='LRTGSFLastCaloExtension',
            GFFTrkPartContainerName='LRT'+flags.Egamma.Keys.Output.GSFTrackParticles)
        )

    # Add calo seeded central algorithms
    if flags.Egamma.doCaloSeeded:

        from egammaAlgs.egammaRecBuilderConfig import (
            egammaRecBuilderCfg)
        from egammaTools.EMTrackMatchBuilderConfig import (
            EMTrackMatchBuilderCfg)
        emextLRTCache = acc.popToolsAndMerge(
            EMExtrapolationToolsLRTCacheCfg(flags))
        lrtemtrackmatch = acc.popToolsAndMerge(EMTrackMatchBuilderCfg(
            flags,
            name='LRTEMTrackMatchBuilder',
            TrackParticlesName='LRT'+flags.Egamma.Keys.Output.GSFTrackParticles,
            ExtrapolationTool=emextLRTCache)
        )
        acc.merge(egammaRecBuilderCfg(
            flags,
            name='LRTegammaRecBuilder',
            egammaRecContainer="LRT"+flags.Egamma.Keys.Internal.EgammaRecs,
            TrackMatchBuilderTool=lrtemtrackmatch,
            doConversions=False)
        )

        from egammaAlgs.egammaSuperClusterBuilderConfig import (
            electronSuperClusterBuilderCfg)
        acc.merge(electronSuperClusterBuilderCfg(
            flags,
            name='LRTelectronSuperClusterBuilder',
            InputEgammaRecContainerName='LRT'+flags.Egamma.Keys.Internal.EgammaRecs,
            SuperElectronRecCollectionName='LRT' +
            flags.Egamma.Keys.Internal.ElectronSuperRecs,
            SuperClusterCollectionName='LRTElectronSuperClusters',
            TrackMatchBuilderTool=lrtemtrackmatch)
        )

        from egammaAlgs.topoEgammaBuilderConfig import (
            topoEgammaBuilderCfg)
        from egammaTools.EMClusterToolConfig import (
            EMClusterToolCfg)
        LRTEMClusterTool = acc.popToolsAndMerge(EMClusterToolCfg(
            flags,
            name='LRTEMClusterTool',
            OutputClusterContainerName='LRT'+flags.Egamma.Keys.Output.CaloClusters)
        )
        acc.merge(topoEgammaBuilderCfg(
            flags,
            name='LRTtopoEgammaBuilder',
            InputElectronRecCollectionName='LRT' +
            flags.Egamma.Keys.Internal.ElectronSuperRecs,
            ElectronOutputName='LRT'+flags.Egamma.Keys.Output.Electrons,
            EMClusterTool=LRTEMClusterTool,
            doPhotons=False)
        )

    # Add truth association
    if flags.Egamma.doTruthAssociation:

        from egammaAlgs.egammaTruthAssociationConfig import (
            egammaTruthAssociationCfg)
        acc.merge(egammaTruthAssociationCfg(
            flags,
            name='LRTegammaTruthAssociationAlg',
            ElectronContainerName='LRT'+flags.Egamma.Keys.Output.Electrons,
            EgammaTruthContainerName='LRT'+flags.Egamma.Keys.Output.TruthParticles,
            MatchPhotons=False,
            MatchForwardElectrons=False)
        )

    mlog.info("EGamma LRT reconstruction configured")

    return acc


if __name__ == "__main__":
    from AthenaCommon.Configurable import Configurable
    Configurable.configurableRun3Behavior = True
    from AthenaConfiguration.AllConfigFlags import ConfigFlags as flags
    from AthenaConfiguration.TestDefaults import defaultTestFiles
    from AthenaConfiguration.MainServicesConfig import MainServicesCfg
    flags.Input.Files = defaultTestFiles.RDO
    flags.Output.doWriteAOD = True  # To test the AOD parts

    acc = MainServicesCfg(flags)
    acc.merge(EGammaLRTReconstructionCfg(flags))
    acc.printConfig(withDetails=True,
                    printDefaults=True)

    with open("egammalrtbuilderconfig.pkl", "wb") as f:
        acc.store(f)
