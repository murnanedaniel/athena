# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.Enums import Format

def InDetRecPreProcessingSiliconCfg(flags, **kwargs):
    acc = ComponentAccumulator()
    # ------------------------------------------------------------
    #
    # ----------- Data-Preparation stage
    #
    # ------------------------------------------------------------
    #
    # ----------- PrepRawData creation from Raw Data Objects
    #

    #
    # --- Slim BCM RDOs by zero-suppressing
    #   
    if flags.Detector.EnableBCM:
        from InDetConfig.BCM_ZeroSuppressionConfig import BCM_ZeroSuppressionCfg
        acc.merge(BCM_ZeroSuppressionCfg(flags))
    
    #
    # -- Pixel Clusterization
    #
    if flags.Detector.EnablePixel and (   flags.Input.Format is Format.BS
                                       or 'PixelRDOs' in flags.Input.Collections):

        #
        # --- PixelClusterization algorithm
        #
        from InDetConfig.InDetPrepRawDataFormationConfig import PixelClusterizationCfg
        acc.merge(PixelClusterizationCfg(flags))
        if flags.InDet.doSplitReco :
            from InDetConfig.InDetPrepRawDataFormationConfig import PixelClusterizationPUCfg
            acc.merge(PixelClusterizationPUCfg(flags))
    #
    # --- SCT Clusterization
    #
    if flags.Detector.EnableSCT and (   flags.Input.Format is Format.BS
                                     or 'SCT_RDOs' in flags.Input.Collections):

        #
        # --- SCT_Clusterization algorithm
        #
        from InDetConfig.InDetPrepRawDataFormationConfig import SCTClusterizationCfg
        acc.merge(SCTClusterizationCfg(flags))
        if flags.InDet.doSplitReco :
            from InDetConfig.InDetPrepRawDataFormationConfig import SCTClusterizationPUCfg
            acc.merge(SCTClusterizationPUCfg(flags))

    #
    # ----------- form SpacePoints from clusters in SCT and Pixels
    #
    #
    from InDetConfig.SiSpacePointFormationConfig import InDetSiTrackerSpacePointFinderCfg
    acc.merge(InDetSiTrackerSpacePointFinderCfg(flags))

    # this truth must only be done if you do PRD and SpacePointformation
    # If you only do the latter (== running on ESD) then the needed input (simdata)
    # is not in ESD but the resulting truth (clustertruth) is already there ...
    if (flags.InDet.doTruth
        and (not flags.Detector.EnableSCT   or 'SCT_SDO_Map'  in flags.Input.Collections)
        and (not flags.Detector.EnablePixel or 'PixelSDO_Map' in flags.Input.Collections)) :

        from InDetConfig.InDetTruthAlgsConfig import InDetPRD_MultiTruthMakerSiCfg
        acc.merge(InDetPRD_MultiTruthMakerSiCfg(flags))
        if flags.InDet.doSplitReco:
            from InDetConfig.InDetTruthAlgsConfig import InDetPRD_MultiTruthMakerSiPUCfg
            acc.merge(InDetPRD_MultiTruthMakerSiPUCfg(flags))

    return acc

def ITkRecPreProcessingSiliconCfg(flags, **kwargs):
    acc = ComponentAccumulator()
    # ------------------------------------------------------------
    #
    # ----------- Data-Preparation stage
    #
    # ------------------------------------------------------------
    #
    # ----------- PrepRawData creation from Raw Data Objects
    #

    #
    # --- Slim BCM RDOs by zero-suppressing
    #
    if flags.Detector.EnableBCMPrime:
        from InDetConfig.BCM_ZeroSuppressionConfig import BCM_ZeroSuppressionCfg
        acc.merge(BCM_ZeroSuppressionCfg(flags))

    #
    # --- Deducing flags
    #
    doAthenaClustering = False
    doActsClustering = False
    doAthenaSpacePointFormation = False
    doActsSpacePointFormation = False

    from InDetConfig.ITkConfigFlags import TrackingComponent
    if TrackingComponent.AthenaChain in flags.ITk.Tracking.recoChain:
        doAthenaClustering = True
        doAthenaSpacePointFormation = True
    if TrackingComponent.ActsChain in flags.ITk.Tracking.recoChain:
        doActsClustering = True
        doActsSpacePointFormation = True
    if TrackingComponent.ValidateActsClusters in flags.ITk.Tracking.recoChain:
        doActsClustering = True
        doAthenaSpacePointFormation = True
    if TrackingComponent.ValidateActsSpacePoints in flags.ITk.Tracking.recoChain:
        doAthenaClustering = True
        doActsSpacePointFormation = True
    if TrackingComponent.ValidateActsSeeds in flags.ITk.Tracking.recoChain:
        doAthenaClustering = True
        doActsSpacePointFormation = True

    convertInDetClusters = doAthenaClustering and not doActsClustering and doActsSpacePointFormation
    convertXAODClusters = doActsClustering and not doAthenaClustering and doAthenaSpacePointFormation
    # TO-DO: flags for Space Point convertions [not available right now]

    #
    # -- Configuration check (dependencies) is done by the scheduler
    # -- We may want to put it here in the future as well
    #
         
    #
    # -- Clusterization Algorithms
    #
    if doAthenaClustering:
        from InDetConfig.InDetPrepRawDataFormationConfig import AthenaTrkClusterizationCfg
        acc.merge(AthenaTrkClusterizationCfg(flags))
        
    if doActsClustering:
        from ActsTrkClusterization.ActsTrkClusterizationConfig import ActsTrkClusterizationCfg
        acc.merge(ActsTrkClusterizationCfg(flags))

    #
    # ---  Cluster EDM converters
    #
    if convertInDetClusters:
        if not flags.Detector.EnableITkPixel or not flags.Detector.EnableITkStrip:
            raise RuntimeError("Cluster EDM converter (InDet -> xAOD) must be activated for both Pixel and Strips")
        #
        # --- InDet -> xAOD Cluster EDM converter
        #
        from InDetConfig.InDetPrepRawDataFormationConfig import ITkInDetToXAODClusterConversionCfg
        acc.merge(ITkInDetToXAODClusterConversionCfg(flags))

    if convertXAODClusters:
        if not flags.Detector.EnableITkPixel or not flags.Detector.EnableITkStrip:
            raise RuntimeError("Cluster EDM converter (xAOD -> InDet) must be activated for both Pixel and Strips")
        #
        # --- xAOD -> InDet Cluster EDM converter
        #
        from InDetConfig.InDetPrepRawDataFormationConfig import ITkXAODToInDetClusterConversionCfg
        acc.merge(ITkXAODToInDetClusterConversionCfg(flags))
        

    #
    # ----------- form SpacePoints from clusters in SCT and Pixels
    #
    if doAthenaSpacePointFormation:
        from InDetConfig.SiSpacePointFormationConfig import ITkSiTrackerSpacePointFinderCfg
        acc.merge(ITkSiTrackerSpacePointFinderCfg(flags))

    if doActsSpacePointFormation:
        from TrkConfig.ActsTrkSpacePointFormationConfig import ActsTrkSpacePointFormationCfg
        acc.merge(ActsTrkSpacePointFormationCfg(flags))


    #
    # --- Space Point EDM converters
    #

    # this truth must only be done if you do PRD and SpacePointformation
    # If you only do the latter (== running on ESD) then the needed input (simdata)
    # is not in ESD but the resulting truth (clustertruth) is already there ...
    if flags.ITk.Tracking.doTruth:
        from InDetConfig.InDetTruthAlgsConfig import ITkPRD_MultiTruthMakerSiCfg
        acc.merge(ITkPRD_MultiTruthMakerSiCfg(flags))

    return acc


if __name__ == "__main__":
    from AthenaConfiguration.AllConfigFlags import initConfigFlags
    flags = initConfigFlags()

    from AthenaConfiguration.TestDefaults import defaultTestFiles
    flags.Input.Files = defaultTestFiles.RDO_RUN2

    flags.InDet.Tracking.doPixelClusterSplitting = True

    numThreads=1
    flags.Concurrency.NumThreads=numThreads
    flags.Concurrency.NumConcurrentEvents=numThreads


    flags.lock()
    flags.dump()
    
    from AthenaConfiguration.MainServicesConfig import MainServicesCfg
    top_acc = MainServicesCfg(flags)

    from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
    top_acc.merge(PoolReadCfg(flags))

    top_acc.merge(InDetRecPreProcessingSiliconCfg(flags))

    iovsvc = top_acc.getService('IOVDbSvc')
    iovsvc.OutputLevel=5

    top_acc.printConfig()
    top_acc.run(25)
    top_acc.store(open("test_SiliconPreProcessing.pkl", "wb"))
