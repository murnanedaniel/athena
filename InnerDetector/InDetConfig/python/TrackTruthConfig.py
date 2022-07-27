# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory

# -------------------------------------------------------------------------
#
# ------- fragment to handle track truth association
#
# -------------------------------------------------------------------------

def InDetTruthMatchToolCfg(flags, name='InDetTruthMatchTool', **kwargs) :
    acc = ComponentAccumulator()

    kwargs.setdefault("WeightPixel", 10.)
    kwargs.setdefault("WeightSCT", 5.)
    kwargs.setdefault("WeightTRT", 1.)

    InDetTruthMatchSimilarityTool = CompFactory.Trk.TruthMatchRatio(name = name, **kwargs)
    acc.setPrivateTools(InDetTruthMatchSimilarityTool)
    return acc

def TrackTruthSimilaritySelectorCfg(flags, DetailedTruth, TracksTruth, name='Selector', **kwargs) :
    acc = ComponentAccumulator()

    InDetTruthMatchSimilarityTool = acc.popToolsAndMerge(InDetTruthMatchToolCfg(flags))

    kwargs.setdefault("DetailedTrackTruthName", DetailedTruth)
    kwargs.setdefault("OutputName", TracksTruth)
    kwargs.setdefault("TrackTruthSimilarityTool", InDetTruthMatchSimilarityTool)

    acc.addEventAlgo(CompFactory.TrackTruthSimilaritySelector(name = TracksTruth+name, **kwargs))
    return acc

def InDetTrackTruthCfg(flags, Tracks="CombinedInDetTracks", DetailedTruth="CombinedInDetTracksDetailedTruth", TracksTruth="CombinedInDetTracksTruthCollection"):
    acc = ComponentAccumulator()
    #
    # --- Enable the detailed track truth
    #
    from InDetConfig.InDetTruthAlgsConfig import InDetDetailedTrackTruthMakerCfg
    acc.merge(InDetDetailedTrackTruthMakerCfg(flags, Tracks, DetailedTruth))
    #
    # --- Detailed to old TrackTruth
    #
    acc.merge(TrackTruthSimilaritySelectorCfg(flags, DetailedTruth, TracksTruth))

    return acc


if __name__ == "__main__":
    from AthenaConfiguration.AllConfigFlags import ConfigFlags

    numThreads=1
    ConfigFlags.Concurrency.NumThreads=numThreads
    ConfigFlags.Concurrency.NumConcurrentEvents=numThreads # Might change this later, but good enough for the moment.

    ConfigFlags.Detector.GeometryPixel = True
    ConfigFlags.Detector.GeometrySCT = True
    ConfigFlags.Detector.GeometryTRT   = True

    ConfigFlags.InDet.Tracking.doPixelClusterSplitting = True

    from AthenaConfiguration.TestDefaults import defaultTestFiles
    ConfigFlags.Input.Files = defaultTestFiles.RDO_RUN2
    ConfigFlags.lock()
    ConfigFlags.dump()

    from AthenaConfiguration.MainServicesConfig import MainServicesCfg
    top_acc = MainServicesCfg(ConfigFlags)

    from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
    top_acc.merge(PoolReadCfg(ConfigFlags))

    ################## SiliconPreProcessing Configurations ###################
    from InDetConfig.SiliconPreProcessing import InDetRecPreProcessingSiliconCfg
    top_acc.merge(InDetRecPreProcessingSiliconCfg(ConfigFlags))
    #################### TRTPreProcessing Configurations #####################
    from InDetConfig.TRTPreProcessing import TRTPreProcessingCfg
    top_acc.merge(TRTPreProcessingCfg(ConfigFlags))
    
    #//// TrackingSiPatternConfig configurations from Temporary location /////
    ################# SiSPSeededTrackFinder Configurations ###################

    InputCollections = []

    SiSPSeededTrackCollectionKey = 'SiSPSeededPixelTracks'
    ResolvedTrackCollectionKey = 'ResolvedPixelTracks'
    from InDetConfig.TrackingSiPatternConfig import SiSPSeededTrackFinderCfg
    top_acc.merge(SiSPSeededTrackFinderCfg( ConfigFlags,
                                            InputCollections = InputCollections, 
                                            SiSPSeededTrackCollectionKey = SiSPSeededTrackCollectionKey))
    ##########################################################################
    #################### InDetTrackTruth Configurations ######################

    InputTrackCollection = 'SiSPSeededPixelTracks'
    InputDetailedTrackTruth = 'DetailedTrackTruth'
    InputTrackCollectionTruth = 'TrackTruthCollection'
    
    top_acc.merge(InDetTrackTruthCfg(flags=ConfigFlags,
                                     Tracks = InputTrackCollection,
                                     DetailedTruth = InputDetailedTrackTruth,
                                     TracksTruth = InputTrackCollectionTruth))
    #################################################################
    top_acc.printConfig()
    top_acc.run(25)
    top_acc.store(open("test_TrackTruthConfig.pkl", "wb"))
