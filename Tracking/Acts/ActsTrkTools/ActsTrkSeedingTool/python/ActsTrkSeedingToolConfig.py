#  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration

from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory
from ActsInterop import UnitConstants
from ActsInterop.ActsConfigFlags import SeedingStrategy 

def ActsTrkITkPixelSeedingToolCfg(flags,
                                  **kwargs) -> ComponentAccumulator:
    acc = ComponentAccumulator()
    ## For ITkPixel
    kwargs.setdefault("numSeedIncrement" , float("inf"))
    kwargs.setdefault("deltaZMax" , float("inf"))
    kwargs.setdefault("maxPtScattering", float("inf"))
    acc.setPrivateTools(CompFactory.ActsTrk.SeedingTool(name = "ActsSeedingTool_ITkPixel", **kwargs))
    return acc

def ActsTrkITkStripSeedingToolCfg(flags,
                                  **kwargs) -> ComponentAccumulator:
    acc = ComponentAccumulator()
    ## For ITkStrip, change properties that have to be modified w.r.t. the default values
    # For SpacePointGridConfig
    kwargs.setdefault("gridRMax" , 1000. * UnitConstants.mm)
    kwargs.setdefault("deltaRMax" , 600. * UnitConstants.mm)
    kwargs.setdefault("impactMax" , 20. * UnitConstants.mm)
    # For SeedfinderConfig
    kwargs.setdefault("rMax" , 1200. * UnitConstants.mm)
    kwargs.setdefault("deltaRMinTopSP" , 20. * UnitConstants.mm)
    kwargs.setdefault("deltaRMaxTopSP" , 300. * UnitConstants.mm)
    kwargs.setdefault("deltaRMinBottomSP" , 20. * UnitConstants.mm)
    kwargs.setdefault("deltaRMaxBottomSP" , 300. * UnitConstants.mm)
    kwargs.setdefault("deltaZMax" , 900. * UnitConstants.mm)
    kwargs.setdefault("interactionPointCut" , False)
    kwargs.setdefault("arithmeticAverageCotTheta" , True)
    kwargs.setdefault("zBinsCustomLooping" , [6, 7, 5, 8, 4, 9, 3, 10, 2, 11, 1])
    kwargs.setdefault("skipPreviousTopSP" , False)
    kwargs.setdefault("deltaRMiddleMinSPRange" , 30 * UnitConstants.mm)
    kwargs.setdefault("deltaRMiddleMaxSPRange" , 150 * UnitConstants.mm)
    kwargs.setdefault("useDetailedDoubleMeasurementInfo" , True)
    kwargs.setdefault("maxPtScattering", float("inf"))
    # For SeedFilterConfig
    kwargs.setdefault("useDeltaRorTopRadius" , False)
    kwargs.setdefault("seedConfirmationInFilter" , False)
    kwargs.setdefault("impactWeightFactor" , 1.)
    kwargs.setdefault("compatSeedLimit" , 4)
    kwargs.setdefault("numSeedIncrement" , 1.)
    kwargs.setdefault("seedWeightIncrement" , 10100.)
    kwargs.setdefault("maxSeedsPerSpMConf" , 100)
    kwargs.setdefault("maxQualitySeedsPerSpMConf" , 100)
    # For seeding algorithm
    kwargs.setdefault("zBinNeighborsBottom" , [(0,1),(0,1),(0,1),(0,2),(0,1),(0,0),(-1,0),(-2,0),(-1,0),(-1,0),(-1,0)])

    acc.setPrivateTools(CompFactory.ActsTrk.SeedingTool(name = "ActsSeedingTool_ITkStrip", **kwargs))
    return acc

def ActsTrkITkPixelOrthogonalSeedingToolCfg(flags,
                                            **kwargs) -> ComponentAccumulator:
    acc = ComponentAccumulator()
    ## For ITkPixel, use default values for ActsTrk::OrthogonalSeedingTool 
    acc.setPrivateTools(CompFactory.ActsTrk.OrthogonalSeedingTool(name = "OrthogonalSeedingTool_ITkPixel", **kwargs))
    return acc

def ActsTrkITkStripOrthogonalSeedingToolCfg(flags,
                                            **kwargs) -> ComponentAccumulator:
    acc = ComponentAccumulator()
    ## For ITkStrip, change properties that have to be modified w.r.t. the default values 
    kwargs.setdefault("impactMax" , 20. * UnitConstants.mm)
    kwargs.setdefault('rMax', 1200. * UnitConstants.mm)
    kwargs.setdefault("deltaRMinTopSP" , 20. * UnitConstants.mm)
    kwargs.setdefault("deltaRMaxTopSP" , 300. * UnitConstants.mm)
    kwargs.setdefault("deltaRMinBottomSP" , 20. * UnitConstants.mm)
    kwargs.setdefault("deltaRMaxBottomSP" , 300. * UnitConstants.mm)
    kwargs.setdefault("deltaZMax" , 900. * UnitConstants.mm)
    kwargs.setdefault("interactionPointCut" , False)
    kwargs.setdefault("skipPreviousTopSP", False)
    kwargs.setdefault("impactWeightFactor" , 1.)
    kwargs.setdefault("compatSeedLimit" , 4)
    kwargs.setdefault("seedWeightIncrement" , 10100.)
    kwargs.setdefault("numSeedIncrement" , 1.)
    kwargs.setdefault("seedConfirmationInFilter" , False)
    kwargs.setdefault("maxSeedsPerSpMConf" , 100)
    kwargs.setdefault("maxQualitySeedsPerSpMConf" , 100)
    kwargs.setdefault("useDeltaRorTopRadius" , False)
    kwargs.setdefault("rMinMiddle", 33. * UnitConstants.mm)
    kwargs.setdefault("rMaxMiddle", 1200. * UnitConstants.mm)

    acc.setPrivateTools(CompFactory.ActsTrk.OrthogonalSeedingTool(name = "OrthogonalSeedingTool_ITkStrip", **kwargs))
    return acc

def  ActsTrkSiSpacePointsSeedMakerCfg(flags,
                                      name: str = 'ActsTrkSiSpacePointsSeedMaker',
                                      **kwargs) -> ComponentAccumulator:
    assert isinstance(name, str)

    acc = ComponentAccumulator()

    kwargs['name'] = name

    # Main properties
    from ActsInterop.TrackingComponentConfigurer import TrackingComponentConfigurer
    configuration_settings = TrackingComponentConfigurer(flags) 

    kwargs.setdefault('SpacePointsPixelName', 'ITkPixelSpacePoints')
    kwargs.setdefault('SpacePointsStripName', 'ITkStripSpacePoints')
    kwargs.setdefault('SpacePointsOverlapName', 'ITkOverlapSpacePoints')
    kwargs.setdefault('usePixel', flags.ITk.Tracking.ActiveConfig.useITkPixel and flags.ITk.Tracking.ActiveConfig.useITkPixelSeeding)
    kwargs.setdefault('useStrip', flags.ITk.Tracking.ActiveConfig.useITkStrip and flags.ITk.Tracking.ActiveConfig.useITkStripSeeding)
    kwargs.setdefault('useOverlapSpCollection', flags.ITk.Tracking.ActiveConfig.useITkStrip and flags.ITk.Tracking.ActiveConfig.useITkStripSeeding)
    kwargs.setdefault('doSpacePointConversion', not (configuration_settings.doActsSpacePoint and configuration_settings.doAthenaToActsCluster))
    kwargs.setdefault('ActsTrkSpacePointsPixelName'    , "ITkPixelSpacePoints")
    kwargs.setdefault('ActsTrkSpacePointsStripName'    , "ITkStripSpacePoints")
    kwargs.setdefault('ActsTrkSpacePointsOverlapName'  , "ITkStripOverlapSpacePoints")
    kwargs.setdefault('PixelClusterContainerKey', "ITkPixelClusters")
    kwargs.setdefault('StripClusterContainerKey', "ITkStripClusters")

    if flags.ITk.Tracking.ActiveConfig.usePrdAssociationTool:
        # not all classes have that property !!!
        kwargs.setdefault('PRDtoTrackMap', 'ITkPRDtoTrackMap'+ flags.ITk.Tracking.ActiveConfig.extension)

    # Acts Seed Tools
    # Do not overwrite if already present in `kwargs`
    seedTool_pixel = None
    if 'SeedToolPixel' not in kwargs:
        if flags.Acts.SeedingStrategy is SeedingStrategy.Orthogonal:
            seedTool_pixel = acc.popToolsAndMerge(ActsTrkITkPixelOrthogonalSeedingToolCfg(flags))
        else:
            seedTool_pixel = acc.popToolsAndMerge(ActsTrkITkPixelSeedingToolCfg(flags))

    seedTool_strip = None
    if 'SeedToolStrip' not in kwargs:
        if flags.Acts.SeedingStrategy is SeedingStrategy.Orthogonal:
            seedTool_strip = acc.popToolsAndMerge(ActsTrkITkStripOrthogonalSeedingToolCfg(flags))
        else:
            seedTool_strip = acc.popToolsAndMerge(ActsTrkITkStripSeedingToolCfg(flags))

    kwargs.setdefault('SeedToolPixel', seedTool_pixel)
    kwargs.setdefault('SeedToolStrip', seedTool_strip)

    # Validation
    if flags.ITk.Tracking.writeSeedValNtuple:
        kwargs.setdefault('WriteNtuple', True)
        HistService = CompFactory.THistSvc(Output = ["valNtuples DATAFILE='SeedMakerValidation.root' OPT='RECREATE'"])
        acc.addService(HistService)

    acc.setPrivateTools(CompFactory.ActsTrk.SiSpacePointsSeedMaker(**kwargs))
    return acc
