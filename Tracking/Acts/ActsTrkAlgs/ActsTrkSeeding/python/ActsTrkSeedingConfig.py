# Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration

from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory
from ActsTrkSeedingTool.ActsTrkSeedingToolConfig import ActsTrkSeedingToolCfg
from ActsTrkTrackParamsEstimationTool.ActsTrkTrackParamsEstimationToolConfig import TrackParamsEstimationToolCfg
from ActsGeometry.ActsGeometryConfig import ActsTrackingGeometryToolCfg
        
# ACTS algorithm using Athena objects upstream
def ActsTrkITkPixelSeedingFromAthenaCfg(flags,
                                        name: str = 'ActsTrkPixelSeedingAlgorithm',
                                        **kwargs):
    acc = ComponentAccumulator()

    # Need To add additional tool(s)
    # Tracking Geometry Tool
    acc_geo, geoTool = ActsTrackingGeometryToolCfg(flags)
    acc.merge(acc_geo)

    # ATLAS Converter Tool
    converterTool = CompFactory.ActsATLASConverterTool(TrackingGeometryTool=geoTool)
    
    # Track Param Estimation Tool
    trackEstimationTool = acc.popToolsAndMerge(TrackParamsEstimationToolCfg(flags))

    seedTool = None
    if "SeedTool" not in kwargs:
        seedTool = acc.popToolsAndMerge(ActsTrkSeedingToolCfg(flags,
                                                              configuration = 'Pixel'))

    kwargs.setdefault('InputSpacePoints', ['ITkPixelSpacePoints'])
    kwargs.setdefault('OutputSeeds', 'ITkPixelSeeds')
    kwargs.setdefault('SeedTool', seedTool)
    kwargs.setdefault('TrackingGeometryTool', geoTool)
    kwargs.setdefault('ATLASConverterTool', converterTool)
    kwargs.setdefault('TrackParamsEstimationTool', trackEstimationTool)
    kwargs.setdefault('OutputEstimatedTrackParameters', 'ITkPixelEstimatedTrackParams')
    kwargs.setdefault('UsePixel', True)
    kwargs.setdefault('DetectorElements', 'ITkPixelDetectorElementCollection')

    acc.addEventAlgo(CompFactory.ActsTrk.SeedingFromAthenaAlgorithm(name, **kwargs))
    return acc


def ActsTrkITkStripSeedingFromAthenaCfg(flags,
                                        name: str = 'ActsTrkStripSeedingAlgorithm',
                                        **kwargs):
    acc = ComponentAccumulator()

    # Need To add additional tool(s)
    # Tracking Geometry Tool
    acc_geo, geoTool = ActsTrackingGeometryToolCfg(flags)
    acc.merge(acc_geo)

    # ATLAS Converter Tool
    converterTool = CompFactory.ActsATLASConverterTool(TrackingGeometryTool=geoTool)
    
    # Track Param Estimation Tool
    trackEstimationTool = acc.popToolsAndMerge(TrackParamsEstimationToolCfg(flags))

    seedTool = None
    if "SeedTool" not in kwargs:
        seedTool = acc.popToolsAndMerge(ActsTrkSeedingToolCfg(flags,
                                                              configuration = 'Strip'))

    kwargs.setdefault('InputSpacePoints', ['ITkStripSpacePoints', 'ITkStripOverlapSpacePoints'])
    kwargs.setdefault('OutputSeeds', 'ITkStripSeeds')
    kwargs.setdefault('SeedTool', seedTool)
    kwargs.setdefault('TrackingGeometryTool', geoTool)
    kwargs.setdefault('ATLASConverterTool', converterTool)
    kwargs.setdefault('TrackParamsEstimationTool', trackEstimationTool)
    kwargs.setdefault('OutputEstimatedTrackParameters', 'ITkStripEstimatedTrackParams')
    kwargs.setdefault('UsePixel', False)
    kwargs.setdefault('DetectorElements', 'ITkStripDetectorElementCollection')

    acc.addEventAlgo(CompFactory.ActsTrk.SeedingFromAthenaAlgorithm(name, **kwargs))
    return acc


def ActsTrkSeedingFromAthenaCfg(flags):
    acc = ComponentAccumulator()
    if flags.Detector.EnableITkPixel:
        acc.merge(ActsTrkITkPixelSeedingFromAthenaCfg(flags))
    if flags.Detector.EnableITkStrip:
        acc.merge(ActsTrkITkStripSeedingFromAthenaCfg(flags))
    return acc


# ACTS only algorithm
def ActsTrkITkPixelSeedingCfg(flags,
                              name: str = 'ActsTrkPixelSeedingAlgorithm',
                              **kwargs):
    acc = ComponentAccumulator()

    # tools
    # check if options have already the seed tool defined
    # and make sure it is not a None
    seedTool = None
    if "SeedTool" not in kwargs:
        seedTool = acc.popToolsAndMerge(ActsTrkSeedingToolCfg(flags,
                                                              configuration = 'Pixel'))
        
    kwargs.setdefault('InputSpacePoints', ['ITkPixelSpacePoints'])
    kwargs.setdefault('OutputSeeds', 'ITkPixelSeeds')
    kwargs.setdefault('SeedTool', seedTool)

    acc.addEventAlgo(CompFactory.ActsTrk.SeedingAlgorithm(name, **kwargs))
    return acc

def ActsTrkITkStripSeedingCfg(flags,
                              name: str = 'ActsTrkStripSeedingAlgorithm',
                              **kwargs):
    acc = ComponentAccumulator()

    # tools
    # check if options have already the seed tool defined
    # and make sure it is not a None
    seedTool = None
    if "SeedTool" not in kwargs:
        seedTool = acc.popToolsAndMerge(ActsTrkSeedingToolCfg(flags,
                                                              configuration = 'Strip'))
        
    kwargs.setdefault('InputSpacePoints', ['ITkStripSpacePoints', 'ITkStripOverlapSpacePoints'])
    kwargs.setdefault('OutputSeeds', 'ITkStripSeeds')
    kwargs.setdefault('SeedTool', seedTool)

    acc.addEventAlgo(CompFactory.ActsTrk.SeedingAlgorithm(name, **kwargs))
    return acc

def ActsTrkSeedingCfg(flags):
    acc = ComponentAccumulator()
    if flags.Detector.EnableITkPixel:
        acc.merge(ActsTrkITkPixelSeedingCfg(flags))
    if flags.Detector.EnableITkStrip:
        acc.merge(ActsTrkITkStripSeedingCfg(flags))
    return acc

