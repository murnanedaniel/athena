# Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
# Define method to construct configures Sec Vtx Finder alg
# attempted by N Ribaric (@LancasterUNI) neza.ribaric@cern.ch


from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory
import AthenaCommon.Constants as Lvl


def InDetIterativeSecVtxFinderToolCfg(flags, name="InDetIterativeSecVtxFinderTool", **kwargs):


  acc = ComponentAccumulator()
  
  from TrkConfig.TrkVertexSeedFinderToolsConfig import IndexedCrossDistancesSeedFinderCfg
  kwargs.setdefault("SeedFinder",acc.popToolsAndMerge(IndexedCrossDistancesSeedFinderCfg(flags,
                                 name                = "IndexedCrossDistancesSeedFinder",
                                 trackdistcutoff     = 0.01 ,
                                 trackdistexppower   = 2,
                                 maximumTracksNoCut  = 30,
                                 maximumDistanceCut  = 7.5,
                                 useweights          = True)))

  
  from TrkConfig.TrkVertexFitterUtilsConfig import KalmanVertexUpdatorCfg
  #kwargs['VertexUpdator'] = acc.getPrimaryAndMerge(KalmanVertexUpdatorCfg(flags))
 
  from TrkConfig.TrkVertexFittersConfig import AdaptiveVertexFitterCfg
  kwargs.setdefault("VertexFitterTool",acc.popToolsAndMerge(AdaptiveVertexFitterCfg(flags,
                                       name                = "AdaptiveVxFitterToolIncSecVtx",
                                       VertexUpdator       = acc.getPrimaryAndMerge(KalmanVertexUpdatorCfg(flags)),
                                       MaxIterations       = 8000,
                                       MaxDistToLinPoint   = 0.2,
                                       InitialError        = 0.2,
                                       DoSmoothing         = True)))

  from InDetConfig.InDetTrackSelectionToolConfig import InDetTrackSelectionTool_TrackTools_Cfg
  kwargs.setdefault("BaseTrackSelector",acc.popToolsAndMerge(InDetTrackSelectionTool_TrackTools_Cfg(flags,
                                        name = "InDetTrackSelectionTool",
                                        CutLevel = "NoCut",
                                        minPt = 1000.,
                                        maxD0 = 500.0,
                                        maxZ0 = 1500.,
                                        maxSigmaD0 = -1.0,
                                        maxSigmaZ0SinTheta = -1.0,
                                        maxChiSqperNdf = 5.0,
                                        maxAbsEta = 2.5,
                                        minNInnermostLayerHits = 0,
                                        minNPixelHits = 0,
                                        maxNPixelHoles = 1,
                                        minNSctHits = 2,
                                        minNTrtHits = 0,
                                        minNSiHits = 0,
                                        maxNSiSharedHits = 6)))

  from InDetConfig.InDetTrackSelectionToolConfig import InDetSecVtxTrackSelectionToolCfg
  kwargs.setdefault("SecVtxTrackSelector",acc.popToolsAndMerge(InDetSecVtxTrackSelectionToolCfg(flags,
                                          name = "SecVtxTrackSelector",
                                          minNPixelHitsAtZeroTRT = 2,
                                          minTotalHits = 0,
                                          minD0 = 0.1)))


  from TrkConfig.TrkVertexFitterUtilsConfig import AtlasImpactPoint3dEstimatorCfg
  kwargs.setdefault("ImpactPoint3dEstimator",acc.popToolsAndMerge(AtlasImpactPoint3dEstimatorCfg(flags)))

  from TrkConfig.TrkVertexFitterUtilsConfig import FullLinearizedTrackFactoryCfg
  kwargs.setdefault("LinearizedTrackFactory",acc.popToolsAndMerge(FullLinearizedTrackFactoryCfg(flags)))

  kwargs.setdefault("significanceCutSeeding",9.)
  kwargs.setdefault("maxCompatibilityCutSeeding",18.)
  kwargs.setdefault("minTrackWeightAtVtx",0.02)
  kwargs.setdefault("maxVertices",20)
  kwargs.setdefault("TrackInnerOuterFraction",0.95)
  kwargs.setdefault("MomentumProjectionOnDirection",-999.9)
  kwargs.setdefault("SeedsMinimumDistance",0.1)

  vtxFlags = flags.Tracking.PriVertex
  print("doMaxTracksCut ",vtxFlags.doMaxTracksCut)
  print("MaxTracks ",vtxFlags.maxTracks)
  kwargs.setdefault("doMaxTracksCut",vtxFlags.doMaxTracksCut)
  kwargs.setdefault("MaxTracks",vtxFlags.maxTracks)

  kwargs["VertexFilterLevel"] = 0
  kwargs.setdefault("OutputLevel",Lvl.VERBOSE)

  acc.setPrivateTools(CompFactory.InDet.InDetIterativeSecVtxFinderTool(name, **kwargs))
  return acc
