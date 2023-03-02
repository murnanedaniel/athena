# Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
# Define method to construct configures Sec Vtx Finder alg
# attempted by N Ribaric (@LancasterUNI) neza.ribaric@cern.ch


from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory

def InDetSecVtxFinderAlgCfg(flags, name="InDetSecVtxFinderAlg", **kwargs):
  """ return a CA with configured SecVtxFinderAlg """
  acc = ComponentAccumulator()

  from InDetIncSecVxFinderTool.InDetIncSecVxFinderToolConfig import InDetIterativeSecVtxFinderToolCfg
  kwargs.setdefault("InclusiveVertexFinderTool", acc.getPrimaryAndMerge(InDetIterativeSecVtxFinderToolCfg(flags)))

  from InDetAdaptiveMultiSecVtxFinderTool.InDetAdaptiveMultiSecVtxFinderToolConfig import InDetAdaptiveMultiSecVtxFinderToolCfg
  kwargs.setdefault("AdaptiveMultiVertexFinderTool", acc.getPrimaryAndMerge(InDetAdaptiveMultiSecVtxFinderToolCfg(flags)))

  from TrkConfig.TrkVertexToolsConfig import SecVertexMergingToolCfg
  kwargs.setdefault("VertexMergingTool", acc.getPrimaryAndMerge(SecVertexMergingToolCfg(flags,
                                         MininumDistance = 5.0,
                                         CompatibilityDimension = 2)))

  acc.addEventAlgo(CompFactory.InDet.InDetSecVtxFinder(name, **kwargs))
  return acc

