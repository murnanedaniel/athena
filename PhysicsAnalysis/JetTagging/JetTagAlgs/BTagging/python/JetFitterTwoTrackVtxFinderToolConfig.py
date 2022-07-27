# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory

from BTagging.JetFitterSequentialVertexFitterConfig import JetFitterSequentialVertexFitterCfg

def InDetJetFitterTwoTrackVtxFinderToolCfg(name, suffix="", useBTagFlagsDefaults = True, **options):

    """Sets up a InDetJetFitterTwoTrackVtxFinderTool tool and returns it.

    The following options have BTaggingFlags defaults:
    
    ID_maxR                             default: 1150.
    ID_maxZ                             default: 2727.
    twoVertexProbabilityCut             default: 3.4e-2
    revertFromPositiveToNegativeTags    default: False

    input:             name: The name of the tool (should be unique).
      useBTagFlagsDefaults : Whether to use BTaggingFlags defaults for options that are not specified.
                  **options: Python dictionary with options for the tool.
    output: The actual tool, which can then by added to ToolSvc via ToolSvc += output."""
    acc = ComponentAccumulator()
    if useBTagFlagsDefaults:
        jetFitterSequentialVertexFitter = acc.popToolsAndMerge(JetFitterSequentialVertexFitterCfg('JFSeqVxFitter'+suffix))
        defaults = { 'ID_maxR' : 1150.,
                     'ID_maxZ' : 2727.,
                     'twoVertexProbabilityCut' : 0.034,
                     'revertFromPositiveToNegativeTags' : True if (suffix=="FLIP_SIGN") else False,
#                     'CrossDistancesSeedFinder' : ,
                     'SequentialVertexFitter' : jetFitterSequentialVertexFitter }
        for option in defaults:
            options.setdefault(option, defaults[option])

    options['name'] = name
    acc.setPrivateTools( CompFactory.InDet.JetFitterTwoTrackVtxFinderTool(**options) )
    return acc




