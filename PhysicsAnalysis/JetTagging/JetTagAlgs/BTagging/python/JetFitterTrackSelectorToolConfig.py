# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory

from BTagging.InDetJetFitterUtilsConfig import InDetJetFitterUtilsCfg
from BTagging.InDetImprovedJetFitterTrackSelectorToolConfig import InDetImprovedJetFitterTrackSelectorToolCfg
from TrkConfig.AtlasExtrapolatorConfig import AtlasExtrapolatorCfg

def InDetJetFitterTrackSelectorToolCfg(flags, name, suffix="", useBTagFlagsDefaults = True, **options):
    """Sets up a InDetJetFitterTrackSelectorTool  tool and returns it.

    The following options have BTaggingFlags defaults:
    
    revertFromPositiveToNegativeTags          default: False
    cutCompPrimaryVertexForPosLifetimeTracks  default: 1e-1
    cutCompPrimaryVertexForNegLifetimeTracks  default: 5e-2

    input:             name: The name of the tool (should be unique).
      useBTagFlagsDefaults : Whether to use BTaggingFlags defaults for options that are not specified.
                  **options: Python dictionary with options for the tool.
    output: The actual tool, which can then by added to ToolSvc via ToolSvc += output."""
    acc = ComponentAccumulator()
    if useBTagFlagsDefaults:
        jetFitterExtrapolator= acc.popToolsAndMerge(AtlasExtrapolatorCfg(flags,'JFExtrapolator'+suffix))
        inDetJetFitterUtils = acc.popToolsAndMerge(InDetJetFitterUtilsCfg(flags,'InDetJFUtils'+suffix))
        inDetImprovedJetFitterTrackSelectorTool = acc.popToolsAndMerge(InDetImprovedJetFitterTrackSelectorToolCfg(flags, 'InDetImprovedJFTrackSelTool'+suffix))
        #revert the signed track IP sign if the suffix is 'FLIP_SIGN'. 
        defaults = { 'revertFromPositiveToNegativeTags' : True if (suffix=="FLIP_SIGN") else False,
                     'cutCompPrimaryVertexForPosLifetimeTracks' : 0.1,
                     'cutCompPrimaryVertexForNegLifetimeTracks' : 0.05,
                     'Extrapolator' : jetFitterExtrapolator ,
                     'InDetJetFitterUtils' : inDetJetFitterUtils,
                     'TrackSelector' : inDetImprovedJetFitterTrackSelectorTool }
        for option in defaults:
            options.setdefault(option, defaults[option])
            
    options['name'] = name
    acc.setPrivateTools( CompFactory.InDet.JetFitterTrackSelectorTool(**options) )
    return acc



