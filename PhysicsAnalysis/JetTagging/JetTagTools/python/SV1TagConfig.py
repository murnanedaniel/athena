# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory
from JetTagTools.NewLikelihoodToolConfig import NewLikelihoodToolCfg

def SV1TagCfg( flags, name = 'SV1Tag', scheme = '', useBTagFlagsDefaults = True, **options ):
    """Sets up a SV1Tag tool and returns it.

    The following options have BTaggingFlags defaults:

    Runmodus                            default: BTagging.RunModus
    referenceType                       default: BTagging.ReferenceType
    SVAlgType                           default: "SV1"
    jetCollectionList                   default: BTaggingFlags.Jets
    LikelihoodTool                      default: None
    SecVxFinderName                     default: "SV1"
    UseCHypo                            default: True

    input:             name: The name of the tool (should be unique).
          useBTagFlagsDefaults : Whether to use BTaggingFlags defaults for options that are not specified.
                  **options: Python dictionary with options for the tool.
    output: The actual tool."""    
    acc = ComponentAccumulator()
    options['name'] = name
    if 'Flip' in name:
        options['xAODBaseName'] = 'SV1Flip'
    else:
        options['xAODBaseName'] = 'SV1'
    if useBTagFlagsDefaults:
        likelihood = acc.popToolsAndMerge(NewLikelihoodToolCfg(flags, 'SV1NewLikelihoodTool', 'SV1', scheme))
        defaults = { 'Runmodus'                         : flags.BTagging.RunModus,
                     'referenceType'                    : flags.BTagging.ReferenceType,
                     'jetPtMinRef'                      : flags.BTagging.JetPtMinRef,
                     'SaveProbabilities'                : flags.BTagging.SaveSV1Probabilities,
                     'SVAlgType'                        : 'SV1',
                     'jetCollectionList'                : [], #used only in reference mode
                     'SecVxFinderName'                  : 'SV1',
                     'UseCHypo'                         : True,
                     'LikelihoodTool'                   : likelihood }
        for option in defaults:
            options.setdefault(option, defaults[option])
    acc.setPrivateTools(CompFactory.Analysis.SVTag(**options))

    return acc

