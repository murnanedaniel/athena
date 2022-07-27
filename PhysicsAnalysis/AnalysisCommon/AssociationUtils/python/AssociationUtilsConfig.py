# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration                                                                                                                                                                                                            
from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory

def OverlapRemovalToolCfg(ConfigFlags, 
                          masterName='OverlapRemovalTool',
                          inputLabel='selected', outputLabel='overlaps',
                          bJetLabel='', maxElePtForBJetAwareOR = 100. * 1000,
                          boostedLeptons=False,
                          outputPassValue=False,
                          linkOverlapObjects=False,
                          doEleEleOR=False,
                          doElectrons=True, doMuons=True, doJets=True,
                          doTaus=True, doPhotons=True, doFatJets=False, doMuPFJetOR=False,
                          **kwargs):

    """
    Provides the pre-configured overlap removal recommendations.
    All overlap tools will be (private) added to the master tool
    which is then returned by this function.

    Arguments:
      masterName         - set the name of the master tool.
      inputLabel         - set the InputLabel property for all tools.
      outputLabel        - set the OutputLabel property for all tools.
      bJetLabel          - set user bjet decoration name. Leave blank to
                           disable btag-aware overlap removal.
      maxElePtForBJetAwareOR  - set the maximum electron pT for which b-tag
                           aware overlap removal is done. Set to negative
                           value to use for all electrons.
      boostedLeptons     - enable sliding dR cones for boosted lepton
                           analyses.
      outputPassValue    - set the OutputPassValue property for all tools
                           which determines whether passing objects are
                           marked with true or false.
      linkOverlapObjects - enable ElementLinks to overlap objects.
      doEleEleOR         - enable electron-electron overlap removal.
      doMuPFJetOR        - enable the pflow jet removal for muons
      doXXXX             - these flags enable/disable object types to
                           configure tools for: doElectrons, doMuons,
                           doJets, doTaus, doPhotons, doFatJets.
      kwargs             - additional properties to be applied to all tools.
                           For example: OutputLevel.
    """

    # These properties can be applied to all tools
    common_args = {
        'InputLabel' : inputLabel,
        'OutputLabel' : outputLabel,
        'OutputPassValue' : outputPassValue
    }
    # Extend with additional user-defined global properties
    common_args.update(kwargs)

    # Configure the master tool
    orTool = CompFactory.ORUtils.OverlapRemovalTool(masterName, **common_args)

    # Overlap tools share an additional common property for object linking
    common_args['LinkOverlapObjects'] = linkOverlapObjects

    # Muon-PFlow fake-jet
    if doMuPFJetOR:
        orTool.MuPFJetORT = CompFactory.ORUtils.MuPFJetOverlapTool('MuPFJetORT', BJetLabel=bJetLabel, **common_args)

    # Electron-electron
    if doElectrons and doEleEleOR:
        orTool.EleEleORT = CompFactory.ORUtils.EleEleOverlapTool('EleEleORT', **common_args)

    # Electron-muon
    if doElectrons and doMuons:
        orTool.EleMuORT = CompFactory.ORUtils.EleMuSharedTrkOverlapTool('EleMuORT', **common_args)
    # Electron-jet
    if doElectrons and doJets:
        orTool.EleJetORT = CompFactory.ORUtils.EleJetOverlapTool('EleJetORT',
                                                                 BJetLabel=bJetLabel,
                                                                 MaxElePtForBJetAwareOR=maxElePtForBJetAwareOR,
                                                                 UseSlidingDR=boostedLeptons,
                                                                 **common_args)
    # Muon-jet
    if doMuons and doJets:
        orTool.MuJetORT = CompFactory.ORUtils.MuJetOverlapTool('MuJetORT',
                                                               BJetLabel=bJetLabel,
                                                               UseSlidingDR=boostedLeptons,
                                                               **common_args)

    # Tau-electron
    if doTaus and doElectrons:
        orTool.TauEleORT = CompFactory.ORUtils.DeltaROverlapTool('TauEleORT', DR=0.2, **common_args)
    # Tau-muon
    if doTaus and doMuons:
        orTool.TauMuORT = CompFactory.ORUtils.DeltaROverlapTool('TauMuORT', DR=0.2, **common_args)
    # Tau-jet
    if doTaus and doJets:
        orTool.TauJetORT = CompFactory.ORUtils.DeltaROverlapTool('TauJetORT', DR=0.2, **common_args)

    # Photon-electron
    if doPhotons and doElectrons:
        orTool.PhoEleORT = CompFactory.ORUtils.DeltaROverlapTool('PhoEleORT', **common_args)
    # Photon-muon
    if doPhotons and doMuons:
        orTool.PhoMuORT = CompFactory.ORUtils.DeltaROverlapTool('PhoMuORT', **common_args)
    # Photon-jet
    if doPhotons and doJets:
        orTool.PhoJetORT = CompFactory.ORUtils.DeltaROverlapTool('PhoJetORT', **common_args)

    # Electron-fatjet
    if doElectrons and doFatJets:
        orTool.EleFatJetORT = CompFactory.ORUtils.DeltaROverlapTool('EleFatJetORT', DR=1.0, **common_args)
    # Jet-fatjet
    if doJets and doFatJets:
        orTool.JetFatJetORT = CompFactory.ORUtils.DeltaROverlapTool('JetFatJetORT', DR=1.0, **common_args)

    acc = ComponentAccumulator()
    acc.setPrivateTools(orTool)
    return acc


