# Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration

from TrigUpgradeTest.ElectronMenuConfig import l2CaloRecoCfg, l2CaloHypoCfg
from TriggerMenuMT.HLTMenuConfig.Menu.MenuComponents import MenuSequence, \
    ChainStep, Chain, getChainStepName

from TrigEgammaHypo.TrigL2CaloHypoTool import TrigL2CaloHypoToolFromDict
from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaCommon.CFElements import createStepView


# TODO remove once full tracking is in place
def fakeHypoAlgCfg(flags, name="FakeHypoForElectron"):
    from TrigUpgradeTest.TrigUpgradeTestConf import HLTTest__TestHypoAlg
    return HLTTest__TestHypoAlg( name, Input="" )


def generateChains( flags,  chainDict ):
    import pprint
    pprint.pprint( chainDict )

    firstStepName = getChainStepName('Electron', 1)
    stepReco, stepView = createStepView(firstStepName)

    accCalo = ComponentAccumulator()
    accCalo.addSequence(stepView)

    l2CaloReco = l2CaloRecoCfg(flags)
    accCalo.merge(l2CaloReco, sequenceName=stepReco.getName())

    l2CaloHypo =  l2CaloHypoCfg( flags, name = 'L2ElectronCaloHypo',
                                 CaloClusters = 'L2CaloEMClusters')
    l2CaloHypo.HypoTools=[ TrigL2CaloHypoToolFromDict( chainDict ) ]

    accCalo.addEventAlgo(l2CaloHypo, sequenceName=stepView.getName())

    fastCaloSequence = MenuSequence( Sequence    = l2CaloReco.sequence(),
                                     Maker       = l2CaloReco.inputMaker(),
                                     Hypo        = l2CaloHypo,
                                     HypoToolGen = None, 
                                     CA = accCalo)

    accCalo.printConfig()

    fastCaloStep = ChainStep(firstStepName, [fastCaloSequence])


    secondStepName = getChainStepName('Electron', 2)
    stepReco, stepView = createStepView(secondStepName)

    accTrk = ComponentAccumulator()
    accTrk.addSequence(stepView)

    # # # fast ID
    from TrigUpgradeTest.InDetConfig import indetInViewRecoCfg
    fastInDetReco = indetInViewRecoCfg(flags, viewMakerName='ElectronInDet')
    accTrk.merge( fastInDetReco, sequenceName=stepReco.getName() )
    # TODO once tracking fully works remove fake hypos

    from TrigUpgradeTest.TrigUpgradeTestConf import HLTTest__TestHypoTool

    # TODO remove once full tracking is in place
    fakeHypoAlg = fakeHypoAlgCfg(flags, name='FakeHypoForElectron')


    def makeFakeHypoTool(name, cfg):
        return HLTTest__TestHypoTool(name)

    fakeHypoAlg.HypoTools = [ makeFakeHypoTool(chainDict['chainName'], None) ]

    accTrk.addEventAlgo(fakeHypoAlg, sequenceName=stepView.getName())

    fastInDetSequence = MenuSequence( Sequence    = fastInDetReco.sequence(),
                                      Maker       = fastInDetReco.inputMaker(),
                                      Hypo        = fakeHypoAlg,
                                      HypoToolGen = None,
                                      CA = accTrk)

    fastInDetStep = ChainStep( secondStepName, [fastInDetSequence] )
    
    # # # EF calo

    # # # EF ID
    
    # # # offline egamma

    chain = Chain( chainDict['chainName'], chainDict['L1item'], [fastCaloStep, fastInDetStep] )
    return chain
