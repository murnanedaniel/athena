#
#  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
#

# menu components   
from TriggerMenuMT.HLT.Config.MenuComponents import MenuSequence, RecoFragmentsPool, algorithmCAToGlobalWrapper
from AthenaCommon.CFElements import seqAND
from ViewAlgs.ViewAlgsConf import EventViewCreatorAlgorithm
from DecisionHandling.DecisionHandlingConf import ViewCreatorPreviousROITool


def tag(ion):
    return 'precision' + ('HI' if ion is True else '') + 'CaloElectron'


def precisionCaloSequence(flags, ion=False, variant=''):
    """ Creates PrecisionCalo sequence """

    from TriggerMenuMT.HLT.Egamma.TrigEgammaKeys import  getTrigEgammaKeys
    TrigEgammaKeys = getTrigEgammaKeys(variant, ion=ion)

    # EV creator
    InViewRoIs="PrecisionCaloRoIs" + variant
    precisionCaloViewsMaker = EventViewCreatorAlgorithm('IM' + tag(ion) + variant)
    precisionCaloViewsMaker.ViewFallThrough = True
    precisionCaloViewsMaker.RoIsLink = "initialRoI" # Merge inputs based on their initial L1 ROI
    roiTool = ViewCreatorPreviousROITool()
    # Note: This step processes Decision Objects which have followed either Electron reco, Photon reco, or both.
    # For Decision Object which have followed both, there is an ambiguity about which ROI should be used in this
    # merged step. In such cases we break the ambiguity by specifying that the Electron ROI is to be used.
    roiTool.RoISGKey = TrigEgammaKeys.fastTrackingRoIContainer
    precisionCaloViewsMaker.RoITool = roiTool
    precisionCaloViewsMaker.InViewRoIs = InViewRoIs
    precisionCaloViewsMaker.Views = tag(ion) + 'Views' + variant
    precisionCaloViewsMaker.RequireParentView = True

    # reco sequence
    from TriggerMenuMT.HLT.Electron.PrecisionCaloRecoSequences import precisionCaloRecoSequence
    (precisionCaloInViewSequence, sequenceOut) = precisionCaloRecoSequence(flags, InViewRoIs, ion, variant)
        
    precisionCaloViewsMaker.ViewNodeName = precisionCaloInViewSequence.name()

    theSequence = seqAND(tag(ion) + 'Sequence' + variant, [])

    if ion is True:
        # add UE subtraction for heavy ion e/gamma triggers
        # NOTE: UE subtraction requires an average pedestal to be calculated
        # using the full event (FS info), and has to be done outside of the
        # event views in this sequence. the egammaFSRecoSequence is thus placed
        # before the precisionCaloInViewSequence.
        from TriggerMenuMT.HLT.Egamma.TrigEgammaFactories import egammaFSCaloRecoSequence
        egammaFSRecoSequence = egammaFSCaloRecoSequence(flags)
        theSequence += egammaFSRecoSequence

    from TrigGenericAlgs.TrigGenericAlgsConfig import ROBPrefetchingAlgCfg_Calo
    robPrefetchAlg = algorithmCAToGlobalWrapper(ROBPrefetchingAlgCfg_Calo, flags, nameSuffix=precisionCaloViewsMaker.name())[0]

    # connect EVC and reco
    theSequence += [precisionCaloViewsMaker, robPrefetchAlg, precisionCaloInViewSequence]
    return (theSequence, precisionCaloViewsMaker, sequenceOut)



def precisionCaloMenuSequence(flags, name, is_probe_leg=False, ion=False, variant=''):
    """ Creates precisionCalo MENU sequence """

    (sequence, precisionCaloViewsMaker, sequenceOut) = RecoFragmentsPool.retrieve(precisionCaloSequence, flags, ion=ion, variant=variant)

    #Hypo
    from TrigEgammaHypo.TrigEgammaHypoConf import TrigEgammaPrecisionCaloHypoAlg
    from TrigEgammaHypo.TrigEgammaPrecisionCaloHypoTool import TrigEgammaPrecisionCaloHypoToolFromDict

    thePrecisionCaloHypo = TrigEgammaPrecisionCaloHypoAlg(name + tag(ion) + 'Hypo'+ variant)
    thePrecisionCaloHypo.CaloClusters = sequenceOut

    return MenuSequence( flags,
                         Sequence    = sequence,
                         Maker       = precisionCaloViewsMaker, 
                         Hypo        = thePrecisionCaloHypo,
                         HypoToolGen = TrigEgammaPrecisionCaloHypoToolFromDict,
                         IsProbe     = is_probe_leg)


def precisionCaloMenuSequence_LRT(flags, name, is_probe_leg=False):
    return precisionCaloMenuSequence(flags, name, is_probe_leg=is_probe_leg, ion=False, variant='_LRT')
