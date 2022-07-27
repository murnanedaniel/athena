# Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration

from TriggerMenuMT.HLT.Config.MenuComponents import RecoFragmentsPool, MenuSequence
from AthenaCommon.CFElements import seqAND
from AthenaConfiguration.ComponentFactory import CompFactory
from TrigEDMConfig.TriggerEDMRun3 import recordable
from TrigGenericAlgs.TrigGenericAlgsConfig import TrigEventInfoRecorderAlgCfg

#from GaudiKernel.Constants import (VERBOSE, DEBUG, INFO, WARNING, ERROR, FATAL)



def TLAJetSequence (flags, jetsIn):
    
    ## add the InputMaker (event context)    
    tlaJetInputMakerAlg = CompFactory.InputMakerForRoI("IMTLAJets_"+jetsIn)#,RoIsLink="initialRoI")
    tlaJetInputMakerAlg.RoITool = CompFactory.ViewCreatorPreviousROITool()
    tlaJetInputMakerAlg.mergeUsingFeature = True
    
    # configure an instance of TrigEventInfoRecorderAlg
    
    eventInfoRecorderAlgCfg = TrigEventInfoRecorderAlgCfg("TrigEventInfoRecorderAlg_TLA")
    eventInfoRecorderAlgCfg.decorateTLA = True
    eventInfoRecorderAlgCfg.trigEventInfoKey = recordable("HLT_TCEventInfo_TLA")
    eventInfoRecorderAlgCfg.primaryVertexInputName = "HLT_IDVertex_FS"



    tlaJetAthSequence = seqAND( "TLAJetAthSequence_"+jetsIn, [tlaJetInputMakerAlg, eventInfoRecorderAlgCfg] )
    jetsOut = recordable(jetsIn+"_TLA")
    return (tlaJetAthSequence, tlaJetInputMakerAlg, jetsOut)


def TLAJetMenuSequence( flags, jetsIn, attachBtag=True ):
    
    # retrieves the sequence
    (tlaJetAthSequence, tlaJetInputMakerAlg, jetsOut) = RecoFragmentsPool.retrieve(TLAJetSequence, flags, jetsIn=jetsIn)  
    #  add the hypo
    from TrigHLTJetHypo.TrigHLTJetHypoConf import TrigJetTLAHypoAlg
    from TrigHLTJetHypo.TrigJetHypoToolConfig import trigJetTLAHypoToolFromDict

    hypo = TrigJetTLAHypoAlg("TrigJetTLAHypoAlg_"+jetsIn) 
    btagJetTool = CompFactory.TrigBtagTLATool("BtagTLATool_"+jetsIn)
    hypo.AttachBtag = attachBtag  # avoid attaching btag if creating EMTopo Jets with no tracking
    if hypo.AttachBtag:
        btagJetTool.TLAOutputBTaggingCollection = recordable(jetsOut+"_BTagging")
    hypo.BtagJetTool = btagJetTool

    hypo.TLAOutputName = jetsOut

    return MenuSequence( Sequence    = tlaJetAthSequence,
                         Maker       = tlaJetInputMakerAlg,
                         Hypo        = hypo,
                         HypoToolGen = trigJetTLAHypoToolFromDict
                         )

