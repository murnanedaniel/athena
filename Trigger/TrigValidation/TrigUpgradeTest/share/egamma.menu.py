#
#  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
#

include("TrigUpgradeTest/testHLT_MT.py")

#workaround to prevent online trigger folders to be enabled
from InDetTrigRecExample.InDetTrigFlags import InDetTrigFlags
InDetTrigFlags.useConditionsClasses.set_Value_and_Lock(False)


from InDetRecExample.InDetJobProperties import InDetFlags
InDetFlags.doCaloSeededBrem = False


from InDetRecExample.InDetJobProperties import InDetFlags
InDetFlags.InDet25nsec = True 
InDetFlags.doPrimaryVertex3DFinding = False 
InDetFlags.doPrintConfigurables = False
InDetFlags.doResolveBackTracks = True 
InDetFlags.doSiSPSeededTrackFinder = True
InDetFlags.doTRTPhaseCalculation = True
InDetFlags.doTRTSeededTrackFinder = True
InDetFlags.doTruth = False
InDetFlags.init()

# PixelLorentzAngleSvc and SCTLorentzAngleSvc
include("InDetRecExample/InDetRecConditionsAccess.py")

from AthenaCommon.AlgSequence import AlgSequence
topSequence = AlgSequence()


from InDetRecExample.InDetKeys import InDetKeys

# provide a minimal menu information
if globalflags.InputFormat.is_bytestream():
   topSequence.L1DecoderTest.ctpUnpacker.OutputLevel=DEBUG
   topSequence.L1DecoderTest.roiUnpackers[0].OutputLevel=DEBUG
## testChains = ["HLT_e3_etcut", "HLT_e5_etcut", "HLT_e7_etcut", "HLT_2e3_etcut", "HLT_e3e5_etcut"]
## topSequence.L1DecoderTest.prescaler.Prescales = ["HLT_e3_etcut:2", "HLT_2e3_etcut:2.5"]

from TrigUpgradeTest.HLTCFConfig import decisionTree_From_Chains
from TrigUpgradeTest.MenuComponents import NodeSequence, MenuSequence, Chain, ChainStep2


from TrigT2CaloEgamma.TrigT2CaloEgammaConfig import T2CaloEgamma_FastAlgo
theFastCaloAlgo=T2CaloEgamma_FastAlgo("FastCaloAlgo" )
theFastCaloAlgo.OutputLevel=VERBOSE
theFastCaloAlgo.ClustersName="L2CaloClusters"
svcMgr.ToolSvc.TrigDataAccess.ApplyOffsetCorrection=False

 
from TrigMultiVarHypo.TrigL2CaloRingerFexMTInit import init_ringer
trigL2CaloRingerFexMT = init_ringer()
trigL2CaloRingerFexMT.OutputLevel = DEBUG    



from AthenaCommon.CFElements import parOR, seqOR, seqAND, stepSeq

from DecisionHandling.DecisionHandlingConf import RoRSeqFilter, DumpDecisions

from ViewAlgs.ViewAlgsConf import TestEventViewCreatorAlgorithm

fastCaloInViewAlgs = seqAND("fastCaloInViewAlgs", [theFastCaloAlgo, trigL2CaloRingerFexMT])

## filterL1RoIsAlg = RoRSeqFilter("filterL1RoIsAlg")
## filterL1RoIsAlg.Input = ["EMRoIDecisions"]
## filterL1RoIsAlg.Output = ["FilteredEMRoIDecisions"]
## filterL1RoIsAlg.Chains = testChains
## filterL1RoIsAlg.OutputLevel = DEBUG


fastCaloViewsMaker = TestEventViewCreatorAlgorithm("fastCaloViewsMaker", OutputLevel=DEBUG)
fastCaloViewsMaker.ViewFallThrough = True
#fastCaloViewsMaker.Decisions = "FilteredEMRoIDecisions" # from EMRoIsUnpackingTool
fastCaloViewsMaker.RoIsLink = "initialRoI" # -||-
fastCaloViewsMaker.InViewRoIs = "EMCaloRoIs" # contract with the fastCalo
fastCaloViewsMaker.Views = "EMCaloViews"
#fastCaloViewsMaker.AlgorithmNameSequence = [theFastCaloAlgo.getName(), trigL2CaloRingerFexMT.getName()]
fastCaloViewsMaker.ViewNodeName = "fastCaloInViewAlgs"
theFastCaloAlgo.RoIs = fastCaloViewsMaker.InViewRoIs

CaloViewVerify = CfgMgr.AthViews__ViewDataVerifier("FastCaloViewDataVerifier")
CaloViewVerify.DataObjects = [('TrigRoiDescriptorCollection' , 'StoreGateSvc+fastCaloViewsMaker_InViewRoIs_out')]



#from TrigEgammaHypo.TrigEgammaHypoConf import TrigL2CaloHypoAlg
#from TrigEgammaHypo.TrigL2CaloHypoTool import TrigL2CaloHypoToolFromName
from TrigEgammaHypo.TrigEgammaHypoConf import TestTrigL2CaloHypoAlg
theFastCaloHypo = TestTrigL2CaloHypoAlg("L2CaloHypo")
theFastCaloHypo.OutputLevel = DEBUG
theFastCaloHypo.RunInView=True
## theFastCaloHypo.L1Decisions = "EMRoIDecisions"
## theFastCaloHypo.Views = fastCaloViewsMaker.Views
theFastCaloHypo.CaloClusters = theFastCaloAlgo.ClustersName
#theFastCaloHypo.RoIs = fastCaloViewsMaker.InViewRoIs
#theFastCaloHypo.Decisions = "EgammaCaloDecisions"
#theFastCaloHypo.HypoTools =  [ TrigL2CaloHypoToolFromName( c ) for c in testChains ]
#[ TrigL2CaloHypoToolFromName("HLT_e5_etcut"),   TrigL2CaloHypoToolFromName("HLT_e7_etcut") , TrigL2CaloHypoToolFromName("HLT_2e3_etcut"), TrigL2CaloHypoToolFromName("HLT_e3e5_etcut") ]

## for t in theFastCaloHypo.HypoTools:
##   t.OutputLevel = DEBUG
  
# topSequence += theFastCaloHypo

#caloDecisionsDumper = DumpDecisions("caloDecisionsDumper", OutputLevel=DEBUG, Decisions = theFastCaloHypo.Output )  

fastCaloViewSequence = seqAND("fastCaloViewSequence", [fastCaloViewsMaker, fastCaloInViewAlgs ])

fastCaloSequence =  seqAND("fastCaloSequence", [fastCaloViewSequence])


#egammaCaloStep = stepSeq("egammaCaloStep", filterL1RoIsAlg, [ fastCaloSequence,  caloDecisionsDumper ])
#NodeSequence("muNodeSeq1",  Maker=muIM, Sequence=mustep1_sequence, Hypo=muHypo, Seed="L1MU", HypoToolClassName=MuStep1HypoTool())

fastCalo_NodeSequence = NodeSequence("fastCalo_NodeSequence",
                                         Sequence=fastCaloSequence,
                                         Maker=fastCaloViewsMaker,
                                         Hypo=theFastCaloHypo,
                                         HypoToolClassName="TrigL2CaloHypoToolConf",
                                         Seed="L1EM")

fastCaloSequence = MenuSequence("egammaCaloStep", nodeSeqList=[fastCalo_NodeSequence])

#########################################
# second step:  tracking.....
#########################################
#


from TrigUpgradeTest.InDetSetup import makeInDetAlgs

(viewAlgs, eventAlgs) = makeInDetAlgs()
from TrigFastTrackFinder.TrigFastTrackFinder_Config import TrigFastTrackFinder_eGamma

theFTF = TrigFastTrackFinder_eGamma()
theFTF.isRoI_Seeded = True
viewAlgs.append(theFTF)


# A simple algorithm to confirm that data has been inherited from parent view
# Required to satisfy data dependencies
ViewVerify = CfgMgr.AthViews__ViewDataVerifier("electronViewDataVerifier")
ViewVerify.DataObjects = [('xAOD::TrigEMClusterContainer','StoreGateSvc+L2CaloClusters')]
ViewVerify.OutputLevel = DEBUG
viewAlgs.append(ViewVerify)

TrackParticlesName = ""
for viewAlg in viewAlgs:
  if viewAlg.name() == "InDetTrigTrackParticleCreatorAlg":
    TrackParticlesName = viewAlg.TrackParticlesName
    

from TrigEgammaHypo.TrigL2ElectronFexMTConfig import L2ElectronFex_1
theElectronFex= L2ElectronFex_1()
theElectronFex.TrigEMClusterName = theFastCaloAlgo.ClustersName
theElectronFex.TrackParticlesName = TrackParticlesName
theElectronFex.ElectronsName="Electrons"
theElectronFex.OutputLevel=VERBOSE



## filterCaloRoIsAlg = RoRSeqFilter("filterCaloRoIsAlg")
## filterCaloRoIsAlg.Input = [theFastCaloHypo.Decisions]
## filterCaloRoIsAlg.Output = ["Filtered"+theFastCaloHypo.Decisions]
## filterCaloRoIsAlg.Chains = testChains
## filterCaloRoIsAlg.OutputLevel = DEBUG


l2ElectronViewsMaker = TestEventViewCreatorAlgorithm("l2ElectronViewsMaker", OutputLevel=DEBUG)
# topSequence += l2ElectronViewsMaker
#l2ElectronViewsMaker.Decisions = filterCaloRoIsAlg.Output[0] # output of L2CaloHypo
l2ElectronViewsMaker.RoIsLink = "roi" # -||-
l2ElectronViewsMaker.InViewRoIs = "EMIDRoIs" # contract with the fastCalo
l2ElectronViewsMaker.Views = "EMElectronViews"
l2ElectronViewsMaker.ViewFallThrough = True


for viewAlg in viewAlgs:
  if viewAlg.properties().has_key("RoIs"):
    viewAlg.RoIs = l2ElectronViewsMaker.InViewRoIs
  if viewAlg.properties().has_key("roiCollectionName"):
    viewAlg.roiCollectionName = l2ElectronViewsMaker.InViewRoIs
theElectronFex.RoIs = l2ElectronViewsMaker.InViewRoIs    

electronInViewAlgs = parOR("electronInViewAlgs", viewAlgs + [ theElectronFex ])

l2ElectronViewsMaker.ViewNodeName = "electronInViewAlgs"


from TrigEgammaHypo.TrigEgammaHypoConf import TestTrigL2ElectronHypoAlg
#from TrigEgammaHypo.TrigL2ElectronHypoTool import TrigL2ElectronHypoToolFromName
theElectronHypo = TestTrigL2ElectronHypoAlg()
#theElectronHypo.Views = l2ElectronViewsMaker.Views
theElectronHypo.Electrons = theElectronFex.ElectronsName
theElectronHypo.RunInView=True
#theElectronHypo.ClusterDecisions = theFastCaloHypo.Decisions 
#theElectronHypo.ElectronDecisions = "ElectronL2Decisions"
theElectronHypo.OutputLevel = VERBOSE
#theElectronHypo.HypoTools = [ TrigL2ElectronHypoToolFromName( c ) for c in testChains ]

#for t in theElectronHypo.HypoTools:
#  t.OutputLevel = VERBOSE
# topSequence += theElectronHypo
# InDetCacheCreatorTrigViews,

#electronDecisionsDumper = DumpDecisions("electronDecisionsDumper", OutputLevel=DEBUG, Decisions = theElectronHypo.Output )    

electronViewSequence = seqAND("electronViewSequence", eventAlgs + [l2ElectronViewsMaker, electronInViewAlgs ] )
electronSequence = seqAND("electronSequence", [electronViewSequence] )

electron_NodeSequence = NodeSequence("electron_NodeSequence",
                                        Maker=l2ElectronViewsMaker,                                        
                                        Sequence=electronSequence,
                                        Hypo=theElectronHypo,
                                        HypoToolClassName="TrigL2ElectronHypoToolConf",
                                        Seed="L1EM")


#egammaIDStep = stepSeq("egammaIDStep", filterCaloRoIsAlg, [ electronSequence,  electronDecisionsDumper ] )
electronSequence = MenuSequence("electronStep", nodeSeqList=[electron_NodeSequence])

##########################################
# menu
##########################################
#testChains = ["HLT_e3_etcut", "HLT_e5_etcut", "HLT_e7_etcut", "HLT_2e3_etcut", "HLT_e3e5_etcut"]

# map L1 decisions for menu
for unpack in topSequence.L1DecoderTest.roiUnpackers:
    if unpack.name() is "EMRoIsUnpackingTool":
        unpack.Decisions="L1EM"
    if unpack.name() is "MURoIsUnpackingTool":
        unpack.Decisions="L1MU"


testChains  = [
   Chain(name='HLT_e3_etcut', Seed="L1_EM3",   \
             ChainSteps=[ ChainStep2("Step1_e3_etcut", [fastCaloSequence]),
                          ChainStep2("Step2_e3_etcut", [electronSequence])]  ),

    ## Chain(name='HLT_e3_etcut', Seed="L1_EM3",   \
    ##         ChainSteps=[ ChainStep("Step1_e3_etcut", [SequenceHypoTool(fastCaloSequence,step1_e3_etcut() )])]  ),
 
    Chain(name='HLT_e5_etcut', Seed="L1_EM3",   \
              ChainSteps=[ChainStep2("Step1_e5_etcut", [fastCaloSequence])]),
    Chain(name='HLT_e7_etcut', Seed="L1_EM3",   \
              ChainSteps=[ChainStep2("Step1_e7_etcut", [fastCaloSequence])]),
    #Chain(name='HLT_2e3_etcut', Seed="L1_EM3",   \
    #          ChainSteps=[ChainStep("Step1_2e3_etcut", [SequenceHypoTool(fastCaloSequence, 2e3_etcut() )])])
    ]
    
topSequence.L1DecoderTest.prescaler.Prescales = ["HLT_e3_etcut:2", "HLT_2e3_etcut:2.5"]
 


##########################################
# CF construction
##########################################
from AthenaCommon.CFElements import parOR, seqOR, seqAND, stepSeq

from DecisionHandling.DecisionHandlingConf import TriggerSummaryAlg
summary = TriggerSummaryAlg( "TriggerSummaryAlg" )
summary.InputDecision = "HLTChains"
summary.FinalDecisions = [ "ElectronL2Decisions", "MuonL2Decisions" ]
from TrigOutputHandling.TrigOutputHandlingConf import HLTEDMCreator
edmCreator = HLTEDMCreator()
edmCreator.TrigCompositeContainer = [ "EgammaCaloDecisions", "ElectronL2Decisions", "MuonL2Decisions", "EMRoIDecisions", "METRoIDecisions", "MURoIDecisions", "HLTChainsResult" ]
summary.OutputTools = [ edmCreator ]
summary.OutputLevel = DEBUG

steps = seqAND("EgammaMenu_HLTSteps"  )
decisionTree_From_Chains(steps, testChains)
steps += summary

from TrigSteerMonitor.TrigSteerMonitorConf import TrigSignatureMoniMT
mon = TrigSignatureMoniMT()
mon.FinalDecisions = [ "ElectronL2Decisions", "MuonL2Decisions", "WhateverElse" ]
from TrigUpgradeTest.TestUtils import MenuTest
mon.ChainsList = [ x.split(":")[1] for x in  MenuTest.CTPToChainMapping ]
mon.OutputLevel = DEBUG

import AthenaPoolCnvSvc.WriteAthenaPool
from OutputStreamAthenaPool.OutputStreamAthenaPool import  createOutputStream
StreamESD=createOutputStream("StreamESD","myESD.pool.root",True)
StreamESD.OutputLevel=VERBOSE
topSequence.remove( StreamESD )
def addTC(name):   
   StreamESD.ItemList += [ "xAOD::TrigCompositeContainer#"+name, "xAOD::TrigCompositeAuxContainer#"+name+"Aux." ]

for tc in edmCreator.TrigCompositeContainer:
   addTC( tc )

addTC("HLTSummary")

print "ESD file content " 
print StreamESD.ItemList


hltTop = seqOR( "hltTop", [ steps, mon, summary, StreamESD ] )
topSequence += hltTop


##########################################  


  
from AthenaCommon.AlgSequence import dumpSequence
dumpSequence(topSequence)

print "Now some debug"
print theElectronFex
#print fastCaloSequence
print theFastCaloHypo
print fastCaloViewsMaker
print l2ElectronViewsMaker
print ViewVerify


