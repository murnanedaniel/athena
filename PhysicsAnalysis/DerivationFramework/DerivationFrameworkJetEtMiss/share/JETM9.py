#====================================================================
# JETM9.py 
# reductionConf flag JETM9 in Reco_tf.py   
#====================================================================

from DerivationFrameworkCore.DerivationFrameworkMaster import *
from DerivationFrameworkJetEtMiss.JetCommon import *
from DerivationFrameworkJetEtMiss.ExtendedJetCommon import *
#from DerivationFrameworkJetEtMiss.METCommon import *

#====================================================================
# SKIMMING TOOL 
#====================================================================
from DerivationFrameworkJetEtMiss.TriggerLists import *
triggers = jetTriggers

# NOTE: need to be able to OR isSimulated as an OR with the trigger
orstr =' || '
trigger = '('+orstr.join(triggers)+')'
expression = trigger+' || (EventInfo.eventTypeBitmask==1)'

from DerivationFrameworkTools.DerivationFrameworkToolsConf import DerivationFramework__xAODStringSkimmingTool
JETM9SkimmingTool = DerivationFramework__xAODStringSkimmingTool(name = "JETM9SkimmingTool1",
                                                                    expression = expression)
ToolSvc += JETM9SkimmingTool

#=======================================
# CREATE PRIVATE SEQUENCE
#=======================================

jetm9Seq = CfgMgr.AthSequencer("JETM9Sequence")
DerivationFrameworkJob += jetm9Seq
#jetm9Seq = DerivationFrameworkJob

#====================================================================
# SET UP STREAM   
#====================================================================
streamName = derivationFlags.WriteDAOD_JETM9Stream.StreamName
fileName   = buildFileName( derivationFlags.WriteDAOD_JETM9Stream )
JETM9Stream = MSMgr.NewPoolRootStream( streamName, fileName )
JETM9Stream.AcceptAlgs(["JETM9Kernel"])
# for thinning
from AthenaServices.Configurables import ThinningSvc, createThinningSvc
augStream = MSMgr.GetStream( streamName )
evtStream = augStream.GetEventStream()
svcMgr += createThinningSvc( svcName="JETM9ThinningSvc", outStreams=[evtStream] )

# Truth particle thinning
thinningTools = []
from AthenaCommon.GlobalFlags import globalflags
if globalflags.DataSource()=='geant4':
    from DerivationFrameworkMCTruth.DerivationFrameworkMCTruthConf import DerivationFramework__MenuTruthThinning
    JETM9TruthThinning = DerivationFramework__MenuTruthThinning(name                  = "JETM9TruthThinning",
                                                                ThinningService            = "JETM9ThinningSvc",
                                                                WriteStatus3               = True,
                                                                PreserveAncestors          = True,
                                                                WriteFirstN                = 10)
    # from DerivationFrameworkMCTruth.DerivationFrameworkMCTruthConf import DerivationFramework__GenericTruthThinning
    # JETM9TruthParticleThinning = DerivationFramework__GenericTruthThinning(name                    = "JETM9TruthThinning",
    #                                                                        ThinningService         = "JETM9TruthThinningSvc",
    #                                                                        ParticlesKey            = "TruthParticles",  
    #                                                                        ParticleSelectionString = "")
    ToolSvc += JETM9TruthThinning
    thinningTools.append(JETM9TruthThinning)

#=======================================
# CREATE THE DERIVATION KERNEL ALGORITHM   
#=======================================

from DerivationFrameworkCore.DerivationFrameworkCoreConf import DerivationFramework__DerivationKernel
jetm9Seq += CfgMgr.DerivationFramework__DerivationKernel(	name = "JETM9Kernel", 
                                                                SkimmingTools = [JETM9SkimmingTool],
                                                                ThinningTools = thinningTools)

#====================================================================
# Special jets
#====================================================================

#AntiKt4PV0TrackJets
#addStandardJets("AntiKt", 0.2, "PV0Track", 2000, mods="track_ungroomed", algseq=jetm9Seq, outputGroup="JETM9")
#addStandardJets("AntiKt", 0.4, "PV0Track", 2000, mods="track_ungroomed", algseq=jetm9Seq, outputGroup="JETM9")
#AntiKt4PV0TrackJets
addAntiKt2PV0TrackJets(jetm9Seq, "JETM9")
addAntiKt4PV0TrackJets(jetm9Seq, "JETM9")

OutputJets["JETM9"] = ["AntiKt4EMTopoJets","AntiKt4LCTopoJets"]

if globalflags.DataSource()=='geant4':
#    addStandardJets("AntiKt", 0.4, "Truth", mods="truth_ungroomed", ptmin=5000, algseq=jetm9Seq, outputGroup="JETM9")
#    addStandardJets("AntiKt", 0.4, "TruthWZ", mods="truth_ungroomed", ptmin=5000, algseq=jetm9Seq, outputGroup="JETM9")
     addAntiKt4TruthJets(jetm9Seq, "JETM9")
     addAntiKt4TruthWZJets(jetm9Seq, "JETM9")

#====================================================================
# Jets for R-scan 
#====================================================================
from JetRec.JetRecStandard import jtm
from JetRec.JetRecConf import JetAlgorithm

def addRscanJets(jetalg,radius,inputtype,sequence,outputlist):
    jetname = "{0}{1}{2}Jets".format(jetalg,int(radius*10),inputtype)
    algname = "jetalg"+jetname

    if not hasattr(sequence,algname):
        if inputtype == "Truth":
            addStandardJets(jetalg, radius, "Truth", mods="truth_ungroomed", ptmin=5000, algseq=sequence, outputGroup=outputlist)
        if inputtype == "TruthWZ":
            addStandardJets(jetalg, radius, "TruthWZ", mods="truth_ungroomed", ptmin=5000, algseq=sequence, outputGroup=outputlist)
        elif inputtype == "LCTopo":
            addStandardJets(jetalg, radius, "LCTopo", mods="lctopo_ungroomed",
                            ghostArea=0.01, ptmin=2000, ptminFilter=7000, calibOpt="none", algseq=sequence, outputGroup=outputlist)

for radius in [0.2, 0.3, 0.5, 0.6, 0.7, 0.8]:
    if jetFlags.useTruth:
        addRscanJets("AntiKt",radius,"Truth",jetm9Seq,"JETM9")
        addRscanJets("AntiKt",radius,"TruthWZ",jetm9Seq,"JETM9")
    addRscanJets("AntiKt",radius,"LCTopo",jetm9Seq,"JETM9")

#====================================================================
# Add the containers to the output stream - slimming done here
#====================================================================
from DerivationFrameworkCore.SlimmingHelper import SlimmingHelper
JETM9SlimmingHelper = SlimmingHelper("JETM9SlimmingHelper")
JETM9SlimmingHelper.SmartCollections = ["AntiKt4EMTopoJets","AntiKt4LCTopoJets","PrimaryVertices"]
JETM9SlimmingHelper.AllVariables = ["TruthEvents","MuonSegments","Kt4EMTopoOriginEventShape","Kt4LCTopoOriginEventShape","Kt4EMPFlowEventShape"]
JETM9SlimmingHelper.ExtraVariables = ["TruthVertices.z"]

# Trigger content
JETM9SlimmingHelper.IncludeJetTriggerContent = True

# Add the jet containers to the stream
addJetOutputs(JETM9SlimmingHelper,["JETM9"])
JETM9SlimmingHelper.AppendContentToStream(JETM9Stream)
JETM9Stream.RemoveItem("xAOD::TrigNavigation#*")
JETM9Stream.RemoveItem("xAOD::TrigNavigationAuxInfo#*")
