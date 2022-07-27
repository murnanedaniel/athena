
#====================================================================
# JETM1.py
# reductionConf flag JETM1 in Reco_tf.py
#====================================================================

from DerivationFrameworkCore.DerivationFrameworkMaster import DerivationFrameworkIsMonteCarlo, DerivationFrameworkJob, buildFileName

from DerivationFrameworkPhys import PhysCommon

from JetRecConfig.StandardSmallRJets import AntiKt4PV0Track, AntiKt4EMPFlow, AntiKt4EMPFlowNoPtCut, AntiKt4EMTopoLowPt, AntiKt4EMPFlowCSSK

from DerivationFrameworkJetEtMiss.JetCommon import addJetOutputs, addDAODJets, OutputJets

#====================================================================
# SKIMMING TOOL
#====================================================================
from DerivationFrameworkJetEtMiss import TriggerLists
triggers = TriggerLists.jetTrig()

JETM1SkimmingTools = []

if not DerivationFrameworkIsMonteCarlo:

    from DerivationFrameworkTools.DerivationFrameworkToolsConf import DerivationFramework__TriggerSkimmingTool
    JETM1TrigSkimmingTool = DerivationFramework__TriggerSkimmingTool( name                   = "JETM1TrigSkimmingTool1",
                                                                      TriggerListOR          = triggers )
    ToolSvc += JETM1TrigSkimmingTool

    expression = 'HLT_xe120_pufit_L1XE50'
    from DerivationFrameworkTools.DerivationFrameworkToolsConf import DerivationFramework__xAODStringSkimmingTool
    JETM1OfflineSkimmingTool = DerivationFramework__xAODStringSkimmingTool(name       = "JETM1OfflineSkimmingTool1",
                                                                           expression = expression)
    ToolSvc += JETM1OfflineSkimmingTool

    # OR of the above two selections
    from DerivationFrameworkTools.DerivationFrameworkToolsConf import DerivationFramework__FilterCombinationOR
    JETM1ORTool = DerivationFramework__FilterCombinationOR(name="JETM1ORTool", FilterList=[JETM1TrigSkimmingTool,JETM1OfflineSkimmingTool] )
    ToolSvc+=JETM1ORTool

    JETM1SkimmingTools += [JETM1ORTool]

#=======================================
# CREATE PRIVATE SEQUENCE
#=======================================

jetm1Seq = CfgMgr.AthSequencer("JETM1Sequence")
DerivationFrameworkJob += jetm1Seq

#====================================================================
# SET UP STREAM
#====================================================================
streamName = derivationFlags.WriteDAOD_JETM1Stream.StreamName
fileName   = buildFileName( derivationFlags.WriteDAOD_JETM1Stream )
JETM1Stream = MSMgr.NewPoolRootStream( streamName, fileName )
JETM1Stream.AcceptAlgs(["JETM1Kernel"])

#====================================================================
# THINNING TOOLS
#====================================================================
thinningTools = []

# TrackParticles associated with Muons
from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__MuonTrackParticleThinning
JETM1MuonTPThinningTool = DerivationFramework__MuonTrackParticleThinning(name     = "JETM1MuonTPThinningTool",
                                                                    StreamName              = JETM1Stream.Name,
                                                                    MuonKey                 = "Muons",
                                                                    InDetTrackParticlesKey  = "InDetTrackParticles")
ToolSvc += JETM1MuonTPThinningTool
thinningTools.append(JETM1MuonTPThinningTool)

# TrackParticles associated with electrons
from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__EgammaTrackParticleThinning
JETM1ElectronTPThinningTool = DerivationFramework__EgammaTrackParticleThinning(name                    = "JETM1ElectronTPThinningTool",
                                                                               StreamName              = JETM1Stream.Name,
                                                                               SGKey                   = "Electrons",
                                                                               InDetTrackParticlesKey  = "InDetTrackParticles")
ToolSvc += JETM1ElectronTPThinningTool
thinningTools.append(JETM1ElectronTPThinningTool)

thinning_expression = "InDetTrackParticles.JETM1DFLoose && ( abs(InDetTrackParticles.d0) < 5.0*mm ) && ( abs(DFCommonInDetTrackZ0AtPV*sin(InDetTrackParticles.theta)) < 5.0*mm )"

# TrackParticles associated with small-R jets
from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__JetTrackParticleThinning
JETM1Akt4JetTPThinningTool = DerivationFramework__JetTrackParticleThinning( name          = "JETM1Akt4JetTPThinningTool",
                                                                            StreamName              = streamName,
                                                                            JetKey                  = "AntiKt4EMTopoJets",
                                                                            SelectionString         = "AntiKt4EMTopoJets.pt > 18*GeV",
                                                                            TrackSelectionString    = thinning_expression, 
                                                                            InDetTrackParticlesKey  = "InDetTrackParticles")
ToolSvc += JETM1Akt4JetTPThinningTool
thinningTools.append(JETM1Akt4JetTPThinningTool)

JETM1Akt4PFlowJetTPThinningTool = DerivationFramework__JetTrackParticleThinning( name          = "JETM1Akt4PFlowJetTPThinningTool",
                                                                                 StreamName              = streamName,
                                                                                 JetKey                  = "AntiKt4EMPFlowJets",
                                                                                 SelectionString         = "AntiKt4EMPFlowJets.pt > 18*GeV",
                                                                                 TrackSelectionString    = thinning_expression,
                                                                                 InDetTrackParticlesKey  = "InDetTrackParticles")
ToolSvc += JETM1Akt4PFlowJetTPThinningTool
thinningTools.append(JETM1Akt4PFlowJetTPThinningTool)

#=======================================
# Augmentation tools
#=======================================

augmentationTools = []

from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__InDetTrackSelectionToolWrapper
JETM1TrackSelectionTool = DerivationFramework__InDetTrackSelectionToolWrapper(name = "JETM1TrackSelectionTool",
                                                                              ContainerName = "InDetTrackParticles",
                                                                              DecorationName = "JETM1DFLoose" )

JETM1TrackSelectionTool.TrackSelectionTool.CutLevel = "Loose"
ToolSvc += JETM1TrackSelectionTool
augmentationTools.append(JETM1TrackSelectionTool)

#=======================================
# CREATE THE DERIVATION KERNEL ALGORITHM
#=======================================

from DerivationFrameworkCore.DerivationFrameworkCoreConf import DerivationFramework__DerivationKernel
jetm1Seq += CfgMgr.DerivationFramework__DerivationKernel("JETM1Kernel" ,
                                                         AugmentationTools = augmentationTools,
                                                         SkimmingTools = JETM1SkimmingTools,
                                                         ThinningTools = thinningTools)

#=======================================
# SCHEDULE ADDITIONAL JET DECORATIONS
#=======================================

from JetRecConfig.StandardJetMods import stdJetModifiers
from JetRecConfig.JetRecConfig import getModifier

bJVTTool = getModifier(AntiKt4EMPFlow, stdJetModifiers['bJVT'], stdJetModifiers['bJVT'].modspec)
jetm1Seq += CfgMgr.JetDecorationAlg('bJVTAlg', JetContainer='AntiKt4EMPFlowJets', Decorators=[bJVTTool])

#=======================================
# R = 0.4 track-jets (needed for Rtrk)
#=======================================
jetList = [AntiKt4PV0Track]

#=======================================
# SCHEDULE SMALL-R JETS WITH NO PT CUT
#=======================================
if DerivationFrameworkIsMonteCarlo:
    jetList += [AntiKt4EMPFlowNoPtCut, AntiKt4EMTopoLowPt]

#=======================================
# CSSK R = 0.4 EMPFlow jets
#=======================================
jetList += [AntiKt4EMPFlowCSSK]

addDAODJets(jetList,DerivationFrameworkJob)

OutputJets["JETM1"] = ["AntiKt4PV0TrackJets","AntiKt4EMPFlowCSSKJets","AntiKt4EMPFlowNoPtCutJets","AntiKt4EMTopoLowPtJets"]

#====================================================================
# Add the containers to the output stream - slimming done here
#====================================================================
from DerivationFrameworkCore.SlimmingHelper import SlimmingHelper
JETM1SlimmingHelper = SlimmingHelper("JETM1SlimmingHelper")

JETM1SlimmingHelper.SmartCollections = ["Electrons", "Photons", "Muons", "PrimaryVertices",
                                        "InDetTrackParticles",
                                        "AntiKt4EMTopoJets","AntiKt4EMPFlowJets",
                                        "AntiKt10UFOCSSKJets",
                                        "AntiKt10LCTopoTrimmedPtFrac5SmallR20Jets",
                                        "AntiKt10UFOCSSKSoftDropBeta100Zcut10Jets",
                                        "BTagging_AntiKt4EMPFlow"]

JETM1SlimmingHelper.ExtraVariables  = ["AntiKt4EMTopoJets.DFCommonJets_QGTagger_NTracks.DFCommonJets_QGTagger_TracksWidth.DFCommonJets_QGTagger_TracksC1",
                                       "AntiKt4EMPFlowJets.DFCommonJets_QGTagger_NTracks.DFCommonJets_QGTagger_TracksWidth.DFCommonJets_QGTagger_TracksC1",
                                       "AntiKt4EMPFlowJets.GhostTower",
                                       "AntiKt4EMPFlowJets.passOnlyBJVT","AntiKt4EMPFlowJets.DFCommonJets_bJvt",
                                       "InDetTrackParticles.truthMatchProbability", 
                                       "AntiKt10UFOCSSKSoftDropBeta100Zcut10Jets.zg.rg.NumTrkPt1000.TrackWidthPt1000.GhostMuonSegmentCount.EnergyPerSampling.GhostTrack",
                                       "AntiKt10LCTopoTrimmedPtFrac5SmallR20Jets.zg.rg",
                                       "AntiKt10UFOCSSKJets.NumTrkPt1000.TrackWidthPt1000.GhostMuonSegmentCount.EnergyPerSampling.GhostTrack",
                                       "TruthVertices.barcode.z"]

JETM1SlimmingHelper.AllVariables = [ "MuonSegments", "EventInfo",
                                     "Kt4EMTopoOriginEventShape","Kt4EMPFlowEventShape","Kt4EMPFlowPUSBEventShape","Kt4EMPFlowNeutEventShape",
                                     "CaloCalFwdTopoTowers"]


if DerivationFrameworkIsMonteCarlo:
    from DerivationFrameworkMCTruth import MCTruthCommon
    MCTruthCommon.addBosonsAndDownstreamParticles(generations=4,rejectHadronChildren=True)
    MCTruthCommon.addTopQuarkAndDownstreamParticles(generations=4,rejectHadronChildren=True)

    JETM1SlimmingHelper.AppendToDictionary = {
        'TruthElectrons':'xAOD::TruthParticleContainer','TruthElectronsAux':'xAOD::TruthParticleAuxContainer',
        'TruthMuons':'xAOD::TruthParticleContainer','TruthMuonsAux':'xAOD::TruthParticleAuxContainer',
        'TruthPhotons':'xAOD::TruthParticleContainer','TruthPhotonsAux':'xAOD::TruthParticleAuxContainer',
        'TruthBosonsWithDecayParticles':'xAOD::TruthParticleContainer','TruthBosonsWithDecayParticlesAux':'xAOD::TruthParticleAuxContainer',
        'TruthBosonsWithDecayVertices':'xAOD::TruthVertexContainer','TruthBosonsWithDecayVerticesAux':'xAOD::TruthVertexAuxContainer',
        'TruthTopQuarkWithDecayParticles':'xAOD::TruthParticleContainer','TruthTopQuarkWithDecayParticlesAux':'xAOD::TruthParticleAuxContainer',
        'TruthTopQuarkWithDecayVertices':'xAOD::TruthVertexContainer','TruthTopQuarkWithDecayVerticesAux':'xAOD::TruthVertexAuxContainer',
        'TruthParticles':'xAOD::TruthParticleContainer','TruthParticlesAux':'xAOD::TruthParticleAuxContainer',
        'Kt4EMPFlowNeutEventShape':'xAOD::EventShape','Kt4EMPFlowNeutEventShapeAux':'xAOD::EventShapeAuxInfo'
    }
    JETM1SlimmingHelper.AllVariables += ["TruthMuons", "TruthElectrons", "TruthPhotons", "TruthTopQuarkWithDecayParticles", "TruthBosonsWithDecayParticles"]
    JETM1SlimmingHelper.AllVariables += ["AntiKt4TruthJets", "InTimeAntiKt4TruthJets", "OutOfTimeAntiKt4TruthJets", "TruthParticles"]
    JETM1SlimmingHelper.SmartCollections += ["AntiKt4TruthWZJets"]
else:
    JETM1SlimmingHelper.AppendToDictionary = {
        'Kt4EMPFlowNeutEventShape':'xAOD::EventShape','Kt4EMPFlowNeutEventShapeAux':'xAOD::EventShapeAuxInfo'
    }

# Trigger content
JETM1SlimmingHelper.IncludeJetTriggerContent = True

# Add the jet containers to the stream
addJetOutputs(JETM1SlimmingHelper,
              ["JETM1"],
              JETM1SlimmingHelper.SmartCollections,
              ["AntiKt10EMPFlowJets", "AntiKt10EMPFlowCSSKJets"]# veto list
)

JETM1SlimmingHelper.AppendContentToStream(JETM1Stream)
JETM1Stream.RemoveItem("xAOD::TrigNavigation#*")
JETM1Stream.RemoveItem("xAOD::TrigNavigationAuxInfo#*")
