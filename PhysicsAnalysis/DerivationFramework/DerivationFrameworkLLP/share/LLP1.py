# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

#====================================================================
# LLP1.py
# This defines DAOD_LLP1, a DAOD format for Run 3.
# It contains the variables and objects needed for the large majority
# of physics analyses in ATLAS.
# It requires the reductionConf flag LLP1 in Reco_tf.py
#====================================================================
from AthenaCommon import Logging
from AthenaConfiguration.Enums import LHCPeriod
nanolog = Logging.logging.getLogger('LLP1')

from DerivationFrameworkCore.DerivationFrameworkMaster import buildFileName
from DerivationFrameworkCore.DerivationFrameworkMaster import DerivationFrameworkIsMonteCarlo, DerivationFrameworkJob
from DerivationFrameworkPhys import PhysCommon
from DerivationFrameworkPhys import PhysCommonTrigger
from DerivationFrameworkEGamma import EGammaLRT
from DerivationFrameworkMuons import MuonsLRT
EGammaLRT.makeLRTEGammaDF()
MuonsLRT.makeLRTMuonsDF()
from DerivationFrameworkJetEtMiss import METCommon
from DerivationFrameworkJetEtMiss.METCommon import scheduleMETAssocAlg
from DerivationFrameworkCore import LHE3WeightMetadata

from TriggerMenuMT.TriggerAPI.TriggerAPI import TriggerAPI
from TriggerMenuMT.TriggerAPI.TriggerEnums import TriggerPeriod, TriggerType


#====================================================================
# Set up sequence for this format and add to the top sequence
#====================================================================
SeqLLP1 = CfgMgr.AthSequencer("SeqLLP1")
DerivationFrameworkJob += SeqLLP1

#====================================================================
# SET UP STREAM
#====================================================================
streamName = derivationFlags.WriteDAOD_LLP1Stream.StreamName
fileName   = buildFileName( derivationFlags.WriteDAOD_LLP1Stream )
LLP1Stream = MSMgr.NewPoolRootStream( streamName, fileName )
LLP1Stream.AcceptAlgs(["LLP1Kernel"])

### Thinning, skimming and augmentation tools lists
thinningTools       = []
skimmingTools       = []
augmentationTools   = []

# Special sequence
SeqLLP1 = CfgMgr.AthSequencer("SeqLLP1")

#====================================================================
# Turn on LRT
#====================================================================
from InDetRecExample.InDetJobProperties import InDetFlags
InDetFlags.doR3LargeD0.set_Value_and_Lock(True)

#====================================================================
# Run the LRT merger
#====================================================================
MergedElectronContainer = "StdWithLRTElectrons"
MergedMuonContainer = "StdWithLRTMuons"
MergedTrackCollection = "InDetWithLRTTrackParticles"
from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__TrackParticleMerger
LRTAndStandardTrackParticleMerger = DerivationFramework__TrackParticleMerger(name                        = "LRTAndStandardTrackParticleMerger",
                                                                             TrackParticleLocation       = ["InDetTrackParticles","InDetLargeD0TrackParticles"],
                                                                             OutputTrackParticleLocation = MergedTrackCollection,
                                                                             CreateViewColllection       = True)

ToolSvc += LRTAndStandardTrackParticleMerger
SeqLLP1 += CfgMgr.DerivationFramework__CommonAugmentation("InDetWithLRTLRTMerge",
                                                                         AugmentationTools = [LRTAndStandardTrackParticleMerger])


#====================================================================
# Run the LRT Electron merger
#====================================================================
SeqLLP1 += CfgMgr.CP__ElectronLRTMergingAlg(name="LLP1_ElectronLRTMergingAlg",
                                            PromptElectronLocation  = "Electrons",
                                            LRTElectronLocation     = "LRTElectrons",
                                            OutputCollectionName    = MergedElectronContainer,
                                            isDAOD                  = False,
                                            CreateViewCollection    = True)

#====================================================================
# Run the Muon merger
#====================================================================
SeqLLP1 += CfgMgr.CP__MuonLRTMergingAlg(name="LLP1_MuonLRTMergingAlg",
                                        PromptMuonLocation    = "Muons",
                                        LRTMuonLocation       = "MuonsLRT",
                                        OutputMuonLocation    = MergedMuonContainer,
                                        CreateViewCollection  = True)

#====================================================================
# Run VSI
#====================================================================
from VrtSecInclusive.VrtSecInclusive import VrtSecInclusive
from VrtSecInclusive.VrtSecInclusive_Configuration import setupVSI
from TrkExTools.AtlasExtrapolator import AtlasExtrapolator

ToolSvc += AtlasExtrapolator(name     = "AtlasExtrapolator")

ToolSvc += TrackingCommon.getTrackToVertexIPEstimator(name = "LLP1IPETool")
ToolSvc += TrackingCommon.getInDetTrackToVertexTool(name = "LLP1T2VTool")

# set options related to the vertex fitter
from TrkVKalVrtFitter.TrkVKalVrtFitterConf import Trk__TrkVKalVrtFitter
InclusiveVxFitterTool = Trk__TrkVKalVrtFitter(name                = "InclusiveVxFitter",
                                              Extrapolator        = ToolSvc.AtlasExtrapolator,
                                              IterationNumber     = 30,
                                              )
ToolSvc +=  InclusiveVxFitterTool
InclusiveVxFitterTool.OutputLevel = INFO

from InDetRecExample.TrackingCommon import getInDetPixelConditionsSummaryTool
InDetPixelConditionsSummaryTool = getInDetPixelConditionsSummaryTool()
InDetPixelConditionsSummaryTool.UseByteStreamFEI4 = False
InDetPixelConditionsSummaryTool.UseByteStreamFEI3 = False

VrtSecInclusive_InDet = setupVSI( "" )

VrtSecInclusive_InDet.VertexFitterTool             = InclusiveVxFitterTool
VrtSecInclusive_InDet.Extrapolator                 = ToolSvc.AtlasExtrapolator
VrtSecInclusive_InDet.TrackToVertexIPEstimatorTool = ToolSvc.LLP1IPETool
VrtSecInclusive_InDet.TrackToVertexTool            = ToolSvc.LLP1T2VTool
VrtSecInclusive_InDet.FillIntermediateVertices     = False
VrtSecInclusive_InDet.TrackLocation                = MergedTrackCollection
VrtSecInclusive_InDet.PixelConditionsSummaryTool   = InDetPixelConditionsSummaryTool
VrtSecInclusive_InDet.doAugmentDVimpactParametersToMuons = False
VrtSecInclusive_InDet.doAugmentDVimpactParametersToElectrons = False
SeqLLP1 += VrtSecInclusive_InDet

# leptons-only VSI
LeptonsModSuffix = "_LeptonsMod_LRTR3_1p0"
vsi_lepMod = setupVSI( "InDet"+LeptonsModSuffix, AugSuffix=LeptonsModSuffix )

vsi_lepMod.twoTrkVtxFormingD0Cut = 1.0 # loosen d0 cut to 1 mm
vsi_lepMod.doSelectTracksWithLRTCuts = True # apply addtional track cuts inspired by LRT Run 3 optimizations
vsi_lepMod.doSelectTracksFromMuons    = True # do leptons-only vertexing
vsi_lepMod.doRemoveCaloTaggedMuons    = True # do remove calo-tagged muons from track selection
vsi_lepMod.doSelectTracksFromElectrons  = True # do leptons-only vertexing
vsi_lepMod.VertexFitterTool             = InclusiveVxFitterTool
vsi_lepMod.Extrapolator                 = ToolSvc.AtlasExtrapolator
vsi_lepMod.TrackToVertexIPEstimatorTool = ToolSvc.LLP1IPETool
vsi_lepMod.TrackToVertexTool            = ToolSvc.LLP1T2VTool
vsi_lepMod.FillIntermediateVertices     = False
vsi_lepMod.TrackLocation                = MergedTrackCollection
vsi_lepMod.MuonLocation                 = MergedMuonContainer
vsi_lepMod.ElectronLocation                 = MergedElectronContainer
vsi_lepMod.PixelConditionsSummaryTool   = InDetPixelConditionsSummaryTool
vsi_lepMod.doAugmentDVimpactParametersToMuons = False
vsi_lepMod.doAugmentDVimpactParametersToElectrons = False
SeqLLP1 += vsi_lepMod

LLP1VrtSecInclusiveSuffixes = ["",LeptonsModSuffix]

#====================================================================
# SKIMMING
#====================================================================
allperiods = TriggerPeriod.y2015 | TriggerPeriod.y2016 | TriggerPeriod.y2017 | TriggerPeriod.y2018 | TriggerPeriod.future2e34
trig_el  = TriggerAPI.getLowestUnprescaledAnyPeriod(allperiods, triggerType=TriggerType.el,  livefraction=0.8)
trig_mu  = TriggerAPI.getLowestUnprescaledAnyPeriod(allperiods, triggerType=TriggerType.mu,  livefraction=0.8)
trig_g   = TriggerAPI.getLowestUnprescaledAnyPeriod(allperiods, triggerType=TriggerType.g,   livefraction=0.8)
trig_elmu = TriggerAPI.getLowestUnprescaledAnyPeriod(allperiods, triggerType=TriggerType.el, additionalTriggerType=TriggerType.mu,  livefraction=0.8)
trig_mug = TriggerAPI.getLowestUnprescaledAnyPeriod(allperiods, triggerType=TriggerType.mu, additionalTriggerType=TriggerType.g,  livefraction=0.8)

trig_VBF_2018 =["HLT_j55_gsc80_bmv2c1070_split_j45_gsc60_bmv2c1085_split_j45_320eta490", "HLT_j45_gsc55_bmv2c1070_split_2j45_320eta490_L1J25.0ETA23_2J15.31ETA49", "HLT_j80_0eta240_j60_j45_320eta490_AND_2j35_gsc45_bmv2c1070_split", "HLT_ht300_2j40_0eta490_invm700_L1HT150-J20s5.ETA31_MJJ-400-CF_AND_2j35_gsc45_bmv2c1070_split", "HLT_j70_j50_0eta490_invm1100j70_dphi20_deta40_L1MJJ-500-NFF"]


triggers = trig_el + trig_mu + trig_g + trig_elmu + trig_mug + trig_VBF_2018

#remove duplicates
triggers = sorted(list(set(triggers)))

#trigger
from DerivationFrameworkTools.DerivationFrameworkToolsConf import DerivationFramework__TriggerSkimmingTool
LLP1TriggerSkimmingTool = DerivationFramework__TriggerSkimmingTool(name = "LLP1TriggerPreSkimmingTool",
                                                              TriggerListAND = [],
                                                              TriggerListOR  = triggers)
ToolSvc += LLP1TriggerSkimmingTool

print('LLP1 list of triggers used for skimming:')
for trig in triggers: print(trig)

skimmingTools.append(LLP1TriggerSkimmingTool)


#====================================================================
# THINNING
#====================================================================
# ID tracks: See recommedations here:
# https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/DaodRecommendations

# Inner detector group recommendations for indet tracks in analysis
LLP1_thinning_expression = "InDetTrackParticles.DFCommonTightPrimary && abs(DFCommonInDetTrackZ0AtPV)*sin(InDetTrackParticles.theta) < 3.0*mm && InDetTrackParticles.pt > 10*GeV"
from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__TrackParticleThinning
LLP1TrackParticleThinningTool = DerivationFramework__TrackParticleThinning(name                    = "LLP1TrackParticleThinningTool",
                                                                           StreamName              = LLP1Stream.Name,
                                                                           SelectionString         = LLP1_thinning_expression,
                                                                           InDetTrackParticlesKey  = "InDetTrackParticles")

ToolSvc += LLP1TrackParticleThinningTool
thinningTools.append(LLP1TrackParticleThinningTool)

# Include inner detector tracks associated with muons
from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__MuonTrackParticleThinning
LLP1MuonTPThinningTool = DerivationFramework__MuonTrackParticleThinning(name                    = "LLP1MuonTPThinningTool",
                                                                        StreamName              = LLP1Stream.Name,
                                                                        MuonKey                 = "Muons",
                                                                        InDetTrackParticlesKey  = "InDetTrackParticles")

ToolSvc += LLP1MuonTPThinningTool
thinningTools.append(LLP1MuonTPThinningTool)

# Include LRT inner detector tracks associated with LRT muons
from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__MuonTrackParticleThinning
LLP1LRTMuonTPThinningTool = DerivationFramework__MuonTrackParticleThinning(name                    = "LLP1LRRMuonTPThinningTool",
                                                                              StreamName              = LLP1Stream.Name,
                                                                              MuonKey                 = "MuonsLRT",
                                                                              InDetTrackParticlesKey  = "InDetLargeD0TrackParticles")

ToolSvc += LLP1LRTMuonTPThinningTool
thinningTools.append(LLP1LRTMuonTPThinningTool)

# TauJets thinning
# Disabled for 1st production in 2021, to allow use by tau CP group
#tau_thinning_expression = "(TauJets.ptFinalCalib >= 13.*GeV) && (TauJets.nTracks>=1) && (TauJets.nTracks<=3) && (TauJets.RNNJetScoreSigTrans>0.01)"
#from DerivationFrameworkTools.DerivationFrameworkToolsConf import DerivationFramework__GenericObjectThinning
#LLP1TauJetsThinningTool = DerivationFramework__GenericObjectThinning(name            = "LLP1TauJetsThinningTool",
#                                                                     StreamName      = LLP1Stream.Name,
#                                                                     ContainerName   = "TauJets",
#                                                                     SelectionString = tau_thinning_expression)
#ToolSvc += LLP1TauJetsThinningTool
#thinningTools.append(LLP1TauJetsThinningTool)

# Only keep tau tracks (and associated ID tracks) classified as charged tracks
from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__TauTrackParticleThinning
LLP1TauTPThinningTool = DerivationFramework__TauTrackParticleThinning(name                   = "LLP1TauTPThinningTool",
                                                                      StreamName             = LLP1Stream.Name,
                                                                      TauKey                 = "TauJets",
                                                                      InDetTrackParticlesKey = "InDetTrackParticles",
#                                                                      SelectionString        = tau_thinning_expression,
                                                                      DoTauTracksThinning    = True,
                                                                      TauTracksKey           = "TauTracks")
ToolSvc += LLP1TauTPThinningTool
thinningTools.append(LLP1TauTPThinningTool)

# ID tracks associated with high-pt di-tau
from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__DiTauTrackParticleThinning
LLP1DiTauTPThinningTool = DerivationFramework__DiTauTrackParticleThinning(name                    = "LLP1DiTauTPThinningTool",
                                                                          StreamName              = LLP1Stream.Name,
                                                                          DiTauKey                = "DiTauJets",
                                                                          InDetTrackParticlesKey  = "InDetTrackParticles")
ToolSvc += LLP1DiTauTPThinningTool
thinningTools.append(LLP1DiTauTPThinningTool)

# Low-pt di-tau thinning
from DerivationFrameworkTools.DerivationFrameworkToolsConf import DerivationFramework__GenericObjectThinning
LLP1DiTauLowPtThinningTool = DerivationFramework__GenericObjectThinning(name            = "LLP1DiTauLowPtThinningTool",
                                                                        StreamName      = LLP1Stream.Name,
                                                                        ContainerName   = "DiTauJetsLowPt",
                                                                        SelectionString = "DiTauJetsLowPt.nSubjets > 1")
ToolSvc += LLP1DiTauLowPtThinningTool
thinningTools.append(LLP1DiTauLowPtThinningTool)

# ID tracks associated with low-pt ditau
LLP1DiTauLowPtTPThinningTool = DerivationFramework__DiTauTrackParticleThinning(name                    = "LLP1DiTauLowPtTPThinningTool",
                                                                               StreamName              = LLP1Stream.Name,
                                                                               DiTauKey                = "DiTauJetsLowPt",
                                                                               InDetTrackParticlesKey  = "InDetTrackParticles",
                                                                               SelectionString         = "DiTauJetsLowPt.nSubjets > 1")
ToolSvc += LLP1DiTauLowPtTPThinningTool
thinningTools.append(LLP1DiTauLowPtTPThinningTool)

# ID Tracks associated with jets
from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__JetTrackParticleThinning
LLP1JetTPThinningTool = DerivationFramework__JetTrackParticleThinning( name                    = "LLP1JetTPThinningTool",
                                                                          StreamName              = LLP1Stream.Name,
                                                                          JetKey                  = "AntiKt4EMTopoJets",
                                                                          SelectionString         = "(AntiKt4EMTopoJets.pt > 20.*GeV) && (abs(AntiKt4EMTopoJets.eta) < 2.5)",
                                                                          InDetTrackParticlesKey  = "InDetTrackParticles")
ToolSvc += LLP1JetTPThinningTool
thinningTools.append(LLP1JetTPThinningTool)

# LRT Tracks associated with jets
if InDetFlags.doR3LargeD0():
  from DerivationFrameworkLLP.DerivationFrameworkLLPConf import DerivationFramework__JetLargeD0TrackParticleThinning
  LLP1JetLD0TPThinningTool = DerivationFramework__JetLargeD0TrackParticleThinning( name                    = "LLP1JetLD0TPThinningTool",
                                                                                      StreamName              = LLP1Stream.Name,
                                                                                      JetKey                  = "AntiKt4EMTopoJets",
                                                                                      SelectionString         = "(AntiKt4EMTopoJets.pt > 20.*GeV) && (abs(AntiKt4EMTopoJets.eta) < 2.5)",
                                                                                      InDetTrackParticlesKey  = "InDetLargeD0TrackParticles")
  ToolSvc += LLP1JetLD0TPThinningTool
  thinningTools.append(LLP1JetLD0TPThinningTool)

# ID Tracks associated with secondary vertices
from DerivationFrameworkLLP.DerivationFrameworkLLPConf import DerivationFramework__VSITrackParticleThinning
LLP1VSITPThinningTool = DerivationFramework__VSITrackParticleThinning( name                  = "LLP1VSITPThinningTool",
                                                                          StreamName              = LLP1Stream.Name,
                                                                          InDetTrackParticlesKey  = "InDetTrackParticles",
                                                                          AugVerStrings = LLP1VrtSecInclusiveSuffixes)

ToolSvc += LLP1VSITPThinningTool
thinningTools.append(LLP1VSITPThinningTool)

from DerivationFrameworkLLP.DerivationFrameworkLLPConf import DerivationFramework__VSITrackParticleThinning
LLP1LD0VSITPThinningTool = DerivationFramework__VSITrackParticleThinning( name                  = "LLP1LD0VSITPThinningTool",
                                                                            StreamName              = LLP1Stream.Name,
                                                                            InDetTrackParticlesKey  = "InDetLargeD0TrackParticles",
                                                                            AugVerStrings = LLP1VrtSecInclusiveSuffixes)
ToolSvc += LLP1LD0VSITPThinningTool
thinningTools.append(LLP1LD0VSITPThinningTool)

#====================================================================
# Max Cell sum decoration tool
#====================================================================

from LArCabling.LArCablingAccess import LArOnOffIdMapping
#LArOnOffIdMapping()

from DerivationFrameworkCalo.DerivationFrameworkCaloConf import DerivationFramework__MaxCellDecorator

LLP1_MaxCellDecoratorTool = DerivationFramework__MaxCellDecorator( name = "LLP1_MaxCellDecoratorTool",
                                                                      SGKey_electrons = "Electrons",
                                                                      SGKey_photons   = "Photons"
                                                                      )
ToolSvc += LLP1_MaxCellDecoratorTool
augmentationTools.append(LLP1_MaxCellDecoratorTool)

if ConfigFlags.GeoModel.Run == LHCPeriod.Run3:
    LLP1_LRTMaxCellDecoratorTool = DerivationFramework__MaxCellDecorator( name = "LLP1_LRTMaxCellDecoratorTool",
                                                                      SGKey_electrons = "LRTElectrons",
                                                                      )
else:
    LLP1_LRTMaxCellDecoratorTool = DerivationFramework__MaxCellDecorator( name = "LLP1_LRTMaxCellDecoratorTool",
                                                                      SGKey_electrons = "LRTElectrons",
                                                                      SGKey_egammaClusters = "egammaClusters",
                                                                      )

ToolSvc += LLP1_LRTMaxCellDecoratorTool
augmentationTools.append(LLP1_LRTMaxCellDecoratorTool)

#====================================================================
# SUSY Track Particle Calo Cell Decorator 
#====================================================================
from DerivationFrameworkSUSY.DerivationFrameworkSUSYConf import DerivationFramework__TrackParticleCaloCellDecorator
LLP1_TrackParticleCaloCellDecorator = DerivationFramework__TrackParticleCaloCellDecorator(
        name = "LLP1_TrackParticleCaloCellDecorator",
        DecorationPrefix = "LLP1",
        ContainerName    = "InDetTrackParticles")
ToolSvc += LLP1_TrackParticleCaloCellDecorator
augmentationTools.append(LLP1_TrackParticleCaloCellDecorator)
    
## no InDetLargeD0TrackParticlesClusterAssociations for LRT tracks... 
LLP1_TrackParticleCaloCellDecorator_LargeD0 = DerivationFramework__TrackParticleCaloCellDecorator(
        name = "LLP1_TrackParticleCaloCellDecorator_LargeD0",
        DecorationPrefix = "LLP1",
        ContainerName    = "InDetLargeD0TrackParticles")
ToolSvc += LLP1_TrackParticleCaloCellDecorator_LargeD0
augmentationTools.append(LLP1_TrackParticleCaloCellDecorator_LargeD0)

#====================================================================
# TRIGGER MATCHING WITH LRT LEPTONS (does not include tau chains)
#===================================================================
from DerivationFrameworkTrigger.TriggerMatchingHelper import TriggerMatchingHelper
lrt_tm_helper = TriggerMatchingHelper(
    PhysCommonTrigger.trigger_names_notau,
    name="LRTDFTriggerMatchingTool",
    OutputContainerPrefix="LRTTrigMatch_",
    InputElectrons=MergedElectronContainer,
    InputMuons=MergedMuonContainer
)

ToolSvc += lrt_tm_helper.tool
augmentationTools.append(lrt_tm_helper.tool)

#====================================================================
# CREATE THE DERIVATION KERNEL ALGORITHM
#====================================================================
# Add the kernel for thinning (requires the objects be defined)
from DerivationFrameworkCore.DerivationFrameworkCoreConf import DerivationFramework__DerivationKernel
SeqLLP1 += CfgMgr.DerivationFramework__DerivationKernel("LLP1Kernel",
                                                        SkimmingTools = skimmingTools,
                                                        AugmentationTools = augmentationTools,
                                                        ThinningTools = thinningTools)

#====================================================================
# FLAVOUR TAGGING
#====================================================================

from DerivationFrameworkFlavourTag.FtagDerivationConfig import FtagJetCollections
FtagJetCollections(['AntiKt4EMTopoJets'],SeqLLP1)


#====================================================================
# CONTENTS
#====================================================================
from DerivationFrameworkCore.SlimmingHelper import SlimmingHelper
LLP1SlimmingHelper = SlimmingHelper("LLP1SlimmingHelper")

LLP1SlimmingHelper.SmartCollections = ["EventInfo",
                                       "Electrons",
                                       "LRTElectrons",
                                       "Photons",
                                       "Muons",
                                       "MuonsLRT",
                                       "PrimaryVertices",
                                       "InDetTrackParticles",
                                       "InDetLargeD0TrackParticles",
                                       "AntiKt4EMTopoJets",
                                       "AntiKt4EMPFlowJets",
                                       "BTagging_AntiKt4EMTopo",
                                       "BTagging_AntiKt4EMPFlow",
                                       "BTagging_AntiKtVR30Rmax4Rmin02Track",
                                       "MET_Baseline_AntiKt4EMTopo",
                                       "MET_Baseline_AntiKt4EMPFlow",
                                       "TauJets",
                                       "TauJets_MuonRM",
                                       "DiTauJets",
                                       "DiTauJetsLowPt",
                                       "AntiKt10LCTopoTrimmedPtFrac5SmallR20Jets",
                                       "AntiKt10UFOCSSKSoftDropBeta100Zcut10Jets",
                                       "AntiKtVR30Rmax4Rmin02PV0TrackJets",
                                      ]

LLP1SlimmingHelper.AllVariables = ["MSDisplacedVertex",
                                   "MuonSpectrometerTrackParticles",
                                   "MuonSegments",
                                   "MSonlyTracklets",
                                   "CombinedMuonTrackParticles",
                                   "ExtrapolatedMuonTrackParticles",
                                   "CombinedMuonsLRTTrackParticles",
                                   "ExtraPolatedMuonsLRTTrackParticles",
                                   "MSOnlyExtraPolatedMuonsLRTTrackParticles",
                                  ]


excludedVertexAuxData = "-vxTrackAtVertex.-MvfFitInfo.-isInitialized.-VTAV"
StaticContent = []
StaticContent += ["xAOD::VertexContainer#SoftBVrtClusterTool_Tight_Vertices"]
StaticContent += ["xAOD::VertexAuxContainer#SoftBVrtClusterTool_Tight_VerticesAux." + excludedVertexAuxData]
StaticContent += ["xAOD::VertexContainer#SoftBVrtClusterTool_Medium_Vertices"]
StaticContent += ["xAOD::VertexAuxContainer#SoftBVrtClusterTool_Medium_VerticesAux." + excludedVertexAuxData]
StaticContent += ["xAOD::VertexContainer#SoftBVrtClusterTool_Loose_Vertices"]
StaticContent += ["xAOD::VertexAuxContainer#SoftBVrtClusterTool_Loose_VerticesAux." + excludedVertexAuxData]

for wp in LLP1VrtSecInclusiveSuffixes:
    StaticContent += ["xAOD::VertexContainer#VrtSecInclusive_SecondaryVertices" + wp]
    StaticContent += ["xAOD::VertexAuxContainer#VrtSecInclusive_SecondaryVertices" + wp + "Aux."]

LLP1SlimmingHelper.StaticContent = StaticContent

# Trigger content
LLP1SlimmingHelper.IncludeTriggerNavigation = False
LLP1SlimmingHelper.IncludeJetTriggerContent = False
LLP1SlimmingHelper.IncludeMuonTriggerContent = False
LLP1SlimmingHelper.IncludeEGammaTriggerContent = False
LLP1SlimmingHelper.IncludeJetTauEtMissTriggerContent = False
LLP1SlimmingHelper.IncludeTauTriggerContent = False
LLP1SlimmingHelper.IncludeEtMissTriggerContent = False
LLP1SlimmingHelper.IncludeBJetTriggerContent = False
LLP1SlimmingHelper.IncludeBPhysTriggerContent = False
LLP1SlimmingHelper.IncludeMinBiasTriggerContent = False

# Truth containers
if DerivationFrameworkIsMonteCarlo:
   LLP1SlimmingHelper.AppendToDictionary = {'TruthEvents':'xAOD::TruthEventContainer','TruthEventsAux':'xAOD::TruthEventAuxContainer',
                                            'MET_Truth':'xAOD::MissingETContainer','MET_TruthAux':'xAOD::MissingETAuxContainer',
                                            'TruthElectrons':'xAOD::TruthParticleContainer','TruthElectronsAux':'xAOD::TruthParticleAuxContainer',
                                            'TruthMuons':'xAOD::TruthParticleContainer','TruthMuonsAux':'xAOD::TruthParticleAuxContainer',
                                            'TruthPhotons':'xAOD::TruthParticleContainer','TruthPhotonsAux':'xAOD::TruthParticleAuxContainer',
                                            'TruthTaus':'xAOD::TruthParticleContainer','TruthTausAux':'xAOD::TruthParticleAuxContainer',
                                            'TruthNeutrinos':'xAOD::TruthParticleContainer','TruthNeutrinosAux':'xAOD::TruthParticleAuxContainer',
                                            'TruthBSM':'xAOD::TruthParticleContainer','TruthBSMAux':'xAOD::TruthParticleAuxContainer',
                                            'TruthBoson':'xAOD::TruthParticleContainer','TruthBosonAux':'xAOD::TruthParticleAuxContainer',
                                            'TruthTop':'xAOD::TruthParticleContainer','TruthTopAux':'xAOD::TruthParticleAuxContainer',
                                            'TruthForwardProtons':'xAOD::TruthParticleContainer','TruthForwardProtonsAux':'xAOD::TruthParticleAuxContainer',
                                            'BornLeptons':'xAOD::TruthParticleContainer','BornLeptonsAux':'xAOD::TruthParticleAuxContainer',
                                            'TruthBosonsWithDecayParticles':'xAOD::TruthParticleContainer','TruthBosonsWithDecayParticlesAux':'xAOD::TruthParticleAuxContainer',
                                            'TruthBosonsWithDecayVertices':'xAOD::TruthVertexContainer','TruthBosonsWithDecayVerticesAux':'xAOD::TruthVertexAuxContainer',
                                            'TruthBSMWithDecayParticles':'xAOD::TruthParticleContainer','TruthBSMWithDecayParticlesAux':'xAOD::TruthParticleAuxContainer',
                                            'TruthBSMWithDecayVertices':'xAOD::TruthVertexContainer','TruthBSMWithDecayVerticesAux':'xAOD::TruthVertexAuxContainer',
                                            'HardScatterParticles':'xAOD::TruthParticleContainer','HardScatterParticlesAux':'xAOD::TruthParticleAuxContainer',
                                            'HardScatterVertices':'xAOD::TruthVertexContainer','HardScatterVerticesAux':'xAOD::TruthVertexAuxContainer',
                                            'TruthHFWithDecayParticles':'xAOD::TruthParticleContainer','TruthHFWithDecayParticlesAux':'xAOD::TruthParticleAuxContainer',
                                            'TruthHFWithDecayVertices':'xAOD::TruthVertexContainer','TruthHFWithDecayVerticesAux':'xAOD::TruthVertexAuxContainer',
                                            'TruthCharm':'xAOD::TruthParticleContainer','TruthCharmAux':'xAOD::TruthParticleAuxContainer',
                                            'TruthPrimaryVertices':'xAOD::TruthVertexContainer','TruthPrimaryVerticesAux':'xAOD::TruthVertexAuxContainer',
                                            'AntiKt10TruthTrimmedPtFrac5SmallR20Jets':'xAOD::JetContainer', 'AntiKt10TruthTrimmedPtFrac5SmallR20JetsAux':'xAOD::JetAuxContainer'
                                           }

   from DerivationFrameworkMCTruth.MCTruthCommon import addTruth3ContentToSlimmerTool
   addTruth3ContentToSlimmerTool(LLP1SlimmingHelper)
   LLP1SlimmingHelper.AllVariables += ['TruthHFWithDecayParticles','TruthHFWithDecayVertices','TruthCharm','TruthPileupParticles','InTimeAntiKt4TruthJets','OutOfTimeAntiKt4TruthJets']

LLP1SlimmingHelper.ExtraVariables += ["AntiKt10TruthTrimmedPtFrac5SmallR20Jets.Tau1_wta.Tau2_wta.Tau3_wta.D2.GhostBHadronsFinalCount",
                                      "Electrons.TruthLink","LRTElectrons.TruthLink",
                                      "Electrons.LHValue.DFCommonElectronsLHVeryLooseNoPixResult.maxEcell_time.maxEcell_energy.maxEcell_gain.maxEcell_onlId.maxEcell_x.maxEcell_y.maxEcell_z.f3",
                                      "LRTElectrons.LHValue.DFCommonElectronsLHVeryLooseNoPixResult.maxEcell_time.maxEcell_energy.maxEcell_gain.maxEcell_onlId.maxEcell_x.maxEcell_y.maxEcell_z.f3",
                                      "Photons.maxEcell_time.maxEcell_energy.maxEcell_gain.maxEcell_onlId.maxEcell_x.maxEcell_y.maxEcell_z.f3",
                                      "egammaClusters.phi_sampl.eta0.phi0",
                                      "LRTegammaClusters.phi_sampl.eta0.phi0",
                                      "Muons.TruthLink", "MuonsLRT.TruthLink",
                                      "Photons.TruthLink",
                                      "AntiKt4EMTopoJets.DFCommonJets_QGTagger_truthjet_nCharged.DFCommonJets_QGTagger_truthjet_pt.DFCommonJets_QGTagger_truthjet_eta.DFCommonJets_QGTagger_NTracks.DFCommonJets_QGTagger_TracksWidth.DFCommonJets_QGTagger_TracksC1.PartonTruthLabelID.ConeExclBHadronsFinal.ConeExclCHadronsFinal.GhostBHadronsFinal.GhostCHadronsFinal.GhostBHadronsFinalCount.GhostBHadronsFinalPt.GhostCHadronsFinalCount.GhostCHadronsFinalPt.GhostBHadronsFinal.GhostCHadronsFinal.GhostTrack.GhostTrackCount.GhostTrackLRT.GhostTrackLRTCount",
                                      "AntiKt4EMPFlowJets.DFCommonJets_QGTagger_truthjet_nCharged.DFCommonJets_QGTagger_truthjet_pt.DFCommonJets_QGTagger_truthjet_eta.DFCommonJets_QGTagger_NTracks.DFCommonJets_QGTagger_TracksWidth.DFCommonJets_QGTagger_TracksC1.PartonTruthLabelID.DFCommonJets_fJvt.ConeExclBHadronsFinal.ConeExclCHadronsFinal.GhostBHadronsFinal.GhostCHadronsFinal.GhostBHadronsFinalCount.GhostBHadronsFinalPt.GhostCHadronsFinalCount.GhostCHadronsFinalPt.GhostBHadronsFinal.GhostCHadronsFinal",
                                      "AntiKtVR30Rmax4Rmin02TrackJets_BTagging201903.GhostBHadronsFinal.GhostCHadronsFinal.GhostBHadronsFinalCount.GhostBHadronsFinalPt.GhostCHadronsFinalCount.GhostCHadronsFinalPt.GhostTausFinal.GhostTausFinalCount",
                                      "AntiKtVR30Rmax4Rmin02TrackJets_BTagging201810.GhostBHadronsFinal.GhostCHadronsFinal.GhostBHadronsFinalCount.GhostBHadronsFinalPt.GhostCHadronsFinalCount.GhostCHadronsFinalPt.GhostTausFinal.GhostTausFinalCount",
                                      "TruthPrimaryVertices.t.x.y.z",
                                      "PrimaryVertices.t.x.y.z",
                                      "InDetTrackParticles.d0.z0.vz.TTVA_AMVFVertices.TTVA_AMVFWeights.eProbabilityHT.truthParticleLink.truthMatchProbability.radiusO.fFirstHit.hitPattern.LLP1_CaloCelldEta.LLP1_CaloCelldPhi.LLP1_CaloCelldR.LLP1_CaloCelldX.LLP1_CaloCelldY.LLP1_CaloCelldZ.LLP1_CaloCellE.LLP1_CaloCellEta.LLP1_CaloCellPhi.LLP1_CaloCellR.LLP1_CaloCellSampling.LLP1_CaloCellTime.LLP1_CaloCellX.LLP1_CaloCellY.LLP1_CaloCellZ",
                                      "InDetLargeD0TrackParticles.d0.z0.vz.TTVA_AMVFVertices.TTVA_AMVFWeights.eProbabilityHT.truthParticleLink.truthMatchProbability.radiusOfFirstHit.hitPattern.LLP1_CaloCelldEta.LLP1_CaloCelldPhi.LLP1_CaloCelldR.LLP1_CaloCelldX.LLP1_CaloCelldY.LLP1_CaloCelldZ.LLP1_CaloCellE.LLP1_CaloCellEta.LLP1_CaloCellPhi.LLP1_CaloCellR.LLP1_CaloCellSampling.LLP1_CaloCellTime.LLP1_CaloCellX.LLP1_CaloCellY.LLP1_CaloCellZ",
                                      "GSFTrackParticles.d0.z0.vz.TTVA_AMVFVertices.TTVA_AMVFWeights.eProbabilityHT.truthParticleLink.truthMatchProbability.radiusOfFirstHit.numberOfPixelHoles.numberOfSCTHoles.numberDoF.chiSquared.hitPattern",
                                      "LRTGSFTrackParticles.d0.z0.vz.TTVA_AMVFVertices.TTVA_AMVFWeights.eProbabilityHT.truthParticleLink.truthMatchProbability.radiusOfFirstHit.numberOfPixelHoles.numberOfSCTHoles.numberDoF.chiSquared.hitPattern",
                                      "EventInfo.hardScatterVertexLink.timeStampNSOffset",
                                      "TauJets.dRmax.etOverPtLeadTrk",
                                      "HLT_xAOD__TrigMissingETContainer_TrigEFMissingET.ex.ey",
                                      "HLT_xAOD__TrigMissingETContainer_TrigEFMissingET_mht.ex.ey"]


VSITrackAuxVars = [
    "is_selected", "is_associated", "is_svtrk_final", "pt_wrtSV", "eta_wrtSV",
    "phi_wrtSV", "d0_wrtSV", "z0_wrtSV", "errP_wrtSV", "errd0_wrtSV",
    "errz0_wrtSV", "chi2_toSV"
]

for suffix in LLP1VrtSecInclusiveSuffixes:
    LLP1SlimmingHelper.ExtraVariables += [ "InDetTrackParticles." + '.'.join( [ var + suffix for var in VSITrackAuxVars] ) ]
    LLP1SlimmingHelper.ExtraVariables += [ "InDetLargeD0TrackParticles." + '.'.join( [ var + suffix for var in VSITrackAuxVars] ) ]
    LLP1SlimmingHelper.ExtraVariables += [ "GSFTrackParticles." + '.'.join( [ var + suffix for var in VSITrackAuxVars] ) ]
    LLP1SlimmingHelper.ExtraVariables += [ "LRTGSFTrackParticles." + '.'.join( [ var + suffix for var in VSITrackAuxVars] ) ]

# Add trigger matching
# Run 2
PhysCommonTrigger.trigmatching_helper_notau.add_to_slimming(LLP1SlimmingHelper)
PhysCommonTrigger.trigmatching_helper_tau.add_to_slimming(LLP1SlimmingHelper)

# LRT Lepton trigger matching variables
lrt_tm_helper.add_to_slimming(LLP1SlimmingHelper)


# Final construction of output stream
LLP1SlimmingHelper.AppendContentToStream(LLP1Stream)


if not hasattr(ServiceMgr, 'THistSvc'):
    from GaudiSvc.GaudiSvcConf import THistSvc
    ServiceMgr += THistSvc()

ServiceMgr.THistSvc.Output += [ "AANT  DATAFILE='VrtSecInclusive.root' OPT='RECREATE'" ]
