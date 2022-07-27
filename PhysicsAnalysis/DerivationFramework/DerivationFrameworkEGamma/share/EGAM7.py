#********************************************************************
# EGAM7.py - keep events passing OR of electron triggers, or inclusive
#            electron selection, to retain fake electron candidates
# reductionConf flag EGAM7 in Reco_tf.py
# author: giovanni.marchiori@cern.ch
#********************************************************************

from DerivationFrameworkJetEtMiss.JetCommon import addDAODJets
from JetRecConfig.StandardSmallRJets import AntiKt4Truth,AntiKt4TruthDressedWZ,AntiKt4PV0Track
from DerivationFrameworkCore.DerivationFrameworkMaster import buildFileName
from DerivationFrameworkCore.DerivationFrameworkMaster import (
    DerivationFrameworkIsMonteCarlo, DerivationFrameworkJob)
from DerivationFrameworkPhys import PhysCommon
from DerivationFrameworkEGamma.EGammaCommon import *
from DerivationFrameworkEGamma.EGAM7ExtraContent import *


#====================================================================
# read common DFEGamma settings from egammaDFFlags
#====================================================================
from DerivationFrameworkEGamma.egammaDFFlags import jobproperties
jobproperties.egammaDFFlags.print_JobProperties("full")

# additional settings for this derivation
thinCells = True

#====================================================================
# check if we run on data or MC
#====================================================================
print("DerivationFrameworkIsMonteCarlo: ", DerivationFrameworkIsMonteCarlo)
if not DerivationFrameworkIsMonteCarlo:
    ExtraContainersTrigger += ExtraContainersTriggerDataOnly


#====================================================================
# Set up sequence for this format and add to the top sequence
#====================================================================
EGAM7Sequence = CfgMgr.AthSequencer("EGAM7Sequence")
DerivationFrameworkJob += EGAM7Sequence


#====================================================================
# SET UP STREAM
#====================================================================
streamName = derivationFlags.WriteDAOD_EGAM7Stream.StreamName
fileName   = buildFileName( derivationFlags.WriteDAOD_EGAM7Stream )
EGAM7Stream = MSMgr.NewPoolRootStream( streamName, fileName )
# Only events that pass the filters listed below are written out.
# Name must match that of the kernel above
# AcceptAlgs  = logical OR of filters
# RequireAlgs = logical AND of filters
EGAM7Stream.AcceptAlgs(["EGAM7Kernel"])


### Thinning and augmentation tools lists
augmentationTools = []
thinningTools=[]


#====================================================================
# SET UP AUGMENTATIONS
#====================================================================


#====================================================================
# Max Cell sum decoration tool
#====================================================================
from DerivationFrameworkCalo.DerivationFrameworkCaloConf import DerivationFramework__MaxCellDecorator
EGAM7_MaxCellDecoratorTool = DerivationFramework__MaxCellDecorator( name                    = "EGAM7_MaxCellDecoratorTool",
                                                                    SGKey_electrons         = "Electrons",
                                                                    SGKey_photons           = "Photons",
                                                                    )
ToolSvc += EGAM7_MaxCellDecoratorTool
augmentationTools += [EGAM7_MaxCellDecoratorTool]


#====================================================================
# Gain and cluster energies per layer decoration tool
#====================================================================
from DerivationFrameworkCalo.DerivationFrameworkCaloFactories import GainDecorator, getGainDecorations, getClusterEnergyPerLayerDecorator, getClusterEnergyPerLayerDecorations
EGAM7_GainDecoratorTool = GainDecorator()
ToolSvc += EGAM7_GainDecoratorTool
augmentationTools += [EGAM7_GainDecoratorTool]

cluster_sizes = (3,7), (5,5), (7,11)
EGAM7_ClusterEnergyPerLayerDecorators = [getClusterEnergyPerLayerDecorator(neta, nphi)() for neta, nphi in cluster_sizes]
augmentationTools += EGAM7_ClusterEnergyPerLayerDecorators


#====================================================================
# SET UP THINNING
#====================================================================

print('WARNING, Thinning of trigger navigation has to be properly implemented in R22')
#from DerivationFrameworkCore.ThinningHelper import ThinningHelper
#EGAM7ThinningHelper = ThinningHelper( "EGAM7ThinningHelper" )
#EGAM7ThinningHelper.TriggerChains = '(^(?!.*_[0-9]*(mu|j|xe|tau|ht|xs|te))(?!HLT_[eg].*_[0-9]*[eg][0-9].*)(?!HLT_eb.*)(?!.*larpeb.*)(?!HLT_.*_AFP_.*)(HLT_[eg].*))|HLT_e.*_Zee.*'
#EGAM7ThinningHelper.AppendToStream( EGAM7Stream, ExtraContainersTrigger )


# Track thinning
if jobproperties.egammaDFFlags.doEGammaDAODTrackThinning:

    TrackThinningKeepElectronTracks = True
    TrackThinningKeepPhotonTracks = True
    TrackThinningKeepJetTracks = False
    TrackThinningKeepMuonTracks = False
    TrackThinningKeepTauTracks = False
    TrackThinningKeepPVTracks = False

    # Tracks associated with Electrons
    if (TrackThinningKeepElectronTracks) :
        from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__EgammaTrackParticleThinning
        EGAM7ElectronTPThinningTool = DerivationFramework__EgammaTrackParticleThinning( name                    = "EGAM7ElectronTPThinningTool",
                                                                                        StreamName              = streamName,
                                                                                        SGKey                   = "Electrons",
                                                                                        GSFTrackParticlesKey    = "GSFTrackParticles",
                                                                                        InDetTrackParticlesKey  = "InDetTrackParticles",
                                                                                        SelectionString         = "Electrons.pt > 0*GeV",
                                                                                        BestMatchOnly = True,
                                                                                        ConeSize = 0.3)
        ToolSvc += EGAM7ElectronTPThinningTool
        print(EGAM7ElectronTPThinningTool)
        thinningTools.append(EGAM7ElectronTPThinningTool)

    # Tracks associated with Photons
    if (TrackThinningKeepPhotonTracks) :
        from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__EgammaTrackParticleThinning
        EGAM7PhotonTPThinningTool = DerivationFramework__EgammaTrackParticleThinning( name                    = "EGAM7PhotonTPThinningTool",
                                                                                      StreamName              = streamName,
                                                                                      SGKey                   = "Photons",
                                                                                      GSFTrackParticlesKey    = "GSFTrackParticles",
                                                                                      InDetTrackParticlesKey  = "InDetTrackParticles",
                                                                                      SelectionString         = "Photons.pt > 0*GeV",
                                                                                      BestMatchOnly = True,
                                                                                      ConeSize = 0.3)
        ToolSvc += EGAM7PhotonTPThinningTool
        print(EGAM7PhotonTPThinningTool)
        thinningTools.append(EGAM7PhotonTPThinningTool)

    # Tracks associated with Jets
    if (TrackThinningKeepJetTracks) :
        from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__JetTrackParticleThinning
        EGAM7JetTPThinningTool = DerivationFramework__JetTrackParticleThinning( name                    = "EGAM7JetTPThinningTool",
                                                                                StreamName              = streamName,
                                                                                JetKey                  = "AntiKt4EMPFlowJets",
                                                                                InDetTrackParticlesKey  = "InDetTrackParticles")
        ToolSvc += EGAM7JetTPThinningTool
        print(EGAM7JetTPThinningTool)
        thinningTools.append(EGAM7JetTPThinningTool)

    # Tracks associated with Muons
    if (TrackThinningKeepMuonTracks) :
        from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__MuonTrackParticleThinning
        EGAM7MuonTPThinningTool = DerivationFramework__MuonTrackParticleThinning( name                    = "EGAM7MuonTPThinningTool",
                                                                                  StreamName              = streamName,
                                                                                  MuonKey                 = "Muons",
                                                                                  InDetTrackParticlesKey  = "InDetTrackParticles")
        ToolSvc += EGAM7MuonTPThinningTool
        print(EGAM7MuonTPThinningTool)
        thinningTools.append(EGAM7MuonTPThinningTool)

    # Tracks associated with Taus
    if (TrackThinningKeepTauTracks) :
        from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__TauTrackParticleThinning
        EGAM7TauTPThinningTool = DerivationFramework__TauTrackParticleThinning( name                    = "EGAM7TauTPThinningTool",
                                                                                StreamName              = streamName,
                                                                                TauKey                  = "TauJets",
                                                                                ConeSize                = 0.6,
                                                                                InDetTrackParticlesKey  = "InDetTrackParticles")
        ToolSvc += EGAM7TauTPThinningTool
        print(EGAM7TauTPThinningTool)
        thinningTools.append(EGAM7TauTPThinningTool)

    # Tracks from primary vertex
    thinning_expression = "InDetTrackParticles.DFCommonTightPrimary && abs(DFCommonInDetTrackZ0AtPV)*sin(InDetTrackParticles.theta) < 3.0*mm && InDetTrackParticles.pt > 10*GeV"
    if (TrackThinningKeepPVTracks) :
        from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__TrackParticleThinning
        EGAM7TPThinningTool = DerivationFramework__TrackParticleThinning( name                    = "EGAM7TPThinningTool",
                                                                          StreamName              = streamName,
                                                                          SelectionString         = thinning_expression,
                                                                          InDetTrackParticlesKey  = "InDetTrackParticles")
        ToolSvc += EGAM7TPThinningTool
        print(EGAM7TPThinningTool)
        thinningTools.append(EGAM7TPThinningTool)

# Truth thinning
if DerivationFrameworkIsMonteCarlo:
    truth_cond_WZH = "((abs(TruthParticles.pdgId) >= 23) && (abs(TruthParticles.pdgId) <= 25))" # W, Z and Higgs
    truth_cond_lep = "((abs(TruthParticles.pdgId) >= 11) && (abs(TruthParticles.pdgId) <= 16))" # Leptons
    truth_cond_top = "((abs(TruthParticles.pdgId) ==  6))"                                     # Top quark
    truth_cond_gam = "((abs(TruthParticles.pdgId) == 22) && (TruthParticles.pt > 1*GeV))"       # Photon
    truth_cond_finalState = '(TruthParticles.status == 1 && TruthParticles.barcode < 200000)'   # stable particles
    truth_expression = '(' + truth_cond_WZH + ' ||  ' + truth_cond_lep +' || '+truth_cond_top +' || '+truth_cond_gam + ') || (' + truth_cond_finalState+')'

    from DerivationFrameworkMCTruth.DerivationFrameworkMCTruthConf import DerivationFramework__GenericTruthThinning
    EGAM7TruthThinningTool = DerivationFramework__GenericTruthThinning(name                    = "EGAM7TruthThinningTool",
                                                                       StreamName              = streamName,
                                                                       ParticleSelectionString = truth_expression,
                                                                       PreserveDescendants     = False,
                                                                       PreserveGeneratorDescendants     = True,
                                                                       PreserveAncestors      = True)
    ToolSvc += EGAM7TruthThinningTool
    thinningTools.append(EGAM7TruthThinningTool)


#============ Thin cells for EGAM7 =================================
# Keep only calo cells associated with the egammaClusters collection
#====================================================================
if thinCells:
    from DerivationFrameworkCalo.CaloCellDFGetter import thinCaloCellsForDF
    thinCaloCellsForDF (inputClusterKeys=["egammaClusters"],
                        streamName = EGAM7Stream.Name,
                        outputCellKey = "DFEGAMCellContainer")


#====================================================================
# Setup the skimming criteria
#====================================================================

#====================================================================
# offline-based selection (1 central electron with pT>4.5 GeV)
#====================================================================
requirement_object = 'Electrons.pt > 4.5*GeV'
objectSelection = 'count('+requirement_object+') >= 1'
from DerivationFrameworkTools.DerivationFrameworkToolsConf import DerivationFramework__xAODStringSkimmingTool
EGAM7_OfflineSkimmingTool = DerivationFramework__xAODStringSkimmingTool( name = "EGAM7_OfflineSkimmingTool",
                                                                         expression = objectSelection)
ToolSvc += EGAM7_OfflineSkimmingTool
print("EGAM7 offline skimming tool:", EGAM7_OfflineSkimmingTool)


#====================================================================
# trigger-based selection
# prescaled _etcut triggers
# prescaled _loose triggers
# prescaled _lhloose triggers
#====================================================================
triggers =  ['HLT_e4_etcut'        ]
triggers += ['HLT_e5_etcut'        ]
triggers += ['HLT_e9_etcut'        ]            
triggers += ['HLT_e10_etcut_L1EM7' ]            
triggers += ['HLT_e14_etcut'       ]            
triggers += ['HLT_e15_etcut_L1EM7' ]
triggers += ['HLT_e17_etcut_L1EM15']            
triggers += ['HLT_e20_etcut_L1EM12']            
triggers += ['HLT_e25_etcut_L1EM15']            
triggers += ['HLT_e30_etcut_L1EM15']            
triggers += ['HLT_e40_etcut_L1EM15']            
triggers += ['HLT_e50_etcut_L1EM15']            
triggers += ['HLT_e60_etcut'       ]            
triggers += ['HLT_e80_etcut'       ]            
triggers += ['HLT_e100_etcut'      ]            
triggers += ['HLT_e120_etcut'      ]            
triggers += ['HLT_g10_etcut'       ]            
triggers += ['HLT_g20_etcut_L1EM12']            
triggers += ['HLT_g200_etcut'      ]            

triggers += ['HLT_e5_lhloose'                      ]
triggers += ['HLT_e5_lhvloose'                     ]
triggers += ['HLT_e5_loose'                        ]
triggers += ['HLT_e5_vloose'                       ]
triggers += ['HLT_e10_lhvloose_L1EM7'              ]
triggers += ['HLT_e10_vloose_L1EM7'                ]
triggers += ['HLT_e12_lhloose'                     ]
triggers += ['HLT_e12_lhloose_L1EM10VH'            ]
triggers += ['HLT_e12_lhvloose_L1EM10VH'           ]
triggers += ['HLT_e12_loose'                       ]
triggers += ['HLT_e12_loose_L1EM10VH'              ]
triggers += ['HLT_e12_vloose_L1EM10VH'             ]
triggers += ['HLT_e15_lhloose_L1EM13VH'            ]
triggers += ['HLT_e15_lhvloose_L1EM13VH'           ]
triggers += ['HLT_e15_lhvloose_L1EM7'              ]
triggers += ['HLT_e15_loose_L1EM13VH'              ]
triggers += ['HLT_e15_vloose_L1EM13VH'             ]
triggers += ['HLT_e15_vloose_L1EM7'                ]
triggers += ['HLT_e17_lhloose'                     ]
triggers += ['HLT_e17_lhloose_L1EM15'              ]
triggers += ['HLT_e17_lhloose_cutd0dphideta_L1EM15']
triggers += ['HLT_e17_lhloose_nod0_L1EM15'         ]
triggers += ['HLT_e17_lhloose_nodeta_L1EM15'       ]
triggers += ['HLT_e17_lhloose_nodphires_L1EM15'    ]
triggers += ['HLT_e17_lhloose_L1EM15VHJJ1523ETA49' ]
triggers += ['HLT_e17_lhvloose'                    ]
triggers += ['HLT_e17_loose'                       ]
triggers += ['HLT_e17_loose_L1EM15'                ]
triggers += ['HLT_e17_loose_L1EM15VHJJ1523ETA49'   ]
triggers += ['HLT_e17_vloose'                      ]
triggers += ['HLT_e20_lhvloose'                    ]
triggers += ['HLT_e20_lhvloose_L1EM12'             ]
triggers += ['HLT_e20_vloose'                      ]
triggers += ['HLT_e20_vloose_L1EM12'               ]
triggers += ['HLT_e25_lhvloose_L1EM15'             ]
triggers += ['HLT_e25_vloose_L1EM15'               ]
triggers += ['HLT_e30_lhvloose_L1EM15'             ]
triggers += ['HLT_e30_vloose_L1EM15'               ]
triggers += ['HLT_e40_lhvloose'                    ]
triggers += ['HLT_e40_lhvloose_L1EM15'             ]
triggers += ['HLT_e40_vloose_L1EM15'               ]
triggers += ['HLT_e50_lhvloose_L1EM15'             ]
triggers += ['HLT_e50_vloose_L1EM15'               ]
triggers += ['HLT_e60_loose'                       ]
triggers += ['HLT_e60_vloose'                      ]
triggers += ['HLT_e60_lhvloose'                    ]
triggers += ['HLT_e70_etcut'                       ]
triggers += ['HLT_e70_lhloose'                     ]
triggers += ['HLT_e70_lhvloose'                    ]
triggers += ['HLT_e70_loose'                       ]
triggers += ['HLT_e70_vloose'                      ]
triggers += ['HLT_e80_lhvloose'                    ]
triggers += ['HLT_e80_vloose'                      ]
triggers += ['HLT_e100_lhvloose'                   ]
triggers += ['HLT_e100_vloose'                     ]
triggers += ['HLT_e120_lhvloose'                   ]
triggers += ['HLT_e120_lhloose'                    ]
triggers += ['HLT_e120_loose'                      ]
triggers += ['HLT_e120_vloose'                     ]
triggers += ['HLT_e140_etcut'                      ]
triggers += ['HLT_e160_etcut'                      ]
triggers += ['HLT_e180_etcut'                      ]
triggers += ['HLT_e200_etcut'                      ]
triggers += ['HLT_e250_etcut'                      ]
triggers += ['HLT_e300_etcut'                      ]
triggers += ['HLT_g250_etcut'                      ]
triggers += ['HLT_g300_etcut'                      ]

from DerivationFrameworkTools.DerivationFrameworkToolsConf import DerivationFramework__TriggerSkimmingTool
EGAM7_TriggerSkimmingTool = DerivationFramework__TriggerSkimmingTool(   name = "EGAM7_TriggerSkimmingTool", TriggerListOR = triggers)
ToolSvc += EGAM7_TriggerSkimmingTool
print("EGAM7 trigger skimming tool:", EGAM7_TriggerSkimmingTool)

#====================================================================
# make the AND of offline and trigger-based selections
#====================================================================
from DerivationFrameworkTools.DerivationFrameworkToolsConf import DerivationFramework__FilterCombinationAND
EGAM7_SkimmingTool = DerivationFramework__FilterCombinationAND(name="EGAM7SkimmingTool", FilterList=[EGAM7_OfflineSkimmingTool,EGAM7_TriggerSkimmingTool] )
ToolSvc+=EGAM7_SkimmingTool



#=======================================
# CREATE THE DERIVATION KERNEL ALGORITHM
#=======================================
from DerivationFrameworkCore.DerivationFrameworkCoreConf import DerivationFramework__DerivationKernel
print("EGAM7 skimming tools: ", [EGAM7_SkimmingTool])
print("EGAM7 thinning tools: ", thinningTools)
print("EGAM7 augmentation tools: ", augmentationTools)
EGAM7Sequence += CfgMgr.DerivationFramework__DerivationKernel("EGAM7Kernel",
                                                              AugmentationTools = augmentationTools,
                                                              SkimmingTools = [EGAM7_SkimmingTool],
                                                              ThinningTools = thinningTools
                                                              )


#====================================================================
# JET/MET
#====================================================================
jetList = [AntiKt4PV0Track]
if DerivationFrameworkIsMonteCarlo:
    jetList += [AntiKt4Truth, AntiKt4TruthDressedWZ]
addDAODJets(jetList, EGAM7Sequence)


#====================================================================
# FLAVOUR TAGGING
#====================================================================
#not available yet in r22
#from DerivationFrameworkFlavourTag.FtagRun3DerivationConfig import FtagJetCollections
#FtagJetCollections(['AntiKt4PV0TrackJets'], EGAM7Sequence)


#====================================================================
# SET UP SLIMMING
#====================================================================
from DerivationFrameworkCore.SlimmingHelper import SlimmingHelper
EGAM7SlimmingHelper = SlimmingHelper("EGAM7SlimmingHelper")
EGAM7SlimmingHelper.SmartCollections = ["Electrons",
                                        "Photons",
                                        "Muons",
                                        "TauJets",
                                        "MET_Baseline_AntiKt4EMPFlow",
                                        "AntiKt4EMPFlowJets",
                                        "BTagging_AntiKt4EMPFlow",
#                                        "BTagging_AntiKt4Track",
                                        "InDetTrackParticles",
                                        "PrimaryVertices" ]
if DerivationFrameworkIsMonteCarlo:
    EGAM7SlimmingHelper.SmartCollections += ["AntiKt4TruthJets",
                                             "AntiKt4TruthDressedWZJets"]

# Add egamma trigger objects
EGAM7SlimmingHelper.IncludeEGammaTriggerContent = True

# Extra variables
EGAM7SlimmingHelper.ExtraVariables = ExtraContentAll
EGAM7SlimmingHelper.AllVariables = ExtraContainersElectrons
EGAM7SlimmingHelper.AllVariables += ExtraContainersTrigger

# not available yet in R22
#from DerivationFrameworkFlavourTag.BTaggingContent import BTaggingStandardContent
#EGAM7SlimmingHelper.ExtraVariables.extend(BTaggingStandardContent("AntiKt4PV0TrackJets"))


if DerivationFrameworkIsMonteCarlo:
    EGAM7SlimmingHelper.ExtraVariables += ExtraContentAllTruth
    EGAM7SlimmingHelper.AllVariables += ExtraContainersTruth
else:
    EGAM7SlimmingHelper.ExtraVariables += ExtraContainersTriggerDataOnly

for tool in EGAM7_ClusterEnergyPerLayerDecorators:
    EGAM7SlimmingHelper.ExtraVariables.extend( getClusterEnergyPerLayerDecorations( tool ) )

# Add event info
if jobproperties.egammaDFFlags.doEGammaEventInfoSlimming:
    EGAM7SlimmingHelper.SmartCollections.append("EventInfo")
else:
    EGAM7SlimmingHelper.AllVariables += ["EventInfo"]

# Add detailed shower shape variables
from DerivationFrameworkEGamma.ElectronsCPDetailedContent import *
EGAM7SlimmingHelper.ExtraVariables += ElectronsCPDetailedContent
EGAM7SlimmingHelper.ExtraVariables += GSFTracksCPDetailedContent
from DerivationFrameworkEGamma.PhotonsCPDetailedContent import *
EGAM7SlimmingHelper.ExtraVariables += PhotonsCPDetailedContent

# This line must come after we have finished configuring EGAM7SlimmingHelper
EGAM7SlimmingHelper.AppendContentToStream(EGAM7Stream)

#Add CellContainer
if thinCells:
    EGAM7Stream.AddItem("CaloCellContainer#DFEGAMCellContainer")
else:
    EGAM7Stream.AddItem("CaloCellContainer#AllCalo")
# Add the cluster->cells links
EGAM7Stream.AddItem("CaloClusterCellLinkContainer#egammaClusters_links")
