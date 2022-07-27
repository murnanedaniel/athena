#********************************************************************
# EGAM6.py
# reductionConf flag EGAM6 in Reco_tf.py
# identical to EGAM1 but with loose electron cuts for start of run with poor ID/ECAL alignment
# author: giovanni.marchiori@cern.ch
#********************************************************************

from DerivationFrameworkCore.DerivationFrameworkMaster import buildFileName
from DerivationFrameworkCore.DerivationFrameworkMaster import (
    DerivationFrameworkIsMonteCarlo, DerivationFrameworkJob)
from DerivationFrameworkPhys import PhysCommon
from DerivationFrameworkEGamma.EGammaCommon import *
from DerivationFrameworkEGamma.EGAM1ExtraContent import *
from DerivationFrameworkJetEtMiss.JetCommon import addDAODJets
from DerivationFrameworkCore.DerivationFrameworkCoreConf import (
    DerivationFramework__DerivationKernel)
from DerivationFrameworkTools.DerivationFrameworkToolsConf import (
    DerivationFramework__xAODStringSkimmingTool)
from DerivationFrameworkCalo.DerivationFrameworkCaloFactories import (
    GainDecorator, getGainDecorations, getClusterEnergyPerLayerDecorator,
    getClusterEnergyPerLayerDecorations)
from DerivationFrameworkCalo.DerivationFrameworkCaloConf import (
    DerivationFramework__MaxCellDecorator)
from DerivationFrameworkEGamma.DerivationFrameworkEGammaConf import (
    DerivationFramework__EGInvariantMassTool)


#====================================================================
# read common DFEGamma settings from egammaDFFlags
#====================================================================
from DerivationFrameworkEGamma.egammaDFFlags import jobproperties
jobproperties.egammaDFFlags.print_JobProperties("full")

RecomputeElectronSelectors = True
#RecomputeElectronSelectors = False


#====================================================================
# check if we run on data or MC
#====================================================================
print("DerivationFrameworkIsMonteCarlo: ", DerivationFrameworkIsMonteCarlo)


#====================================================================
# Set up sequence for this format and add to the top sequence
#====================================================================
EGAM6Sequence = CfgMgr.AthSequencer("EGAM6Sequence")
DerivationFrameworkJob += EGAM6Sequence


#====================================================================
# SET UP STREAM (to be done early in the game to set up thinning Svc
#====================================================================
streamName = derivationFlags.WriteDAOD_EGAM6Stream.StreamName
fileName   = buildFileName( derivationFlags.WriteDAOD_EGAM6Stream )
EGAM6Stream = MSMgr.NewPoolRootStream( streamName, fileName )
# Only events that pass the filters listed below are written out.
# Name must match that of the kernel above
# AcceptAlgs  = logical OR of filters
# RequireAlgs = logical AND of filters
EGAM6Stream.AcceptAlgs(["EGAM6Kernel"])


### Thinning and augmentation tools lists
augmentationTools = []
thinningTools=[]


#====================================================================
# SET UP SKIMMING
#====================================================================


# SELECTION FOR CALIBRATION

#====================================================================
# Z->ee selection based on single e trigger:
# 1 tight e, central, pT>25 GeV
# 1 medium e, pT>20 GeV
# OS, mee>60 GeV
#====================================================================

# switch to likelihood selectors only as soon as they're commissioned (and used in trigger)

if RecomputeElectronSelectors :
    requirement_tag = '(Electrons.DFCommonElectronsLHTight) && Electrons.pt > 24.5*GeV'
    requirement_probe = '(Electrons.DFCommonElectronsLHMedium) && Electrons.pt > 19.5*GeV'
else :
    requirement_tag = '(Electrons.LHTight) && Electrons.pt > 24.5*GeV'
    requirement_probe = '(Electrons.LHMedium) && Electrons.pt > 19.5*GeV'

EGAM6_ZEEMassTool1 = DerivationFramework__EGInvariantMassTool(
    name = "EGAM6_ZEEMassTool1",
    Object1Requirements = requirement_tag,
    Object2Requirements = requirement_probe,
    StoreGateEntryName = "EGAM6_DiElectronMass",
    Mass1Hypothesis = 0.511*MeV,
    Mass2Hypothesis = 0.511*MeV,
    Container1Name = "Electrons",
    Container2Name = "Electrons",
    CheckCharge = True,
    DoTransverseMass = False,
    MinDeltaR = 0.0)

ToolSvc += EGAM6_ZEEMassTool1
augmentationTools += [EGAM6_ZEEMassTool1]
print(EGAM6_ZEEMassTool1)


#====================================================================
# Z->ee selection based on di-electron trigger
# 2 medium e, central, pT>20 GeV
# OS, mee>60 GeV
#====================================================================

if RecomputeElectronSelectors:
    requirement = '(Electrons.DFCommonElectronsLHLoose || Electrons.DFCommonElectronsLHMedium) && Electrons.pt > 19.5*GeV'
else:
    requirement = '(Electrons.LHLoose || Electrons.LHMedium) && Electrons.pt > 19.5*GeV'

EGAM6_ZEEMassTool2 = DerivationFramework__EGInvariantMassTool(
    name = "EGAM6_ZEEMassTool2",
    Object1Requirements = requirement,
    Object2Requirements = requirement,
    StoreGateEntryName = "EGAM6_DiElectronMass2",
    Mass1Hypothesis = 0.511*MeV,
    Mass2Hypothesis = 0.511*MeV,
    Container1Name = "Electrons",
    Container2Name = "Electrons",
    CheckCharge = True,
    DoTransverseMass = False,
    MinDeltaR = 0.0)

ToolSvc += EGAM6_ZEEMassTool2
augmentationTools += [EGAM6_ZEEMassTool2]
print(EGAM6_ZEEMassTool2)


# SELECTION FOR T&P

#====================================================================
# Z->ee selection based on single e trigger, for reco (central) and ID SF(central+fwd)
# 1 tight e, central, pT>25 GeV
# 1 e, pT>15 GeV if central, >20 GeV if forward
# OS+SS, mee>60 GeV
#====================================================================

# switch to likelihood selectors only as soon as they're commissioned (and used in trigger)
if RecomputeElectronSelectors :
    requirement_tag = '(Electrons.DFCommonElectronsLHLoose || Electrons.DFCommonElectronsLHMedium) && Electrons.pt > 24.5*GeV'
else :
    requirement_tag = '(Electrons.LHLoose || Electrons.LHMedium) && Electrons.pt > 24.5*GeV'

# central electrons: collection = Electrons, pt>14.5 GeV
requirement_probe = 'Electrons.pt > 14.5*GeV'

EGAM6_ZEEMassTool3 = DerivationFramework__EGInvariantMassTool(
    name = "EGAM6_ZEEMassTool3",
    Object1Requirements = requirement_tag,
    Object2Requirements = requirement_probe,
    StoreGateEntryName = "EGAM6_DiElectronMass3",
    Mass1Hypothesis = 0.511*MeV,
    Mass2Hypothesis = 0.511*MeV,
    Container1Name = "Electrons",
    Container2Name = "Electrons",
    CheckCharge = False,
    DoTransverseMass = False,
    MinDeltaR = 0.0)

ToolSvc += EGAM6_ZEEMassTool3
augmentationTools += [EGAM6_ZEEMassTool3]
print(EGAM6_ZEEMassTool3)


#====================================================================
# Z->eg selection based on single e trigger, for reco SF (central)
# 1 tight e, central, pT>25 GeV
# 1 gamma, pT>15 GeV, central
# OS+SS, mee>60 GeV
#====================================================================

# switch to likelihood selectors only as soon as they're commissioned (and used in trigger)
if RecomputeElectronSelectors:
    requirement_tag = '(Electrons.DFCommonElectronsLHLoose || Electrons.DFCommonElectronsLHMedium) && Electrons.pt > 24.5*GeV'
else:
    requirement_tag = '(Electrons.LHLoose || Electrons.LHMedium) && Electrons.pt > 24.5*GeV'
requirement_probe = 'DFCommonPhotons_et > 14.5*GeV'

EGAM6_ZEGMassTool = DerivationFramework__EGInvariantMassTool(
    name = "EGAM6_ZEGMassTool",
    Object1Requirements = requirement_tag,
    Object2Requirements = requirement_probe,
    StoreGateEntryName = "EGAM6_ElectronPhotonMass",
    Mass1Hypothesis = 0.511*MeV,
    Mass2Hypothesis = 0.511*MeV,
    Container1Name = "Electrons",
    Container2Name = "Photons",
    Pt2BranchName = "DFCommonPhotons_et",
    Eta2BranchName = "DFCommonPhotons_eta",
    Phi2BranchName = "DFCommonPhotons_phi",
    CheckCharge = False,
    DoTransverseMass = False,
    MinDeltaR = 0.0)

ToolSvc += EGAM6_ZEGMassTool
augmentationTools += [EGAM6_ZEGMassTool]
print(EGAM6_ZEGMassTool)


# Skimming criteria
expression = 'count(EGAM6_DiElectronMass > 60.0*GeV)>=1 || count(EGAM6_DiElectronMass2 > 60.0*GeV)>=1 || count(EGAM6_DiElectronMass3 > 60.0*GeV)>=1 ||  count (EGAM6_ElectronPhotonMass > 60.0*GeV)>=1'
EGAM6_SkimmingTool = DerivationFramework__xAODStringSkimmingTool( name = "EGAM6_SkimmingTool",
                                                                 expression = expression)
ToolSvc += EGAM6_SkimmingTool
print("EGAM6 skimming tool:", EGAM6_SkimmingTool)



#====================================================================
# DECORATION TOOLS
#====================================================================


#====================================================================
# Gain and cluster energies per layer decoration tool
#====================================================================
EGAM6_GainDecoratorTool = GainDecorator()
ToolSvc += EGAM6_GainDecoratorTool
augmentationTools += [EGAM6_GainDecoratorTool]

cluster_sizes = (3,5), (5,7), (7,7), (7,11)
EGAM6_ClusterEnergyPerLayerDecorators = [getClusterEnergyPerLayerDecorator(neta, nphi)() for neta, nphi in cluster_sizes]
augmentationTools += EGAM6_ClusterEnergyPerLayerDecorators


#====================================================================
# Max Cell sum decoration tool
#====================================================================                                                        
EGAM6_MaxCellDecoratorTool = DerivationFramework__MaxCellDecorator(
    name                    = "EGAM6_MaxCellDecoratorTool",
    SGKey_electrons         = "Electrons",
    SGKey_photons           = "Photons",
)
ToolSvc += EGAM6_MaxCellDecoratorTool
augmentationTools += [EGAM6_MaxCellDecoratorTool]


#====================================================================
# SET UP THINNING
#====================================================================

print('WARNING, Thinning of trigger navigation has to be properly implemented in R22')
#from DerivationFrameworkCore.ThinningHelper import ThinningHelper
#EGAM6ThinningHelper = ThinningHelper( "EGAM6ThinningHelper" )
#EGAM6ThinningHelper.TriggerChains = '(^(?!.*_[0-9]*(mu|j|xe|tau|ht|xs|te))(?!HLT_[eg].*_[0-9]*[eg][0-9].*)(?!HLT_eb.*)(?!.*larpeb.*)(?!HLT_.*_AFP_.*)(HLT_[eg].*))|HLT_e.*_Zee.*'
#EGAM6ThinningHelper.AppendToStream( EGAM6Stream, ExtraContainersTrigger )

if jobproperties.egammaDFFlags.doEGammaDAODTrackThinning:

    TrackThinningKeepElectronTracks = True
    TrackThinningKeepPhotonTracks = True
    TrackThinningKeepAllElectronTracks = True
    TrackThinningKeepJetTracks = False
    TrackThinningKeepMuonTracks = False
    TrackThinningKeepTauTracks = False
    TrackThinningKeepPVTracks = True

    # Tracks associated with Electrons
    if (TrackThinningKeepElectronTracks) : 
        from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__EgammaTrackParticleThinning
        EGAM6ElectronTPThinningTool = DerivationFramework__EgammaTrackParticleThinning( name                    = "EGAM6ElectronTPThinningTool",
                                                                                        StreamName              = streamName,
                                                                                        SGKey                   = "Electrons",
                                                                                        GSFTrackParticlesKey    = "GSFTrackParticles",        
                                                                                        InDetTrackParticlesKey  = "InDetTrackParticles",
                                                                                        SelectionString         = "Electrons.pt > 0*GeV",
                                                                                        BestMatchOnly = True,
                                                                                        ConeSize = 0.3)
        ToolSvc += EGAM6ElectronTPThinningTool
        print(EGAM6ElectronTPThinningTool)
        thinningTools.append(EGAM6ElectronTPThinningTool)

    # Tracks associated with Electrons (all tracks, large cone, for track isolation studies of the selected electrons)
    if (TrackThinningKeepAllElectronTracks) :
        EGAM6ElectronTPThinningTool2 = DerivationFramework__EgammaTrackParticleThinning( name                    = "EGAM6ElectronTPThinningTool2",
                                                                                         StreamName              = streamName,
                                                                                         SGKey                   = "Electrons",
                                                                                         GSFTrackParticlesKey    = "GSFTrackParticles",        
                                                                                         InDetTrackParticlesKey  = "InDetTrackParticles",
                                                                                         SelectionString         = "Electrons.pt > 4*GeV",
                                                                                         BestMatchOnly = False,
                                                                                         ConeSize = 0.6)

        ToolSvc += EGAM6ElectronTPThinningTool2
        print(EGAM6ElectronTPThinningTool2)
        thinningTools.append(EGAM6ElectronTPThinningTool2)

    # Tracks associated with Photons
    if (TrackThinningKeepPhotonTracks) : 
        from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__EgammaTrackParticleThinning
        EGAM6PhotonTPThinningTool = DerivationFramework__EgammaTrackParticleThinning( name                    = "EGAM6PhotonTPThinningTool",
                                                                                      StreamName              = streamName,
                                                                                      SGKey                   = "Photons",
                                                                                     GSFTrackParticlesKey    = "GSFTrackParticles",        
                                                                                      InDetTrackParticlesKey  = "InDetTrackParticles",
                                                                                      SelectionString         = "Photons.pt > 0*GeV",
                                                                                      BestMatchOnly = True,
                                                                                      ConeSize = 0.3)
        ToolSvc += EGAM6PhotonTPThinningTool
        print(EGAM6PhotonTPThinningTool)
        thinningTools.append(EGAM6PhotonTPThinningTool)

    # Tracks associated with Jets
    if (TrackThinningKeepJetTracks) : 
        from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__JetTrackParticleThinning
        EGAM6JetTPThinningTool = DerivationFramework__JetTrackParticleThinning( name                    = "EGAM6JetTPThinningTool",
                                                                                StreamName              = streamName,
                                                                                JetKey                  = "AntiKt4EMPFlowJets",
                                                                                InDetTrackParticlesKey  = "InDetTrackParticles")
        ToolSvc += EGAM6JetTPThinningTool
        print(EGAM6JetTPThinningTool)
        thinningTools.append(EGAM6JetTPThinningTool)


    # Tracks associated with Muons
    if (TrackThinningKeepMuonTracks) : 
        from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__MuonTrackParticleThinning
        EGAM6MuonTPThinningTool = DerivationFramework__MuonTrackParticleThinning( name                    = "EGAM6MuonTPThinningTool",
                                                                                  StreamName              = streamName,
                                                                                  MuonKey                 = "Muons",
                                                                                  InDetTrackParticlesKey  = "InDetTrackParticles")
        ToolSvc += EGAM6MuonTPThinningTool
        print(EGAM6MuonTPThinningTool)
        thinningTools.append(EGAM6MuonTPThinningTool)
        

    # Tracks associated with Taus
    if (TrackThinningKeepTauTracks) : 
        from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__TauTrackParticleThinning
        EGAM6TauTPThinningTool = DerivationFramework__TauTrackParticleThinning( name                    = "EGAM6TauTPThinningTool",
                                                                                StreamName              = streamName,
                                                                                TauKey                  = "TauJets",
                                                                                ConeSize                = 0.6,
                                                                                InDetTrackParticlesKey  = "InDetTrackParticles")
        ToolSvc += EGAM6TauTPThinningTool
        print(EGAM6TauTPThinningTool)
        thinningTools.append(EGAM6TauTPThinningTool)

    # Tracks from primary vertex
    if (TrackThinningKeepPVTracks) : 
        from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__TrackParticleThinning
        EGAM6TPThinningTool = DerivationFramework__TrackParticleThinning( name                    = "EGAM6TPThinningTool",
                                                                          StreamName              = streamName,
                                                                          SelectionString         = "InDetTrackParticles.DFCommonTightPrimary && abs( DFCommonInDetTrackZ0AtPV * sin(InDetTrackParticles.theta)) < 3.0*mm",
                                                                          InDetTrackParticlesKey  = "InDetTrackParticles")
        ToolSvc += EGAM6TPThinningTool
        print(EGAM6TPThinningTool)
        thinningTools.append(EGAM6TPThinningTool)


# Truth thinning
if DerivationFrameworkIsMonteCarlo:
    truth_cond_WZH = "((abs(TruthParticles.pdgId) >= 23) && (abs(TruthParticles.pdgId) <= 25))" # W, Z and Higgs
    truth_cond_lep = "((abs(TruthParticles.pdgId) >= 11) && (abs(TruthParticles.pdgId) <= 16))" # Leptons
    truth_cond_top = "((abs(TruthParticles.pdgId) ==  6))"                                     # Top quark
    truth_cond_gam = "((abs(TruthParticles.pdgId) == 22) && (TruthParticles.pt > 1*GeV))"       # Photon
    truth_cond_finalState = '(TruthParticles.status == 1 && TruthParticles.barcode < 200000)'   # stable particles
    truth_expression = '(' + truth_cond_WZH + ' ||  ' + truth_cond_lep +' || '+truth_cond_top +' || '+truth_cond_gam + ') || (' + truth_cond_finalState+')'
    
    from DerivationFrameworkMCTruth.DerivationFrameworkMCTruthConf import DerivationFramework__GenericTruthThinning
    EGAM6TruthThinningTool = DerivationFramework__GenericTruthThinning(name                    = "EGAM6TruthThinningTool",
                                                                       StreamName              = streamName,
                                                                       ParticleSelectionString = truth_expression,
                                                                       PreserveDescendants     = False,
                                                                       PreserveGeneratorDescendants     = True,
                                                                       PreserveAncestors      = True)


    ToolSvc += EGAM6TruthThinningTool
    thinningTools.append(EGAM6TruthThinningTool)
    
print("EGAM6 thinningTools: ", thinningTools)



#=======================================
# CREATE THE DERIVATION KERNEL ALGORITHM
#=======================================
print("EGAM6 skimming tools: ", [EGAM6_SkimmingTool])
print("EGAM6 thinning tools: ", thinningTools)
print("EGAM6 augmentation tools: ", augmentationTools)
EGAM6Sequence += CfgMgr.DerivationFramework__DerivationKernel(
    "EGAM6Kernel",
    AugmentationTools = augmentationTools,
    SkimmingTools = [EGAM6_SkimmingTool],
    ThinningTools = thinningTools
)


#====================================================================
# JET/MET
#====================================================================
from JetRecConfig.StandardSmallRJets import AntiKt4Truth,AntiKt4TruthDressedWZ
jetList = []
if DerivationFrameworkIsMonteCarlo:
    jetList += [AntiKt4Truth, AntiKt4TruthDressedWZ]
addDAODJets(jetList, EGAM6Sequence)


#====================================================================
# FLAVOUR TAGGING   
#====================================================================
from DerivationFrameworkFlavourTag.FtagRun3DerivationConfig import FtagJetCollection
FtagJetCollection('AntiKt4EMPFlowJets',EGAM6Sequence)


#====================================================================
# SET UP SLIMMING
#====================================================================
from DerivationFrameworkCore.SlimmingHelper import SlimmingHelper
EGAM6SlimmingHelper = SlimmingHelper("EGAM6SlimmingHelper")

EGAM6SlimmingHelper.SmartCollections = [
                                        "Electrons",
                                        "Photons",
                                        "Muons",
                                        "TauJets",
                                        "MET_Baseline_AntiKt4EMPFlow",
                                        "AntiKt4EMPFlowJets",
                                        "BTagging_AntiKt4EMPFlow",
                                        "InDetTrackParticles",
                                        "PrimaryVertices"
                                        ]

if DerivationFrameworkIsMonteCarlo:
    EGAM6SlimmingHelper.SmartCollections += ["AntiKt4TruthJets",
                                             "AntiKt4TruthDressedWZJets"]
	
# Add egamma trigger objects
EGAM6SlimmingHelper.IncludeEGammaTriggerContent = True

# read list of extra content from EGAM1 file (output of EGAM6 and EGAM1 is the same)
EGAM6SlimmingHelper.ExtraVariables = ExtraContentAll
EGAM6SlimmingHelper.AllVariables = ExtraContainersElectrons
EGAM6SlimmingHelper.AllVariables += ExtraContainersTrigger
if not DerivationFrameworkIsMonteCarlo:
    EGAM6SlimmingHelper.AllVariables += ExtraContainersTriggerDataOnly

if DerivationFrameworkIsMonteCarlo:
    EGAM6SlimmingHelper.ExtraVariables += ExtraContentAllTruth
    EGAM6SlimmingHelper.AllVariables += ExtraContainersTruth


for tool in EGAM6_ClusterEnergyPerLayerDecorators:
    EGAM6SlimmingHelper.ExtraVariables.extend( getClusterEnergyPerLayerDecorations( tool ) )

# Add event info
if jobproperties.egammaDFFlags.doEGammaEventInfoSlimming:
    EGAM6SlimmingHelper.SmartCollections.append("EventInfo")
else:
    EGAM6SlimmingHelper.AllVariables += ["EventInfo"]

# This line must come after we have finished configuring EGAM6SlimmingHelper
EGAM6SlimmingHelper.AppendContentToStream(EGAM6Stream)

# Add Derived Egamma CellContainer
# from DerivationFrameworkEGamma.EGammaCellCommon import CellCommonThinning
# CellCommonThinning(EGAM6Stream)

#Add full CellContainer
EGAM6Stream.AddItem("CaloCellContainer#AllCalo")
EGAM6Stream.AddItem("CaloClusterCellLinkContainer#egammaClusters_links")
