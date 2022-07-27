# ********************************************************************
# EGAM2.py
# J/psi->ee derivation for calibration
# reductionConf flag EGAM2 in Reco_tf.py
# author: giovanni.marchiori@cern.ch
# ********************************************************************

from DerivationFrameworkEGamma.PhotonsCPDetailedContent import *
from DerivationFrameworkEGamma.ElectronsCPDetailedContent import *
from DerivationFrameworkCore.SlimmingHelper import SlimmingHelper
from DerivationFrameworkJetEtMiss.JetCommon import addDAODJets
from JetRecConfig.StandardSmallRJets import AntiKt4Truth, AntiKt4TruthDressedWZ
from DerivationFrameworkCore.DerivationFrameworkCoreConf import (
    DerivationFramework__DerivationKernel)
from DerivationFrameworkTools.DerivationFrameworkToolsConf import (
    DerivationFramework__FilterCombinationOR)
from DerivationFrameworkTools.DerivationFrameworkToolsConf import (
    DerivationFramework__TriggerSkimmingTool)
from DerivationFrameworkTools.DerivationFrameworkToolsConf import (
    DerivationFramework__xAODStringSkimmingTool)
from DerivationFrameworkCalo.DerivationFrameworkCaloFactories import (
    GainDecorator, getGainDecorations,
    getClusterEnergyPerLayerDecorator, getClusterEnergyPerLayerDecorations)
from DerivationFrameworkEGamma.DerivationFrameworkEGammaConf import (
    DerivationFramework__EGInvariantMassTool)
from DerivationFrameworkCore.DerivationFrameworkMaster import buildFileName
from DerivationFrameworkCore.DerivationFrameworkMaster import (
    DerivationFrameworkIsMonteCarlo, DerivationFrameworkJob)
from DerivationFrameworkPhys import PhysCommon
from DerivationFrameworkEGamma.EGammaCommon import *
from DerivationFrameworkEGamma.EGAM2ExtraContent import *


# ====================================================================
# read common DFEGamma settings from egammaDFFlags
# ====================================================================
from DerivationFrameworkEGamma.egammaDFFlags import jobproperties
jobproperties.egammaDFFlags.print_JobProperties("full")

# this could also go in egammaDFFlags
RecomputeElectronSelectors = True
#RecomputeElectronSelectors = False


# ====================================================================
# check if we run on data or MC
# ====================================================================
print("DerivationFrameworkIsMonteCarlo: ", DerivationFrameworkIsMonteCarlo)


# ====================================================================
# Set up sequence for this format and add to the top sequence
# ====================================================================
EGAM2Sequence = CfgMgr.AthSequencer("EGAM2Sequence")
DerivationFrameworkJob += EGAM2Sequence


# ====================================================================
# SET UP STREAM
# ====================================================================
streamName = derivationFlags.WriteDAOD_EGAM2Stream.StreamName
fileName = buildFileName(derivationFlags.WriteDAOD_EGAM2Stream)
EGAM2Stream = MSMgr.NewPoolRootStream(streamName, fileName)
# Only events that pass the filters listed below are written out.
# Name must match that of the kernel above
# AcceptAlgs  = logical OR of filters
# RequireAlgs = logical AND of filters
EGAM2Stream.AcceptAlgs(["EGAM2Kernel"])


# Thinning and augmentation tools lists
augmentationTools = []
thinningTools = []


# ====================================================================
# SET UP AUGMENTATIONS
# ====================================================================


# ====================================================================
# 1. di-electron invariant mass for events passing the J/psi->ee
#    selection for e/gamma calibration
#
#    2 tight or medium e (depends on Run2 triggers..), pT>4.5 GeV, OS
#    1<Mee<5 GeV (applied in skimming step later)
# ====================================================================
electronPtRequirement = '(Electrons.pt > 4.5*GeV)'
if RecomputeElectronSelectors:
    electronQualityRequirement = '(Electrons.DFCommonElectronsLHMedium)'
else:
    electronQualityRequirement = '(Electrons.LHMedium)'
requirement_el = '(' + electronQualityRequirement + \
    '&&' + electronPtRequirement + ')'


EGAM2_JPSIEEMassTool = DerivationFramework__EGInvariantMassTool(
    name="EGAM2_JPSIEEMassTool",
    Object1Requirements=requirement_el,
    Object2Requirements=requirement_el,
    StoreGateEntryName="EGAM2_DiElectronMass",
    Mass1Hypothesis=0.511*MeV,
    Mass2Hypothesis=0.511*MeV,
    Container1Name="Electrons",
    Container2Name="Electrons",
    CheckCharge=True,
    DoTransverseMass=False,
    MinDeltaR=0.0)
ToolSvc += EGAM2_JPSIEEMassTool
augmentationTools += [EGAM2_JPSIEEMassTool]
print(EGAM2_JPSIEEMassTool)

# ====================================================================
# 2. di-electron invariant mass for events passing the J/psi->ee
#    selection based on Jpsi->e+cluster trigger, for low pT (7-20 GeV)
#    central photon efficiencies with tag and probe
#
# Tag: 1 tight e, central, pT>4.5 GeV
# Probe: 1 e, central, pT>4.5 GeV
# OS+SS
# dR>0.15
# 1<mee<6 GeV (applied in skimming step later)
# ====================================================================
if RecomputeElectronSelectors:
    requirement_el_tag = '(Electrons.DFCommonElectronsLHTight) && (Electrons.pt > 4.5*GeV)'
else:
    requirement_el_tag = '(Electrons.LHTight) && (Electrons.pt > 4.5*GeV)'
requirement_el_probe = 'Electrons.pt > 4.5*GeV'

EGAM2_JPSIEEMassTool2 = DerivationFramework__EGInvariantMassTool(
    name="EGAM2_JPSIEEMassTool2",
    Object1Requirements=requirement_el_tag,
    Object2Requirements=requirement_el_probe,
    StoreGateEntryName="EGAM2_DiElectronMass2",
    Mass1Hypothesis=0.511*MeV,
    Mass2Hypothesis=0.511*MeV,
    Container1Name="Electrons",
    Container2Name="Electrons",
    CheckCharge=False,
    DoTransverseMass=False,
    MinDeltaR=0.15)
ToolSvc += EGAM2_JPSIEEMassTool2
augmentationTools += [EGAM2_JPSIEEMassTool2]
print(EGAM2_JPSIEEMassTool2)


# ====================================================================
# Gain and cluster energies per layer decoration tool
# ====================================================================
EGAM2_GainDecoratorTool = GainDecorator()
ToolSvc += EGAM2_GainDecoratorTool
augmentationTools += [EGAM2_GainDecoratorTool]

cluster_sizes = (3, 7), (5, 5), (7, 11)
EGAM2_ClusterEnergyPerLayerDecorators = [getClusterEnergyPerLayerDecorator(
    neta, nphi)() for neta, nphi in cluster_sizes]
augmentationTools += EGAM2_ClusterEnergyPerLayerDecorators


# ====================================================================
# SET UP THINNING
# ====================================================================

print('WARNING, Thinning of trigger navigation has to be properly implemented in R22')
#from DerivationFrameworkCore.ThinningHelper import ThinningHelper
#EGAM2ThinningHelper = ThinningHelper( "EGAM2ThinningHelper" )
#EGAM2ThinningHelper.TriggerChains = '(^(?!.*_[0-9]*(mu|j|xe|tau|ht|xs|te))(?!HLT_[eg].*_[0-9]*[eg][0-9].*)(?!HLT_eb.*)(?!.*larpeb.*)(?!HLT_.*_AFP_.*)(HLT_[eg].*))|HLT_e.*_Jpsiee.*'
#EGAM2ThinningHelper.AppendToStream( EGAM2Stream, ExtraContainersTrigger )


# Track thinning
if jobproperties.egammaDFFlags.doEGammaDAODTrackThinning:

    TrackThinningKeepElectronTracks = True
    TrackThinningKeepPhotonTracks = True
    TrackThinningKeepJetTracks = False
    TrackThinningKeepMuonTracks = False
    TrackThinningKeepTauTracks = False
    TrackThinningKeepPVTracks = False

    # Tracks associated with Jets
    if (TrackThinningKeepJetTracks):
        from DerivationFrameworkInDet.DerivationFrameworkInDetConf import (
            DerivationFramework__JetTrackParticleThinning)
        EGAM2JetTPThinningTool = DerivationFramework__JetTrackParticleThinning(
            name="EGAM2JetTPThinningTool",
            StreamName=streamName,
            JetKey="AntiKt4EMPFlowJets",
            InDetTrackParticlesKey="InDetTrackParticles")
        ToolSvc += EGAM2JetTPThinningTool
        print(EGAM2JetTPThinningTool)
        thinningTools.append(EGAM2JetTPThinningTool)

    # Tracks associated with Muons
    if (TrackThinningKeepMuonTracks):
        from DerivationFrameworkInDet.DerivationFrameworkInDetConf import (
            DerivationFramework__MuonTrackParticleThinning)
        EGAM2MuonTPThinningTool = DerivationFramework__MuonTrackParticleThinning(
            name="EGAM2MuonTPThinningTool",
            StreamName=streamName,
            MuonKey="Muons",
            InDetTrackParticlesKey="InDetTrackParticles")
        ToolSvc += EGAM2MuonTPThinningTool
        print(EGAM2MuonTPThinningTool)
        thinningTools.append(EGAM2MuonTPThinningTool)

    # Tracks associated with Electrons
    if (TrackThinningKeepElectronTracks):
        from DerivationFrameworkInDet.DerivationFrameworkInDetConf import (
            DerivationFramework__EgammaTrackParticleThinning)
        EGAM2ElectronTPThinningTool = DerivationFramework__EgammaTrackParticleThinning(
            name="EGAM2ElectronTPThinningTool",
            StreamName=streamName,
            SGKey="Electrons",
            GSFTrackParticlesKey="GSFTrackParticles",
            InDetTrackParticlesKey="InDetTrackParticles",
            SelectionString="Electrons.pt > 0*GeV",
            BestMatchOnly=True,
            ConeSize=0.3)
        ToolSvc += EGAM2ElectronTPThinningTool
        print(EGAM2ElectronTPThinningTool)
        thinningTools.append(EGAM2ElectronTPThinningTool)

    # Tracks associated with Photons
    if (TrackThinningKeepPhotonTracks):
        from DerivationFrameworkInDet.DerivationFrameworkInDetConf import (
            DerivationFramework__EgammaTrackParticleThinning)
        EGAM2PhotonTPThinningTool = DerivationFramework__EgammaTrackParticleThinning(
            name="EGAM2PhotonTPThinningTool",
            StreamName=streamName,
            SGKey="Photons",
            GSFTrackParticlesKey="GSFTrackParticles",
            InDetTrackParticlesKey="InDetTrackParticles",
            SelectionString="Photons.pt > 0*GeV",
            BestMatchOnly=True,
            ConeSize=0.3)

        ToolSvc += EGAM2PhotonTPThinningTool
        print(EGAM2PhotonTPThinningTool)
        thinningTools.append(EGAM2PhotonTPThinningTool)

    # Tracks associated with Taus
    if (TrackThinningKeepTauTracks):
        from DerivationFrameworkInDet.DerivationFrameworkInDetConf import (
            DerivationFramework__TauTrackParticleThinning)
        EGAM2TauTPThinningTool = DerivationFramework__TauTrackParticleThinning(
            name="EGAM2TauTPThinningTool",
            StreamName=streamName,
            TauKey="TauJets",
            ConeSize=0.6,
            InDetTrackParticlesKey="InDetTrackParticles")
        ToolSvc += EGAM2TauTPThinningTool
        print(EGAM2TauTPThinningTool)
        thinningTools.append(EGAM2TauTPThinningTool)

    # Tracks from primary vertex
    thinning_expression = "InDetTrackParticles.DFCommonTightPrimary && abs(DFCommonInDetTrackZ0AtPV)*sin(InDetTrackParticles.theta) < 3.0*mm && InDetTrackParticles.pt > 10*GeV"
    if (TrackThinningKeepPVTracks):
        from DerivationFrameworkInDet.DerivationFrameworkInDetConf import (
            DerivationFramework__TrackParticleThinning)
        EGAM2TPThinningTool = DerivationFramework__TrackParticleThinning(
            name="EGAM2TPThinningTool",
            StreamName=streamName,
            SelectionString=thinning_expression,
            InDetTrackParticlesKey="InDetTrackParticles")
        ToolSvc += EGAM2TPThinningTool
        print(EGAM2TPThinningTool)
        thinningTools.append(EGAM2TPThinningTool)


# ====================================================================
# Setup the skimming criteria
# ====================================================================

# off-line based selection
expression_calib = '(count(EGAM2_DiElectronMass > 1.0*GeV && EGAM2_DiElectronMass < 5.0*GeV)>=1)'
expression_TP = '(count(EGAM2_DiElectronMass2 > 1.0*GeV && EGAM2_DiElectronMass2 < 6.0*GeV)>=1)'
expression = expression_calib + ' || ' + expression_TP
EGAM2_OfflineSkimmingTool = DerivationFramework__xAODStringSkimmingTool(
    name="EGAM2_OfflineSkimmingTool",
    expression=expression)
ToolSvc += EGAM2_OfflineSkimmingTool
print("EGAM2 offline skimming tool:", EGAM2_OfflineSkimmingTool)

# trigger-based selection
triggers = ['HLT_e5_lhtight_e4_etcut_Jpsiee']
triggers += ['HLT_e5_lhtight_nod0_e4_etcut_Jpsiee']
triggers += ['HLT_e5_lhtight_e4_etcut']
triggers += ['HLT_e5_lhtight_nod0_e4_etcut']

triggers += ['HLT_e9_lhtight_e4_etcut_Jpsiee']
triggers += ['HLT_e9_lhtight_nod0_e4_etcut_Jpsiee']
triggers += ['HLT_e9_etcut_e5_lhtight_nod0_Jpsiee']
triggers += ['HLT_e9_etcut_e5_lhtight_Jpsiee']

triggers += ['HLT_e14_etcut_e5_lhtight_Jpsiee']
triggers += ['HLT_e14_etcut_e5_lhtight_nod0_Jpsiee']
triggers += ['HLT_e14_lhtight_e4_etcut_Jpsiee']
triggers += ['HLT_e14_lhtight_nod0_e4_etcut_Jpsiee']

triggers += ['HLT_e5_lhtight_nod0_e4_etcut_Jpsiee_L1RD0_FILLED']
triggers += ['HLT_e5_lhtight_nod0_e9_etcut_Jpsiee']
triggers += ['HLT_e5_lhtight_nod0_e14_etcut_Jpsiee']
triggers += ['HLT_e5_lhtight_nod0_e9_etcut_Jpsiee_L1JPSI-1M5-EM7']
triggers += ['HLT_e9_lhtight_nod0_e4_etcut_Jpsiee_L1JPSI-1M5-EM7']
triggers += ['HLT_e5_lhtight_nod0_e14_etcut_Jpsiee_L1JPSI-1M5-EM12']
triggers += ['HLT_e14_lhtight_nod0_e4_etcut_Jpsiee_L1JPSI-1M5-EM12']

EGAM2_TriggerSkimmingTool = DerivationFramework__TriggerSkimmingTool(
    name="EGAM2_TriggerSkimmingTool",
    TriggerListOR=triggers)
ToolSvc += EGAM2_TriggerSkimmingTool
print("EGAM2 trigger skimming tool:", EGAM2_TriggerSkimmingTool)

# do the OR of trigger-based and offline-based selection
EGAM2_SkimmingTool = DerivationFramework__FilterCombinationOR(
    name="EGAM2SkimmingTool",
    FilterList=[EGAM2_OfflineSkimmingTool, EGAM2_TriggerSkimmingTool])
ToolSvc += EGAM2_SkimmingTool


# =======================================
# CREATE THE DERIVATION KERNEL ALGORITHM
# =======================================

print("EGAM2 skimming tools: ", [EGAM2_SkimmingTool])
print("EGAM2 thinning tools: ", thinningTools)
print("EGAM2 augmentation tools: ", augmentationTools)
EGAM2Sequence += CfgMgr.DerivationFramework__DerivationKernel(
    "EGAM2Kernel",
    AugmentationTools=augmentationTools,
    SkimmingTools=[EGAM2_SkimmingTool],
    ThinningTools=thinningTools
)

# ====================================================================
# JET/MET
# ====================================================================
jetList=[]
if DerivationFrameworkIsMonteCarlo:
    jetList += [AntiKt4Truth, AntiKt4TruthDressedWZ]
addDAODJets(jetList, EGAM2Sequence)


# ====================================================================
# SET UP SLIMMING
# ====================================================================
EGAM2SlimmingHelper = SlimmingHelper("EGAM2SlimmingHelper")
EGAM2SlimmingHelper.SmartCollections = ["Electrons",
                                        "Photons",
                                        "Muons",
                                        "TauJets",
                                        "MET_Baseline_AntiKt4EMPFlow",
                                        "AntiKt4EMPFlowJets",
                                        "BTagging_AntiKt4EMPFlow",
                                        "InDetTrackParticles",
                                        "PrimaryVertices"]
# muons, tau, MET, b-tagging could be switched off if not needed and use too much space
if DerivationFrameworkIsMonteCarlo:
    EGAM2SlimmingHelper.SmartCollections += ["AntiKt4TruthJets",
                                             "AntiKt4TruthDressedWZJets"]

# Add egamma trigger objects
EGAM2SlimmingHelper.IncludeEGammaTriggerContent = True

# Extra variables
EGAM2SlimmingHelper.ExtraVariables = ExtraContentAll
EGAM2SlimmingHelper.AllVariables = ExtraContainersElectrons
EGAM2SlimmingHelper.AllVariables += ExtraContainersTrigger

if DerivationFrameworkIsMonteCarlo:
    EGAM2SlimmingHelper.ExtraVariables += ExtraContentAllTruth
    EGAM2SlimmingHelper.AllVariables += ExtraContainersTruth
else:
    EGAM2SlimmingHelper.ExtraVariables += ExtraContainersTriggerDataOnly

for tool in EGAM2_ClusterEnergyPerLayerDecorators:
    EGAM2SlimmingHelper.ExtraVariables.extend(
        getClusterEnergyPerLayerDecorations(tool))


# Add event info
if jobproperties.egammaDFFlags.doEGammaEventInfoSlimming:
    EGAM2SlimmingHelper.SmartCollections.append("EventInfo")
else:
    EGAM2SlimmingHelper.AllVariables += ["EventInfo"]

# Add detailed shower shape variables
EGAM2SlimmingHelper.ExtraVariables += ElectronsCPDetailedContent
EGAM2SlimmingHelper.ExtraVariables += GSFTracksCPDetailedContent
EGAM2SlimmingHelper.ExtraVariables += PhotonsCPDetailedContent

# This line must come after we have finished configuring EGAM2SlimmingHelper
EGAM2SlimmingHelper.AppendContentToStream(EGAM2Stream)
