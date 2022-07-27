#********************************************************************
# EGAM5.py
# W->enu reduction
# reductionConf flag EGAM5 in Reco_tf.py
# author: giovanni.marchiori@cern.ch
#********************************************************************

from DerivationFrameworkCore.DerivationFrameworkMaster import buildFileName
from DerivationFrameworkCore.DerivationFrameworkMaster import (
    DerivationFrameworkIsMonteCarlo, DerivationFrameworkJob)
from DerivationFrameworkPhys import PhysCommon
from DerivationFrameworkEGamma.EGammaCommon import *
from DerivationFrameworkEGamma.EGAM5ExtraContent import *


#====================================================================
# read common DFEGamma settings from egammaDFFlags
#====================================================================
from DerivationFrameworkEGamma.egammaDFFlags import jobproperties
jobproperties.egammaDFFlags.print_JobProperties("full")

# this could also go in egammaDFFlags
RecomputeElectronSelectors = True
#RecomputeElectronSelectors = False


#====================================================================
# check if we run on data or MC
#====================================================================
print("DerivationFrameworkIsMonteCarlo: ", DerivationFrameworkIsMonteCarlo)
if not DerivationFrameworkIsMonteCarlo:
    ExtraContainersTrigger += ExtraContainersTriggerDataOnly


#====================================================================
# Set up sequence for this format and add to the top sequence
#====================================================================
EGAM5Sequence = CfgMgr.AthSequencer("EGAM5Sequence")
DerivationFrameworkJob += EGAM5Sequence


#====================================================================
# SET UP STREAM
#====================================================================
streamName = derivationFlags.WriteDAOD_EGAM5Stream.StreamName
fileName   = buildFileName( derivationFlags.WriteDAOD_EGAM5Stream )
EGAM5Stream = MSMgr.NewPoolRootStream( streamName, fileName )
# Only events that pass the filters listed below are written out.
# Name must match that of the kernel above
# AcceptAlgs  = logical OR of filters
# RequireAlgs = logical AND of filters
EGAM5Stream.AcceptAlgs(["EGAM5Kernel"])


### Thinning and augmentation tools lists
augmentationTools = []
thinningTools=[]


#====================================================================
# SET UP AUGMENTATIONS
#====================================================================

#====================================================================
# transverse mass of W->enu candidates for electron calibration
#====================================================================

# could add track isolation (if included in single electron trigger..)
if RecomputeElectronSelectors :
    requirement_el = '(Electrons.DFCommonElectronsLHTight) && Electrons.pt > 24.5*GeV'
else :
    requirement_el = '(Electrons.LHTight) && Electrons.pt > 24.5*GeV'

from DerivationFrameworkEGamma.DerivationFrameworkEGammaConf import DerivationFramework__EGTransverseMassTool
EGAM5_MTTool = DerivationFramework__EGTransverseMassTool(
    name = "EGAM5_MTTool",
    ObjectRequirements = requirement_el,
    METmin = 25*GeV,
    StoreGateEntryName = "WENU_TransverseMass",
    ObjectMassHypothesis = 0.511*MeV,
    ObjectContainerName = "Electrons",
    METContainerName = "MET_Core_AntiKt4EMPFlow",
)
ToolSvc += EGAM5_MTTool
print(EGAM5_MTTool)
augmentationTools += [EGAM5_MTTool]


#====================================================================
# Max Cell sum decoration tool
#====================================================================
from DerivationFrameworkCalo.DerivationFrameworkCaloConf import DerivationFramework__MaxCellDecorator
EGAM5_MaxCellDecoratorTool = DerivationFramework__MaxCellDecorator( name                    = "EGAM5_MaxCellDecoratorTool",
                                                                    SGKey_electrons         = "Electrons",
                                                                    SGKey_photons           = "Photons")
ToolSvc += EGAM5_MaxCellDecoratorTool
augmentationTools += [EGAM5_MaxCellDecoratorTool]


#====================================================================
# Gain and cluster energies per layer decoration tool
#====================================================================
from DerivationFrameworkCalo.DerivationFrameworkCaloFactories import GainDecorator, getGainDecorations, getClusterEnergyPerLayerDecorator, getClusterEnergyPerLayerDecorations
EGAM5_GainDecoratorTool = GainDecorator()
ToolSvc += EGAM5_GainDecoratorTool
augmentationTools += [EGAM5_GainDecoratorTool]

cluster_sizes = (3,7), (5,5), (7,11)
EGAM5_ClusterEnergyPerLayerDecorators = [getClusterEnergyPerLayerDecorator(neta, nphi)() for neta, nphi in cluster_sizes]
augmentationTools += EGAM5_ClusterEnergyPerLayerDecorators


#====================================================================
# SET UP THINNING
#====================================================================
print('WARNING, Thinning of trigger navigation has to be properly implemented in R22')
#from DerivationFrameworkCore.ThinningHelper import ThinningHelper
#EGAM5ThinningHelper = ThinningHelper( "EGAM5ThinningHelper" )
#EGAM5ThinningHelper.TriggerChains = '(^(?!.*_[0-9]*(mu|j|xe|tau|ht|xs|te))(?!HLT_[eg].*_[0-9]*[eg][0-9].*)(?!HLT_eb.*)(?!.*larpeb.*)(?!HLT_.*_AFP_.*)(HLT_[eg].*))'
#EGAM5ThinningHelper.AppendToStream( EGAM5Stream, ExtraContainersTrigger )


# Track thinning
# See recommedations here:
# https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/DaodRecommendations
if jobproperties.egammaDFFlags.doEGammaDAODTrackThinning:

    TrackThinningKeepElectronTracks = True
    TrackThinningKeepPhotonTracks = True
    TrackThinningKeepJetTracks = False
    TrackThinningKeepMuonTracks = False
    TrackThinningKeepTauTracks = False
    TrackThinningKeepPVTracks = True

    # Tracks associated with Jets
    if (TrackThinningKeepJetTracks) : 
        from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__JetTrackParticleThinning
        EGAM5JetTPThinningTool = DerivationFramework__JetTrackParticleThinning(
            name                    = "EGAM5JetTPThinningTool",
            StreamName              = streamName,
            JetKey                  = "AntiKt4EMPFlowJets",
            InDetTrackParticlesKey  = "InDetTrackParticles")
        ToolSvc += EGAM5JetTPThinningTool
        print(EGAM5JetTPThinningTool)
        thinningTools.append(EGAM5JetTPThinningTool)
    
    # Tracks associated with Muons
    if (TrackThinningKeepMuonTracks) :
        from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__MuonTrackParticleThinning
        EGAM5MuonTPThinningTool = DerivationFramework__MuonTrackParticleThinning(
            name                    = "EGAM5MuonTPThinningTool",
            StreamName              = streamName,
            MuonKey                 = "Muons",
            InDetTrackParticlesKey  = "InDetTrackParticles")
        ToolSvc += EGAM5MuonTPThinningTool
        print(EGAM5MuonTPThinningTool)
        thinningTools.append(EGAM5MuonTPThinningTool)
    
    # Tracks associated with Electrons
    if (TrackThinningKeepElectronTracks) : 
        from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__EgammaTrackParticleThinning
        EGAM5ElectronTPThinningTool = DerivationFramework__EgammaTrackParticleThinning(
            name                    = "EGAM5ElectronTPThinningTool",
            StreamName              = streamName,
            SGKey                   = "Electrons",
            GSFTrackParticlesKey    = "GSFTrackParticles",        
            InDetTrackParticlesKey  = "InDetTrackParticles",
            SelectionString         = "Electrons.pt > 0*GeV",
            BestMatchOnly = True,
            ConeSize = 0.3)
        ToolSvc += EGAM5ElectronTPThinningTool
        print(EGAM5ElectronTPThinningTool)
        thinningTools.append(EGAM5ElectronTPThinningTool)
        
    # Tracks associated with Photons
    if (TrackThinningKeepPhotonTracks) : 
        from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__EgammaTrackParticleThinning
        EGAM5PhotonTPThinningTool = DerivationFramework__EgammaTrackParticleThinning(
            name                    = "EGAM5PhotonTPThinningTool",
            StreamName              = streamName,
            SGKey                   = "Photons",
            GSFTrackParticlesKey    = "GSFTrackParticles",        
            InDetTrackParticlesKey  = "InDetTrackParticles",
            SelectionString         = "Photons.pt > 0*GeV",
            BestMatchOnly = True,
            ConeSize = 0.3)
        
        ToolSvc += EGAM5PhotonTPThinningTool
        print(EGAM5PhotonTPThinningTool)
        thinningTools.append(EGAM5PhotonTPThinningTool)

    # Tracks associated with Taus
    if (TrackThinningKeepTauTracks) : 
        from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__TauTrackParticleThinning
        EGAM5TauTPThinningTool = DerivationFramework__TauTrackParticleThinning(
            name                    = "EGAM5TauTPThinningTool",
            StreamName              = streamName,
            TauKey                  = "TauJets",
            ConeSize                = 0.6,
            InDetTrackParticlesKey  = "InDetTrackParticles")
        ToolSvc += EGAM5TauTPThinningTool
        print(EGAM5TauTPThinningTool)
        thinningTools.append(EGAM5TauTPThinningTool)
        
    # Tracks from primary vertex
    thinning_expression = "InDetTrackParticles.DFCommonTightPrimary && abs(DFCommonInDetTrackZ0AtPV)*sin(InDetTrackParticles.theta) < 3.0*mm && InDetTrackParticles.pt > 10*GeV"
    if (TrackThinningKeepPVTracks) :
        from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__TrackParticleThinning
        EGAM5TPThinningTool = DerivationFramework__TrackParticleThinning(
            name                    = "EGAM5TPThinningTool",
            StreamName              = streamName,
            SelectionString         = thinning_expression,
            InDetTrackParticlesKey  = "InDetTrackParticles")
        ToolSvc += EGAM5TPThinningTool
        print(EGAM5TPThinningTool)
        thinningTools.append(EGAM5TPThinningTool)


#====================================================================
# Setup the skimming criteria
#====================================================================

#====================================================================
# 1st selection: trigger-based (WTP triggers)
#====================================================================

# L1Topo W T&P 
triggers =  ['HLT_e13_etcut_trkcut' ]
triggers +=  ['HLT_e18_etcut_trkcut' ]
## # Non-L1Topo W TP commissioning triggers ==> in MC, in 50 ns data
triggers +=  ['HLT_e13_etcut_trkcut_xs15' ]
triggers +=  ['HLT_e18_etcut_trkcut_xs20' ]
## W T&P triggers ==> not in MC, in 50 ns data
triggers +=  ['HLT_e13_etcut_trkcut_xs15_mt25' ]
triggers +=  ['HLT_e18_etcut_trkcut_xs20_mt35' ]
###W T&P triggers ==> not in MC, not in 50 ns data, will be in 25 ns data
triggers +=  ['HLT_e13_etcut_trkcut_xs15_j20_perf_xe15_2dphi05' ]
triggers +=  ['HLT_e13_etcut_trkcut_xs15_j20_perf_xe15_2dphi05_mt25' ]
triggers +=  ['HLT_e13_etcut_trkcut_j20_perf_xe15_2dphi05_mt25' ]
triggers +=  ['HLT_e13_etcut_trkcut_j20_perf_xe15_2dphi05' ]
triggers +=  ['HLT_e13_etcut_trkcut_xs15_j20_perf_xe15_6dphi05' ]
triggers +=  ['HLT_e13_etcut_trkcut_xs15_j20_perf_xe15_6dphi05_mt25' ]
triggers +=  ['HLT_e13_etcut_trkcut_j20_perf_xe15_6dphi05_mt25' ]
triggers +=  ['HLT_e13_etcut_trkcut_j20_perf_xe15_6dphi05' ]
triggers +=  ['HLT_e18_etcut_trkcut_xs20_j20_perf_xe20_6dphi15' ]
triggers +=  ['HLT_e18_etcut_trkcut_xs20_j20_perf_xe20_6dphi15_mt35' ]
triggers +=  ['HLT_e18_etcut_trkcut_j20_perf_xe20_6dphi15_mt35' ]
triggers +=  ['HLT_e18_etcut_trkcut_j20_perf_xe20_6dphi15' ]

# others
triggers +=  ['HLT_e5_etcut_L1W-05DPHI-JXE-0']
triggers +=  ['HLT_e5_etcut_L1W-10DPHI-JXE-0']
triggers +=  ['HLT_e5_etcut_L1W-15DPHI-JXE-0']
triggers +=  ['HLT_e5_etcut_L1W-10DPHI-EMXE-0']
triggers +=  ['HLT_e5_etcut_L1W-15DPHI-EMXE-0']
triggers +=  ['HLT_e5_etcut_L1W-05DPHI-EMXE-1']
triggers +=  ['HLT_e5_etcut_L1W-05RO-XEHT-0']
triggers +=  ['HLT_e5_etcut_L1W-90RO2-XEHT-0']
triggers +=  ['HLT_e5_etcut_L1W-250RO2-XEHT-0']
triggers +=  ['HLT_e5_etcut_L1W-HT20-JJ15.ETA49']

triggers +=  ['HLT_e13_etcut_L1W-NOMATCH']
triggers +=  ['HLT_e13_etcut_L1W-NOMATCH_W-05RO-XEEMHT']
triggers +=  ['HLT_e13_etcut_L1EM10_W-MT25']
triggers +=  ['HLT_e13_etcut_L1EM10_W-MT30']
triggers +=  ['HLT_e13_etcut_trkcut_L1EM12']
triggers +=  ['HLT_e13_etcut_trkcut_L1EM10_W-MT25_W-15DPHI-JXE-0_W-15DPHI-EMXE']
triggers +=  ['HLT_e13_etcut_trkcut_j20_perf_xe15_6dphi15_mt25']
triggers +=  ['HLT_e13_etcut_trkcut_j20_perf_xe15_6dphi15_mt25_L1EM12_W-MT25_W-15DPHI-JXE-0_W-15DPHI-EMXE_XS20']
triggers +=  ['HLT_e13_etcut_trkcut_j20_perf_xe15_6dphi15_mt25_L1EM12_W-MT25_W-15DPHI-JXE-0_W-15DPHI-EMXE_W-90RO2-XEHT-0']
triggers +=  ['HLT_e13_etcut_trkcut_xs30_xe30_mt35']
triggers +=  ['HLT_e13_etcut_trkcut_xs30_j15_perf_xe30_6dphi05_mt35']
triggers +=  ['HLT_e13_etcut_trkcut_xs30_j15_perf_xe30_6dphi15_mt35']
triggers +=  ['HLT_e13_etcut_trkcut_xs30_j15_perf_xe30_2dphi05_mt35']
triggers +=  ['HLT_e13_etcut_trkcut_xs30_j15_perf_xe30_2dphi15_mt35']
triggers +=  ['HLT_e13_etcut_trkcut_xs30_j15_perf_xe30_2dphi15_mt35_L1EM12_W-MT25_W-15DPHI-JXE-0_W-15DPHI-EMXE_XS20']
triggers +=  ['HLT_e13_etcut_trkcut_xs30_j15_perf_xe30_6dphi15_mt35_L1EM12_W-MT25_W-15DPHI-JXE-0_W-15DPHI-EMXE_W-90RO2-XEHT-0']

triggers +=  ['HLT_e18_etcut_L1EM15_W-MT35']
triggers +=  ['HLT_e18_etcut_trkcut_L1EM15']
triggers +=  ['HLT_e18_etcut_trkcut_L1EM15_W-MT35_W-05DPHI-JXE-0_W-05DPHI-EMXE']
triggers +=  ['HLT_e18_etcut_trkcut_xs30_xe30_mt35']
triggers +=  ['HLT_e18_etcut_trkcut_xs30_j15_perf_xe30_6dphi05_mt35']
triggers +=  ['HLT_e18_etcut_trkcut_xs30_j15_perf_xe30_6dphi15_mt35']
triggers +=  ['HLT_e18_etcut_trkcut_xs30_j15_perf_xe30_2dphi05_mt35']
triggers +=  ['HLT_e18_etcut_trkcut_xs30_j15_perf_xe30_2dphi15_mt35']
triggers +=  ['HLT_e18_etcut_trkcut_xs30_j15_perf_xe30_6dphi05_mt35_L1EM15_W-MT35_W-05DPHI-JXE-0_W-05DPHI-EM15XE_XS30'] 
triggers +=  ['HLT_e18_etcut_trkcut_xs30_j15_perf_xe30_6dphi15_mt35_L1EM15_W-MT35_W-05DPHI-JXE-0_W-05DPHI-EM15XE_XS30'] 
triggers +=  ['HLT_e18_etcut_trkcut_xs30_j15_perf_xe30_2dphi05_mt35_L1EM15_W-MT35_W-250RO2-XEHT-0_W-05DPHI-JXE-0_W-05DPHI-EM15XE'] 
triggers +=  ['HLT_e18_etcut_trkcut_xs30_j15_perf_xe30_2dphi15_mt35_L1EM15_W-MT35_W-250RO2-XEHT-0_W-05DPHI-JXE-0_W-05DPHI-EM15XE'] 
triggers +=  ['HLT_e18_etcut_trkcut_xs30_j15_perf_xe30_2dphi15_mt35_L1EM15_W-MT35_W-05DPHI-JXE-0_W-05DPHI-EM15XE_XS30']
triggers +=  ['HLT_e18_etcut_trkcut_xs30_j15_perf_xe30_6dphi05_mt35_L1EM15_W-MT35_W-250RO2-XEHT-0_W-05DPHI-JXE-0_W-05DPHI-EM15XE']
triggers +=  ['HLT_e18_etcut_trkcut_xs30_j15_perf_xe30_6dphi15_mt35_L1EM15_W-MT35_W-250RO2-XEHT-0_W-05DPHI-JXE-0_W-05DPHI-EM15XE']
triggers +=  ['HLT_e18_etcut_trkcut_xs30_j15_perf_xe30_6dphi15_mt35_L1EM15_W-MT35_W-250RO2-XEHT-0_W-15DPHI-JXE-0_W-15DPHI-EM15XE']
triggers +=  ['HLT_e18_etcut_trkcut_xs30_xe30_mt35_L1EM15_W-MT35_W-05DPHI-JXE-0_W-05DPHI-EM15XE_XS30'] 
triggers +=  ['HLT_e18_etcut_trkcut_xs30_xe30_mt35_L1EM15_W-MT35_W-250RO2-XEHT-0_W-05DPHI-JXE-0_W-05DPHI-EM15XE'] 
triggers +=  ['HLT_e18_etcut_trkcut_xs30_xe30_mt35_L1EM15_W-MT35_W-250RO2-XEHT-0_W-15DPHI-JXE-0_W-15DPHI-EM15XE']
triggers +=  ['HLT_e18_etcut_trkcut_j20_perf_xe20_6dphi15_mt35_L1EM15_W-MT35_W-05DPHI-JXE-0_W-05DPHI-EM15XE_XS30'] 
triggers +=  ['HLT_e18_etcut_trkcut_j20_perf_xe20_6dphi15_mt35_L1EM15_W-MT35_W-250RO2-XEHT-0_W-05DPHI-JXE-0_W-05DPHI-EM15XE']

# added for 2017
triggers += ['HLT_e60_etcut']
triggers += ['HLT_e60_etcut_L1EM24VHIM']
triggers += ['HLT_e60_etcut_trkcut_L1EM24VHIM_j15_perf_xe60_6dphi15_mt35']
triggers += ['HLT_e60_etcut_trkcut_L1EM24VHIM_xe60_mt35']
triggers += ['HLT_e60_etcut_trkcut_L1EM24VHIM_xs30_j15_perf_xe30_6dphi15_mt35']
triggers += ['HLT_e60_etcut_trkcut_L1EM24VHIM_xs30_xe30_mt35']
triggers += ['HLT_e60_lhmedium_nod0']
triggers += ['HLT_e60_lhmedium_nod0_L1EM24VHI']
triggers += ['HLT_e60_lhmedium_nod0_L1EM24VHIM']
triggers += ['HLT_e60_lhvloose_nod0']
triggers += ['HLT_e60_etcut_trkcut_j15_perf_xe60_6dphi05_mt35']
triggers += ['HLT_e60_etcut_trkcut_xs30_j15_perf_xe30_6dphi05_mt35']
triggers += ['HLT_e70_etcut']
triggers += ['HLT_e70_etcut_L1EM24VHIM']
triggers += ['HLT_e70_lhloose_nod0_L1EM24VHIM_xe70noL1']
triggers += ['HLT_e70_lhloose_nod0_xe70noL1']
triggers += ['HLT_noalg_l1topo_L1EM15']
triggers += ['HLT_noalg_l1topo_L1EM7']
triggers += ['HLT_j80_xe80']
triggers += ['HLT_xe80_tc_lcw_L1XE50']
triggers += ['HLT_xe90_mht_L1XE50']
triggers += ['HLT_xe90_tc_lcw_wEFMu_L1XE50']
triggers += ['HLT_xe90_mht_wEFMu_L1XE50']
triggers += ['HLT_xe110_mht_L1XE50']
triggers += ['HLT_xe110_pufit_L1XE50']

#added for low-mu data analysis, 2017 and 2018 data
triggers += ['HLT_e15_lhloose_nod0_L1EM12']
#added for low-mu data analysis, 2018 data
triggers += ['HLT_xe35']
triggers += ['HLT_e15_etcut_trkcut_xe30noL1']

from DerivationFrameworkTools.DerivationFrameworkToolsConf import DerivationFramework__TriggerSkimmingTool
EGAM5_TriggerSkimmingTool = DerivationFramework__TriggerSkimmingTool( name                   = "EGAM5_TriggerSkimmingTool",
                                                                      TriggerListOR          = triggers)
ToolSvc += EGAM5_TriggerSkimmingTool
print("EGAM5 trigger skimming tool:", EGAM5_TriggerSkimmingTool)

#====================================================
# Second selection: offline, based on mT(enu)
#====================================================
expression = 'count(WENU_TransverseMass>40*GeV)>=1'
from DerivationFrameworkTools.DerivationFrameworkToolsConf import DerivationFramework__xAODStringSkimmingTool
EGAM5_OfflineSkimmingTool = DerivationFramework__xAODStringSkimmingTool( name = "EGAM5_OfflineSkimmingTool",
                                                                         expression = expression)
ToolSvc += EGAM5_OfflineSkimmingTool
print("EGAM5 offline skimming tool:", EGAM5_OfflineSkimmingTool)

#====================================================
# Third selection: mix of offline and online criteria
#====================================================
expression = '(HLT_e60_lhloose_xe60noL1 || HLT_e120_lhloose || HLT_j80_xe80 || HLT_xe70 || HLT_xe80_tc_lcw_L1XE50  || HLT_xe90_mht_L1XE50 ||  HLT_xe90_tc_lcw_wEFMu_L1XE50 || HLT_xe90_mht_wEFMu_L1XE50) && count(Electrons.pt>14.5*GeV)>=1'
EGAM5_ThirdSkimmingTool = DerivationFramework__xAODStringSkimmingTool( name = "EGAM5_ThirdSkimmingTool",
                                                                       expression = expression)
ToolSvc += EGAM5_ThirdSkimmingTool
print("EGAM5 offline skimming tool:", EGAM5_ThirdSkimmingTool)

#====================================================
# Make OR of previous selections
#====================================================
from DerivationFrameworkTools.DerivationFrameworkToolsConf import DerivationFramework__FilterCombinationOR
EGAM5_SkimmingTool = DerivationFramework__FilterCombinationOR(name="EGAM5_SkimmingTool", FilterList=[EGAM5_TriggerSkimmingTool,EGAM5_OfflineSkimmingTool,EGAM5_ThirdSkimmingTool] )
ToolSvc+=EGAM5_SkimmingTool


#=======================================
# CREATE THE DERIVATION KERNEL ALGORITHM
#=======================================
from DerivationFrameworkCore.DerivationFrameworkCoreConf import DerivationFramework__DerivationKernel
print("EGAM5 skimming tools: ", [EGAM5_SkimmingTool])
print("EGAM5 thinning tools: ", thinningTools)
print("EGAM5 augmentation tools: ", augmentationTools)
EGAM5Sequence += CfgMgr.DerivationFramework__DerivationKernel("EGAM5Kernel",
                                                              AugmentationTools = augmentationTools,
                                                              SkimmingTools = [EGAM5_SkimmingTool],
                                                              ThinningTools = thinningTools
                                                              )


#====================================================================
# JET/MET
#====================================================================
from DerivationFrameworkJetEtMiss.JetCommon import addDAODJets
from JetRecConfig.StandardSmallRJets import AntiKt4Truth, AntiKt4TruthDressedWZ
jetList=[]
if DerivationFrameworkIsMonteCarlo:
    jetList += [AntiKt4Truth, AntiKt4TruthDressedWZ]
addDAODJets(jetList, EGAM5Sequence)


#====================================================================
# SET UP SLIMMING
#====================================================================
from DerivationFrameworkCore.SlimmingHelper import SlimmingHelper
EGAM5SlimmingHelper = SlimmingHelper("EGAM5SlimmingHelper")
EGAM5SlimmingHelper.SmartCollections = ["Electrons",
                                        "Photons",
                                        "Muons",
                                        "TauJets",
                                        "MET_Baseline_AntiKt4EMPFlow",
                                        "AntiKt4EMPFlowJets",
                                        "BTagging_AntiKt4EMPFlow",
                                        "InDetTrackParticles",
                                        "PrimaryVertices" ]
if DerivationFrameworkIsMonteCarlo:
    EGAM5SlimmingHelper.SmartCollections += ["AntiKt4TruthJets",
                                             "AntiKt4TruthDressedWZJets"]

# Add egamma trigger objects
EGAM5SlimmingHelper.IncludeEGammaTriggerContent = True

# Extra variables
EGAM5SlimmingHelper.ExtraVariables = ExtraContentAll
EGAM5SlimmingHelper.AllVariables = ExtraContainersElectrons
EGAM5SlimmingHelper.AllVariables += ExtraContainersTrigger

if DerivationFrameworkIsMonteCarlo:
    EGAM5SlimmingHelper.ExtraVariables += ExtraContentAllTruth
    EGAM5SlimmingHelper.AllVariables += ExtraContainersTruth
else:
    EGAM5SlimmingHelper.ExtraVariables += ExtraContainersTriggerDataOnly

for tool in EGAM5_ClusterEnergyPerLayerDecorators:
    EGAM5SlimmingHelper.ExtraVariables.extend( getClusterEnergyPerLayerDecorations( tool ) )

# Add event info
if jobproperties.egammaDFFlags.doEGammaEventInfoSlimming:
    EGAM5SlimmingHelper.SmartCollections.append("EventInfo")
else:
    EGAM5SlimmingHelper.AllVariables += ["EventInfo"]

# Add detailed shower shape variables
from DerivationFrameworkEGamma.ElectronsCPDetailedContent import *
EGAM5SlimmingHelper.ExtraVariables += ElectronsCPDetailedContent
EGAM5SlimmingHelper.ExtraVariables += GSFTracksCPDetailedContent
from DerivationFrameworkEGamma.PhotonsCPDetailedContent import *
EGAM5SlimmingHelper.ExtraVariables += PhotonsCPDetailedContent

# This line must come after we have finished configuring EGAM5SlimmingHelper
EGAM5SlimmingHelper.AppendContentToStream(EGAM5Stream)

