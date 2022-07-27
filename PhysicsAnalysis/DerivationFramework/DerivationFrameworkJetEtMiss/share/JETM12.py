#====================================================================
# JETM12.py 
# reductionConf flag JETM12 in Reco_tf.py   
#====================================================================

from DerivationFrameworkCore.DerivationFrameworkMaster import DerivationFrameworkIsMonteCarlo, DerivationFrameworkJob, buildFileName
from DerivationFrameworkJetEtMiss.JetCommon import OutputJets
from DerivationFrameworkPhys import PhysCommon

#====================================================================
# SKIMMING TOOL 
#====================================================================

from DerivationFrameworkJetEtMiss import TriggerLists
metTriggers = TriggerLists.MET_Trig()
muTriggers = TriggerLists.single_mu_Trig()
orstr  = ' || '
andstr = ' && '
trackRequirements = '(InDetTrackParticles.pt > 10.*GeV && InDetTrackParticles.TrkIsoPt1000_ptcone20 < 0.12*InDetTrackParticles.pt && InDetTrackParticles.DFCommonTightPrimary && abs(DFCommonInDetTrackZ0AtPV) < 5.0*mm )'
trackRequirementsMu = '(InDetTrackParticles.pt > 70.*GeV && InDetTrackParticles.TrkIsoPt1000_ptcone20 < 0.12*InDetTrackParticles.pt && InDetTrackParticles.DFCommonTightPrimary && abs(DFCommonInDetTrackZ0AtPV) < 5.0*mm )'
trackRequirementsTtbar = '(InDetTrackParticles.pt > 25.*GeV && InDetTrackParticles.TrkIsoPt1000_ptcone20 < 0.12*InDetTrackParticles.pt && InDetTrackParticles.DFCommonTightPrimary && abs(DFCommonInDetTrackZ0AtPV) < 5.0*mm )'
jetRequirementsTtbar = '( AntiKt4EMTopoJets.pt > 20*GeV && BTagging_AntiKt4EMTopo.MV2c10_discriminant > 0.11 )'
expressionW = '( (' + orstr.join(metTriggers) + ' )' + andstr + '( count('+trackRequirements+') >=1 ) )'
expressionMu = '( (' + orstr.join(muTriggers) + ' )' + andstr + '( count('+trackRequirementsMu+') >=1 ) )'
expressionTtbar = '( (' + orstr.join(muTriggers) + ' )' + andstr + '( count('+trackRequirementsTtbar+') >=1 )' + andstr + '( count('+trackRequirements+') >=2 )' + andstr + '( count('+jetRequirementsTtbar+') >=1 ) )'
expression = '( '+expressionW+' || '+expressionMu+' || '+expressionTtbar+' )'

from DerivationFrameworkTools.DerivationFrameworkToolsConf import DerivationFramework__xAODStringSkimmingTool
JETM12SkimmingTool = DerivationFramework__xAODStringSkimmingTool( name = "JETM12SkimmingTool1",
                                                                    expression = expression)
ToolSvc += JETM12SkimmingTool

#====================================================================
# SET UP STREAM   
#====================================================================
streamName = derivationFlags.WriteDAOD_JETM12Stream.StreamName
fileName   = buildFileName( derivationFlags.WriteDAOD_JETM12Stream )
JETM12Stream = MSMgr.NewPoolRootStream( streamName, fileName )
JETM12Stream.AcceptAlgs(["JETM12Kernel"])

#=======================================
# ESTABLISH THE THINNING HELPER
#=======================================

from DerivationFrameworkCore.ThinningHelper import ThinningHelper
JETM12ThinningHelper = ThinningHelper( "JETM12ThinningHelper" )
JETM12ThinningHelper.AppendToStream( JETM12Stream )

thinningTools = []
AugmentationToolsSkim   = []
AugmentationTools   = []

#====================================================================
# ISOLATION TOOL 
#====================================================================
#Track selection
from IsolationTool.IsolationToolConf import xAOD__TrackIsolationTool
TrackIsoTool = xAOD__TrackIsolationTool("TrackIsoTool")
TrackIsoTool.TrackSelectionTool.maxZ0SinTheta= 3.
TrackIsoTool.TrackSelectionTool.minPt= 1000.
TrackIsoTool.TrackSelectionTool.CutLevel= "Loose"
ToolSvc += TrackIsoTool

TrackIsoTool500 = xAOD__TrackIsolationTool("TrackIsoTool500")
TrackIsoTool500.TrackSelectionTool.maxZ0SinTheta= 3.
TrackIsoTool500.TrackSelectionTool.minPt= 500.
TrackIsoTool500.TrackSelectionTool.CutLevel= "Loose"
ToolSvc += TrackIsoTool500

from DerivationFrameworkSUSY.DerivationFrameworkSUSYConf import DerivationFramework__trackIsolationDecorator
import ROOT, cppyy
cppyy.load_library('libxAODCoreRflxDict')
cppyy.load_library('libxAODPrimitivesDict')
isoPar = ROOT.xAOD.Iso
Pt1000IsoTrackDecorator = DerivationFramework__trackIsolationDecorator(name = "Pt1000IsoTrackDecorator",
                                                                TrackIsolationTool = TrackIsoTool,
                                                                TargetContainer = "InDetTrackParticles",
                                                                ptcones = [isoPar.ptcone40,isoPar.ptcone30,isoPar.ptcone20],
                                                                Prefix = 'TrkIsoPt1000_'
                                                               )
Pt500IsoTrackDecorator = DerivationFramework__trackIsolationDecorator(name = "Pt500IsoTrackDecorator",
                                                                TrackIsolationTool = TrackIsoTool500,
                                                                TargetContainer = "InDetTrackParticles",
                                                                ptcones = [isoPar.ptcone40,isoPar.ptcone30,isoPar.ptcone20],
                                                                Prefix = 'TrkIsoPt500_'
                                                               )
ToolSvc += Pt1000IsoTrackDecorator
ToolSvc += Pt500IsoTrackDecorator

AugmentationToolsSkim.append(Pt1000IsoTrackDecorator)
AugmentationTools.append(Pt500IsoTrackDecorator)

#====================================================================
# THINNING TOOLS 
#====================================================================
# TrackParticles directly
from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__TrackParticleThinning
JETM12TPThinningTool = DerivationFramework__TrackParticleThinning(name = "JETM12TPThinningTool",
                                                                 StreamName              = streamName,
                                                                 SelectionString         = "( InDetTrackParticles.pt > 10*GeV && InDetTrackParticles.DFCommonTightPrimary && abs(DFCommonInDetTrackZ0AtPV) < 5.0*mm )",
                                                                 InDetTrackParticlesKey  = "InDetTrackParticles")
ToolSvc += JETM12TPThinningTool
thinningTools.append(JETM12TPThinningTool)

# TrackParticles associated with Muons
from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__MuonTrackParticleThinning
JETM12MuonTPThinningTool = DerivationFramework__MuonTrackParticleThinning(name     = "JETM12MuonTPThinningTool",
                                                                    StreamName              = streamName,
                                                                    MuonKey                 = "Muons",
                                                                    InDetTrackParticlesKey  = "InDetTrackParticles")
ToolSvc += JETM12MuonTPThinningTool
thinningTools.append(JETM12MuonTPThinningTool)

# TrackParticles associated with electrons
from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__EgammaTrackParticleThinning
JETM12ElectronTPThinningTool = DerivationFramework__EgammaTrackParticleThinning(name                    = "JETM12ElectronTPThinningTool",
                                                                               StreamName              = streamName,
                                                                               SGKey                   = "Electrons",
                                                                               InDetTrackParticlesKey  = "InDetTrackParticles")
ToolSvc += JETM12ElectronTPThinningTool
thinningTools.append(JETM12ElectronTPThinningTool)

# TrackParticles associated with photons
JETM12PhotonTPThinningTool = DerivationFramework__EgammaTrackParticleThinning(name                    = "JETM12PhotonTPThinningTool",
                                                                             StreamName              = streamName,
                                                                             SGKey                   = "Photons",
                                                                             InDetTrackParticlesKey  = "InDetTrackParticles")
ToolSvc += JETM12PhotonTPThinningTool
thinningTools.append(JETM12PhotonTPThinningTool)

# TrackParticles associated with taus
from DerivationFrameworkInDet.DerivationFrameworkInDetConf import DerivationFramework__TauTrackParticleThinning
JETM12TauTPThinningTool = DerivationFramework__TauTrackParticleThinning( name            = "JETM12TauTPThinningTool",
                                                                        StreamName      = streamName,
                                                                        TauKey          = "TauJets",
                                                                        InDetTrackParticlesKey  = "InDetTrackParticles")
ToolSvc += JETM12TauTPThinningTool
thinningTools.append(JETM12TauTPThinningTool)

##################################
# CaloCalTopoCluster thinning
from DerivationFrameworkCalo.DerivationFrameworkCaloConf import DerivationFramework__CaloClusterThinning
JETM12CaloClusterThinning  = DerivationFramework__CaloClusterThinning(name                  = "JETM12CaloClusterThinning",
                                                                       StreamName            = streamName,
                                                                       SGKey                 = "InDetTrackParticles",
                                                                       TopoClCollectionSGKey = "CaloCalTopoClusters",
                                                                       SelectionString = "( InDetTrackParticles.pt > 10*GeV && InDetTrackParticles.DFCommonTightPrimary && abs(DFCommonInDetTrackZ0AtPV) < 5.0*mm )",
                                                                       ConeSize = 0.6,
                                                                       )
ToolSvc += JETM12CaloClusterThinning
thinningTools.append(JETM12CaloClusterThinning)
###################################

# Truth particle thinning
doTruthThinning = True
from AthenaCommon.GlobalFlags import globalflags
if doTruthThinning and DerivationFrameworkIsMonteCarlo:
    truth_cond_status    = "( (TruthParticles.status == 1) && (TruthParticles.barcode < 200000) && (TruthParticles.pt > 8*GeV) )"            # high pt pions for E/p
    truth_cond_Lepton = "((abs(TruthParticles.pdgId) >= 11) && (abs(TruthParticles.pdgId) <= 16) && (TruthParticles.barcode < 200000))" # Leptons
    truth_expression = '('+truth_cond_status+' || '+truth_cond_Lepton +')'
  
    from DerivationFrameworkMCTruth.DerivationFrameworkMCTruthConf import DerivationFramework__GenericTruthThinning
    JETM12TruthThinningTool = DerivationFramework__GenericTruthThinning( name = "JETM12TruthThinningTool",
                                                                        StreamName              = streamName,
                                                                        ParticleSelectionString = truth_expression,
                                                                        PreserveDescendants     = False,
                                                                        PreserveGeneratorDescendants = True,
                                                                        PreserveAncestors = False)
  
    ToolSvc += JETM12TruthThinningTool
    thinningTools.append(JETM12TruthThinningTool)

#=======================================
# CREATE PRIVATE SEQUENCE
#=======================================

jetm12Seq = CfgMgr.AthSequencer("JETM12Sequence")
DerivationFrameworkJob += jetm12Seq

#=======================================
# CREATE THE DERIVATION KERNEL ALGORITHM   
#=======================================
from DerivationFrameworkCore.DerivationFrameworkCoreConf import DerivationFramework__DerivationKernel
jetm12Seq += CfgMgr.DerivationFramework__DerivationKernel("JETM12KernelSkim",
                                                         SkimmingTools = [JETM12SkimmingTool],
                                                         AugmentationTools = AugmentationToolsSkim)
jetm12Seq += CfgMgr.DerivationFramework__DerivationKernel("JETM12Kernel",
                                                         ThinningTools = thinningTools,
                                                         AugmentationTools = AugmentationTools)


OutputJets["JETM12"] = []

if DerivationFrameworkIsMonteCarlo:
  from DerivationFrameworkMCTruth.MCTruthCommon import addStandardTruthContents
  addStandardTruthContents()

#====================================================================
# Add the containers to the output stream - slimming done here
#====================================================================
from DerivationFrameworkCore.SlimmingHelper import SlimmingHelper
JETM12SlimmingHelper = SlimmingHelper("JETM12SlimmingHelper")
JETM12SlimmingHelper.SmartCollections = ["EventInfo",
                                         "Electrons", "Photons", "Muons", "TauJets",
                                         "InDetTrackParticles", "PrimaryVertices",
                                         "MET_Baseline_AntiKt4EMTopo",
                                         "MET_Baseline_AntiKt4EMPFlow",
                                         "AntiKt4EMTopoJets","AntiKt4EMPFlowJets","AntiKt4TruthJets",
                                         "BTagging_AntiKt4EMPFlow"]

JETM12SlimmingHelper.AllVariables = ["MuonTruthParticles","TruthParticles", "TruthVertices",
                                     "MuonSegments","InDetTrackParticles",
                                     "Kt4EMTopoOriginEventShape","Kt4EMPFlowEventShape","MET_Truth","CaloCalTopoClusters",
                                     "TruthMuons","TruthElectrons","TruthPhotons","TruthTaus","TruthNeutrinos",
                                    ]

JETM12SlimmingHelper.ExtraVariables = ["InDetTrackParticles.TrkIsoPt1000_ptcone40.TrkIsoPt1000_ptcone30.TrkIsoPt1000_ptcone20.TrkIsoPt500_ptcone40.TrkIsoPt500_ptcone30.TrkIsoPt500_ptcone20"]

JETM12SlimmingHelper.AppendContentToStream(JETM12Stream)
