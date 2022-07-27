# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

# MCTruthCommonConfig
# Contains the configuration for the common truth containers/decorations used in analysis DAODs
# including PHYS(LITE)

from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory
from DerivationFrameworkMCTruth.TruthDerivationToolsConfig import TruthDecayCollectionMakerCfg

def TruthMetaDataWriterCfg(flags, name):
    acc = ComponentAccumulator()
    theTruthMetaDataWriter =  CompFactory.DerivationFramework.TruthMetaDataWriter(name)
    acc.addPublicTool(theTruthMetaDataWriter)
    CommonAugmentation = CompFactory.DerivationFramework.CommonAugmentation
    acc.addEventAlgo(CommonAugmentation(f"{name}Kernel", AugmentationTools = [theTruthMetaDataWriter]))
    return acc 

def HepMCtoXAODTruthCfg(ConfigFlags):
    """Conversion of HepMC to xAOD truth"""
    acc = ComponentAccumulator()

    # Only run for MC input
    if ConfigFlags.Input.isMC is False:
        raise RuntimeError("Common MC truth building requested for non-MC input")

    # Local steering flag to identify EVNT input
    # Commented because the block it is needed for isn't working (TruthMetaData)
    isEVNT = False

    # Ensure EventInfoCnvAlg is scheduled
    if ("EventInfo#EventInfo" and "xAOD::EventInfo#EventInfo") not in ConfigFlags.Input.TypedCollections:
        from xAODEventInfoCnv.xAODEventInfoCnvConfig import EventInfoCnvAlgCfg
        acc.merge(EventInfoCnvAlgCfg(ConfigFlags, inputKey="McEventInfo", outputKey="EventInfo", disableBeamSpot=True))

    # Build truth collection if input is HepMC. Must be scheduled first to allow slimming.
    # Input file is event generator output (EVNT)
    from xAODTruthCnv.xAODTruthCnvConfigNew import GEN_EVNT2xAODCfg 
    if "McEventCollection#GEN_EVENT" in ConfigFlags.Input.TypedCollections:                  
        acc.merge(GEN_EVNT2xAODCfg(ConfigFlags,name="GEN_EVNT2xAOD",AODContainerName="GEN_EVENT"))
        isEVNT = True 
    # Input file is simulation output (HITS)
    elif "McEventCollection#TruthEvent" in ConfigFlags.Input.TypedCollections:
        acc.merge(GEN_EVNT2xAODCfg(name="GEN_EVNT2xAOD",AODContainerName="TruthEvent"))
    # Input file already has xAOD truth. Don't do anything.
    elif "xAOD::TruthEventContainer#TruthEvents" in ConfigFlags.Input.TypedCollections:
        pass
    else:
        raise RuntimeError("No recognised HepMC truth information found in the input")

    # If it isn't available, make a truth meta data object (will hold MC Event Weights)
    if "TruthMetaDataContainer#TruthMetaData" not in ConfigFlags.Input.TypedCollections and not isEVNT:
        # If we are going to be making the truth collection (isEVNT) then this will be made elsewhere
        acc.merge(TruthMetaDataWriterCfg(ConfigFlags, name = 'DFCommonTruthMetaDataWriter'))

    return acc

# TODO: this should be in the JetCommon.py, but it isn't yet CA compatible
def AddDAODJetsCACfg(flags,jetlist):

    from JetRecConfig.JetRecConfig import getJetAlgs, reOrderAlgs
    from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator

    acc = ComponentAccumulator()

    for jd in jetlist:
        algs, jetdef_i = getJetAlgs(flags, jd, True)
        algs = reOrderAlgs( [a for a in algs if a is not None])
        for a in algs:
            acc.addEventAlgo(a)

    return(acc)



# Helper for adding truth jet collections via new jet config
def AddTruthJetsCfg(ConfigFlags):
    acc = ComponentAccumulator()

    from JetRecConfig.StandardSmallRJets import AntiKt4Truth,AntiKt4TruthWZ,AntiKt4TruthDressedWZ,AntiKtVRTruthCharged
    from JetRecConfig.StandardLargeRJets import AntiKt10TruthTrimmed,AntiKt10TruthSoftDrop

    jetList = [AntiKt4Truth,AntiKt4TruthWZ,AntiKt4TruthDressedWZ,AntiKtVRTruthCharged,
               AntiKt10TruthTrimmed,AntiKt10TruthSoftDrop]

    acc.merge(AddDAODJetsCACfg(ConfigFlags,jetList))

    return acc

# Helper for scheduling the truth MET collection
def AddTruthMETCfg(ConfigFlags):

    acc = ComponentAccumulator()

    # Only do this if the truth MET is not present
    # This should handle EVNT correctly without an explicit check
    if ( "MissingETContainer#MET_Truth") not in ConfigFlags.Input.TypedCollections:
        from METReconstruction.METTruth_Cfg import METTruth_Cfg
        acc.merge(METTruth_Cfg(ConfigFlags))

    return acc

def PreJetMCTruthAugmentationsCfg(ConfigFlags, **kwargs):

    acc = ComponentAccumulator()

    augmentationToolsList = []

    # These augmentations do *not* require truth jets at all
    # If requested, add a decoration to photons that were used in the dressing

    from DerivationFrameworkMCTruth.TruthDerivationToolsConfig import ( DFCommonTruthElectronDressingToolCfg, 
    DFCommonTruthMuonDressingToolCfg, DFCommonTruthClassificationToolCfg, DFCommonTruthMuonToolCfg, DFCommonTruthElectronToolCfg,
    DFCommonTruthPhotonToolSimCfg, DFCommonTruthNeutrinoToolCfg, DFCommonTruthBottomToolCfg, DFCommonTruthTopToolCfg, 
    DFCommonTruthBosonToolCfg, DFCommonTruthBSMToolCfg, DFCommonTruthForwardProtonToolCfg, DFCommonTruthElectronIsolationTool1Cfg,
    DFCommonTruthElectronIsolationTool2Cfg, DFCommonTruthMuonIsolationTool1Cfg, DFCommonTruthMuonIsolationTool2Cfg, 
    DFCommonTruthPhotonIsolationTool1Cfg, DFCommonTruthPhotonIsolationTool2Cfg, DFCommonTruthPhotonIsolationTool3Cfg ) 
    
    # schedule the special truth building tools and add them to a common augmentation; note taus are handled separately below
    for item in [ DFCommonTruthClassificationToolCfg, DFCommonTruthMuonToolCfg, DFCommonTruthElectronToolCfg,
    DFCommonTruthPhotonToolSimCfg, DFCommonTruthNeutrinoToolCfg, DFCommonTruthBottomToolCfg, DFCommonTruthTopToolCfg,
    DFCommonTruthBosonToolCfg, DFCommonTruthBSMToolCfg, DFCommonTruthElectronIsolationTool1Cfg,
    DFCommonTruthElectronIsolationTool2Cfg, DFCommonTruthMuonIsolationTool1Cfg, DFCommonTruthMuonIsolationTool2Cfg,
    DFCommonTruthPhotonIsolationTool1Cfg, DFCommonTruthPhotonIsolationTool2Cfg, DFCommonTruthPhotonIsolationTool3Cfg]:    
        augmentationToolsList.append(acc.getPrimaryAndMerge(item()))
    augmentationToolsList.append(acc.getPrimaryAndMerge(DFCommonTruthForwardProtonToolCfg(ConfigFlags)))

    if 'decorationDressing' in kwargs:
        augmentationToolsList.append(acc.getPrimaryAndMerge(DFCommonTruthElectronDressingToolCfg(decorationName = kwargs['decorationDressing'])))
        augmentationToolsList.append(acc.getPrimaryAndMerge(DFCommonTruthMuonDressingToolCfg(decorationName = kwargs['decorationDressing'])))

    CommonAugmentation = CompFactory.DerivationFramework.CommonAugmentation
    acc.addEventAlgo(CommonAugmentation(name = "MCTruthCommonPreJetKernel", AugmentationTools = augmentationToolsList))

    return(acc)


def PostJetMCTruthAugmentationsCfg(ConfigFlags, **kwargs):

    acc = ComponentAccumulator()

    # Tau collections are built separately
    # truth tau matching needs truth jets, truth electrons and truth muons
    from DerivationFrameworkTau.TauTruthCommonConfig import TauTruthToolsCfg
    acc.merge(TauTruthToolsCfg(ConfigFlags))
    from DerivationFrameworkMCTruth.TruthDerivationToolsConfig import DFCommonTruthTauDressingToolCfg
    augmentationToolsList = [ acc.getPrimaryAndMerge(DFCommonTruthTauDressingToolCfg()) ]

    #Save the post-shower HT and MET filter values that will make combining filtered samples easier (adds to the EventInfo)
    from DerivationFrameworkMCTruth.GenFilterToolConfig import GenFilterToolCfg
    # schedule the special truth building tools and add them to a common augmentation; note taus are handled separately below
    from DerivationFrameworkMCTruth.TruthDerivationToolsConfig import DFCommonTruthDressedWZQGLabelToolCfg
    augmentationToolsList += [ acc.getPrimaryAndMerge(GenFilterToolCfg()) ,
                               acc.getPrimaryAndMerge(DFCommonTruthDressedWZQGLabelToolCfg())]

    # SUSY signal decorations
    from DerivationFrameworkSUSY.DecorateSUSYProcessConfig import IsSUSYSignalRun3
    if IsSUSYSignalRun3(ConfigFlags):
        from DerivationFrameworkSUSY.DecorateSUSYProcessConfig import DecorateSUSYProcessCfg
        augmentationToolsList += DecorateSUSYProcessCfg('MCTruthCommon')

    CommonAugmentation = CompFactory.DerivationFramework.CommonAugmentation
    acc.addEventAlgo(CommonAugmentation(name              = "MCTruthCommonPostJetKernel", 
                                        AugmentationTools = augmentationToolsList))

    # add SoW of individual SUSY final states, relies on augmentation from DecorateSUSYProcess()
    if IsSUSYSignalRun3(ConfigFlags):
        from DerivationFrameworkSUSY.SUSYWeightMetadataConfig import AddSUSYWeightsCfg
        acc.merge(AddSUSYWeightsCfg(ConfigFlags))

    return(acc)

# This adds the entirety of TRUTH3
def AddStandardTruthContentsCfg(ConfigFlags,
                                decorationDressing='dressedPhoton',
                                includeTausInDressingPhotonRemoval=False,
                                prefix=''):

    acc = ComponentAccumulator()

    # Schedule HepMC->xAOD truth conversion
    acc.merge(HepMCtoXAODTruthCfg(ConfigFlags))

    # Local flag
    isEVNT = False
    if "McEventCollection#GEN_EVENT" in ConfigFlags.Input.TypedCollections: isEVNT = True
    # Tools that must come before jets
    acc.merge(PreJetMCTruthAugmentationsCfg(ConfigFlags,decorationDressing = decorationDressing))
    # Should photons that are dressed onto taus also be removed from truth jets?
    if includeTausInDressingPhotonRemoval:
        acc.getPublicTool("DFCommonTruthTauDressingTool").decorationName=decorationDressing
    # Jets and MET
    acc.merge(AddTruthJetsCfg(ConfigFlags))
    acc.merge(AddTruthMETCfg(ConfigFlags))
    # Tools that must come after jets
    acc.merge(PostJetMCTruthAugmentationsCfg(ConfigFlags, decorationDressing = decorationDressing))
    # Add back the navigation contect for the collections we want
    acc.merge(AddTruthCollectionNavigationDecorationsCfg(["TruthElectrons", "TruthMuons", "TruthPhotons", "TruthTaus", "TruthNeutrinos", "TruthBSM", "TruthBottom", "TruthTop", "TruthBoson"], prefix=prefix))
    # Some more additions for standard TRUTH3
    acc.merge(AddBosonsAndDownstreamParticlesCfg())
    if isEVNT: acc.merge(AddLargeRJetD2Cfg())
    # Special collection for BSM particles
    acc.merge(AddBSMAndDownstreamParticlesCfg())
    # Special collection for Born leptons
    acc.merge(AddBornLeptonCollectionCfg())
    # Special collection for hard scatter (matrix element) - save TWO extra generations of particles
    acc.merge(AddHardScatterCollectionCfg(2))
    # Energy density for isolation corrections
    if isEVNT: acc.merge(AddTruthEnergyDensityCfg())

    return acc

def AddParentAndDownstreamParticlesCfg(generations=1,
                                       parents=[6],
                                       prefix='TopQuark',
                                       collection_prefix=None,
                                       rejectHadronChildren=False):
    """Configure tools for adding immediate parents and descendants"""
    acc = ComponentAccumulator()
    collection_name=collection_prefix+'WithDecay' if collection_prefix is not None else 'Truth'+prefix+'WithDecay'
    # Set up a tool to keep the W/Z/H bosons and all downstream particles
    collection_maker = acc.getPrimaryAndMerge(TruthDecayCollectionMakerCfg(name                 ='DFCommon'+prefix+'AndDecaysTool',
                                                                           NewCollectionName    = collection_name,
                                                                           PDGIDsToKeep         = parents,
                                                                           Generations          = generations,
                                                                           RejectHadronChildren = rejectHadronChildren))
    CommonAugmentation = CompFactory.DerivationFramework.CommonAugmentation
    kernel_name = 'MCTruthCommon'+prefix+'AndDecaysKernel'
    acc.addEventAlgo(CommonAugmentation(kernel_name, AugmentationTools = [collection_maker] ))
    return acc

# Next two don't seem to be used for anything...
## Add taus and their downstream particles (immediate and further decay products) in a special collection
#def addTausAndDownstreamParticles(kernel=None, generations=1):
#    return addParentAndDownstreamParticles(kernel=kernel,
#                                    generations=generations,
#                                    parents=[15],
#                                    prefix='Tau')
#
## Add W bosons and their downstream particles
#def addWbosonsAndDownstreamParticles(kernel=None, generations=1,
#                                     rejectHadronChildren=False):
#    return addParentAndDownstreamParticles(kernel=kernel,
#                                           generations=generations,
#                                           parents=[24],
#                                           prefix='Wboson',
#                                           rejectHadronChildren=rejectHadronChildren)

# Add W/Z/H bosons and their downstream particles (notice "boson" here does not include photons and gluons)
def AddBosonsAndDownstreamParticlesCfg(generations=1,
                                       rejectHadronChildren=False):
    """Add bosons and downstream particles (not photons/gluons)"""
    return AddParentAndDownstreamParticlesCfg(generations          = generations,
                                              parents              = [23,24,25],
                                              prefix               = 'Bosons',
                                              rejectHadronChildren = rejectHadronChildren)

# Following commented methods don't seem to be used for anything...

#def addBottomQuarkAndDownstreamParticles(kernel=None, generations=1, rejectHadronChildren=False):
#   return addParentAndDownstreamParticles(kernel=kernel,
#                                          generations=generations,
#                                          parents=[5],
#                                          prefix='BottomQuark',
#                                          rejectHadronChildren=rejectHadronChildren)
#
#def addTopQuarkAndDownstreamParticles(kernel=None, generations=1,
#                                      rejectHadronChildren=False):
#   return addParentAndDownstreamParticles(kernel=kernel,
#                                          generations=generations,
#                                          parents=[6],
#                                          prefix='TopQuark',
#                                          rejectHadronChildren=rejectHadronChildren)
#
## Add electrons, photons, and their downstream particles in a special collection
#def addEgammaAndDownstreamParticles(kernel=None, generations=1):
#    return addParentAndDownstreamParticles(kernel=kernel,
#                                           generations=generations,
#                                           parents=[11,22],
#                                           prefix='Egamma')
#

# Add b/c-hadrons and their downstream particles (immediate and further decay products) in a special collection
def AddHFAndDownstreamParticlesCfg(ConfigFlags, **kwargs):
    """Add b/c-hadrons and their downstream particles"""
    kwargs.setdefault("addB",True)
    kwargs.setdefault("addC",True)
    kwargs.setdefault("generations",-1)
    kwargs.setdefault("prefix",'')
    acc = ComponentAccumulator()
    # Set up a tool to keep b- and c-quarks and all downstream particles
    from DerivationFrameworkMCTruth.TruthDerivationToolsConfig import TruthDecayCollectionMakerCfg
    DFCommonHFAndDecaysTool = acc.getPrimaryAndMerge(TruthDecayCollectionMakerCfg( 
        name=kwargs['prefix']+"DFCommonHFAndDecaysTool",
        NewCollectionName=kwargs['prefix']+"TruthHFWithDecay",
        KeepBHadrons=kwargs['addB'],
        KeepCHadrons=kwargs['addC'],
        Generations=kwargs['generations']))
    acc.addEventAlgo(CompFactory.DerivationFramework.CommonAugmentation(
        kwargs['prefix']+"MCTruthCommonHFAndDecaysKernel",
        AugmentationTools = [DFCommonHFAndDecaysTool] ))
    return acc


# Add a one-vertex-per event "primary vertex" container
def AddPVCollectionCfg(ConfigFlags):
    """Add a one-vertex-per event "primary vertex" container"""
    acc = ComponentAccumulator()
    # Set up a tool to keep the primary vertices
    from DerivationFrameworkMCTruth.TruthDerivationToolsConfig import TruthPVCollectionMakerCfg
    DFCommonTruthPVCollTool = acc.getPrimaryAndMerge(TruthPVCollectionMakerCfg( 
        ConfigFlags,
        name="DFCommonTruthPVCollTool",
        NewCollectionName="TruthPrimaryVertices"))
    acc.addEventAlgo(CompFactory.DerivationFramework.CommonAugmentation(
        "MCTruthCommonTruthPVCollKernel",
        AugmentationTools = [DFCommonTruthPVCollTool] ))
    return acc

# Add a mini-collection for the hard scatter and N subsequent generations
def AddHardScatterCollectionCfg(generations=1):
    """Add a mini-collection for the hard scatter and N subsequent generations"""
    # Set up a tool to keep the taus and all downstream particles
    acc = ComponentAccumulator()
    from DerivationFrameworkMCTruth.TruthDerivationToolsConfig import HardScatterCollectionMakerCfg
    DFCommonHSCollectionTool = acc.getPrimaryAndMerge(HardScatterCollectionMakerCfg(name = "DFCommonHSCollectionTool",
                                                                                    NewCollectionName = "HardScatter",
                                                                                    Generations        = generations))
    CommonAugmentation = CompFactory.DerivationFramework.CommonAugmentation
    acc.addEventAlgo(CommonAugmentation(name              = "MCTruthCommonHSCollectionKernel",
                                        AugmentationTools = [DFCommonHSCollectionTool] ))
    return acc

# Add navigation decorations on the truth collections
def AddTruthCollectionNavigationDecorationsCfg(TruthCollections=[], prefix=''):
    """Tool to add navigation decorations on the truth collections"""
    acc = ComponentAccumulator() 
    if len(TruthCollections)==0: return
    # Set up a tool to add the navigation decorations
    from DerivationFrameworkMCTruth.TruthDerivationToolsConfig import TruthNavigationDecoratorCfg
    DFCommonTruthNavigationDecorator = acc.getPrimaryAndMerge(TruthNavigationDecoratorCfg(name             = prefix+'DFCommonTruthNavigationDecorator',
                                                                                          InputCollections = TruthCollections))
    CommonAugmentation = CompFactory.DerivationFramework.CommonAugmentation
    acc.addEventAlgo(CommonAugmentation(prefix+"MCTruthNavigationDecoratorKernel",
                                        AugmentationTools = [DFCommonTruthNavigationDecorator] ))
    return acc

# Add BSM particles and their downstream particles (immediate and further decay products) in a special collection
def AddBSMAndDownstreamParticlesCfg(generations=-1):
    """Add BSM particles and their downstream particles in a special collection"""
    acc = ComponentAccumulator()
    # Set up a tool to keep the taus and all downstream particles
     
    DFCommonBSMAndDecaysTool = acc.getPrimaryAndMerge(TruthDecayCollectionMakerCfg(name              = "DFCommonBSMAndDecaysTool",
                                                                                   NewCollectionName = "TruthBSMWithDecay",
                                                                                   KeepBSM           = True,
                                                                                   Generations       = generations))
    CommonAugmentation = CompFactory.DerivationFramework.CommonAugmentation
    acc.addEventAlgo(CommonAugmentation(name              = "MCTruthCommonBSMAndDecaysKernel",
                                        AugmentationTools = [DFCommonBSMAndDecaysTool] ))
    return acc

# Add a mini-collection for the born leptons
def AddBornLeptonCollectionCfg():
    """Add born leptons as a mini collection"""
    acc = ComponentAccumulator()
    # Set up a tool to keep the taus and all downstream particles
    from DerivationFrameworkMCTruth.TruthDerivationToolsConfig import TruthBornLeptonCollectionMakerCfg
    DFCommonBornLeptonCollTool = acc.getPrimaryAndMerge(TruthBornLeptonCollectionMakerCfg(name              = "DFCommonBornLeptonCollTool",
                                                                                          NewCollectionName ="BornLeptons"))
    CommonAugmentation = CompFactory.DerivationFramework.CommonAugmentation
    acc.addEventAlgo(CommonAugmentation("MCTruthCommonBornLeptonsKernel", AugmentationTools = [DFCommonBornLeptonCollTool] ))
    return acc

def AddLargeRJetD2Cfg():
    """Add large-R jet D2 variable"""
    #Extra classifier for D2 variable
    acc = ComponentAccumulator()
    from DerivationFrameworkMCTruth.TruthDerivationToolsConfig import TruthD2DecoratorCfg
    theTruthD2Decorator = acc.getPrimaryAndMerge(TruthD2DecoratorCfg(name            = "TruthD2Decorator",
                                                                     JetContainerKey = "AntiKt10TruthSoftDropBeta100Zcut10Jets",
                                                                     DecorationName  = "D2"))
    DerivationKernel = CompFactory.DerivationFramework.DerivationKernel
    acc.addEventAlgo(DerivationKernel("TRUTHD2Kernel", AugmentationTools = [theTruthD2Decorator] ))
    return acc

# Truth energy density tools
def AddTruthEnergyDensityCfg():
    """Truth energy density tools"""
    acc = ComponentAccumulator()
    from EventShapeTools.EventDensityConfig import configEventDensityTool
    from JetRecConfig.StandardJetConstits import stdConstitDic as cst
    EventDensityAthAlg = CompFactory.EventDensityAthAlg 
    # Algorithms for the energy density - needed only if e/gamma hasn't set things up already
    DFCommonTruthCentralEDTool = configEventDensityTool("DFCommonTruthCentralEDTool",
                                                        cst.Truth,
                                                        0.5,
                                                        AbsRapidityMax      = 1.5,
                                                        OutputContainer     = "TruthIsoCentralEventShape",
                                                       )
    acc.addPublicTool(DFCommonTruthCentralEDTool, primary = True)
    acc.addEventAlgo(EventDensityAthAlg("DFCommonTruthCentralEDAlg", EventDensityTool = DFCommonTruthCentralEDTool ))
    DFCommonTruthForwardEDTool = configEventDensityTool("DFCommonTruthForwardEDTool",
                                                        cst.Truth,
                                                        0.5,
                                                        AbsRapidityMin      = 1.5,
                                                        AbsRapidityMax      = 3.0,
                                                        OutputContainer     = "TruthIsoForwardEventShape",
                                                       )
    acc.addPublicTool(DFCommonTruthForwardEDTool, primary = True)
    acc.addEventAlgo(EventDensityAthAlg("DFCommonTruthForwardEDAlg", EventDensityTool = DFCommonTruthForwardEDTool ))

    # Now add the tool to do the decoration
    DFCommonTruthEDDecorator = CompFactory.DerivationFramework.TruthEDDecorator("DFCommonTruthEDDecorator",
                                                                                EventInfoName="EventInfo",
                                                                                EnergyDensityKeys=["TruthIsoCentralEventShape","TruthIsoForwardEventShape"],
                                                                                DecorationSuffix="_rho"
                                                                               )
    acc.addPublicTool(DFCommonTruthEDDecorator, primary = True)
    
    DerivationKernel = CompFactory.DerivationFramework.DerivationKernel
    acc.addEventAlgo(DerivationKernel("DFCommonTruthEDKernel", AugmentationTools = [DFCommonTruthEDDecorator] ))
    return acc


# Sets up modifiers to move pointers to old truth collections to new mini truth collections
def AddMiniTruthCollectionLinksCfg(ConfigFlags, **kwargs):
    """Tool to move pointers to new mini truth collections"""
    acc = ComponentAccumulator()
    kwargs.setdefault("doElectrons",True)
    kwargs.setdefault("doPhotons",True)
    kwargs.setdefault("doMuons",True) 
    aug_tools = []
    from DerivationFrameworkMCTruth.TruthDerivationToolsConfig import TruthLinkRepointToolCfg
    if kwargs['doElectrons']:
        electron_relink = acc.getPrimaryAndMerge(TruthLinkRepointToolCfg(
            ConfigFlags,
            name="ElMiniCollectionTruthLinkTool",
            RecoCollection="Electrons", 
            TargetCollections=["TruthMuons","TruthPhotons","TruthElectrons"]))
        aug_tools += [ electron_relink ]
    if kwargs['doPhotons']:
        photon_relink = acc.getPrimaryAndMerge(TruthLinkRepointToolCfg(
            ConfigFlags,
            name="PhMiniCollectionTruthLinkTool",
            RecoCollection="Photons", 
            TargetCollections=["TruthMuons","TruthPhotons","TruthElectrons"]))
        aug_tools += [ photon_relink ]
    if kwargs['doMuons']:
        muon_relink = acc.getPrimaryAndMerge(TruthLinkRepointToolCfg(
            ConfigFlags,
            name="MuMiniCollectionTruthLinkTool",
            RecoCollection="Muons", 
            TargetCollections=["TruthMuons","TruthPhotons","TruthElectrons"]))
        aug_tools += [ muon_relink ]
    acc.addEventAlgo(CompFactory.DerivationFramework.DerivationKernel(
        "MiniCollectionTruthLinkKernel",
        AugmentationTools = aug_tools ))
    return acc

def addTruth3ContentToSlimmerTool(slimmer):
    slimmer.AllVariables += [
        "MET_Truth",
        "TruthElectrons",
        "TruthMuons",
        "TruthPhotons",
        "TruthTaus",
        "TruthNeutrinos",
        "TruthBSM",
        "TruthBottom",
        "TruthTop",
        "TruthBoson",
        "TruthForwardProtons",
        "BornLeptons",
        "TruthBosonsWithDecayParticles",
        "TruthBosonsWithDecayVertices",
        "TruthBSMWithDecayParticles",
        "TruthBSMWithDecayVertices",
        "HardScatterParticles",
        "HardScatterVertices",
    ]
    slimmer.ExtraVariables += [
        "AntiKt4TruthDressedWZJets.GhostCHadronsFinalCount.GhostBHadronsFinalCount.pt.HadronConeExclTruthLabelID.ConeTruthLabelID.PartonTruthLabelID.TrueFlavor",
        "AntiKt10TruthSoftDropBeta100Zcut10Jets.pt.Tau1_wta.Tau2_wta.Tau3_wta.D2",
        "TruthEvents.Q.XF1.XF2.PDGID1.PDGID2.PDFID1.PDFID2.X1.X2.crossSection"]
