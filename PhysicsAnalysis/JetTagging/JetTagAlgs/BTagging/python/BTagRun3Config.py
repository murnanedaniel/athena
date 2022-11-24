# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory
from AthenaConfiguration.Enums import BeamType
from BTagging.JetParticleAssociationAlgConfig import JetParticleAssociationAlgCfg
from BTagging.JetBTaggingAlgConfig import JetBTaggingAlgCfg
from BTagging.JetSecVertexingAlgConfig import JetSecVertexingAlgCfg
from BTagging.JetSecVtxFindingAlgConfig import JetSecVtxFindingAlgCfg
from BTagging.BTagTrackAugmenterAlgConfig import BTagTrackAugmenterAlgCfg
from FlavorTagDiscriminants.BTagJetAugmenterAlgConfig import (
    BTagJetAugmenterAlgCfg)
from FlavorTagDiscriminants.BTagMuonAugmenterAlgConfig import (
    BTagMuonAugmenterAlgCfg)
from FlavorTagDiscriminants.FlavorTagNNConfig import FlavorTagNNCfg
from JetTagCalibration.JetTagCalibConfig import JetTagCalibCfg
from BTagging.BTaggingFlags import BTaggingFlags
from BTagging.BTaggingConfiguration import getConfiguration
from OutputStreamAthenaPool.OutputStreamConfig import addToESD, addToAOD

# this is where you add the new trainings!
def GetTaggerTrainingMap(jet_collection_list):
    derivationTrainingMap = {
        "AntiKt4EMPFlow": [
            "BTagging/201903/rnnip/antikt4empflow/network.json",
            "BTagging/201903/dl1r/antikt4empflow/network.json",
            "BTagging/20210519r22/dl1r/antikt4empflow/network.json",
            "BTagging/20210729/dipsLoose/antikt4empflow/network.json",  # old r22 trainings
            "BTagging/20210729/dips/antikt4empflow/network.json",
            "BTagging/20210824r22/dl1dLoose/antikt4empflow/network.json",  # “recommended tagger” which is DL1dLoose20210824r22 named DL1dv00 in EDM
            "BTagging/20210824r22/dl1d/antikt4empflow/network.json",
            "BTagging/20210824r22/dl1r/antikt4empflow/network.json",
            "BTagging/20220314/dipsLoose/antikt4empflow/network.json",  # new r22 training
            "BTagging/20220509/dl1dLoose/antikt4empflow/network.json",  # new "recommended tagger" named DL1dv01 in EDM
            "BTagging/20220509/gn1/antikt4empflow/network.onnx",
        ],
        "AntiKt4EMPFlowCustomVtx": [
            "BTagging/201903/rnnip/antikt4empflow/network.json",
            "BTagging/201903/dl1r/antikt4empflow/network.json",
        ],
        "AntiKt4EMTopo": [
            "BTagging/201903/rnnip/antikt4empflow/network.json",
            "BTagging/201903/dl1r/antikt4empflow/network.json",
            "BTagging/20210519r22/dl1r/antikt4empflow/network.json",
            "BTagging/20210729/dipsLoose/antikt4empflow/network.json",  # old r22 trainings
            "BTagging/20210729/dips/antikt4empflow/network.json",
            "BTagging/20210824r22/dl1dLoose/antikt4empflow/network.json",  # “recommended tagger” which is DL1dLoose20210824r22 named DL1dv00 in EDM
            "BTagging/20210824r22/dl1d/antikt4empflow/network.json",
            "BTagging/20210824r22/dl1r/antikt4empflow/network.json",
            "BTagging/20220314/dipsLoose/antikt4empflow/network.json",  # new r22 training
            "BTagging/20220509/dl1dLoose/antikt4empflow/network.json",  # new "recommended tagger" named DL1dv01 in EDM
        ],
        "AntiKtVR30Rmax4Rmin02Track": [
            "BTagging/201903/rnnip/antiktvr30rmax4rmin02track/network.json",
            "BTagging/201903/dl1r/antiktvr30rmax4rmin02track/network.json",
            "BTagging/20210519r22/dl1r/antikt4empflow/network.json",
            "BTagging/20210729/dipsLoose/antikt4empflow/network.json",  # old r22 trainings
            "BTagging/20210729/dips/antikt4empflow/network.json",
            "BTagging/20210824r22/dl1dLoose/antikt4empflow/network.json",  # “recommended tagger” which is DL1dLoose20210824r22 named DL1dv00 in EDM
            "BTagging/20210824r22/dl1d/antikt4empflow/network.json",
            "BTagging/20210824r22/dl1r/antikt4empflow/network.json",
            "BTagging/20220314/dipsLoose/antikt4empflow/network.json",  # new r22 training
            "BTagging/20220509/dl1dLoose/antikt4empflow/network.json",  # new "recommended tagger" named DL1dv01 in EDM
        ],
    }

    return derivationTrainingMap[jet_collection_list]


def RetagRenameInputContainerCfg(suffix, JetCollectionShort, tracksKey='InDetTrackParticles', addRenameMaps=None):
    
    acc=ComponentAccumulator()
    AddressRemappingSvc, ProxyProviderSvc=CompFactory.getComps("AddressRemappingSvc","ProxyProviderSvc",)
    AddressRemappingSvc = AddressRemappingSvc("AddressRemappingSvc")
    AddressRemappingSvc.TypeKeyRenameMaps += ['xAOD::JetAuxContainer#' + JetCollectionShort + 'Jets.BTagTrackToJetAssociator->' + JetCollectionShort + 'Jets.BTagTrackToJetAssociator_' + suffix]
    AddressRemappingSvc.TypeKeyRenameMaps += ['xAOD::JetAuxContainer#' + JetCollectionShort + 'Jets.JFVtx->' + JetCollectionShort + 'Jets.JFVtx_' + suffix]
    AddressRemappingSvc.TypeKeyRenameMaps += ['xAOD::JetAuxContainer#' + JetCollectionShort + 'Jets.SecVtx->' + JetCollectionShort + 'Jets.SecVtx_' + suffix]
    AddressRemappingSvc.TypeKeyRenameMaps += ['xAOD::JetAuxContainer#' + JetCollectionShort + 'Jets.btaggingLink->' + JetCollectionShort + 'Jets.btaggingLink_' + suffix]
    AddressRemappingSvc.TypeKeyRenameMaps += ['xAOD::BTaggingContainer#BTagging_' + JetCollectionShort + '->BTagging_' + JetCollectionShort + '_' + suffix]
    AddressRemappingSvc.TypeKeyRenameMaps += ['xAOD::BTaggingAuxContainer#BTagging_' + JetCollectionShort + 'Aux.->BTagging_' + JetCollectionShort + '_' + suffix+"Aux."]
    AddressRemappingSvc.TypeKeyRenameMaps += ['xAOD::VertexContainer#BTagging_' + JetCollectionShort + 'SecVtx->BTagging_' + JetCollectionShort + 'SecVtx_' + suffix]
    AddressRemappingSvc.TypeKeyRenameMaps += ['xAOD::VertexAuxContainer#BTagging_' + JetCollectionShort + 'SecVtxAux.->BTagging_' + JetCollectionShort + 'SecVtx_' + suffix+"Aux."]
    AddressRemappingSvc.TypeKeyRenameMaps += ['xAOD::BTagVertexContainer#BTagging_' + JetCollectionShort + 'JFVtx->BTagging_' + JetCollectionShort + 'JFVtx_' + suffix]
    AddressRemappingSvc.TypeKeyRenameMaps += ['xAOD::BTagVertexAuxContainer#BTagging_' + JetCollectionShort + 'JFVtxAux.->BTagging_' + JetCollectionShort + 'JFVtx_' + suffix+"Aux."]
    AddressRemappingSvc.TypeKeyRenameMaps += ['xAOD::TrackParticleAuxContainer#' + tracksKey + '.TrackCompatibility->' + tracksKey + '.TrackCompatibility_' + suffix]
    AddressRemappingSvc.TypeKeyRenameMaps += ['xAOD::TrackParticleAuxContainer#' + tracksKey + '.btagIp_d0->' + tracksKey + '.btagIp_d0_' + suffix]
    AddressRemappingSvc.TypeKeyRenameMaps += ['xAOD::TrackParticleAuxContainer#' + tracksKey + '.btagIp_z0SinTheta->' + tracksKey + '.btagIp_z0SinTheta_' + suffix]
    AddressRemappingSvc.TypeKeyRenameMaps += ['xAOD::TrackParticleAuxContainer#' + tracksKey + '.btagIp_d0Uncertainty->' + tracksKey + '.btagIp_d0Uncertainty_' + suffix]
    AddressRemappingSvc.TypeKeyRenameMaps += ['xAOD::TrackParticleAuxContainer#' + tracksKey + '.btagIp_z0SinThetaUncertainty->' + tracksKey + '.btagIp_z0SinThetaUncertainty_' + suffix]
    AddressRemappingSvc.TypeKeyRenameMaps += ['xAOD::TrackParticleAuxContainer#' + tracksKey + '.btagIp_trackMomentum->' + tracksKey + '.btagIp_trackMomentum_' + suffix]
    AddressRemappingSvc.TypeKeyRenameMaps += ['xAOD::TrackParticleAuxContainer#' + tracksKey + '.btagIp_trackDisplacement->' + tracksKey + '.btagIp_trackDisplacement_' + suffix]
    AddressRemappingSvc.TypeKeyRenameMaps += ['xAOD::TrackParticleAuxContainer#' + tracksKey + '.JetFitter_TrackCompatibility_antikt4empflow->' + tracksKey + '.JetFitter_TrackCompatibility_antikt4empflow_' + suffix]
    AddressRemappingSvc.TypeKeyRenameMaps += ['xAOD::JetAuxContainer#' + JetCollectionShort + 'Jets.TracksForBTagging->' + JetCollectionShort + 'Jets.TracksForBTagging' + suffix]
    AddressRemappingSvc.TypeKeyRenameMaps += ['xAOD::JetAuxContainer#' + JetCollectionShort + 'Jets.TracksForBTaggingOverPtThreshold->' + JetCollectionShort + 'Jets.TracksForBTaggingOverPtThreshold' + suffix]
    AddressRemappingSvc.TypeKeyRenameMaps += ['xAOD::JetAuxContainer#' + JetCollectionShort + 'Jets.MuonsForBTagging->' + JetCollectionShort + 'Jets.MuonsForBTagging' + suffix]
    AddressRemappingSvc.TypeKeyRenameMaps += ['xAOD::JetAuxContainer#' + JetCollectionShort + 'Jets.MuonsForBTaggingOverPtThreshold->' + JetCollectionShort + 'Jets.MuonsForBTaggingOverPtThreshold' + suffix]
    
    # add extra mappings if present
    if addRenameMaps:
        AddressRemappingSvc.TypeKeyRenameMaps += addRenameMaps
    
    acc.addService(AddressRemappingSvc)
    acc.addService(ProxyProviderSvc(ProviderNames = [ "AddressRemappingSvc" ]))

    return acc


def BTagRecoSplitCfg(inputFlags, JetCollection=['AntiKt4EMTopo','AntiKt4EMPFlow']):

    result=ComponentAccumulator()

    # Can only configure b-tagging for collisions; not cosmics, etc.
    if inputFlags.Beam.Type is not BeamType.Collisions:
        return result

    result.merge(JetTagCalibCfg(inputFlags))

    #Track Augmenter
    result.merge(BTagTrackAugmenterAlgCfg(inputFlags))

    for jc in JetCollection:
        result.merge(
            BTagAlgsCfg(
                inputFlags,
                JetCollection=jc,
                nnList=GetTaggerTrainingMap(jc),
                muons='', # muon augmentation isn't thread safe, disable
                renameTrackJets=True,
                AddedJetSuffix='Jets'
            )
        )

    # By default, in Run3 we don't write out BTagging containers in AOD or ESD
    # following allows to write them out when using Reco_tf.py --CA run 3 style configuration
    
    if inputFlags.Output.doWriteAOD and inputFlags.Jet.WriteToAOD:
     result.merge(addBTagToOutput(inputFlags, JetCollection, toAOD=True, toESD=False))     
     
    if inputFlags.Output.doWriteESD:
     result.merge(addBTagToOutput(inputFlags, JetCollection, toAOD=False, toESD=True))
    
    # Invoking the alhorithm saving hits in the vicinity of jets, with proper flags
    if inputFlags.BTagging.Trackless:
        BTaggingFlags.DoJetHitAssociation=True
        from JetHitAssociation.JetHitAssociationConfig import JetHitAssociationCfg
        result.merge(JetHitAssociationCfg(inputFlags))
        BTaggingAODList = ['xAOD::TrackMeasurementValidationContainer#JetAssociatedPixelClusters',
                           'xAOD::TrackMeasurementValidationAuxContainer#JetAssociatedPixelClustersAux.']
        BTaggingAODList += ['xAOD::TrackMeasurementValidationContainer#JetAssociatedSCTClusters',
                            'xAOD::TrackMeasurementValidationAuxContainer#JetAssociatedSCTClustersAux.']
        result.merge(addToAOD(inputFlags, BTaggingAODList))

    return result


def BTagAlgsCfg(inputFlags,
                JetCollection,
                nnList=[],
                TaggerList=None,
                SecVertexers=None,
                trackCollection='InDetTrackParticles',
                primaryVertices='PrimaryVertices',
                muons='Muons',
                BTagCollection=None,
                renameTrackJets=False,
                AddedJetSuffix=''):

    # If things aren't specified in the arguments, we'll read them
    # from the config flags
    if TaggerList is None:
        TaggerList = inputFlags.BTagging.taggerList
    if SecVertexers is None:
        SecVertexers = ['JetFitter', 'SV1']
        if inputFlags.BTagging.RunFlipTaggers:
            SecVertexers += ['JetFitterFlip','SV1Flip']
    jet = JetCollection
    jetcol_no_suffix = JetCollection
    jetcol = JetCollection + AddedJetSuffix
    if renameTrackJets is True:
        jetcol_no_suffix = jet.replace("Track", "PV0Track")
        jetcol = jetcol.replace("Track", "PV0Track")

    if BTagCollection is None:
        BTagCollection = inputFlags.BTagging.OutputFiles.Prefix + jet

    # Names of element link vectors that are stored on the jet and
    # BTagging object. These are added and read out by the packages
    # that are configured below: in principal you should be able to
    # change these without changing the final b-tagging output.
    JetTrackAssociator = 'TracksForBTagging'
    BTagTrackAssociator = 'BTagTrackToJetAssociator'
    JetMuonAssociator = 'MuonsForBTagging'
    BTagMuonAssociator = 'Muons'

    # List of input VxSecVertexInfo containers
    VxSecVertexInfoNameList = []

    #List of secondary vertex finders
    secVtxFinderxAODBaseNameList = []
    result = ComponentAccumulator()

    # Associate tracks to the jet
    result.merge(JetParticleAssociationAlgCfg(
        inputFlags,
        jetcol,
        trackCollection,
        JetTrackAssociator,
    ))

    if muons:
        result.merge(JetParticleAssociationAlgCfg(
            inputFlags, jetcol, muons, JetMuonAssociator))

    # Build secondary vertices
    for sv in SecVertexers:
        BTagVxSecVertexInfoName = sv + 'VxSecVertexInfo_' + jet
        VxSecVertexInfoNameList.append(BTagVxSecVertexInfoName)
        secVtxFinderxAODBaseNameList.append(sv)
        AlgName = (jet + '_' + sv).lower()
        result.merge(JetSecVtxFindingAlgCfg(
            inputFlags,
            BTagVxSecVertexInfoName = BTagVxSecVertexInfoName,
            SVAlgName = AlgName + '_secvtxfinding',
            JetCollection = jetcol,
            PrimaryVertexCollectionName = primaryVertices,
            SVFinder = sv,
            TracksToTag = JetTrackAssociator,
        ))
        result.merge(JetSecVertexingAlgCfg(
            inputFlags,
            BTagVxSecVertexInfoName = BTagVxSecVertexInfoName,
            SVAlgName = AlgName + '_secvtx',
            BTaggingCollection = BTagCollection,
            JetCollection = jetcol,
            TrackCollection = trackCollection,
            PrimaryVertexCollectionName = primaryVertices,
            SVFinder = sv, 
        ))

    # Create the b-tagging object, and run the older b-tagging algorithms
    secVtxFinderTrackNameList = [ BTagTrackAssociator ] * len(SecVertexers)
    result.merge(
        JetBTaggingAlgCfg(
            inputFlags,
            BTaggingCollection=BTagCollection,
            JetCollection=jetcol,
            JetColNoJetsSuffix=jetcol_no_suffix,
            PrimaryVertexCollectionName=primaryVertices,
            TaggerList=TaggerList,
            Tracks=JetTrackAssociator,
            Muons=JetMuonAssociator if muons else '',
            VxSecVertexInfoNameList = VxSecVertexInfoNameList,
            secVtxFinderxAODBaseNameList = secVtxFinderxAODBaseNameList,
            secVtxFinderTrackNameList = secVtxFinderTrackNameList,
            OutgoingTracks=BTagTrackAssociator,
            OutgoingMuons=BTagMuonAssociator,
        )
    )

    # Add some high level information to the b-tagging object we
    # created above
    result.merge(
        BTagJetAugmenterAlgCfg(
            inputFlags,
            BTagCollection=BTagCollection,
            Associator=BTagTrackAssociator,
            TrackCollection=trackCollection,
        )
    )
    
    #add also Flip tagger information
    if inputFlags.BTagging.RunFlipTaggers:
       result.merge(
           BTagJetAugmenterAlgCfg(
               inputFlags,
               BTagCollection=BTagCollection,
               Associator=BTagTrackAssociator,
               TrackCollection=trackCollection,
               doFlipTagger=True,
           )
       ) 

    
    if muons:
        result.merge(
            BTagMuonAugmenterAlgCfg(
                inputFlags,
                BTagCollection=BTagCollection,
                Associator=BTagMuonAssociator,
                MuonCollection=muons,
            )
        )

    # Add the final taggers based on neural networks
    for dl2 in nnList:
        result.merge(
            FlavorTagNNCfg(
                inputFlags,
                BTagCollection,
                TrackCollection=trackCollection,
                NNFile=dl2)
        )
        # add flip taggers, sometimes
        if inputFlags.BTagging.RunFlipTaggers:
            for flip_config in _get_flip_config(dl2):
                result.merge(
                    FlavorTagNNCfg(
                        inputFlags,
                        BTaggingCollection=BTagCollection,
                        TrackCollection=trackCollection,
                        NNFile=dl2,
                        FlipConfig=flip_config,
                    )
                )

    return result


def _get_flip_config(nn_path):
    """
    Schedule NN-based IP 'flip' taggers (rnnipflip and dipsflip) -
    this should for the moment only run on the low-level taggers and
    not on 'dl1x'.

    FlipConfig is "STANDARD" by default - for flip tagger set up with
    option "NEGATIVE_IP_ONLY" (flip sign of d0 and use only (flipped)
    positive d0 values).

    Returns a list of flip configurations, or [] for things we don't flip.
    """
    #flipping of DL1r with 2019 taggers does not work at the moment
    if (('dl1d' in nn_path) or ('dl1r' in nn_path and '201903' not in nn_path)):
        return ['FLIP_SIGN']
    if 'rnnip' in nn_path or 'dips' in nn_path or 'gn1' in nn_path:
        return ['NEGATIVE_IP_ONLY']
    else:
        return []


def addBTagToOutput(inputFlags, JetCollectionList, toAOD=True, toESD=True):
    """Write out the BTagging containers as defined by JetCollectionList
    In Run3 we don't write out BTagging in AOD or ESD : this function is for convenience and testing purpose.
    """
    result = ComponentAccumulator()

    BTaggingAODList =  BTaggingFlags.btaggingAODList

    BTagConf = getConfiguration()
    for coll in JetCollectionList:
      BTagConf.RegisterOutputContainersForJetCollection(coll)

    BTaggingAODList = BTaggingFlags.btaggingAODList if toAOD else []
    BTaggingESDList = BTaggingFlags.btaggingESDList if toESD else []

    if toESD:
        result.merge(addToESD(inputFlags, BTaggingESDList))
    if toAOD:
        result.merge(addToAOD(inputFlags, BTaggingAODList))

    return result
