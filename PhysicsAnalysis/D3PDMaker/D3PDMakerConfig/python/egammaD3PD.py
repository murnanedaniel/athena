# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

#
# @file D3PDMakerConfig/python/egammaD3PD.py
# @author scott snyder <snyder@bnl.gov>
# @date Aug, 2009
# @brief Construct an egamma D3PD.
#

from D3PDMakerConfig.D3PDMakerFlags           import D3PDMakerFlags
from EventCommonD3PDMaker.EventInfoD3PDObject import EventInfoD3PDObject
from egammaD3PDMaker.ElectronD3PDObject       import ElectronD3PDObject
from egammaD3PDMaker.PhotonD3PDObject         import PhotonD3PDObject
from MuonD3PDMaker.MuonD3PDObject             import MuonD3PDObject
# from MuonD3PDMaker.MuonTriggerBitsD3PDObject  import MuonTriggerBitsD3PDObject
from JetD3PDMaker.JetD3PDObject               import JetD3PDObject
from CaloD3PDMaker.MBTSD3PDObject             import MBTSD3PDObject
from CaloD3PDMaker.MBTSTimeD3PDObject         import MBTSTimeD3PDObject
from CaloD3PDMaker.LArCollisionTimeD3PDObject import LArCollisionTimeD3PDObject
from CaloD3PDMaker.CollisionDecisionD3PDObject import CollisionDecisionD3PDObject
# from EventCommonD3PDMaker.LBMetadataConfig    import LBMetadataConfig
from TruthD3PDMaker.TruthEventD3PDObject      import TruthEventD3PDObject
# from TruthD3PDMaker.GenEventD3PDObject        import GenEventD3PDObject
from TruthD3PDMaker.TruthParticleD3PDObject   import TruthParticleD3PDObject
from TrackD3PDMaker.xAODVertexD3PDObject      import PrimaryxAODVertexD3PDObject
from TrackD3PDMaker.xAODTrackD3PDObject       import xAODTrackParticleD3PDObject
#from TrigEgammaD3PDMaker.TrigEgammaD3PD       import TrigEgammaD3PDObjects
from RecExConfig.RecFlags                     import rec
from RecExConfig.ObjKeyStore                  import cfgKeyStore


from AthenaCommon.AlgSequence import AlgSequence
topSequence = AlgSequence()


def _args (level, name, kwin, **kw):
    kw = kw.copy()
    kw['level'] = level
    for (k, v) in kwin.items():
        if k.startswith (name + '_'):
            kw[k[len(name)+1:]] = v
    return kw


def createxAOD (seq):
    from D3PDMakerConfig.makexAOD import makexAOD
    makexAOD (seq, 'xAOD::TruthParticleContainer',
              'TruthParticles', 'GEN_AOD',
              xaod_truth_event_key = 'TruthEvents',
              xaod_truth_vertex_key = 'TruthVertices',
              xaod_truth_links_key = 'xAODTruthLinks')
    makexAOD (seq, 'xAOD::CaloClusterContainer',  'CaloCalTopoCluster')
    makexAOD (seq, 'xAOD::CaloClusterContainer',  'egClusterCollection')
    makexAOD (seq, 'xAOD::CaloClusterContainer',  'LArClusterEMFrwd')
    makexAOD (seq, 'xAOD::CaloClusterContainer',  'LArClusterEMSofte')
    makexAOD (seq, 'xAOD::CaloClusterContainer',  'EMTopoSW35')
    makexAOD (seq, 'xAOD::CaloClusterContainer',  'HLT')
    makexAOD (seq, 'xAOD::CaloClusterContainer',  'HLT_TrigCaloClusterMaker')
    makexAOD (seq, 'xAOD::CaloClusterContainer',  'HLT_TrigCaloClusterMaker_slw')
    makexAOD (seq, 'xAOD::CaloClusterContainer',  'MuonClusterCollection')
    makexAOD (seq, 'xAOD::CaloClusterContainer',  'Tau1P3PCellCluster')
    makexAOD (seq, 'xAOD::CaloClusterContainer',  'Tau1P3PEM012ClusterContainer')
    makexAOD (seq, 'xAOD::CaloClusterContainer',  'Tai1P3PPi0ClusterContainer')
    makexAOD (seq, 'xAOD::JetContainer',          'AntiKt4TopoEMJets')
    makexAOD (seq, 'xAOD::JetContainer',          'AntiKt4TruthJets')

    #makexAOD (seq, 'xAOD::TrackParticleContainer',
    #          'InDetTrackParticles', 'Tracks',
    #          truth_key = 'TrackTruthCollection')
    makexAOD (seq, 'xAOD::TrackParticleContainer',
              'GSFInDetTrackParticles', 'GSFTracks',
              truth_key = 'GSFTrackTruthCollection')

    makexAOD (seq, 'xAOD::TrackParticleContainer',
              'InDetTrackParticles', 'TrackParticleCandidate',
              truth_key = 'TrackParticleTruthCollection')
    makexAOD (seq, 'xAOD::TrackParticleContainer',
              'GSFTrackParticles', 'GSFTrackParticleCandidate',
              truth_key = 'GSFTrackParticleTruthCollection',
              trackmap_key = 'GSFTrackAssociation',
              base_xaod_key = 'InDetTrackParticles',
              base_key = 'TrackParticleCandidate')

    makexAOD (seq, 'xAOD::TrackParticleContainer',
              'MuonboyTrackParticles',
              truth_key = 'TrackParticleTruthCollection')
    makexAOD (seq, 'xAOD::TrackParticleContainer',
              'StacoTrackParticles',
              truth_key = 'TrackParticleTruthCollection')
    makexAOD (seq, 'xAOD::TrackParticleContainer',
              'MuidExtrTrackParticles',
              truth_key = 'TrackParticleTruthCollection')
    makexAOD (seq, 'xAOD::TrackParticleContainer',
              'MuidCombTrackParticles',
              truth_key = 'TrackParticleTruthCollection')
    makexAOD (seq, 'xAOD::TrackParticleContainer',
              'ExtrapolatedMuonSpectrometerParticles',
              truth_key = 'TrackParticleTruthCollection')
    makexAOD (seq, 'xAOD::TrackParticleContainer',
              'CombinedfitMuonparticles',
              truth_key = 'TrackParticleTruthCollection')
    makexAOD (seq, 'xAOD::TrackParticleContainer',
              'HLT_InDetTrigParticleCreation_Electron_EFID',
              truth_key = '')
    makexAOD (seq, 'xAOD::TrackParticleContainer',
              'HLT_InDetTrigParticleCreationTRTOnly_Electron_EFID',
              truth_key = '')

    makexAOD (seq, 'xAOD::TrackParticleContainer',
              'HLT_InDetTrigParticleCreation_Muon_EFID',
              truth_key = '')
    makexAOD (seq, 'xAOD::TrackParticleContainer',
              'HLT_InDetTrigParticleCreation_FullScan_EFID',
              truth_key = '')
    makexAOD (seq, 'xAOD::TrackParticleContainer',
              'HLT_InDetTrigParticleCreation_Bphysics_EFID',
              truth_key = '')
    makexAOD (seq, 'xAOD::TrackParticleContainer',
              'HLT_InDetTrigParticleCreationTRTOnly_Muon_EFID',
              truth_key = '')
    makexAOD (seq, 'xAOD::TrackParticleContainer',
              'HLT_InDetTrigParticleCreationTRTOnly_FullScan_EFID',
              truth_key = '')

    makexAOD (seq, 'xAOD::VertexContainer', 'AllPhotonsVxCandidates',
              TPContainerName = 'GSFTrackParticles')
    makexAOD (seq, 'xAOD::VertexContainer', 'VxPrimaryCandidate',
              key = 'PrimaryVertices',
              TPContainerName = 'InDetTrackParticles')

    makexAOD (seq, 'xAOD::ElectronContainer',
              'Electrons', 'ElectronAODCollection',
              xAODContainerFrwdName = 'ForwardElectrons')
    makexAOD (seq, 'xAOD::ElectronContainer',
              'HLT_egamma_Electrons', 'HLT_egamma_Electrons',
              xAODContainerFrwdName = 'HLT_egamma_Electrons_fwd',
              forTrigger = True)
    makexAOD (seq, 'xAOD::PhotonContainer',
              'Photons', 'PhotonAODCollection',
              topo_cluster_xaod_key = 'EMTopoSW35',
              vertex_xaod_key = 'AllPhotonsVxCandidates')
    makexAOD (seq, 'xAOD::PhotonContainer',
              'HLT_egamma_Photons', 'HLT_egamma_Photons',
              topo_cluster_xaod_key = '',
              vertex_xaod_key = '',
              forTrigger = True)

    makexAOD (seq, 'xAOD::MuonContainer',
              'StacoMuonCollection',
              xaod_tp_key = 'InDetTrackParticles',
              xaod_sa_key = 'MuonboyTrackParticles',
              xaod_cb_key = 'StacoTrackParticles')
    makexAOD (seq, 'xAOD::MuonContainer',
              'MuidMuonCollection',
              xaod_tp_key = 'InDetTrackParticles',
              xaod_sa_key = 'MuidExtrTrackParticles',
              xaod_cb_key = 'MuidCombTrackParticles')
    makexAOD (seq, 'xAOD::MuonContainer',
              'CaloMuonCollection',
              xaod_tp_key = 'InDetTrackParticles',
              xaod_sa_key = '',
              xaod_cb_key = '')
    makexAOD (seq, 'xAOD::MuonContainer',
              'HLT_MuonEFInfo')
    #makexAOD (seq, 'xAOD::MuonContainer',
    #          'HLT_eMuonEFInfo')
    #makexAOD (seq, 'xAOD::MuonContainer',
    #          'Muons',
    #          xaod_tp_key = 'InDetTrackParticles',
    #          xaod_sa_key = 'ExtrapolatedMuonSpectrometerParticles',
    #          xaod_cb_key = 'CombinedfitMuonparticles')
    makexAOD (seq, 'xAOD::TrigEMClusterContainer',
              'HLT_TrigT2CaloEgamma')


    makexAOD (seq, 'xAOD::MissingETContainer',
              'MET_Core_AntiKt4EMTopo',
              'MET_RefFinal')

    makexAOD (seq, 'xAOD::TrigDecision', 'xTrigDecision', 'TrigDecision')

    from xAODTriggerCnv.xAODRoICreator import xAODRoICreator
    if cfgKeyStore.isInInput ('LVL1_ROI', 'LVL1_ROI'):
        xAODRoICreator (seq)
    return


# Make a container merging xAOD central and forward electrons.
import AthenaPython.PyAthena as PyAthena
from AthenaPython.PyAthena import StatusCode
class MergeElectrons (PyAthena.Alg):
    def __init__ (self, name = 'MergeElectrons', **kw):
        super (MergeElectrons, self).__init__ (name = name, **kw)
        return
    def execute (self):
        import ROOT
        sg=PyAthena.py_svc('StoreGateSvc')
        e1 = sg['Electrons']
        e2 = sg['ForwardElectrons']
        enew = ROOT.DataVector(ROOT.xAOD.Electron_v1) (1) #VIEW_ELEMENTS
        for e in e1: enew.push_back(e)
        for e in e2: enew.push_back(e)
        ROOT.SetOwnership (enew, False)
        sg.record (enew, 'AllElectrons')
        cfgKeyStore.addTransient ('xAOD::ElectronContainer', 'AllElectrons')

        # Make sure these aux variables are defined at this point.
        ROOT.xAOD.ElectronAuxContainer()
        ROOT.xAOD.CaloClusterAuxContainer()
               
        return StatusCode.Success
        


def egammaD3PD (alg = None,
                trigalg = None,
                level = 10,
                file = 'egamma.root',
                tuplename = 'egamma',
                **kw):

    preseq = AlgSequence (D3PDMakerFlags.PreD3PDAlgSeqName())
    createxAOD (preseq)

    if not cfgKeyStore.isInInput ('xAOD::ElectronContainer',
                                  'AllElectrons'):
        preseq += MergeElectrons()

    # Core algorithm
    # By default the user of this function should have created an algorithm already.
    # But for backwards compatibility let the function create its own algorithm if
    # needed...
    if not alg:
        from OutputStreamAthenaPool.MultipleStreamManager import MSMgr
        alg = MSMgr.NewRootStream( tuplename, file )

    alg += EventInfoD3PDObject        (**_args (level, 'EventInfo', kw))

    # Segregate trigger decision flags in a separate tree.
    if not trigalg:
        from OutputStreamAthenaPool.MultipleStreamManager import MSMgr
        trigalg = MSMgr.NewRootStream( tuplename + ':' + tuplename + 'TrigDec',
                                       FileName = file,
                                       TreeName = tuplename + 'TrigDec')

    trigalg += EventInfoD3PDObject        (0)
    alg.Stream.trigDecisionTree = trigalg

    alg += LArCollisionTimeD3PDObject (**_args (level, 'LArCollisionTime', kw))
    alg += ElectronD3PDObject         (**_args (level, 'Electron', kw,
                                                EgammaJetSignedIPAndPTRel_target='jet_'))
    alg += PhotonD3PDObject           (**_args (level, 'Photon', kw))
        

    # # Additional muon variables for Zg skims.
    # if not MuonD3PDObject.allBlocknames().has_key("MuonScatteringAngleSignificance"):
    #     from AthenaCommon.AppMgr import ToolSvc
    #     muonScatteringSigTool=None
    #     if hasattr(ToolSvc, "MuonScatteringSigTool"):
    #         muonScatteringSigTool=ToolSvc.MuonScatteringSigTool
    #     from JetTagD3PDMaker import MuonScatteringAngleSignificanceFillerTool
    #     MuonD3PDObject.defineBlock (100, "MuonScatteringAngleSignificance",
    #                                 MuonScatteringAngleSignificanceFillerTool,
    #                                 ScatteringSigTool=muonScatteringSigTool)
    alg += MuonD3PDObject             (**_args (10, 'Muons', kw,
                                                sgkey='StacoMuonCollection,Muons', prefix='mu_',
                                                include = ["EFCBInfoIndex", "EFMGInfoIndex",
                                                           "EFMEInfoIndex",
                                                           "L2CBInfoIndex", "L1InfoIndex",
                                                           "MuonScatteringAngleSignificance"],
                                                exclude = ["EFCBInfo", "EFMGInfo", "EFMEInfo",
                                                           "L2CBInfo", "L1Info"],
                                                allowMissing = True ))
    
    alg += JetD3PDObject              (**_args (0, 'Jet', kw,
                                                include=['JetQual',
                                                         'DQMoments']))
    alg += MBTSD3PDObject             (**_args (level, 'MBTS', kw))
    alg += MBTSTimeD3PDObject         (**_args (level, 'MBTSTime', kw))
    alg += CollisionDecisionD3PDObject(**_args (level, 'CollisionDecision', kw))

    # May be missing in single-beam data.
    alg += PrimaryxAODVertexD3PDObject (**_args (
        1, 'PrimaryVertex', kw,
        #allowMissing = True,
        sgkey = D3PDMakerFlags.VertexSGKey(),
        prefix = 'vxp_',
        storeVertexTrackIndexAssociation = False,
        storeVertexTrackAssociation = True,
        storeDiagonalCovarianceAsErrors = True))

    alg += xAODTrackParticleD3PDObject    (**_args (
        3, 'TrackParticleCandidate', kw,
        trackParametersAtGlobalPerigeeLevelOfDetails = 3,
        storeDiagonalCovarianceAsErrors = True))

    if rec.doTruth():
        from MuonD3PDMaker.TruthMuonD3PDObject        import TruthMuonD3PDObject
        alg += TruthMuonD3PDObject    (**_args (level, 'TruthMuon', kw))
        alg += TruthEventD3PDObject     (**_args (1, 'TruthEvent', kw))
        alg += TruthParticleD3PDObject(**_args (level, 'TruthParticle', kw))

    # alg.MetadataTools += [LBMetadataConfig()]

    # if D3PDMakerFlags.FilterCollCand():
    #     from CaloD3PDMaker.CollisionFilterAlg import CollisionFilterAlg
    #     alg.filterSeq += CollisionFilterAlg (tuplename + '_CollCandFilter')

    return alg
