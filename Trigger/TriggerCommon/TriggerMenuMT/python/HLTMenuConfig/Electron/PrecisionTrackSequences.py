#
#  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
#

from AthenaCommon.CFElements import parOR
from AthenaCommon.GlobalFlags import globalflags
import AthenaCommon.CfgMgr as CfgMgr

#logging
from AthenaCommon.Logging import logging
log = logging.getLogger( 'TriggerMenuMT.HLTMenuConfig.Egamma.PrecisionTracking')
#from TriggerMenuMT.HLTMenuConfig.Egamma.PrecisionCaloSequenceSetup import precisionCaloMenuDefs

from TriggerMenuMT.HLTMenuConfig.Egamma.EgammaDefs import TrigEgammaKeys

def precisionTracking(RoIs, precisionCaloClusters, ion=False):
## Taking Fast Track information computed in 2nd step ##
    from TrigInDetConfig.ConfigSettings import getInDetTrigConfig
    IDTrigConfig = getInDetTrigConfig( 'electron' )

    tag = '_ion' if ion is True else ''

    # TrackCollection="TrigFastTrackFinder_Tracks_Electron"
    ViewVerifyTrk = CfgMgr.AthViews__ViewDataVerifier("FastTrackViewDataVerifier" + tag)
    
    ViewVerifyTrk.DataObjects = [( 'TrigRoiDescriptorCollection' , 'StoreGateSvc+' + RoIs ),
                                 ( 'xAOD::CaloClusterContainer' , 'StoreGateSvc+' + precisionCaloClusters ),
                                 ( 'CaloCellContainer' , 'StoreGateSvc+CaloCells' ),
                                 ( 'SG::AuxElement' , 'StoreGateSvc+EventInfo.averageInteractionsPerCrossing' ),
                                 ( 'IDCInDetBSErrContainer' , 'StoreGateSvc+SCT_FlaggedCondData_TRIG' )]

    # These objects must be loaded from SGIL if not from CondInputLoader
    from AthenaCommon.AlgSequence import AlgSequence
    topSequence = AlgSequence()
    if globalflags.InputFormat.is_bytestream():
      ViewVerifyTrk.DataObjects += [( 'IDCInDetBSErrContainer' , 'StoreGateSvc+PixelByteStreamErrs' ),
                                    ( 'IDCInDetBSErrContainer' , 'StoreGateSvc+SCT_ByteStreamErrs' )]
    else:
      topSequence.SGInputLoader.Load += [( 'TRT_RDO_Container' , 'StoreGateSvc+TRT_RDOs' )]
      ViewVerifyTrk.DataObjects += [( 'TRT_RDO_Container' , 'StoreGateSvc+TRT_RDOs' )]

    """ Precision Track Related Setup.... """
    PTAlgs = []
    PTTracks = []
    PTTrackParticles = []
    
    from TrigInDetConfig.InDetPT import makeInDetPrecisionTracking

    PTTracks, PTTrackParticles, PTAlgs = makeInDetPrecisionTracking( config = IDTrigConfig, verifier = ViewVerifyTrk, rois= RoIs )
    PTSeq = parOR("precisionTrackingInElectrons" + tag, PTAlgs)
    #trackParticles = PTTrackParticles[-1]    
    trackParticles = TrigEgammaKeys.TrigElectronTracksCollectionName

    electronPrecisionTrack = parOR("electronPrecisionTrack" + tag)
    electronPrecisionTrack += ViewVerifyTrk
    electronPrecisionTrack += PTSeq

    return electronPrecisionTrack, trackParticles

