## HeadCalo Stuff
# Defined as function such that the user can change the cut level and minPt


def HeadCaloSetup( cutLevel = "TightPrimary", minPT = 100.0 ):
    try: 
        from TrkExTools.AtlasExtrapolator import AtlasExtrapolator
        from TrackToCalo.TrackToCaloConf import Trk__ParticleCaloExtensionTool
    except:
        mlog.error("could not import TrackToCaloConf.Trk__ParticleCaloExtensionTool")
        mlog.error("could not import TrkExTools.AtlasExtrapolator")
        print traceback.format_exc()
    try:
        from TrackToCalo.TrackToCaloConf import Trk__HeadCaloExtensionBuilder as HeadCaloExtensionBuilder
    except:
        mlog.error("could not import TrackToCaloConf.Trk__HeadCaloExtensionBuilder")
        print traceback.format_exc()
    try:
        from InDetTrackSelectionTool.InDetTrackSelectionToolConf import InDet__InDetTrackSelectionTool
        from InDetTrackSelectorTool.InDetTrackSelectorToolConf import InDet__InDetDetailedTrackSelectorTool
    except:
        mlog.error("could not import InDetTrackSelectionTool.InDet__InDetTrackSelectionTool")
        print traceback.format_exc()
    try:
        from AthenaCommon.AppMgr import ToolSvc
    except:
        mlog.error("could not import ToolSvc")
        print traceback.format_exc()
    try:
        from AthenaCommon.AlgSequence import AlgSequence
    except:
        mlog.error("could not import AlgSequence")
        print traceback.format_exc()

    topSequence=AlgSequence()
    
    theAtlasExtrapolator=AtlasExtrapolator(name = "HeadCaloAtlasExtrapolator")
    theAtlasExtrapolator.DoCaloDynamic = False # this turns off dynamic

    pcExtensionTool = Trk__ParticleCaloExtensionTool(Extrapolator = theAtlasExtrapolator)
    ToolSvc += pcExtensionTool

    HeadTrackCaloExtensionTool = HeadCaloExtensionBuilder(LastCaloExtentionTool = pcExtensionTool)
    TrackSelectionToolHC = InDet__InDetTrackSelectionTool(name            = "HeadCaloTrackSelectionTool",
                                                          minPt           = minPT,
                                                          CutLevel        = cutLevel)#,
                                                        #   maxD0           = 9999.*mm,
                                                        #   maxZ0           = 9999.*mm,                                                                 
                                                        #   minNPixelHits   = 2,  # PixelHits + PixelDeadSensors
                                                        #   minNSctHits     = 0,  # SCTHits + SCTDeadSensors
                                                        #   minNSiHits      = 7,  # PixelHits + SCTHits + PixelDeadSensors + SCTDeadSensors
                                                        #   minNTrtHits     = 0)
    TrackDetailedSelectionToolHC = InDet__InDetDetailedTrackSelectorTool(name = "HeadCaloDetailedTrackSelectionTool",
                                                                         pTMin                = minPT,
                                                                         IPd0Max              = 1.,
                                                                         IPz0Max              = 1.5, 
                                                                         useTrackSummaryInfo  = True,
                                                                         nHitBLayer           = 0, 
                                                                         nHitPix              = 2,  # PixelHits + PixelDeadSensors
                                                                         nHitSct              = 0,  # SCTHits + SCTDeadSensors
                                                                         nHitSi               = 7,  # PixelHits + SCTHits + PixelDeadSensors + SCTDeadSensors
                                                                         nHitTrt              = 0,  # nTRTHits
                                                                         useSharedHitInfo     = False,
                                                                         nSharedBLayer        = 99999,
                                                                         nSharedPix           = 99999,
                                                                         nSharedSct           = 99999,
                                                                         nSharedSi            = 99999,
                                                                         useTrackQualityInfo  = True,
                                                                         fitChi2OnNdfMax      = 99999,
                                                                         TrackSummaryTool     = None,
                                                                         Extrapolator         = theAtlasExtrapolator)

    ToolSvc += TrackSelectionToolHC
    ToolSvc += TrackDetailedSelectionToolHC

    HeadTrackCaloExtensionTool.TrkSelection         = TrackSelectionToolHC
    HeadTrackCaloExtensionTool.TrkDetailedSelection = TrackDetailedSelectionToolHC

    ToolSvc += HeadTrackCaloExtensionTool.LastCaloExtentionTool

    topSequence += HeadTrackCaloExtensionTool

    return True
