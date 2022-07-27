# Blocking the include for after first inclusion
include.block ('InDetRecExample/ConfiguredxAODTrackParticleCreation.py')

# ------------------------------------------------------------
#
# ----------- Setup InDet xAOD::TrackParticle Creation
#
# ------------------------------------------------------------

class ConfiguredxAODTrackParticleCreation:

     def __init__(self, InputTrackCollection = None, InputTrackTruthCollection = None, OutputTrackParticleContainer = None):


         from InDetRecExample.InDetJobProperties import InDetFlags
         from InDetRecExample.InDetKeys import InDetKeys
         #
         # --- get ToolSvc and topSequence
         #
         from AthenaCommon.AppMgr                import ToolSvc
         from AthenaCommon.AlgSequence           import AlgSequence
         from InDetRecExample                    import TrackingCommon
         topSequence = AlgSequence()

         _perigee_expression=InDetFlags.perigeeExpression()
         # need to treat Vertex specifically because at the time of
         # the track particle creation the primary vertex does not yet exist.
         # The problem is solved by first creating track particles wrt. the beam line
         # and correcting the parameters after the vertex finding.
         if _perigee_expression == 'Vertex' :
              _perigee_expression = 'BeamLine'

         #Always the same (so far) so can in principle go in InDetRecLoadTools
         # @TODO shoueld use InDetxAODParticleCreator.getInDetxAODParticleCreatorTool
         from TrkParticleCreator.TrkParticleCreatorConf import Trk__TrackParticleCreatorTool
         InDetxAODParticleCreatorTool = Trk__TrackParticleCreatorTool(name = "InDetxAODParticleCreatorTool"+InputTrackCollection,
                                                                      TrackToVertex           = TrackingCommon.getInDetTrackToVertexTool(),
                                                                      TrackSummaryTool        = TrackingCommon.getInDetTrackSummaryToolSharedHits(),
                                                                      TestPixelLayerTool      = TrackingCommong.getInDetTestPixelLayerToolInner(),
                                                                      ComputeAdditionalInfo   = True,
                                                                      BadClusterID            = InDetFlags.pixelClusterBadClusterID(),
                                                                      KeepParameters          = True,
                                                                      KeepFirstParameters     = InDetFlags.KeepFirstParameters(),
                                                                      PerigeeExpression       = _perigee_expression)

         ToolSvc += InDetxAODParticleCreatorTool
         if (InDetFlags.doPrintConfigurables()):
            printfunc (InDetxAODParticleCreatorTool)

         from xAODTrackingCnv.xAODTrackingCnvConf import xAODMaker__TrackCollectionCnvTool
         InDetTrackCollectionCnvTool = xAODMaker__TrackCollectionCnvTool("InDetTrackCollectionCnvTool"+InputTrackCollection,
                                                                         TrackParticleCreator = InDetxAODParticleCreatorTool)


         from xAODTrackingCnv.xAODTrackingCnvConf import xAODMaker__TrackParticleCnvAlg
         xAODTrackParticleCnvAlg = xAODMaker__TrackParticleCnvAlg(name = "InDetxAODParticleCreatorAlg"+InputTrackCollection,
                                                                  ConvertTracks = True,
                                                                  ConvertTrackParticles = False,
                                                                  TrackContainerName = InputTrackCollection,
                                                                  xAODContainerName = OutputTrackParticleContainer,
                                                                  xAODTrackParticlesFromTracksContainerName = OutputTrackParticleContainer,
                                                                  TrackParticleCreator = InDetxAODParticleCreatorTool,
                                                                  TrackCollectionCnvTool = InDetTrackCollectionCnvTool)

         if (InDetFlags.doTruth() and not InputTrackTruthCollection == ''):
             xAODTrackParticleCnvAlg.AddTruthLink = True
             xAODTrackParticleCnvAlg.TrackTruthContainerName = InputTrackTruthCollection

             from MCTruthClassifier.MCTruthClassifierBase import MCTruthClassifier
             xAODTrackParticleCnvAlg.MCTruthClassifier = MCTruthClassifier

         elif (InDetFlags.doTruth() and InputTrackTruthCollection == ''):
             printfunc ("WARNING: ConfiguredxAODTrackParticleCreation - doTruth = True, but no input Truth collection specified!")
         else:
            xAODTrackParticleCnvAlg.AddTruthLink = False

         if InDetFlags.doTIDE_AmbiTrackMonitoring():
             xAODTrackParticleCnvAlg.AugmentObservedTracks = True
             xAODTrackParticleCnvAlg.TracksMapName = InDetKeys.ObservedTracks()+"Map"

         topSequence += xAODTrackParticleCnvAlg
         if (InDetFlags.doPrintConfigurables()):
            printfunc (xAODTrackParticleCnvAlg)
