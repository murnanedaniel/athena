# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration


from __future__ import print_function

# ------------------------------------------------------------
#
# ----------- Loading of all the Tools needed to configure
#
# ------------------------------------------------------------
#
# common
from AthenaCommon.AppMgr import ToolSvc
from InDetTrigRecExample.InDetTrigFlags import InDetTrigFlags
from InDetTrigRecExample.ConfiguredNewTrackingTrigCuts import EFIDTrackingCuts
InDetTrigCutValues = EFIDTrackingCuts

from AthenaCommon.DetFlags import DetFlags
from AthenaCommon.Logging import logging 
log = logging.getLogger("InDetTrigConfigRecLoadTools.py")

from InDetTrigRecExample.InDetTrigConditionsAccess import PixelConditionsSetup, SCT_ConditionsSetup
from AthenaCommon.CfgGetter import getPublicTool,getPrivateTool
TrigPixelLorentzAngleTool = getPublicTool("PixelLorentzAngleTool")
TrigSCTLorentzAngleTool = getPrivateTool("SCTLorentzAngleTool") 

from SiLorentzAngleTool.SCTLorentzAngleToolSetup import SCTLorentzAngleToolSetup
sctLorentzAngleToolSetup = SCTLorentzAngleToolSetup()


#
# common ClusterMakerTool
#
from SiClusterizationTool.SiClusterizationToolConf import InDet__ClusterMakerTool
InDetTrigClusterMakerTool = \
    InDet__ClusterMakerTool( name = "InDetTrigClusterMakerTool",
                             PixelLorentzAngleTool = TrigPixelLorentzAngleTool,
                             SCTLorentzAngleTool = TrigSCTLorentzAngleTool
                             )
if (InDetTrigFlags.doPrintConfigurables()):
  print (InDetTrigClusterMakerTool)
ToolSvc += InDetTrigClusterMakerTool

from SiSpacePointTool.SiSpacePointToolConf import InDet__SiSpacePointMakerTool
InDetTrigSiSpacePointMakerTool = InDet__SiSpacePointMakerTool(name="InDetTrigSiSpacePointMakerTool")
if (InDetTrigFlags.doPrintConfigurables()):
  print (InDetTrigSiSpacePointMakerTool)
ToolSvc += InDetTrigSiSpacePointMakerTool      
  

#
# ----------- control loading of ROT_creator
#
if InDetTrigFlags.loadRotCreator():

  #4 clusterOnTrack Tools
  #
  from SiClusterOnTrackTool.SiClusterOnTrackToolConf import InDet__SCT_ClusterOnTrackTool
  SCT_ClusterOnTrackTool = InDet__SCT_ClusterOnTrackTool ("SCT_ClusterOnTrackTool",
                                                          CorrectionStrategy = 0,  # do correct position bias
                                                          ErrorStrategy      = 2,  # do use phi dependent errors
                                                          LorentzAngleTool   = TrigSCTLorentzAngleTool)

  ToolSvc += SCT_ClusterOnTrackTool
  if (InDetTrigFlags.doPrintConfigurables()):
    print (SCT_ClusterOnTrackTool)

  # tool to always make conservative pixel cluster errors
  from SiClusterOnTrackTool.SiClusterOnTrackToolConf import InDet__PixelClusterOnTrackTool

  if InDetTrigFlags.doPixelClusterSplitting():
    from TrkNeuralNetworkUtils.TrkNeuralNetworkUtilsConf import Trk__NeuralNetworkToHistoTool
    NeuralNetworkToHistoTool=Trk__NeuralNetworkToHistoTool(name = "NeuralNetworkToHistoTool")
      
    ToolSvc += NeuralNetworkToHistoTool
    if (InDetTrigFlags.doPrintConfigurables()):
      print (NeuralNetworkToHistoTool)
    
    from SiClusterizationTool.SiClusterizationToolConf import InDet__NnClusterizationFactory
    from AtlasGeoModel.CommonGMJobProperties import CommonGeometryFlags as geoFlags
    do_runI = geoFlags.Run() not in ["RUN2", "RUN3"]
    from InDetRecExample.TrackingCommon import createAndAddCondAlg,getPixelClusterNnCondAlg,getPixelClusterNnWithTrackCondAlg
    createAndAddCondAlg( getPixelClusterNnCondAlg,         'PixelClusterNnCondAlg',          GetInputsInfo = do_runI)
    createAndAddCondAlg( getPixelClusterNnWithTrackCondAlg,'PixelClusterNnWithTrackCondAlg', GetInputsInfo = do_runI)
    if do_runI :
      TrigNnClusterizationFactory = InDet__NnClusterizationFactory( name                 = "TrigNnClusterizationFactory",
                                                                    PixelLorentzAngleTool              = TrigPixelLorentzAngleTool,
                                                                    doRunI                             = True,
                                                                    useToT                             = False,
                                                                    useRecenteringNNWithoutTracks      = True,
                                                                    useRecenteringNNWithTracks         = False,
                                                                    correctLorShiftBarrelWithoutTracks = 0,
                                                                    correctLorShiftBarrelWithTracks    = 0.030,
                                                                    NnCollectionReadKey                = 'PixelClusterNN',
                                                                    NnCollectionWithTrackReadKey       = 'PixelClusterNNWithTrack')
    else:
        TrigNnClusterizationFactory = InDet__NnClusterizationFactory( name                         = "TrigNnClusterizationFactory",
                                                                      PixelLorentzAngleTool        = TrigPixelLorentzAngleTool,
                                                                      useToT                       = InDetTrigFlags.doNNToTCalibration(),
                                                                      NnCollectionReadKey          = 'PixelClusterNN',
                                                                      NnCollectionWithTrackReadKey = 'PixelClusterNNWithTrack')

    ToolSvc += TrigNnClusterizationFactory

  else:
    TrigNnClusterizationFactory = None

  if (InDetTrigFlags.doPrintConfigurables()):
    print (TrigNnClusterizationFactory)

  InDetTrigPixelClusterOnTrackTool = InDet__PixelClusterOnTrackTool("InDetTrigPixelClusterOnTrackTool",
                                                                    ErrorStrategy = 2,
                                                                    LorentzAngleTool = TrigPixelLorentzAngleTool,
                                                                    NnClusterizationFactory= TrigNnClusterizationFactory,
  )

  ToolSvc += InDetTrigPixelClusterOnTrackTool

  if (InDetTrigFlags.doPrintConfigurables()):
    print (InDetTrigPixelClusterOnTrackTool)

  # tool to always make conservative sct cluster errors
  from SiClusterOnTrackTool.SiClusterOnTrackToolConf import InDet__SCT_ClusterOnTrackTool
  InDetTrigBroadSCT_ClusterOnTrackTool = InDet__SCT_ClusterOnTrackTool ("InDetTrigBroadSCT_ClusterOnTrackTool",
                                         CorrectionStrategy = 0,  # do correct position bias
                                         ErrorStrategy      = 0,  # do use broad errors
                                         LorentzAngleTool   = TrigSCTLorentzAngleTool)
  ToolSvc += InDetTrigBroadSCT_ClusterOnTrackTool
  if (InDetTrigFlags.doPrintConfigurables()):
    print (InDetTrigBroadSCT_ClusterOnTrackTool)

  #--
  InDetTrigBroadPixelClusterOnTrackTool = InDet__PixelClusterOnTrackTool("InDetTrigBroadPixelClusterOnTrackTool",
                                                                         ErrorStrategy = 0,
                                                                         LorentzAngleTool = TrigPixelLorentzAngleTool,
                                                                         NnClusterizationFactory= TrigNnClusterizationFactory
  )
  ToolSvc += InDetTrigBroadPixelClusterOnTrackTool
  if (InDetTrigFlags.doPrintConfigurables()):
    print (InDetTrigBroadPixelClusterOnTrackTool)

  # load RIO_OnTrackCreator for Inner Detector
  #

  from TrkRIO_OnTrackCreator.TrkRIO_OnTrackCreatorConf import Trk__RIO_OnTrackCreator
  InDetTrigRotCreator = Trk__RIO_OnTrackCreator(name = 'InDetTrigRotCreator',
                                                ToolPixelCluster= InDetTrigPixelClusterOnTrackTool,
                                                ToolSCT_Cluster = SCT_ClusterOnTrackTool,
                                                Mode = 'indet')
  ToolSvc += InDetTrigRotCreator

  if InDetTrigFlags.useBroadClusterErrors():
    InDetTrigRotCreator.ToolPixelCluster = InDetTrigBroadPixelClusterOnTrackTool
    InDetTrigRotCreator.ToolSCT_Cluster  = InDetTrigBroadSCT_ClusterOnTrackTool

  if (InDetTrigFlags.doPrintConfigurables()):
    print (InDetTrigRotCreator)

  #--
  from TrkRIO_OnTrackCreator.TrkRIO_OnTrackCreatorConf import Trk__RIO_OnTrackCreator
  InDetTrigBroadInDetRotCreator = \
      Trk__RIO_OnTrackCreator(name            = 'InDetTrigBroadInDetRotCreator',
                              ToolPixelCluster= InDetTrigBroadPixelClusterOnTrackTool,
                              ToolSCT_Cluster = InDetTrigBroadSCT_ClusterOnTrackTool,
                              Mode            = 'indet')
  ToolSvc += InDetTrigBroadInDetRotCreator
  if (InDetTrigFlags.doPrintConfigurables()):
    print (InDetTrigBroadInDetRotCreator)

  # load error scaling
  #TODO - instanceName?
  from InDetRecExample.TrackingCommon import createAndAddCondAlg, getRIO_OnTrackErrorScalingCondAlg
  createAndAddCondAlg(getRIO_OnTrackErrorScalingCondAlg,'RIO_OnTrackErrorScalingCondAlg')

  #
  # smart ROT creator in case we do the TRT LR in the refit
  #
  if InDetTrigFlags.redoTRT_LR():

    from InDetRecExample.TrackingCommon import getInDetTRT_DriftCircleOnTrackTool
    from TRT_DriftCircleOnTrackTool.TRT_DriftCircleOnTrackToolConf import \
        InDet__TRT_DriftCircleOnTrackUniversalTool
    InDetTrigTRT_RefitRotCreator = \
        InDet__TRT_DriftCircleOnTrackUniversalTool(name  = 'InDetTrigTRT_RefitRotCreator',
                                                   RIOonTrackToolDrift = getInDetTRT_DriftCircleOnTrackTool(), # special settings for trigger needed ?
                                                   ScaleHitUncertainty = 2.5) # fix from Thijs
#    if InDetTrigFlags.doCommissioning():    #introduced for cosmics do not use for collisions
#      InDetTrigTRT_RefitRotCreator.ScaleHitUncertainty = 5.
    ToolSvc += InDetTrigTRT_RefitRotCreator
      
    if (InDetTrigFlags.doPrintConfigurables()):
      print (     InDetTrigTRT_RefitRotCreator)
  
    from TrkRIO_OnTrackCreator.TrkRIO_OnTrackCreatorConf import Trk__RIO_OnTrackCreator
    InDetTrigRefitRotCreator = Trk__RIO_OnTrackCreator(name              = 'InDetTrigRefitRotCreator',
                                                       ToolPixelCluster= InDetTrigPixelClusterOnTrackTool,
                                                       ToolSCT_Cluster     = SCT_ClusterOnTrackTool,
                                                       ToolTRT_DriftCircle = InDetTrigTRT_RefitRotCreator,
                                                       Mode                = 'indet')
    if InDetTrigFlags.useBroadClusterErrors():
      InDetTrigRefitRotCreator.ToolPixelCluster = InDetTrigBroadPixelClusterOnTrackTool
      InDetTrigRefitRotCreator.ToolSCT_Cluster  = InDetTrigBroadSCT_ClusterOnTrackTool

    ToolSvc += InDetTrigRefitRotCreator
    if (InDetTrigFlags.doPrintConfigurables()):
      print (     InDetTrigRefitRotCreator)
         
  else:
    InDetTrigRefitRotCreator = InDetTrigRotCreator

#
# ----------- control loading of the kalman updator
#
if InDetTrigFlags.loadUpdator():
   
  if InDetTrigFlags.kalmanUpdator() == "fast" :
    # fast Kalman updator tool
    from TrkMeasurementUpdator_xk.TrkMeasurementUpdator_xkConf import Trk__KalmanUpdator_xk
    InDetTrigUpdator = Trk__KalmanUpdator_xk(name = 'InDetTrigUpdator')
  elif InDetTrigFlags.kalmanUpdator() == "weight" :
    from TrkMeasurementUpdator.TrkMeasurementUpdatorConf import Trk__KalmanWeightUpdator as ConfiguredWeightUpdator
    InDetTrigUpdator = ConfiguredWeightUpdator(name='InDetTrigUpdator')
  else :
    from TrkMeasurementUpdator.TrkMeasurementUpdatorConf import Trk__KalmanUpdator as ConfiguredKalmanUpdator
    InDetTrigUpdator = ConfiguredKalmanUpdator('InDetTrigUpdator')

  ToolSvc += InDetTrigUpdator
  if (InDetTrigFlags.doPrintConfigurables()):
    print (     InDetTrigUpdator)

#
# ----------- control loading extrapolation
#
if InDetTrigFlags.loadExtrapolator():

  
  #
  # get propagator
  #
  from TrkExSTEP_Propagator.TrkExSTEP_PropagatorConf import Trk__STEP_Propagator
  InDetTrigStepPropagator = Trk__STEP_Propagator(name = 'InDetTrigStepPropagator')
  ToolSvc += InDetTrigStepPropagator
  
  from TrkExRungeKuttaPropagator.TrkExRungeKuttaPropagatorConf import Trk__RungeKuttaPropagator
  InDetTrigRKPropagator = Trk__RungeKuttaPropagator(name = 'InDetTrigRKPropagator')
  ToolSvc += InDetTrigRKPropagator
  
  if InDetTrigFlags.propagatorType() == "STEP":
    InDetTrigPropagator = InDetTrigStepPropagator
  else:
    InDetTrigPropagator = InDetTrigRKPropagator
    
  if (InDetTrigFlags.doPrintConfigurables()):
    print (     InDetTrigPropagator      )

  #
  # Setup the Navigator (default, could be removed)
  #
  from InDetRecExample import TrackingCommon
  InDetTrigNavigator = TrackingCommon.getInDetNavigator('InDetTrigNavigator')

  ToolSvc += InDetTrigNavigator
  if (InDetTrigFlags.doPrintConfigurables()):
    print (     InDetTrigNavigator)

  #
  # Setup the MaterialEffectsUpdator
  #
  from TrkExTools.TrkExToolsConf import Trk__MaterialEffectsUpdator
  InDetTrigMaterialUpdator = Trk__MaterialEffectsUpdator(name = "InDetTrigMaterialEffectsUpdator")

  ToolSvc += InDetTrigMaterialUpdator
  if (InDetTrigFlags.doPrintConfigurables()):
    print (     InDetTrigMaterialUpdator)

  #
  # Set up extrapolator
  #

  InDetTrigSubPropagators = []
  InDetTrigSubUpdators = []
  
  # -------------------- set it depending on the geometry --------------
  # default for ID is (Rk,Mat)
  InDetTrigSubPropagators += [ InDetTrigPropagator.name() ]
  InDetTrigSubUpdators    += [ InDetTrigMaterialUpdator.name() ]
  
  # default for Calo is (Rk,MatLandau)
  InDetTrigSubPropagators += [ InDetTrigPropagator.name() ]
  InDetTrigSubUpdators    += [ InDetTrigMaterialUpdator.name() ]
  
  # default for MS is (STEP,Mat)
  InDetTrigSubPropagators += [ InDetTrigStepPropagator.name() ]
  InDetTrigSubUpdators    += [ InDetTrigMaterialUpdator.name() ]

  from TrkExTools.TrkExToolsConf import Trk__Extrapolator
  InDetTrigExtrapolator = Trk__Extrapolator(name        = 'InDetTrigExtrapolator',
                                            Propagators             = [ InDetTrigRKPropagator, InDetTrigStepPropagator],
                                            MaterialEffectsUpdators = [ InDetTrigMaterialUpdator ],
                                            Navigator               = InDetTrigNavigator,
                                            SubPropagators          = InDetTrigSubPropagators,
                                            SubMEUpdators           = InDetTrigSubUpdators,
                                            #DoCaloDynamic          = False
                                            )

  ToolSvc += InDetTrigExtrapolator
  if (InDetTrigFlags.doPrintConfigurables()):
    print (     InDetTrigExtrapolator                                          )

#
# ----------- control loading of fitters
#
if InDetTrigFlags.loadFitter():
    
  if InDetTrigFlags.trackFitterType() == 'DistributedKalmanFilter' :
   
    from TrkDistributedKalmanFilter.TrkDistributedKalmanFilterConf import Trk__DistributedKalmanFilter
    InDetTrigTrackFitter = Trk__DistributedKalmanFilter(name             = 'InDetTrigTrackFitter',
                                                        ExtrapolatorTool = InDetTrigExtrapolator,
                                                        ROTcreator       = InDetTrigRotCreator
                                                        #sortingReferencePoint = ???
                                                        )
   
  elif InDetTrigFlags.trackFitterType() == 'GlobalChi2Fitter' :
    from InDetRecExample import TrackingCommon
    cond_alg = TrackingCommon.createAndAddCondAlg(TrackingCommon.getTrackingGeometryCondAlg, "AtlasTrackingGeometryCondAlg", name="AtlasTrackingGeometryCondAlg")

    from TrkGlobalChi2Fitter.TrkGlobalChi2FitterConf import Trk__GlobalChi2Fitter
    InDetTrigTrackFitter = Trk__GlobalChi2Fitter(name                  = 'InDetTrigTrackFitter',
                                                 ExtrapolationTool     = InDetTrigExtrapolator,
                                                 NavigatorTool         = InDetTrigNavigator,
                                                 PropagatorTool        = InDetTrigPropagator,		
                                                 RotCreatorTool        = InDetTrigRefitRotCreator,
                                                 BroadRotCreatorTool   = InDetTrigBroadInDetRotCreator,
                                                 MeasurementUpdateTool = InDetTrigUpdator,
                                                 MaterialUpdateTool    = InDetTrigMaterialUpdator,
                                                 StraightLine          = not InDetTrigFlags.solenoidOn(),
                                                 OutlierCut            = 4,
                                                 SignedDriftRadius     = True,
                                                 RecalibrateSilicon    = True,
                                                 RecalibrateTRT        = True,
                                                 ReintegrateOutliers   = True,
                                                 TrackChi2PerNDFCut    = 9,
                                                 TRTExtensionCuts      = True, 
                                                 MaxIterations         = 40,
                                                 Acceleration          = True,
                                                 #Momentum=1000.,
                                                 Momentum=0.,
                                                 TrackingGeometryReadKey=cond_alg.TrackingGeometryWriteKey)
    if InDetTrigFlags.useBroadClusterErrors():
      InDetTrigTrackFitter.RecalibrateSilicon = False

    if InDetTrigFlags.doRefit():
      InDetTrigTrackFitter.BroadRotCreatorTool = None
      InDetTrigTrackFitter.RecalibrateSilicon = False
      InDetTrigTrackFitter.RecalibrateTRT     = False
      InDetTrigTrackFitter.ReintegrateOutliers= False


    if InDetTrigFlags.doRobustReco():
      InDetTrigTrackFitter.OutlierCut         = 10.0
      InDetTrigTrackFitter.TrackChi2PerNDFCut = 20
      InDetTrigTrackFitter.MaxOutliers        = 99
      #only switch off for cosmics InDetTrigTrackFitter.Acceleration       = False
      
    InDetTrigTrackFitterLowPt = Trk__GlobalChi2Fitter(name                  = 'InDetTrigTrackFitterLowPt',
                                                      ExtrapolationTool     = InDetTrigExtrapolator,
                                                      NavigatorTool         = InDetTrigNavigator,
                                                      PropagatorTool        = InDetTrigPropagator,		
                                                      RotCreatorTool        = InDetTrigRefitRotCreator,
                                                      BroadRotCreatorTool   = InDetTrigBroadInDetRotCreator,
                                                      MeasurementUpdateTool = InDetTrigUpdator,
                                                      StraightLine          = not InDetTrigFlags.solenoidOn(),
                                                      OutlierCut            = 5.0,
                                                      SignedDriftRadius     = True,
                                                      RecalibrateSilicon    = True,
                                                      RecalibrateTRT        = True,
                                                      ReintegrateOutliers   = True,
                                                      TrackChi2PerNDFCut    = 10,
                                                      TRTExtensionCuts      = True, 
                                                      MaxIterations         = 40,
                                                      Momentum=0.,
                                                      TrackingGeometryReadKey=cond_alg.TrackingGeometryWriteKey)
    ToolSvc += InDetTrigTrackFitterLowPt



    #for cosmics
    InDetTrigTrackFitterCosmics = \
        Trk__GlobalChi2Fitter(name = 'InDetTrigTrackFitterCosmics',
                              ExtrapolationTool     = InDetTrigExtrapolator,
                              NavigatorTool         = InDetTrigNavigator,
                              PropagatorTool        = InDetTrigPropagator,		
                              RotCreatorTool        = InDetTrigRefitRotCreator,
                              BroadRotCreatorTool   = InDetTrigBroadInDetRotCreator,
                              MeasurementUpdateTool = InDetTrigUpdator,
                              StraightLine          = not InDetTrigFlags.solenoidOn(),
                              OutlierCut            = 10.0,
                              SignedDriftRadius     = True,
                              RecalibrateTRT        = True,
                              ReintegrateOutliers   = True,
                              TrackChi2PerNDFCut    = 20,
                              TRTExtensionCuts      = False, 
                              MaxIterations         = 40,
                              MaxOutliers           = 99,
                              RecalculateDerivatives = True,
                              Momentum=1000,
                              TrackingGeometryReadKey=cond_alg.TrackingGeometryWriteKey)
    InDetTrigTrackFitterCosmics.Acceleration       = False
    ToolSvc += InDetTrigTrackFitterCosmics

    #TRT
    InDetTrigTrackFitterTRT = Trk__GlobalChi2Fitter(name                  = 'InDetTrigTrackFitterTRT',
                                                    ExtrapolationTool     = InDetTrigExtrapolator,
                                                    NavigatorTool         = InDetTrigNavigator,
                                                    PropagatorTool        = InDetTrigPropagator,
                                                    RotCreatorTool        = InDetTrigRefitRotCreator,
                                                    MeasurementUpdateTool = InDetTrigUpdator,
                                                    StraightLine          = not InDetTrigFlags.solenoidOn(),
                                                    ReintegrateOutliers   = False,
                                                    MaxIterations         = 30,
                                                    #MaxOutliers           = 99,
                                                    RecalculateDerivatives= True,
                                                    TrackChi2PerNDFCut    = 999999,
                                                    Momentum=0.,
                                                    TrackingGeometryReadKey=cond_alg.TrackingGeometryWriteKey)
    if InDetTrigFlags.doRobustReco():
      InDetTrigTrackFitterTRT.MaxOutliers=99
            
  elif InDetTrigFlags.trackFitterType() == 'GaussianSumFilter' :
    import egammaRec.EMCommonRefitter
    GSFTrackFitter = egammaRec.EMCommonRefitter.getGSFTrackFitter(name = 'InDetTrigTrackFitter')
    
   # --- end of fitter loading
  ToolSvc += InDetTrigTrackFitter
  if (InDetTrigFlags.doPrintConfigurables()):
    print (     InDetTrigTrackFitter)
    
  if InDetTrigFlags.trackFitterType() != 'GlobalChi2Fitter' :
    InDetTrigTrackFitterTRT=InDetTrigTrackFitter
    InDetTrigTrackFitterLowPt=InDetTrigTrackFitter
    ToolSvc += InDetTrigTrackFitterLowPt
    InDetTrigTrackFitterCosmics=InDetTrigTrackFitter
    
  ToolSvc += InDetTrigTrackFitterTRT
  if (InDetTrigFlags.doPrintConfigurables()):
    print (     InDetTrigTrackFitterTRT)



InDetTrigPixelConditionsSummaryTool = PixelConditionsSetup.summaryTool

if DetFlags.haveRIO.SCT_on():
  from SCT_ConditionsTools.SCT_ConditionsToolsConf import SCT_ConditionsSummaryTool
  InDetTrigSCTConditionsSummaryTool = SCT_ConditionsSummaryTool(SCT_ConditionsSetup.instanceName('InDetSCT_ConditionsSummaryTool'))

  # Fix the conditions tools - please tell me there's a better way
  from AthenaCommon.GlobalFlags import globalflags
  fixedTools = []
  for tool in InDetTrigSCTConditionsSummaryTool.ConditionsTools:
    if hasattr( tool, "SCT_FlaggedCondData" ):
      tool.SCT_FlaggedCondData = "SCT_FlaggedCondData_TRIG"
    if not globalflags.InputFormat.is_bytestream() and "ByteStream" in tool.getName():
      continue
    fixedTools.append( tool )
  InDetTrigSCTConditionsSummaryTool.ConditionsTools = fixedTools

else:
  InDetTrigSCTConditionsSummaryTool = None

#
# ------load association tool from Inner Detector to handle pixel ganged ambiguities
#
if InDetTrigFlags.loadAssoTool():
  from InDetAssociationTools.InDetAssociationToolsConf import InDet__InDetPRD_AssociationToolGangedPixels
  InDetTrigPrdAssociationTool = InDet__InDetPRD_AssociationToolGangedPixels(name = "InDetTrigPrdAssociationTool",
                                                                             PixelClusterAmbiguitiesMapName = "TrigPixelClusterAmbiguitiesMap")
   
  ToolSvc += InDetTrigPrdAssociationTool
  if (InDetTrigFlags.doPrintConfigurables()):
    print (     InDetTrigPrdAssociationTool)

#
# ----------- control loading of Summary Tool
#
if InDetTrigFlags.loadSummaryTool():

  # Load Pixel Layer tool
  # used as private tool only, don't add it to ToolSvc
  from InDetTestPixelLayer.InDetTestPixelLayerConf import InDet__InDetTestPixelLayerTool
  InDetTrigTestPixelLayerTool = InDet__InDetTestPixelLayerTool(name             = "InDetTrigTestPixelLayerTool",
                                                               PixelSummaryTool = InDetTrigPixelConditionsSummaryTool,
                                                               Extrapolator     = InDetTrigExtrapolator,
                                                               CheckActiveAreas = True,
                                                               CheckDeadRegions = True)
  if (InDetTrigFlags.doPrintConfigurables()):
    print ( InDetTrigTestPixelLayerTool)

  from InDetBoundaryCheckTool.InDetBoundaryCheckToolConf import InDet__InDetBoundaryCheckTool
  # used as private tool only, don't add it to ToolSvc
  InDetTrigBoundaryCheckTool = InDet__InDetBoundaryCheckTool(name="InDetTrigBoundaryCheckTool",
                                                             UsePixel=DetFlags.haveRIO.pixel_on(),
                                                             UseSCT=DetFlags.haveRIO.SCT_on(),
                                                             SctSummaryTool = InDetTrigSCTConditionsSummaryTool,
                                                             PixelLayerTool = InDetTrigTestPixelLayerTool
                                                             )

  #
  # Loading Configurable HoleSearchTool
  #
  from InDetTrackHoleSearch.InDetTrackHoleSearchConf import InDet__InDetTrackHoleSearchTool

  InDetTrigHoleSearchTool = InDet__InDetTrackHoleSearchTool(name = "InDetTrigHoleSearchTool",
                                                            Extrapolator = InDetTrigExtrapolator,
                                                            BoundaryCheckTool=InDetTrigBoundaryCheckTool
                                                            )
                                                            #Commissioning = InDetTrigFlags.doCommissioning()) #renamed
  InDetTrigHoleSearchTool.CountDeadModulesAfterLastHit = True  

  ToolSvc += InDetTrigHoleSearchTool
  if (InDetTrigFlags.doPrintConfigurables()):
    print (     InDetTrigHoleSearchTool)

  #Load inner Pixel layer tool
  from InDetRecExample import TrackingCommon
  InDetTrigTestPixelLayerToolInner = TrackingCommon.getInDetTrigTestPixelLayerToolInner()
  ToolSvc += InDetTrigTestPixelLayerToolInner
  if (InDetTrigFlags.doPrintConfigurables()):
    print ( InDetTrigTestPixelLayerToolInner)

  #prevent loading of the pixel dE/dx tool  
  InDetTrigPixelToTPIDTool = None
  if (InDetTrigFlags.doPrintConfigurables()):
    print (    'InDetTrigPixelToTPIDTool ', InDetTrigPixelToTPIDTool)

  #
  # Configrable version of loading the InDetTrackSummaryHelperTool
  #
  from AthenaCommon.GlobalFlags import globalflags
  
  from InDetTrackSummaryHelperTool.InDetTrackSummaryHelperToolConf import InDet__InDetTrackSummaryHelperTool
  from InDetTrigRecExample.InDetTrigConditionsAccess import TRT_ConditionsSetup  # noqa: F401
  from InDetTrigRecExample.InDetTrigCommonTools import InDetTrigTRTStrawStatusSummaryTool
  InDetTrigTrackSummaryHelperTool = InDet__InDetTrackSummaryHelperTool(name          = "InDetTrigSummaryHelper",
                                                                       HoleSearch    = InDetTrigHoleSearchTool,
                                                                       AssoTool      = InDetTrigPrdAssociationTool,
                                                                       DoSharedHits  = False,
                                                                       TRTStrawSummarySvc=InDetTrigTRTStrawStatusSummaryTool,
                                                                       usePixel      = DetFlags.haveRIO.pixel_on(),
                                                                       useSCT        = DetFlags.haveRIO.SCT_on(),
                                                                       useTRT        = DetFlags.haveRIO.TRT_on())

  ToolSvc += InDetTrigTrackSummaryHelperTool
  if (InDetTrigFlags.doPrintConfigurables()):
    print (     InDetTrigTrackSummaryHelperTool)

  InDetTrigTrackSummaryHelperToolSi = InDet__InDetTrackSummaryHelperTool(name          = "InDetTrigSummaryHelperSi",
                                                                         HoleSearch    = InDetTrigHoleSearchTool,
                                                                         AssoTool      = InDetTrigPrdAssociationTool,
                                                                         DoSharedHits  = False,
                                                                         TRTStrawSummarySvc=None,
                                                                         usePixel      = DetFlags.haveRIO.pixel_on(),
                                                                         useSCT        = DetFlags.haveRIO.SCT_on(),
                                                                         useTRT        = False)

  ToolSvc += InDetTrigTrackSummaryHelperToolSi
  if (InDetTrigFlags.doPrintConfigurables()):
    print (     InDetTrigTrackSummaryHelperToolSi)
    
  #
  # Configurable version of TRT_ElectronPidTools
  #
  from IOVDbSvc.CondDB import conddb

  if not (conddb.folderRequested("/TRT/Calib/PID_vector") or \
            conddb.folderRequested("/TRT/Onl/Calib/PID_vector")):
    conddb.addFolderSplitOnline("TRT","/TRT/Onl/Calib/PID_vector","/TRT/Calib/PID_vector",className='CondAttrListVec')
  if not (conddb.folderRequested("/TRT/Calib/ToT/ToTVectors") or \
            conddb.folderRequested("/TRT/Onl/Calib/ToT/ToTVectors")):
    conddb.addFolderSplitOnline("TRT","/TRT/Onl/Calib/ToT/ToTVectors","/TRT/Calib/ToT/ToTVectors",className='CondAttrListVec')
  if not (conddb.folderRequested("/TRT/Calib/ToT/ToTValue") or \
            conddb.folderRequested("/TRT/Onl/Calib/ToT/ToTValue")):
    conddb.addFolderSplitOnline("TRT","/TRT/Onl/Calib/ToT/ToTValue","/TRT/Calib/ToT/ToTValue",className='CondAttrListCollection')
  if InDetTrigFlags.doTRTPIDNN():
    if not (conddb.folderRequested( "/TRT/Calib/PID_NN") or \
           conddb.folderRequested( "/TRT/Onl/Calib/PID_NN")):
      conddb.addFolderSplitOnline( "TRT", "/TRT/Onl/Calib/PID_NN", "/TRT/Calib/PID_NN",className='CondAttrListCollection')
    # FIXME: need to force an override for the online DB until this folder has been added to the latest tag
    conddb.addOverride("/TRT/Onl/Calib/PID_NN", "TRTCalibPID_NN_v2")

  # Calibration DB Tool
  from TRT_ConditionsServices.TRT_ConditionsServicesConf import TRT_CalDbTool
  InDetTRTCalDbTool = TRT_CalDbTool(name = "TRT_CalDbTool")

 
  #
  # Configurable version of TrkTrackSummaryTool
  #
  from TrkTrackSummaryTool.TrkTrackSummaryToolConf import Trk__TrackSummaryTool
  InDetTrigTrackSummaryTool = Trk__TrackSummaryTool(name = "InDetTrigTrackSummaryTool",
                                                    InDetSummaryHelperTool = InDetTrigTrackSummaryHelperTool,
                                                    doSharedHits           = False,
                                                    doHolesInDet           = True,
                                                    #this may be temporary #61512 (and used within egamma later)
                                                    )
  ToolSvc += InDetTrigTrackSummaryTool
  if (InDetTrigFlags.doPrintConfigurables()):
     print (     InDetTrigTrackSummaryTool)

  if InDetTrigFlags.doSharedHits():
    #
    # Configrable version of loading the InDetTrackSummaryHelperTool
    #
    from InDetTrackSummaryHelperTool.InDetTrackSummaryHelperToolConf import InDet__InDetTrackSummaryHelperTool
    InDetTrigTrackSummaryHelperToolSharedHits = InDet__InDetTrackSummaryHelperTool(name         = "InDetTrigSummaryHelperSharedHits",
                                                                                   AssoTool     = InDetTrigPrdAssociationTool,
                                                                                   DoSharedHits = InDetTrigFlags.doSharedHits(),
                                                                                   HoleSearch   = InDetTrigHoleSearchTool,
                                                                                   TRTStrawSummarySvc = InDetTrigTRTStrawStatusSummaryTool)

    ToolSvc += InDetTrigTrackSummaryHelperToolSharedHits
    if (InDetTrigFlags.doPrintConfigurables()):
      print (     InDetTrigTrackSummaryHelperToolSharedHits)
    #
    # Configurable version of TrkTrackSummaryTool
    #
    from TrkTrackSummaryTool.TrkTrackSummaryToolConf import Trk__TrackSummaryTool
    InDetTrigTrackSummaryToolSharedHits = Trk__TrackSummaryTool(name = "InDetTrigTrackSummaryToolSharedHits",
                                                                InDetSummaryHelperTool = InDetTrigTrackSummaryHelperToolSharedHits,
                                                                doSharedHits           = InDetTrigFlags.doSharedHits(),
                                                                doHolesInDet           = True)

    ToolSvc += InDetTrigTrackSummaryToolSharedHits
    if (InDetTrigFlags.doPrintConfigurables()):
      print (     InDetTrigTrackSummaryToolSharedHits)


  else:   
    InDetTrigTrackSummaryToolSharedHits        = InDetTrigTrackSummaryTool   

#
# ----------- control loading of tools which are needed by new tracking and backtracking
#
if InDetTrigFlags.doNewTracking() or InDetTrigFlags.doBackTracking() or InDetTrigFlags.doTrtSegments():
  # Igor's propagator and updator for the pattern
  #
  from TrkExRungeKuttaPropagator.TrkExRungeKuttaPropagatorConf import Trk__RungeKuttaPropagator as Propagator
  InDetTrigPatternPropagator = Propagator(name = 'InDetTrigPatternPropagator')
  
  ToolSvc += InDetTrigPatternPropagator
  if (InDetTrigFlags.doPrintConfigurables()):
     print (     InDetTrigPatternPropagator)

  # fast Kalman updator tool
  #
  from TrkMeasurementUpdator_xk.TrkMeasurementUpdator_xkConf import Trk__KalmanUpdator_xk
  InDetTrigPatternUpdator = Trk__KalmanUpdator_xk(name = 'InDetTrigPatternUpdator')

  ToolSvc += InDetTrigPatternUpdator
  if (InDetTrigFlags.doPrintConfigurables()):
     print (     InDetTrigPatternUpdator)


#
# TRT segment minimum number of drift circles tool
#

from InDetTrackSelectorTool.InDetTrackSelectorToolConf import InDet__InDetTrtDriftCircleCutTool
InDetTrigTRTDriftCircleCut = InDet__InDetTrtDriftCircleCutTool(
  name             = 'InDetTrigTRTDriftCircleCut',
  MinOffsetDCs     = 5,
  UseNewParameterization = True,
  UseActiveFractionSvc   = True #DetFlags.haveRIO.TRT_on()  # Use Thomas's new parameterization by default
  )

ToolSvc += InDetTrigTRTDriftCircleCut
if (InDetTrigFlags.doPrintConfigurables()):
  print (  InDetTrigTRTDriftCircleCut)

InDetTrigTRTDriftCircleCutForPatt = InDet__InDetTrtDriftCircleCutTool(
  name             = 'InDetTrigTRTDriftCircleCutForPatt',
  MinOffsetDCs     = 5,
  UseNewParameterization = InDetTrigCutValues.useNewParameterizationTRT(),
  UseActiveFractionSvc   = True #DetFlags.haveRIO.TRT_on()  # Use Thomas's new parameterization by default
  )

ToolSvc += InDetTrigTRTDriftCircleCut
if (InDetTrigFlags.doPrintConfigurables()):
  print (  InDetTrigTRTDriftCircleCut)




    #default

# from here if InDetTrigFlags.doSiSPSeededTrackFinder():
# was unable to import with this condition from L2 TRTSegMaker

#if InDetTrigFlags.doSiSPSeededTrackFinder():
if InDetTrigFlags.doNewTracking():
  # SCT and Pixel detector elements road builder
  #
  from SiDetElementsRoadTool_xk.SiDetElementsRoadTool_xkConf import InDet__SiDetElementsRoadMaker_xk

  InDetTrigSiDetElementsRoadMaker = \
                                  InDet__SiDetElementsRoadMaker_xk(name = 'InDetTrigSiDetElementsRoadMaker',
                                                                   PropagatorTool = InDetTrigPatternPropagator,
                                                                   usePixel     = DetFlags.haveRIO.pixel_on(), 
                                                                   useSCT       = DetFlags.haveRIO.SCT_on(),
                                                                   RoadWidth    = InDetTrigCutValues.RoadWidth()
                                                                   )
  ToolSvc += InDetTrigSiDetElementsRoadMaker

  # Condition algorithm for InDet__SiDetElementsRoadMaker_xk
  if DetFlags.haveRIO.SCT_on():
    from AthenaCommon.AlgSequence import AthSequencer
    condSeq = AthSequencer("AthCondSeq")
    if not hasattr(condSeq, "InDet__SiDetElementsRoadCondAlg_xk"):
      from SiDetElementsRoadTool_xk.SiDetElementsRoadTool_xkConf import InDet__SiDetElementsRoadCondAlg_xk
      condSeq += InDet__SiDetElementsRoadCondAlg_xk(name = "InDet__SiDetElementsRoadCondAlg_xk")

  # Local combinatorial track finding using space point seed and detector element road
  #
  from SiCombinatorialTrackFinderTool_xk.SiCombinatorialTrackFinderTool_xkConf import InDet__SiCombinatorialTrackFinder_xk
  InDetTrigSiComTrackFinder = \
                            InDet__SiCombinatorialTrackFinder_xk(name = 'InDetTrigSiComTrackFinder',
                                                                 PropagatorTool	= InDetTrigPatternPropagator,
                                                                 UpdatorTool	= InDetTrigPatternUpdator,
                                                                 RIOonTrackTool   = InDetTrigRotCreator,
                                                                 usePixel         = DetFlags.haveRIO.pixel_on(),
                                                                 useSCT           = DetFlags.haveRIO.SCT_on(),   
                                                                 PixelClusterContainer = 'PixelTrigClusters',
                                                                 SCT_ClusterContainer = 'SCT_TrigClusters',
                                                                 PixelSummaryTool = InDetTrigPixelConditionsSummaryTool,
                                                                 SctSummaryTool = InDetTrigSCTConditionsSummaryTool,
                                                                 BoundaryCheckTool = InDetTrigBoundaryCheckTool
                                                                 )															
  if DetFlags.haveRIO.pixel_on():
    # Condition algorithm for SiCombinatorialTrackFinder_xk
    from AthenaCommon.AlgSequence import AthSequencer
    condSeq = AthSequencer("AthCondSeq")
    if not hasattr(condSeq, "InDetSiDetElementBoundaryLinksPixelCondAlg"):
      from SiCombinatorialTrackFinderTool_xk.SiCombinatorialTrackFinderTool_xkConf import InDet__SiDetElementBoundaryLinksCondAlg_xk
      condSeq += InDet__SiDetElementBoundaryLinksCondAlg_xk(name = "InDetSiDetElementBoundaryLinksPixelCondAlg",
                                                            ReadKey = "PixelDetectorElementCollection",
                                                            WriteKey = "PixelDetElementBoundaryLinks_xk")
  if DetFlags.haveRIO.SCT_on():
    # Condition algorithm for SiCombinatorialTrackFinder_xk
    from AthenaCommon.AlgSequence import AthSequencer
    condSeq = AthSequencer("AthCondSeq")
    if not hasattr(condSeq, "InDetSiDetElementBoundaryLinksSCTCondAlg"):
      from SiCombinatorialTrackFinderTool_xk.SiCombinatorialTrackFinderTool_xkConf import InDet__SiDetElementBoundaryLinksCondAlg_xk
      condSeq += InDet__SiDetElementBoundaryLinksCondAlg_xk(name = "InDetSiDetElementBoundaryLinksSCTCondAlg",
                                                            ReadKey = "SCT_DetectorElementCollection",
                                                            WriteKey = "SCT_DetElementBoundaryLinks_xk")
      #to here

import InDetRecExample.TrackingCommon as TrackingCommon
from InDetAmbiTrackSelectionTool.InDetAmbiTrackSelectionToolConf import InDet__InDetAmbiTrackSelectionTool
InDetTrigAmbiTrackSelectionTool = \
    InDet__InDetAmbiTrackSelectionTool(name               = 'InDetTrigAmbiTrackSelectionTool',
                                       DriftCircleCutTool = InDetTrigTRTDriftCircleCut,
                                       AssociationTool = TrackingCommon.getInDetTrigPRDtoTrackMapToolGangedPixels(),
                                       minHits         = InDetTrigCutValues.minClusters(),
                                       minNotShared    = InDetTrigCutValues.minSiNotShared(),
                                       maxShared       = InDetTrigCutValues.maxShared(),
                                       minTRTHits      = 0,  # used for Si only tracking !!!
                                       Cosmics         = False,  #there is a different instance
                                       UseParameterization = False,
                                       # sharedProbCut   = 0.10,
                                       # doPixelSplitting = InDetTrigFlags.doPixelClusterSplitting()
                                       )
 
 
ToolSvc += InDetTrigAmbiTrackSelectionTool
if (InDetTrigFlags.doPrintConfigurables()):
  print (InDetTrigAmbiTrackSelectionTool)

if InDetTrigFlags.doNewTracking():

  #
  # ------ load new track selector (common for all vertexing algorithms, except for the moment VKalVrt
  #
  from InDetTrigRecExample.ConfiguredVertexingTrigCuts import EFIDVertexingCuts
  from InDetTrackSelectionTool.InDetTrackSelectionToolConf import InDet__InDetTrackSelectionTool

  InDetTrigTrackSelectorTool = \
      InDet__InDetTrackSelectionTool(name = "InDetTrigDetailedTrackSelectorTool",
                                     CutLevel                   =  EFIDVertexingCuts.TrackCutLevel(),
                                     minPt                      =  EFIDVertexingCuts.minPT(),
                                     maxD0			=  EFIDVertexingCuts.IPd0Max(),
                                     maxZ0			=  EFIDVertexingCuts.z0Max(),
                                     maxZ0SinTheta              =  EFIDVertexingCuts.IPz0Max(),
                                     maxSigmaD0 = EFIDVertexingCuts.sigIPd0Max(),
                                     maxSigmaZ0SinTheta = EFIDVertexingCuts.sigIPz0Max(),
                                     # maxChiSqperNdf = EFIDVertexingCuts.fitChi2OnNdfMax(), # Seems not to be implemented?
                                     maxAbsEta = EFIDVertexingCuts.etaMax(),
                                     minNInnermostLayerHits = EFIDVertexingCuts.nHitInnermostLayer(),
                                     minNPixelHits = EFIDVertexingCuts.nHitPix(),
                                     maxNPixelHoles = EFIDVertexingCuts.nHolesPix(),
                                     minNSctHits = EFIDVertexingCuts.nHitSct(),
                                     minNTrtHits = EFIDVertexingCuts.nHitTrt(),
                                     minNSiHits = EFIDVertexingCuts.nHitSi(),
                                     TrackSummaryTool =  InDetTrigTrackSummaryTool,
                                     Extrapolator     = InDetTrigExtrapolator,
                                     #TrtDCCutTool     = InDetTrigTRTDriftCircleCut,
                                     )



  ToolSvc += InDetTrigTrackSelectorTool
  if (InDetTrigFlags.doPrintConfigurables()):
    print (     InDetTrigTrackSelectorTool)


# --- set Data/MC flag
isMC = False
if globalflags.DataSource == "geant4" :
    isMC = True

# Calibration DB Service
from TRT_ConditionsServices.TRT_ConditionsServicesConf import TRT_CalDbTool
InDetTRTCalDbTool = TRT_CalDbTool(name = "TRT_CalDbTool")


# TRT_DriftFunctionTool
from TRT_DriftFunctionTool.TRT_DriftFunctionToolConf import TRT_DriftFunctionTool

InDetTrigTRT_DriftFunctionTool = TRT_DriftFunctionTool(name = "InDetTrigTRT_DriftFunctionTool",
                                                       TRTCalDbTool        = InDetTRTCalDbTool,
                                                       AllowDataMCOverride = True,
                                                       ForceData = True,
                                                       IsMC = isMC)

# Second calibration DB Service in case pile-up and physics hits have different calibrations
if DetFlags.overlay.TRT_on() :

    InDetTrigTRTCalDbTool2 = TRT_CalDbTool(name = "TRT_CalDbSvc2")
    InDetTrigTRTCalDbTool2.RtFolderName = "/TRT/Calib/MC/RT"             
    InDetTrigTRTCalDbTool2.T0FolderName = "/TRT/Calib/MC/T0"             
    InDetTrigTRT_DriftFunctionTool.TRTCalDbTool2 = InDetTrigTRTCalDbTool2
    InDetTrigTRT_DriftFunctionTool.IsOverlay = True
    InDetTrigTRT_DriftFunctionTool.IsMC = False

# --- set HT corrections
InDetTrigTRT_DriftFunctionTool.HTCorrectionBarrelXe = 1.5205
InDetTrigTRT_DriftFunctionTool.HTCorrectionEndcapXe = 1.2712
InDetTrigTRT_DriftFunctionTool.HTCorrectionBarrelAr = 1.5205
InDetTrigTRT_DriftFunctionTool.HTCorrectionEndcapAr = 1.2712
         
# --- set ToT corrections
InDetTrigTRT_DriftFunctionTool.ToTCorrectionsBarrelXe = [0., 4.358121, 3.032195, 1.631892, 0.7408397, -0.004113, -0.613288, -0.73758, -0.623346, -0.561229,-0.29828, -0.21344, -0.322892, -0.386718, -0.534751, -0.874178, -1.231799, -1.503689, -1.896464, -2.385958]
InDetTrigTRT_DriftFunctionTool.ToTCorrectionsEndcapXe = [0., 5.514777, 3.342712, 2.056626, 1.08293693, 0.3907979, -0.082819, -0.457485, -0.599706, -0.427493, -0.328962, -0.403399, -0.663656, -1.029428, -1.46008, -1.919092, -2.151582, -2.285481, -2.036822, -2.15805]
InDetTrigTRT_DriftFunctionTool.ToTCorrectionsBarrelAr = [0., 4.358121, 3.032195, 1.631892, 0.7408397, -0.004113, -0.613288, -0.73758, -0.623346, -0.561229, -0.29828, -0.21344, -0.322892, -0.386718, -0.534751, -0.874178, -1.231799, -1.503689, -1.896464, -2.385958]
InDetTrigTRT_DriftFunctionTool.ToTCorrectionsEndcapAr = [0., 5.514777, 3.342712, 2.056626, 1.08293693, 0.3907979, -0.082819, -0.457485, -0.599706, -0.427493, -0.328962, -0.403399, -0.663656, -1.029428, -1.46008, -1.919092, -2.151582, -2.285481, -2.036822, -2.15805]


ToolSvc += InDetTrigTRT_DriftFunctionTool

from AthenaCommon.GlobalFlags import globalflags

# TRT_DriftCircleTool
import AthenaCommon.SystemOfUnits as Units

MinTrailingEdge = 11.0*Units.ns
MaxDriftTime = 60.0*Units.ns
LowGate         = 14.0625*Units.ns # 4.5*3.125 ns
HighGate        = 42.1875*Units.ns # LowGate + 9*3.125 ns
LowGateArgon         = LowGate
HighGateArgon        = HighGate

if globalflags.DataSource == 'data':
    MinTrailingEdge = 11.0*Units.ns
    MaxDriftTime    = 60.0*Units.ns
    LowGate         = 17.1875*Units.ns
    HighGate        = 45.3125*Units.ns
    LowGateArgon    = 18.75*Units.ns
    HighGateArgon   = 43.75*Units.ns



from TRT_DriftCircleTool.TRT_DriftCircleToolConf import InDet__TRT_DriftCircleTool
InDetTrigTRT_DriftCircleTool = InDet__TRT_DriftCircleTool( name = "InDetTrigTRT_DriftCircleTool",
                                                           TRTDriftFunctionTool = InDetTrigTRT_DriftFunctionTool,
                                                           ConditionsSummaryTool           = InDetTrigTRTStrawStatusSummaryTool,
                                                           UseConditionsStatus  = True,
                                                           UseConditionsHTStatus  = True,
                                                           SimpleOutOfTimePileupSupression = False,
                                                           RejectIfFirstBit                = False, # fixes 50 nsec issue 
                                                           MinTrailingEdge                 = MinTrailingEdge,
                                                           MaxDriftTime                    = MaxDriftTime,
                                                           ValidityGateSuppression         = InDetTrigFlags.InDet25nsec(),
                                                           LowGate = LowGate,
                                                           HighGate = HighGate,
                                                           SimpleOutOfTimePileupSupressionArgon = False,# no OOT rejection for argon
                                                           RejectIfFirstBitArgon                = False, # no OOT rejection for argon
                                                           MinTrailingEdgeArgon                 = MinTrailingEdge,
                                                           MaxDriftTimeArgon                    = MaxDriftTime,
                                                           ValidityGateSuppressionArgon         = InDetTrigFlags.InDet25nsec(),
                                                           LowGateArgon                         = LowGateArgon,
                                                           HighGateArgon                        = HighGateArgon,
                                                           useDriftTimeHTCorrection        = True,
                                                           useDriftTimeToTCorrection       = True, # reenable ToT
                                                           )


ToolSvc += InDetTrigTRT_DriftCircleTool
  


# TRT_RodDecoder
from TRT_RawDataByteStreamCnv.TRT_RawDataByteStreamCnvConf import TRT_RodDecoder

InDetTrigTRTRodDecoder = TRT_RodDecoder(name = "InDetTrigTRTRodDecoder",
                                        LoadCompressTableDB = (globalflags.DataSource() != 'geant4'))
ToolSvc += InDetTrigTRTRodDecoder


from TrkTrackSummaryTool.TrkTrackSummaryToolConf import Trk__TrackSummaryTool
InDetTrigFastTrackSummaryTool = Trk__TrackSummaryTool(name = "InDetTrigFastTrackSummaryTool",
                                                      InDetSummaryHelperTool = InDetTrigTrackSummaryHelperToolSi,
                                                      doHolesInDet           = False,
                                                      doSharedHits           = False
                                                      )
ToolSvc += InDetTrigFastTrackSummaryTool
if (InDetTrigFlags.doPrintConfigurables()):
    print      (InDetTrigFastTrackSummaryTool)


from InDetTrigRecExample.InDetTrigConfigRecLoadTools import InDetTrigHoleSearchTool
InDetTrigTrackSummaryToolWithHoleSearch = Trk__TrackSummaryTool(name = "InDetTrigTrackSummaryToolWithHoleSearch",
                                                                InDetSummaryHelperTool = InDetTrigTrackSummaryHelperToolSi,
                                                                doHolesInDet           = True,
                                                                doSharedHits           = False
                                                      )
ToolSvc += InDetTrigTrackSummaryToolWithHoleSearch
if (InDetTrigFlags.doPrintConfigurables()):
    print      (InDetTrigTrackSummaryToolWithHoleSearch)


# HACK to emulate run2 behaviour
from TrkAssociationTools.TrkAssociationToolsConf import Trk__PRDtoTrackMapExchangeTool
InDetTrigPRDtoTrackMapExchangeTool = Trk__PRDtoTrackMapExchangeTool("InDetTrigPRDtoTrackMapExchangeTool")
ToolSvc += InDetTrigPRDtoTrackMapExchangeTool
