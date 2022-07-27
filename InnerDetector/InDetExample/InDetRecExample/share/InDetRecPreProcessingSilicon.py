
# ------------------------------------------------------------
#
# ----------- Data-Preparation stage
#
# ------------------------------------------------------------
#
# ----------- PrepRawData creation from Raw Data Objects
#

if not 'redoPatternRecoAndTracking' in dir():
  redoPatternRecoAndTracking = False

if InDetFlags.doPRDFormation():
   #
   # --- Slim BCM RDOs by zero-suppressing
   #   
   if DetFlags.makeRIO.BCM_on():
      from BCM_ZeroSuppression.BCM_ZeroSuppressionConf import BCM_ZeroSuppression
      InDetBCM_ZeroSuppression = BCM_ZeroSuppression(name             = "InDetBCM_ZeroSuppression",
                                                     BcmContainerName = InDetKeys.BCM_RDOs())
      topSequence += InDetBCM_ZeroSuppression
      if (InDetFlags.doPrintConfigurables()):
         printfunc (InDetBCM_ZeroSuppression)
   
   if DetFlags.makeRIO.pixel_on() or DetFlags.makeRIO.SCT_on():
      #
      # --- SiLorentzAngleTool
      #
      if not hasattr(ToolSvc, "PixelLorentzAngleTool"):
        from SiLorentzAngleTool.PixelLorentzAngleToolSetup import PixelLorentzAngleToolSetup
        pixelLorentzAngleToolSetup = PixelLorentzAngleToolSetup()

      from SiLorentzAngleTool.SCTLorentzAngleToolSetup import SCTLorentzAngleToolSetup
      sctLorentzAngleToolSetup = SCTLorentzAngleToolSetup()

      #
      # --- ClusterMakerTool (public), needed by Pixel and SCT Clusterization
      #
      from SiClusterizationTool.SiClusterizationToolConf import InDet__ClusterMakerTool
      InDetClusterMakerTool = InDet__ClusterMakerTool(name                 = "InDetClusterMakerTool",
                                                      PixelLorentzAngleTool = ToolSvc.PixelLorentzAngleTool,
                                                      SCTLorentzAngleTool = sctLorentzAngleToolSetup.SCTLorentzAngleTool)

      ToolSvc += InDetClusterMakerTool
      if (InDetFlags.doPrintConfigurables()):
        printfunc (InDetClusterMakerTool)
      
   #
   # -- Pixel Clusterization
   #
   if (DetFlags.makeRIO.pixel_on() and InDetFlags.doPixelPRDFormation()) or redoPatternRecoAndTracking:     
      #
      # --- MergedPixelTool (public)
      #
      minSplitProbability  = 0
      minSplitSize         = 3
      clusterSplitProbTool = None
      clusterSplitterTool  = None
      #
      # --- do we use new splittig or not ?
      #
      if InDetFlags.doPixelClusterSplitting():
         #
         # --- Neutral Network version ?
         #
         if InDetFlags.pixelClusterSplittingType() == 'NeuralNet':
            useBeamConstraint = InDetFlags.useBeamConstraint()
            
            # --- new NN prob tool
            MultiplicityContent = [1 , 1 , 1]
            from SiClusterizationTool.SiClusterizationToolConf import InDet__NnPixelClusterSplitProbTool as PixelClusterSplitProbTool
            NnPixelClusterSplitProbTool=PixelClusterSplitProbTool(name                     = "NnPixelClusterSplitProbTool",
                                                                           PriorMultiplicityContent = MultiplicityContent,
                                                                           NnClusterizationFactory  = NnClusterizationFactory,
                                                                           useBeamSpotInfo          = useBeamConstraint)
            ToolSvc += NnPixelClusterSplitProbTool
            if (InDetFlags.doPrintConfigurables()):
              printfunc (NnPixelClusterSplitProbTool)

            # --- remember this prob tool  
            clusterSplitProbTool = NnPixelClusterSplitProbTool
            
            # --- new NN splitter
            from SiClusterizationTool.SiClusterizationToolConf import InDet__NnPixelClusterSplitter as PixelClusterSplitter
            NnPixelClusterSplitter=PixelClusterSplitter(name                                = "NnPixelClusterSplitter",
                                                                 NnClusterizationFactory             = NnClusterizationFactory,
                                                                 ThresholdSplittingIntoTwoClusters   = 0.5, # temp.
                                                                 ThresholdSplittingIntoThreeClusters = 0.25, # temp.
                                                                 SplitOnlyOnBLayer                   = False,
                                                                 useBeamSpotInfo                     = useBeamConstraint)

            
            ToolSvc += NnPixelClusterSplitter
            if (InDetFlags.doPrintConfigurables()):
              printfunc (NnPixelClusterSplitter)

            # remember splitter tool  
            clusterSplitterTool = NnPixelClusterSplitter
            

         #
         # --- Neutral Network version ?
         #
         elif InDetFlags.pixelClusterSplittingType() == 'AnalogClus':      
            # new splitter tool
            from SiClusterizationTool.SiClusterizationToolConf import InDet__TotPixelClusterSplitter
            TotPixelClusterSplitter=InDet__TotPixelClusterSplitter (name="TotPixelClusterSplitter")

            ToolSvc += TotPixelClusterSplitter
            if (InDetFlags.doPrintConfigurables()):
              printfunc (TotPixelClusterSplitter)

            # remember splitter tool    
            clusterSplitterTool = TotPixelClusterSplitter
                        
      #
      # --- now load the framework for the clustering
      #
      from SiClusterizationTool.SiClusterizationToolConf import InDet__MergedPixelsTool
      InDetMergedPixelsTool = InDet__MergedPixelsTool(name                    = "InDetMergedPixelsTool", 
                                                      globalPosAlg            = InDetClusterMakerTool,
                                                      PixelDetElStatus          = "PixelDetectorElementStatus")
      # Enable duplcated RDO check for data15 because duplication mechanism was used.
      from RecExConfig.RecFlags import rec
      if len(rec.projectName())>=6 and rec.projectName()[:6]=="data15":
        InDetMergedPixelsTool.CheckDuplicatedRDO = True
      
      ToolSvc += InDetMergedPixelsTool
      if (InDetFlags.doPrintConfigurables()):
        printfunc (InDetMergedPixelsTool)
                    
      #
      # --- PixelGangedAmbiguitiesFinder tool (public)
      #
      from SiClusterizationTool.SiClusterizationToolConf import InDet__PixelGangedAmbiguitiesFinder
      InDetPixelGangedAmbiguitiesFinder = InDet__PixelGangedAmbiguitiesFinder(name = "InDetPixelGangedAmbiguitiesFinder")
      ToolSvc += InDetPixelGangedAmbiguitiesFinder
      if (InDetFlags.doPrintConfigurables()):
        printfunc (InDetPixelGangedAmbiguitiesFinder)
            
      #
      # --- PixelClusterization algorithm
      #
      from AthenaCommon.AppMgr import ServiceMgr as svcMgr
      if not hasattr(svcMgr, "PixelReadoutManager"):
        from PixelReadoutGeometry.PixelReadoutGeometryConf import InDetDD__PixelReadoutManager
        svcMgr += InDetDD__PixelReadoutManager("PixelReadoutManager")

      from InDetPrepRawDataFormation.InDetPrepRawDataFormationConf import InDet__PixelClusterization
      InDetPixelClusterization = InDet__PixelClusterization(name                    = "InDetPixelClusterization",
                                                            clusteringTool          = InDetMergedPixelsTool,
                                                            gangedAmbiguitiesFinder = InDetPixelGangedAmbiguitiesFinder,
                                                            DataObjectName          = InDetKeys.PixelRDOs(),
                                                            ClustersName            = InDetKeys.PixelClusters())

      from RegionSelector.RegSelToolConfig import makeRegSelTool_Pixel
      InDetPixelClusterization.RegSelTool = makeRegSelTool_Pixel()

      topSequence += InDetPixelClusterization
      if (InDetFlags.doPrintConfigurables()):
         printfunc (InDetPixelClusterization)

      if InDetFlags.doSplitReco() :
        InDetPixelClusterizationPU = InDet__PixelClusterization(name                    = "InDetPixelClusterizationPU",
                                                                clusteringTool          = InDetMergedPixelsTool,
                                                                gangedAmbiguitiesFinder = InDetPixelGangedAmbiguitiesFinder,
                                                                DataObjectName          = InDetKeys.PixelPURDOs(),
                                                                ClustersName            = InDetKeys.PixelPUClusters(),
                                                                AmbiguitiesMap = "PixelClusterAmbiguitiesMapPU")

        from RegionSelector.RegSelToolConfig import makeRegSelTool_Pixel
        InDetPixelClusterizationPU.RegSelTool = makeRegSelTool_Pixel()

        topSequence += InDetPixelClusterizationPU
        if (InDetFlags.doPrintConfigurables()):
          printfunc (InDetPixelClusterizationPU)

   #
   # --- SCT Clusterization
   #
   if DetFlags.makeRIO.SCT_on() and InDetFlags.doSCT_PRDFormation():
   
      #
      # --- SCT_ClusteringTool (public)
      #
      from SiClusterizationTool.SiClusterizationToolConf import InDet__SCT_ClusteringTool
      InDetSCT_ClusteringTool = InDet__SCT_ClusteringTool(name              = "InDetSCT_ClusteringTool",
                                                          globalPosAlg      = InDetClusterMakerTool,
                                                          conditionsTool = InDetSCT_ConditionsSummaryToolWithoutFlagged,
                                                          SCTDetElStatus    = "SCTDetectorElementStatusWithoutFlagged")
      if InDetFlags.selectSCTIntimeHits():
         if InDetFlags.InDet25nsec(): 
            InDetSCT_ClusteringTool.timeBins = "01X" 
         else: 
            InDetSCT_ClusteringTool.timeBins = "X1X" 

      if (InDetFlags.doPrintConfigurables()):
        printfunc (InDetSCT_ClusteringTool)
            

      #
      # --- SCT_Clusterization algorithm
      #
      from InDetPrepRawDataFormation.InDetPrepRawDataFormationConf import InDet__SCT_Clusterization
      InDetSCT_Clusterization = InDet__SCT_Clusterization(name                    = "InDetSCT_Clusterization",
                                                          clusteringTool          = InDetSCT_ClusteringTool,
                                                          # ChannelStatus         = InDetSCT_ChannelStatusAlg,
                                                          DataObjectName          = InDetKeys.SCT_RDOs(),
                                                          ClustersName            = InDetKeys.SCT_Clusters(),
                                                          conditionsTool          = InDetSCT_ConditionsSummaryToolWithoutFlagged,
                                                          SCTDetElStatus            = "SCTDetectorElementStatusWithoutFlagged"
                                                          )
      if InDetFlags.cutSCTOccupancy():
        InDetSCT_Clusterization.maxFiredStrips = 384 #77
      else:
        InDetSCT_Clusterization.maxFiredStrips = 0
      topSequence += InDetSCT_Clusterization
      if (InDetFlags.doPrintConfigurables()):
        printfunc (InDetSCT_Clusterization)

      if InDetFlags.doSplitReco() :
        InDetSCT_ClusterizationPU = InDet__SCT_Clusterization(name                    = "InDetSCT_ClusterizationPU",
                                                              clusteringTool          = InDetSCT_ClusteringTool,
                                                              # ChannelStatus         = InDetSCT_ChannelStatusAlg,
                                                              DataObjectName          = InDetKeys.SCT_PU_RDOs(),
                                                              ClustersName            = InDetKeys.SCT_PU_Clusters(),
                                                              conditionsTool          = InDetSCT_ConditionsSummaryToolWithoutFlagged)
        if InDetFlags.cutSCTOccupancy():
          InDetSCT_ClusterizationPU.maxFiredStrips = 384 #77
        else:
          InDetSCT_ClusterizationPU.maxFiredStrips = 0
        topSequence += InDetSCT_ClusterizationPU
        if (InDetFlags.doPrintConfigurables()):
          printfunc (InDetSCT_ClusterizationPU)

#
# ----------- form SpacePoints from clusters in SCT and Pixels
#
if InDetFlags.doSpacePointFormation():
   #
   # --- SiSpacePointMakerTool (public)
   #
   from SiSpacePointTool.SiSpacePointToolConf import InDet__SiSpacePointMakerTool
   InDetSiSpacePointMakerTool = InDet__SiSpacePointMakerTool(name = "InDetSiSpacePointMakerTool")

   if InDetFlags.doCosmics() or InDetFlags.doBeamHalo():
      InDetSiSpacePointMakerTool.StripLengthTolerance       = 0.05

   ToolSvc += InDetSiSpacePointMakerTool
   if (InDetFlags.doPrintConfigurables()):
     printfunc (InDetSiSpacePointMakerTool)
   
   #
   # SiTrackerSpacePointFinder algorithm
   #
   from SiSpacePointFormation.SiSpacePointFormationConf import InDet__SiTrackerSpacePointFinder
   InDetSiTrackerSpacePointFinder = InDet__SiTrackerSpacePointFinder(name                   = "InDetSiTrackerSpacePointFinder",
                                                                     SiSpacePointMakerTool  = InDetSiSpacePointMakerTool,
                                                                     PixelsClustersName     = InDetKeys.PixelClusters(),
                                                                     SCT_ClustersName       = InDetKeys.SCT_Clusters(),
                                                                     SpacePointsPixelName   = InDetKeys.PixelSpacePoints(),
                                                                     SpacePointsSCTName     = InDetKeys.SCT_SpacePoints(),
                                                                     SpacePointsOverlapName = InDetKeys.OverlapSpacePoints(),
                                                                     ProcessPixels          = DetFlags.haveRIO.pixel_on(),
                                                                     ProcessSCTs            = DetFlags.haveRIO.SCT_on(),
                                                                     ProcessOverlaps        = DetFlags.haveRIO.SCT_on())

   # Condition algorithm for SiTrackerSpacePointFinder
   from AthenaCommon.AlgSequence import AthSequencer
   condSeq = AthSequencer("AthCondSeq")
   if not hasattr(condSeq, "InDetSiElementPropertiesTableCondAlg"):
      from SiSpacePointFormation.SiSpacePointFormationConf import InDet__SiElementPropertiesTableCondAlg
      condSeq += InDet__SiElementPropertiesTableCondAlg(name = "InDetSiElementPropertiesTableCondAlg")

#   if InDetFlags.doDBM():
#     InDetSiTrackerSpacePointFinderDBM = InDet__SiTrackerSpacePointFinder(name                   = "InDetSiTrackerSpacePointFinderDBM",
#                                                                          SiSpacePointMakerTool  = InDetSiSpacePointMakerTool,
#                                                                          PixelsClustersName     = InDetKeys.PixelClusters(),
#                                                                          SCT_ClustersName       = InDetKeys.SCT_Clusters(),
#                                                                          SpacePointsPixelName   = InDetKeys.PixelSpacePoints(),
#                                                                          SpacePointsSCTName     = InDetKeys.SCT_SpacePoints(),
#                                                                          SpacePointsOverlapName = InDetKeys.OverlapSpacePoints(),
#                                                                          ProcessPixels          = DetFlags.haveRIO.pixel_on(),
#                                                                          ProcessSCTs            = DetFlags.haveRIO.SCT_on(),
#                                                                          ProcessOverlaps        = DetFlags.haveRIO.SCT_on(),
#                                                                          OverlapLimitEtaMax     = 5.0,
#                                                                          OverlapLimitEtaMin     = 0
#                                                                          )
#     topSequence += InDetSiTrackerSpacePointFinderDBM

   if InDetFlags.doDBMstandalone():
      InDetSiTrackerSpacePointFinder.OverlapLimitEtaMax = 5.0
      InDetSiTrackerSpacePointFinder.OverlapLimitEtaMin = 0

   if InDetFlags.doCosmics():
      InDetSiTrackerSpacePointFinder.ProcessOverlaps      = False
      InDetSiTrackerSpacePointFinder.OverrideBeamSpot     = True
      InDetSiTrackerSpacePointFinder.VertexZ              = 0
      InDetSiTrackerSpacePointFinder.VertexX              = 0
      InDetSiTrackerSpacePointFinder.VertexY              = 99999999   
      InDetSiTrackerSpacePointFinder.OverlapLimitOpposite = 5

   topSequence += InDetSiTrackerSpacePointFinder
   if (InDetFlags.doPrintConfigurables()):
     printfunc (InDetSiTrackerSpacePointFinder)
     if (InDetFlags.doDBM()):
       printfunc (InDetSiTrackerSpacePointFinderDBM)

# this truth must only be done if you do PRD and SpacePointformation
# If you only do the latter (== running on ESD) then the needed input (simdata)
# is not in ESD but the resulting truth (clustertruth) is already there ...
if InDetFlags.doPRDFormation() and InDetFlags.doSpacePointFormation():
   if InDetFlags.doTruth():
      from InDetTruthAlgs.InDetTruthAlgsConf import InDet__PRD_MultiTruthMaker
      InDetPRD_MultiTruthMakerSi = InDet__PRD_MultiTruthMaker (name                        = 'InDetPRD_MultiTruthMakerSi',
                                                             PixelClusterContainerName   = InDetKeys.PixelClusters(),
                                                             SCTClusterContainerName     = InDetKeys.SCT_Clusters(),
                                                             TRTDriftCircleContainerName = "",
                                                             SimDataMapNamePixel         = InDetKeys.PixelSDOs(),
                                                             SimDataMapNameSCT           = InDetKeys.SCT_SDOs(),
                                                             SimDataMapNameTRT           = "",
                                                             TruthNamePixel              = InDetKeys.PixelClustersTruth(),
                                                             TruthNameSCT                = InDetKeys.SCT_ClustersTruth(),
                                                             TruthNameTRT                = "")
      # a bit complicated, but this is how the truth maker gets to know which detector is on
      if (not DetFlags.haveRIO.pixel_on() or not InDetFlags.doPixelPRDFormation()):
         InDetPRD_MultiTruthMakerSi.PixelClusterContainerName = ""
         InDetPRD_MultiTruthMakerSi.SimDataMapNamePixel       = ""
         InDetPRD_MultiTruthMakerSi.TruthNamePixel            = ""
      if (not DetFlags.haveRIO.SCT_on() or not InDetFlags.doSCT_PRDFormation()):
         InDetPRD_MultiTruthMakerSi.SCTClusterContainerName   = ""
         InDetPRD_MultiTruthMakerSi.SimDataMapNameSCT         = ""
         InDetPRD_MultiTruthMakerSi.TruthNameSCT              = ""

      topSequence += InDetPRD_MultiTruthMakerSi
      if (InDetFlags.doPrintConfigurables()):
        printfunc (InDetPRD_MultiTruthMakerSi)

      if InDetFlags.doSplitReco() :
        InDetPRD_MultiTruthMakerSiPU = InDet__PRD_MultiTruthMaker (name                        = 'InDetPRD_MultiTruthMakerSiPU',
                                                               PixelClusterContainerName   = InDetKeys.PixelPUClusters(),
                                                               SCTClusterContainerName     = InDetKeys.SCT_PU_Clusters(),
                                                               TRTDriftCircleContainerName = "",
                                                               SimDataMapNamePixel         = InDetKeys.PixelPUSDOs(),
                                                               SimDataMapNameSCT           = InDetKeys.SCT_PU_SDOs(),
                                                               SimDataMapNameTRT           = "",
                                                               TruthNamePixel              = InDetKeys.PixelPUClustersTruth(),
                                                               TruthNameSCT                = InDetKeys.SCT_PU_ClustersTruth(),
                                                               TruthNameTRT                = "")
        # a bit complicated, but this is how the truth maker gets to know which detector is on
        if (not DetFlags.haveRIO.pixel_on() or not InDetFlags.doPixelPRDFormation()):
           InDetPRD_MultiTruthMakerSiPU.PixelClusterContainerName = ""
           InDetPRD_MultiTruthMakerSiPU.SimDataMapNamePixel       = ""
           InDetPRD_MultiTruthMakerSiPU.TruthNamePixel            = ""
        if (not DetFlags.haveRIO.SCT_on() or not InDetFlags.doSCT_PRDFormation()):
           InDetPRD_MultiTruthMakerSiPU.SCTClusterContainerName   = ""
           InDetPRD_MultiTruthMakerSiPU.SimDataMapNameSCT         = ""
           InDetPRD_MultiTruthMakerSiPU.TruthNameSCT              = ""

        printfunc (InDetPRD_MultiTruthMakerSiPU)
        InDetPRD_MultiTruthMakerSiPU.OutputLevel=VERBOSE
        topSequence += InDetPRD_MultiTruthMakerSiPU
        if (InDetFlags.doPrintConfigurables()):
          printfunc (InDetPRD_MultiTruthMakerSiPU)
