#
#  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
#

'''
@file InDetMonitoringPixel.py
@brief Top configuration of Pixel Monitoring in Run 3 style but in Run 2 environment
'''
forceOnline = False # for testing of online monitoring and 100LB histograms

doHitMonAlg       = True
doClusterMonAlg   = True
doErrorMonAlg     = True
doMVAMonAlg       = False

from PixelMonitoring.PixelMonitoringConf import PixelAthHitMonAlg
from PixelMonitoring.PixelAthHitMonAlgCfg import PixelAthHitMonAlgCfg

from PixelMonitoring.PixelMonitoringConf import PixelAthClusterMonAlg
from PixelMonitoring.PixelAthClusterMonAlgCfg import PixelAthClusterMonAlgCfg

from PixelMonitoring.PixelMonitoringConf import PixelAthErrorMonAlg
from PixelMonitoring.PixelAthErrorMonAlgCfg import PixelAthErrorMonAlgCfg

from PixelMonitoring.PixelMonitoringConf import PixelAthMVAMonAlg
from PixelMonitoring.PixelAthMVAMonAlgCfg import PixelAthMVAMonAlgCfg

from InDetRecExample.InDetKeys import InDetKeys
from InDetRecExample import TrackingCommon

if forceOnline : 
  isOnline = True
else:
  isOnline = athenaCommonFlags.isOnline()
kwargsHitMonAlg = { 'doOnline'        : isOnline,      #Histograms for online (GlobalMonitoring) running
                     'doLumiBlock'     : not isOnline,     #Turn on/off histograms stored every 1(20) lumi block(s)
                     'doFEPlots'       : True,                                                #Turn on/off per FE-I3 histograms
                     'RDOName'         : InDetKeys.PixelRDOs()
}

kwargsClusMonAlg = { 'doOnline'        : isOnline,     #Histograms for online (GlobalMonitoring) running
                      'doLumiBlock'     : not isOnline,    #Turn on/off histograms stored every 1(20) lumi block(s)
                      'doLowOccupancy'  : InDetFlags.doCosmics(),      #Setting up 1D histogram ranges and binnings, if False, high occupancy i.e. collisions settings will be used
                      'doHeavyIonMon'   : InDetFlags.doHeavyIon(),     #Setting up 1D histogram ranges and binnings for heavy ions
                      'doFEPlots'       : True,                                               #Turn on/off per FE-I3 histograms
                      'ClusterName'     : InDetKeys.PixelClusters(),
                      'TrackName'       : InDetKeys.Tracks()
}

kwargsErrMonAlg = { 'doOnline'        : isOnline,      #Histograms for online (GlobalMonitoring) running
                     'doLumiBlock'     : not isOnline      #Turn on/off histograms stored every 1(20) lumi block(s)
}

kwargsMVAMonAlg = { 'calibFolder'     : '20220503',
                    'RDOName'         : InDetKeys.PixelRDOs(),
                    'ClusterName'     : InDetKeys.PixelClusters(),
                    'TrackParticleContainerName' : InDetKeys.xAODTrackParticleContainer()
}

from AthenaMonitoring.DQMonFlags import DQMonFlags                                                                                                                                      
from AthenaMonitoring import AthMonitorCfgHelperOld
helper = AthMonitorCfgHelperOld(DQMonFlags, "NewPixelMonitoring")

from AthenaMonitoring.FilledBunchFilterTool import GetFilledBunchFilterTool

if doHitMonAlg:
  pixelAthMonAlgHitMonAlg = helper.addAlgorithm(PixelAthHitMonAlg, 'PixelAthHitMonAlg', addFilterTools = [GetFilledBunchFilterTool()])
  for k, v in kwargsHitMonAlg.items():
    setattr(pixelAthMonAlgHitMonAlg, k, v)
  pixelAthMonAlgHitMonAlg.PixelDetElStatus                      = 'PixelDetectorElementStatus'
  pixelAthMonAlgHitMonAlg.PixelDetElStatusActiveOnly            = 'PixelDetectorElementStatusActiveOnly'
  PixelAthHitMonAlgCfg(helper, pixelAthMonAlgHitMonAlg, **kwargsHitMonAlg)

if doClusterMonAlg:
  pixelAthClusterMonAlg = helper.addAlgorithm(PixelAthClusterMonAlg, 'PixelAthClusterMonAlg', addFilterTools = [GetFilledBunchFilterTool()])
  for k, v in kwargsClusMonAlg.items():
    setattr(pixelAthClusterMonAlg, k, v)
  pixelAthClusterMonAlg.HoleSearchTool   = TrackingCommon.getInDetHoleSearchTool()
  pixelAthClusterMonAlg.TrackSelectionTool.UseTrkTrackTools = True
  pixelAthClusterMonAlg.TrackSelectionTool.CutLevel         = "TightPrimary"
  pixelAthClusterMonAlg.TrackSelectionTool.maxNPixelHoles   = 1
  pixelAthClusterMonAlg.TrackSelectionTool.maxD0            = 2
  pixelAthClusterMonAlg.TrackSelectionTool.maxZ0            = 150
  pixelAthClusterMonAlg.TrackSelectionTool.TrackSummaryTool = TrackingCommon.getInDetTrackSummaryTool()
  pixelAthClusterMonAlg.TrackSelectionTool.Extrapolator     = TrackingCommon.getInDetExtrapolator()
  pixelAthClusterMonAlg.PixelDetElStatus                      = 'PixelDetectorElementStatus'
  pixelAthClusterMonAlg.PixelDetElStatusActiveOnly            = 'PixelDetectorElementStatusActiveOnly'
  PixelAthClusterMonAlgCfg(helper, pixelAthClusterMonAlg, **kwargsClusMonAlg)

if doErrorMonAlg:
  pixelAthMonAlgErrorMonAlg = helper.addAlgorithm(PixelAthErrorMonAlg, 'PixelAthErrorMonAlg', addFilterTools = [GetFilledBunchFilterTool()])
  from AthenaCommon.GlobalFlags import globalflags
  isData = (globalflags.DataSource == 'data')

  kwargsErrMonAlg.setdefault(  'PixelByteStreamErrs', 'PixelByteStreamErrs')
  kwargsErrMonAlg.setdefault(  'UseByteStreamFEI4', isData)
  kwargsErrMonAlg.setdefault(  'UseByteStreamFEI3', isData)
  kwargsErrMonAlg.setdefault(  'UseByteStreamRD53', False)
  for k, v in kwargsErrMonAlg.items():
    setattr(pixelAthMonAlgErrorMonAlg, k, v)
  pixelAthMonAlgErrorMonAlg.PixelDetElStatusActiveOnly            = 'PixelDetectorElementStatusActiveOnly'
  PixelAthErrorMonAlgCfg(helper, pixelAthMonAlgErrorMonAlg, **kwargsErrMonAlg)

if doMVAMonAlg:
  pixelAthMVAMonAlg = helper.addAlgorithm(PixelAthMVAMonAlg, 'PixelAthMVAMonAlg', addFilterTools = [GetFilledBunchFilterTool()])
  for k, v in kwargsMVAMonAlg.items():
    setattr(pixelAthMVAMonAlg, k, v)
  pixelAthMVAMonAlg.HoleSearchTool   = TrackingCommon.getInDetHoleSearchTool()
  pixelAthMVAMonAlg.TrackSelectionTool.UseTrkTrackTools = True
  pixelAthMVAMonAlg.TrackSelectionTool.CutLevel         = "TightPrimary"
  #  pixelAthMVAMonAlg.TrackSelectionTool.CutLevel         = "Loose"
  pixelAthMVAMonAlg.TrackSelectionTool.maxNPixelHoles   = 1
  pixelAthMVAMonAlg.TrackSelectionTool.maxD0            = 2
  pixelAthMVAMonAlg.TrackSelectionTool.maxZ0            = 150
  pixelAthMVAMonAlg.TrackSelectionTool.TrackSummaryTool = TrackingCommon.getInDetTrackSummaryTool()
  pixelAthMVAMonAlg.TrackSelectionTool.Extrapolator     = TrackingCommon.getInDetExtrapolator()

  PixelAthMVAMonAlgCfg(helper, pixelAthMVAMonAlg, **kwargsMVAMonAlg)

topSequence += helper.result()



