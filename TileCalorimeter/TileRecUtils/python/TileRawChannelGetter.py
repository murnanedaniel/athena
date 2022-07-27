# Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration

# Author: J. Poveda (Ximo.Poveda@cern.ch)
# TileRawChannel creation from TileDigits 
# TileRawChannelMaker algorithm using

from AthenaCommon.Logging import logging
from AthenaCommon.AthenaCommonFlags import athenaCommonFlags
import traceback

from RecExConfig.Configured import Configured
#from TileRecUtils.TileRawChannelMakerBase import TileRawChannelMakerBase 

class TileRawChannelGetter ( Configured)  :
    """ This getter module creates an instance of the TileRawChannelMaker
    algorithm to obtain TileRawChannel objects from the TileDigits stored
    in a TileDigitsContainer.
    According to the values of the TileRecFlags jobProperties, the
    corresponding AlgTools are added to TileRawChannelMaker algorithm.
    """
    
    _outputType = "TileRawChannelContainer"
    # _output = { _outputType : "TileRawChannelFit" }

    def configure(self):

        mlog = logging.getLogger( 'TileRawChannelGetter::configure:' )
        mlog.info ("entering")

        self._TileRawChannelBuilderFitFilter = None
        self._TileRawChannelBuilderFitFilterCool = None
        self._TileRawChannelBuilderMF = None
        self._TileRawChannelBuilderOF1 = None
        self._TileRawChannelBuilderOpt2Filter = None
        self._TileRawChannelBuilderOptATLAS = None
        self._TileRawChannelBuilderWienerFilter = None

        # Instantiation of the C++ algorithm
        try:        
            from TileRecUtils.TileRecUtilsConf import TileRawChannelMaker               
            theTileRawChannelMaker = TileRawChannelMaker("TileRChMaker")
        except Exception:
            mlog.error("could not import TileRecUtils.TileRawChannelMaker")
            traceback.print_exc()
            return False
    
        self._TileRChMaker = theTileRawChannelMaker

        # Configure TileInfoLoader
        from AthenaCommon.AppMgr import ServiceMgr

        from TileConditions.TileInfoConfigurator import TileInfoConfigurator
        tileInfoConfigurator = TileInfoConfigurator()

        from TileRecUtils.TileRecFlags import jobproperties

        # true for real data, false for MC - GlobalFlags.DataSource.is_data()
        # true for nominal ATLAS configuration - GlobalFlags.DetGeo.is_atlas()
        from AthenaCommon.GlobalFlags import globalflags
        if globalflags.DataSource() == 'data' and jobproperties.TileRecFlags.noiseFilter() < 0:
            # apply noise filter for real data (if this option was not set before)
            jobproperties.TileRecFlags.noiseFilter = 1

        if globalflags.DataSource() == 'data' and not globalflags.isOverlay():
            if jobproperties.TileRecFlags.TileRunType() == 1 :
                tileBeamElemContainer=""
            else:
                tileBeamElemContainer="TileBeamElemCnt"
            if jobproperties.TileRecFlags.readDigits():
                tileDigitsContainer="TileDigitsCnt"
            else:
                tileDigitsContainer=""
            tileRawChannelContainer="TileRawChannelCnt"
        else:
            tileBeamElemContainer=""
            tileDigitsContainer=""
            tileRawChannelContainer=""

        if not globalflags.isOverlay():
            from TileRecUtils.TileDQstatusAlgDefault import TileDQstatusAlgDefault
            TileDQstatusAlgDefault (TileRawChannelContainer = tileRawChannelContainer,
                                    TileDigitsContainer = tileDigitsContainer,
                                    TileBeamElemContainer = tileBeamElemContainer)

        # set time window for amplitude correction if it was not set correctly before
        if jobproperties.TileRecFlags.TimeMaxForAmpCorrection() <= jobproperties.TileRecFlags.TimeMinForAmpCorrection() :
            from AthenaCommon.BeamFlags import jobproperties
            mlog.info("adjusting min/max time of parabolic correction for %s", jobproperties.Beam.bunchSpacing)
            halfBS = jobproperties.Beam.bunchSpacing.get_Value()/2.
            if halfBS > 25.1:
                mlog.info("Bunch spacing is too big, keeping default limits for parabolic correction")
            else:
                jobproperties.TileRecFlags.TimeMinForAmpCorrection = -halfBS
                jobproperties.TileRecFlags.TimeMaxForAmpCorrection = halfBS

        TileFrameLength = ServiceMgr.TileInfoLoader.NSamples
        if TileFrameLength!=7:
            mlog.info("disabling reading of OFC from COOL because Nsamples!=7")
            jobproperties.TileRecFlags.OfcFromCOOL = False

        jobproperties.TileRecFlags.print_JobProperties('tree&value')
        
        # run optimal filter only if readDigits is set
        if jobproperties.TileRecFlags.readDigits():

            from AthenaCommon.AlgSequence import AlgSequence
            topSequence = AlgSequence()

            TilePulseTypes = {0 : 'PHY', 1 : 'PHY', 2 : 'LAS', 4 : 'PHY', 8 : 'CIS'}
            TilePulse = TilePulseTypes[jobproperties.TileRecFlags.TileRunType()]

            from TileConditions.TileCondToolConf import getTileCondToolOfcCool
            toolOfcCool = None
            toolOfcCoolOF1 = None
            toolOfcOnFly = None

            if not jobproperties.TileRecFlags.OfcFromCOOL:
                from TileConditions.TileConditionsConf import TileCondToolOfc
                toolOfcOnFly = TileCondToolOfc()
                toolOfcOnFly.nSamples = TileFrameLength # default = 7
                toolOfcOnFly.OptFilterDeltaCorrelation = False if TileFrameLength == 7 else True  # False - use matrix from DB

            NoiseFilterTools = []
            TileRawChannelContainerDSP = ''
            if jobproperties.TileRecFlags.noiseFilter() == 1:

                if globalflags.DataSource() == 'data':

                    # check if there are DSP thresholds in DB
                    if jobproperties.TileRecFlags.zeroAmplitudeWithoutDigits() and not tileInfoConfigurator.setupCOOLDspThreshold():
                        jobproperties.TileRecFlags.zeroAmplitudeWithoutDigits = False

                    if jobproperties.TileRecFlags.correctPedestalDifference():
                        # check if offline and there are OFCs in DB for OF1 method
                        if not athenaCommonFlags.isOnline():
                            toolOfcCoolOF1 = getTileCondToolOfcCool('COOL', TilePulse, 'OF1', 'TileCondToolOfcCoolOF1')
                        jobproperties.TileRecFlags.correctPedestalDifference = True if toolOfcCoolOF1 else False

                    if jobproperties.TileRecFlags.zeroAmplitudeWithoutDigits() or jobproperties.TileRecFlags.correctPedestalDifference():
                        from TileRecUtils.TileRecUtilsConf import TileRawChannelOF1Corrector
                        theTileRawChannelOF1Corrector = TileRawChannelOF1Corrector()
                        theTileRawChannelOF1Corrector.CorrectPedestalDifference = jobproperties.TileRecFlags.correctPedestalDifference()
                        theTileRawChannelOF1Corrector.ZeroAmplitudeWithoutDigits = jobproperties.TileRecFlags.zeroAmplitudeWithoutDigits()

                        if jobproperties.TileRecFlags.correctPedestalDifference():
                            from TileConditions.TileCondToolConf import getTileCondToolTiming
                            toolOnlineTiming = getTileCondToolTiming('COOL', TilePulse, True, 'TileCondToolOnlineTiming')
                            theTileRawChannelOF1Corrector.TileCondToolTiming = toolOnlineTiming
                            theTileRawChannelOF1Corrector.TileCondToolOfc = toolOfcCoolOF1
                            theTileRawChannelOF1Corrector.TileCondToolNoiseSample.TileOnlineSampleNoise = 'TileOnlineSampleNoise'
                        else:
                            theTileRawChannelOF1Corrector.TileCondToolTiming = None
                            theTileRawChannelOF1Corrector.TileCondToolOfc = None

                        NoiseFilterTools += [theTileRawChannelOF1Corrector]


                from TileRecUtils.TileRecUtilsConf import TileRawChannelNoiseFilter
                theTileRawChannelNoiseFilter = TileRawChannelNoiseFilter()
                NoiseFilterTools += [theTileRawChannelNoiseFilter]

            if globalflags.DataSource() == 'data' and not globalflags.isOverlay():
                if jobproperties.TileRecFlags.correctTimeJumps():
                    from TileRecUtils.TileRecUtilsConf import TileTimeBCOffsetFilter
                    theTileTimeBCOffsetFilter = TileTimeBCOffsetFilter()
                    NoiseFilterTools += [theTileTimeBCOffsetFilter]

                if len(NoiseFilterTools) > 0:
                    from TileRecUtils.TileRecUtilsConf import TileRawChannelCorrectionAlg
                    theTileRawChannelCorrectionAlg = TileRawChannelCorrectionAlg()
                    theTileRawChannelCorrectionAlg.NoiseFilterTools= NoiseFilterTools
                    TileRawChannelContainerDSP = 'TileRawChannelCntCorrected'
                    topSequence += theTileRawChannelCorrectionAlg

            if (jobproperties.TileRecFlags.doTileMF()
                or (not jobproperties.TileRecFlags.OfcFromCOOL()
                    and (jobproperties.TileRecFlags.doTileOF1()
                         or jobproperties.TileRecFlags.doTileWiener()
                         or jobproperties.TileRecFlags.doTileOpt2()
                         or jobproperties.TileRecFlags.doTileOptATLAS()))):

                tileInfoConfigurator.setupCOOLPULSE(type = TilePulse)
                tileInfoConfigurator.setupCOOLAutoCr()

            elif jobproperties.TileRecFlags.doTileFitCool():
                tileInfoConfigurator.setupCOOLPULSE(type = TilePulse)

            if jobproperties.TileRecFlags.OfcFromCOOL():
                if (jobproperties.TileRecFlags.doTileMF()
                    or jobproperties.TileRecFlags.doTileWiener()
                    or jobproperties.TileRecFlags.doTileOpt2()
                    or jobproperties.TileRecFlags.doTileOptATLAS()):

                    tileInfoConfigurator.setupCOOLOFC(type = TilePulse)
                    toolOfcCool = getTileCondToolOfcCool('COOL', TilePulse)

                if jobproperties.TileRecFlags.doTileOF1():
                    tileInfoConfigurator.setupCOOLOFC(type = TilePulse, ofcType = 'OF1')

                    if not toolOfcCoolOF1:
                        toolOfcCoolOF1 = getTileCondToolOfcCool('COOL', TilePulse, 'OF1', 'TileCondToolOfcCoolOF1')


            if jobproperties.TileRecFlags.doTileQIE():
                try:
                    from TileRecUtils.TileRecUtilsConf import TileRawChannelBuilderQIEFilter
                    theTileRawChannelBuilderQIEFilter= TileRawChannelBuilderQIEFilter()
                except Exception:
                    mlog.error("could not get handle to TileRawChannelBuilderQIEFilter Quit")
                    traceback.print_exc()
                    return False
                    
                #TileRawChannelBuilderQIEFilter Options:
                jobproperties.TileRecFlags.TileRawChannelContainer = "TileRawChannelQIE"
                theTileRawChannelBuilderQIEFilter.TileRawChannelContainer = "TileRawChannelQIE"
                theTileRawChannelBuilderQIEFilter.RunType = jobproperties.TileRecFlags.TileRunType()
                theTileRawChannelBuilderQIEFilter.calibrateEnergy = jobproperties.TileRecFlags.calibrateEnergy()
                theTileRawChannelBuilderQIEFilter.correctTime     = jobproperties.TileRecFlags.correctTime()
                theTileRawChannelBuilderQIEFilter.NoiseFilterTools= NoiseFilterTools
                theTileRawChannelBuilderQIEFilter.PedestalMode = 1
                theTileRawChannelBuilderQIEFilter.DSPContainer = TileRawChannelContainerDSP
      
                mlog.info(" adding now TileRawChannelBuilderQIEFilter to the algorithm: %s", theTileRawChannelMaker.name())
      
                theTileRawChannelMaker.TileRawChannelBuilder += [theTileRawChannelBuilderQIEFilter]
            
            # fit with several amplitudes 
            if jobproperties.TileRecFlags.doTileManyAmps():
      
                try:
                    from TileRecUtils.TileRecUtilsConf import TileRawChannelBuilderManyAmps
                    theTileRawChannelBuilderManyAmps= TileRawChannelBuilderManyAmps()
                except Exception:
                    mlog.error("could not get handle to TileRawChannelBuilderManyAmps Quit")
                    traceback.print_exc()
                    return False
      
                #TileRawChannelBuilderManyAmps Options:
                jobproperties.TileRecFlags.TileRawChannelContainer = "TileRawChannelManyAmp"
                theTileRawChannelBuilderManyAmps.TileRawChannelContainer = "TileRawChannelManyAmp"
                theTileRawChannelBuilderManyAmps.RunType = jobproperties.TileRecFlags.TileRunType()
                theTileRawChannelBuilderManyAmps.calibrateEnergy = jobproperties.TileRecFlags.calibrateEnergy()
                theTileRawChannelBuilderManyAmps.correctTime     = jobproperties.TileRecFlags.correctTime()    
                theTileRawChannelBuilderManyAmps.NoiseFilterTools= NoiseFilterTools
                theTileRawChannelBuilderManyAmps.DSPContainer = TileRawChannelContainerDSP
                 
                mlog.info(" adding now TileRawChannelBuilderManyAmps to the algorithm: %s", theTileRawChannelMaker.name())
      
                theTileRawChannelMaker.TileRawChannelBuilder += [theTileRawChannelBuilderManyAmps]
      
            # flat filter - sum of 5 samples
            if jobproperties.TileRecFlags.doTileFlat():
      
                try:
                    from TileRecUtils.TileRecUtilsConf import TileRawChannelBuilderFlatFilter
                    theTileRawChannelBuilderFlatFilter= TileRawChannelBuilderFlatFilter()
                except Exception:
                    mlog.error("could not get handle to TileRawChannelBuilderFlatFilter Quit")
                    traceback.print_exc()
                    return False
      
                #TileRawChannelBuilderFlatFilter Options: 
                jobproperties.TileRecFlags.TileRawChannelContainer = "TileRawChannelFlat"
                theTileRawChannelBuilderFlatFilter.TileRawChannelContainer = "TileRawChannelFlat"
                theTileRawChannelBuilderFlatFilter.RunType = jobproperties.TileRecFlags.TileRunType()
                theTileRawChannelBuilderFlatFilter.calibrateEnergy = jobproperties.TileRecFlags.calibrateEnergy()
                theTileRawChannelBuilderFlatFilter.correctTime     = jobproperties.TileRecFlags.correctTime()    
                theTileRawChannelBuilderFlatFilter.NoiseFilterTools= NoiseFilterTools
                theTileRawChannelBuilderFlatFilter.FrameLength = TileFrameLength
                theTileRawChannelBuilderFlatFilter.SignalLength = TileFrameLength - 1
                theTileRawChannelBuilderFlatFilter.DSPContainer = TileRawChannelContainerDSP
      
                mlog.info(" adding now TileRawChannelBuilderFlatFilter to the algorithm: %s", theTileRawChannelMaker.name())
      
                theTileRawChannelMaker.TileRawChannelBuilder += [theTileRawChannelBuilderFlatFilter]
      
            # Fit method
            if jobproperties.TileRecFlags.doTileFit() or jobproperties.TileRecFlags.doTileOverflowFit():
          
                # configure TileRawChannelMaker here
                try:
                    from TileRecUtils.TileRecUtilsConf import TileRawChannelBuilderFitFilter
                    theTileRawChannelBuilderFitFilter= TileRawChannelBuilderFitFilter()
                except Exception:
                    mlog.error("could not get handle to TileRawChannelBuilderFitFilter Quit")
                    traceback.print_exc()
                    return False
                
                #TileRawChannelBuilderFitFilter Options: 
                theTileRawChannelBuilderFitFilter.RunType = jobproperties.TileRecFlags.TileRunType()
                theTileRawChannelBuilderFitFilter.calibrateEnergy = jobproperties.TileRecFlags.calibrateEnergy()
                theTileRawChannelBuilderFitFilter.correctTime     = jobproperties.TileRecFlags.correctTime()    
                theTileRawChannelBuilderFitFilter.NoiseFilterTools= NoiseFilterTools
                theTileRawChannelBuilderFitFilter.FrameLength = TileFrameLength
                theTileRawChannelBuilderFitFilter.DSPContainer = TileRawChannelContainerDSP
                
                # add the tool to list of tool ( should use ToolHandle eventually)
                mlog.info(" adding now TileRawChannelBuilderFitFilter to the algorithm: %s", theTileRawChannelMaker.name())

                if jobproperties.TileRecFlags.doTileFit():
                    jobproperties.TileRecFlags.TileRawChannelContainer = "TileRawChannelFit"
                    theTileRawChannelBuilderFitFilter.TileRawChannelContainer = "TileRawChannelFit"
                    theTileRawChannelMaker.TileRawChannelBuilder += [theTileRawChannelBuilderFitFilter]
                    self._TileRawChannelBuilderFitFilter = theTileRawChannelMaker.TileRawChannelBuilder['TileRawChannelBuilderFitFilter']
                    
                if jobproperties.TileRecFlags.doTileOverflowFit(): 
                    theTileRawChannelMaker.FitOverflow = True
                    theTileRawChannelMaker.TileRawChannelBuilderFitOverflow = theTileRawChannelBuilderFitFilter
                    mlog.info(" set up TileRawChannelBuilderFitOverflow to TileRawChannelBuilderFitFilter")
                
            # Fit method with reading from COOL
            if jobproperties.TileRecFlags.doTileFitCool():
          
                # configure TileRawChannelMaker here
                try:
                    from TileRecUtils.TileRecUtilsConf import TileRawChannelBuilderFitFilterCool
                    theTileRawChannelBuilderFitFilterCool= TileRawChannelBuilderFitFilterCool()
                except Exception:
                    mlog.error("could not get handle to TileRawChannelBuilderFitFilterCool Quit")
                    traceback.print_exc()
                    return False
                
                #TileRawChannelBuilderFitFilterCool Options: 
                jobproperties.TileRecFlags.TileRawChannelContainer = "TileRawChannelFitCool"
                theTileRawChannelBuilderFitFilterCool.TileRawChannelContainer = "TileRawChannelFitCool"
                theTileRawChannelBuilderFitFilterCool.RunType = jobproperties.TileRecFlags.TileRunType()
                theTileRawChannelBuilderFitFilterCool.calibrateEnergy = jobproperties.TileRecFlags.calibrateEnergy()
                theTileRawChannelBuilderFitFilterCool.correctTime     = jobproperties.TileRecFlags.correctTime()    
                theTileRawChannelBuilderFitFilterCool.NoiseFilterTools= NoiseFilterTools
                theTileRawChannelBuilderFitFilterCool.FrameLength = TileFrameLength
                theTileRawChannelBuilderFitFilterCool.DSPContainer = TileRawChannelContainerDSP
                
                # add the tool to list of tool ( should use ToolHandle eventually)
                mlog.info(" adding now TileRawChannelBuilderFitFilterCool to the algorithm: %s", theTileRawChannelMaker.name())
      
                theTileRawChannelMaker.TileRawChannelBuilder += [theTileRawChannelBuilderFitFilterCool]
                self._TileRawChannelBuilderFitFilterCool = theTileRawChannelMaker.TileRawChannelBuilder['TileRawChannelBuilderFitFilterCool']
                
            # matched filter 
            if jobproperties.TileRecFlags.doTileMF():
      
                try:
                    from TileRecUtils.TileRecUtilsConf import TileRawChannelBuilderMF
                    theTileRawChannelBuilderMF= TileRawChannelBuilderMF()
                except Exception:
                    mlog.error("could not get handle to TileRawChannelBuilderMF Quit")
                    traceback.print_exc()
                    return False
                                    
                # setup COOL to get OFCs, needed for COF to retrieve pulse shape and derivatives
                if jobproperties.TileRecFlags.OfcFromCOOL():
                    theTileRawChannelBuilderMF.TileCondToolOfc = toolOfcCool

                #TileRawChannelBuilderMF Options:
                jobproperties.TileRecFlags.TileRawChannelContainer = "TileRawChannelMF"
                theTileRawChannelBuilderMF.TileRawChannelContainer = "TileRawChannelMF"
                theTileRawChannelBuilderMF.RunType = jobproperties.TileRecFlags.TileRunType()
                theTileRawChannelBuilderMF.calibrateEnergy = jobproperties.TileRecFlags.calibrateEnergy()
                if jobproperties.TileRecFlags.BestPhaseFromCOOL(): # can't correct time and use best phase at the same time
                    theTileRawChannelBuilderMF.correctTime = False
                else:
                    theTileRawChannelBuilderMF.correctTime = jobproperties.TileRecFlags.correctTime()
                theTileRawChannelBuilderMF.BestPhase       = jobproperties.TileRecFlags.BestPhaseFromCOOL()
                theTileRawChannelBuilderMF.NoiseFilterTools = NoiseFilterTools
                theTileRawChannelBuilderMF.MaxIterations = 5 # iterative mode on
                theTileRawChannelBuilderMF.AmplitudeCorrection = False
                theTileRawChannelBuilderMF.TimeFromCOF = False
                theTileRawChannelBuilderMF.AmpMinForAmpCorrection = jobproperties.TileRecFlags.AmpMinForAmpCorrection()
                if jobproperties.TileRecFlags.TimeMaxForAmpCorrection() > jobproperties.TileRecFlags.TimeMinForAmpCorrection():
                    theTileRawChannelBuilderMF.TimeMinForAmpCorrection = jobproperties.TileRecFlags.TimeMinForAmpCorrection()
                    theTileRawChannelBuilderMF.TimeMaxForAmpCorrection = jobproperties.TileRecFlags.TimeMaxForAmpCorrection()

                theTileRawChannelBuilderMF.DSPContainer = TileRawChannelContainerDSP
                mlog.info(" adding now TileRawChannelBuilderMF to the algorithm: %s", theTileRawChannelMaker.name())
      
                theTileRawChannelMaker.TileRawChannelBuilder += [theTileRawChannelBuilderMF]
                self._TileRawChannelBuilderMF = theTileRawChannelMaker.TileRawChannelBuilder['TileRawChannelBuilderMF']

            if jobproperties.TileRecFlags.doTileOF1():
                try:
                    from TileRecUtils.TileRecUtilsConf import TileRawChannelBuilderOpt2Filter
                    theTileRawChannelBuilderOF1 = TileRawChannelBuilderOpt2Filter("TileRawChannelBuilderOF1")
                except Exception:
                    mlog.error("could not get handle to TileRawChannelBuilderOF1 Quit")
                    traceback.print_exc()
                    return False
                
                # setup COOL to get OFCs
                if jobproperties.TileRecFlags.OfcFromCOOL():
                    if toolOfcCoolOF1:
                        theTileRawChannelBuilderOF1.TileCondToolOfc = toolOfcCoolOF1
                    else:
                        # There are no OF1 OFC in the COOL
                        # OFC will be calculated on the fly
                        tileInfoConfigurator.setupCOOLPULSE(type = TilePulse)
                        tileInfoConfigurator.setupCOOLAutoCr()

                #TileRawChannelBuilderOF1 Options:
                jobproperties.TileRecFlags.TileRawChannelContainer = "TileRawChannelOF1"
                theTileRawChannelBuilderOF1.TileRawChannelContainer = "TileRawChannelOF1"
                theTileRawChannelBuilderOF1.RunType = jobproperties.TileRecFlags.TileRunType()
                theTileRawChannelBuilderOF1.calibrateEnergy = jobproperties.TileRecFlags.calibrateEnergy()
                if jobproperties.TileRecFlags.BestPhaseFromCOOL(): # can't correct time and use best phase at the same time
                    theTileRawChannelBuilderOF1.correctTime = False
                else:
                    theTileRawChannelBuilderOF1.correctTime = jobproperties.TileRecFlags.correctTime()
                theTileRawChannelBuilderOF1.BestPhase       = jobproperties.TileRecFlags.BestPhaseFromCOOL()
                theTileRawChannelBuilderOF1.NoiseFilterTools= NoiseFilterTools
                theTileRawChannelBuilderOF1.OF2 = False
                theTileRawChannelBuilderOF1.PedestalMode = -1
                theTileRawChannelBuilderOF1.MaxIterations = 1 # just one iteration
                theTileRawChannelBuilderOF1.Minus1Iteration = False # assume that max sample is at t=0
                theTileRawChannelBuilderOF1.AmplitudeCorrection = jobproperties.TileRecFlags.correctAmplitude()
                theTileRawChannelBuilderOF1.TimeCorrection = False
                theTileRawChannelBuilderOF1.AmpMinForAmpCorrection = jobproperties.TileRecFlags.AmpMinForAmpCorrection()
                if jobproperties.TileRecFlags.TimeMaxForAmpCorrection() > jobproperties.TileRecFlags.TimeMinForAmpCorrection():
                    theTileRawChannelBuilderOF1.TimeMinForAmpCorrection = jobproperties.TileRecFlags.TimeMinForAmpCorrection()
                    theTileRawChannelBuilderOF1.TimeMaxForAmpCorrection = jobproperties.TileRecFlags.TimeMaxForAmpCorrection()
      
                theTileRawChannelBuilderOF1.DSPContainer = TileRawChannelContainerDSP
                mlog.info(" adding now TileRawChannelBuilderOF1 to the algorithm: %s", theTileRawChannelMaker.name())
      
                theTileRawChannelMaker.TileRawChannelBuilder += [theTileRawChannelBuilderOF1]
                self._TileRawChannelBuilderOF1 = theTileRawChannelMaker.TileRawChannelBuilder['TileRawChannelBuilderOF1']

            if jobproperties.TileRecFlags.doTileOpt2():
                try:
                    from TileRecUtils.TileRecUtilsConf import TileRawChannelBuilderOpt2Filter
                    theTileRawChannelBuilderOpt2Filter= TileRawChannelBuilderOpt2Filter()
                except Exception:
                    mlog.error("could not get handle to TileRawChannelBuilderOpt2Filter Quit")
                    traceback.print_exc()
                    return False
                
                # setup COOL to get OFCs
                if jobproperties.TileRecFlags.OfcFromCOOL():
                    theTileRawChannelBuilderOpt2Filter.TileCondToolOfc = toolOfcCool
                else:
                    theTileRawChannelBuilderOpt2Filter.TileCondToolOfc = toolOfcOnFly
                    
                #TileRawChannelBuilderOpt2Filter Options:
                jobproperties.TileRecFlags.TileRawChannelContainer = "TileRawChannelOpt2"
                theTileRawChannelBuilderOpt2Filter.TileRawChannelContainer = "TileRawChannelOpt2"
                theTileRawChannelBuilderOpt2Filter.RunType = jobproperties.TileRecFlags.TileRunType()
                theTileRawChannelBuilderOpt2Filter.calibrateEnergy = jobproperties.TileRecFlags.calibrateEnergy()
                theTileRawChannelBuilderOpt2Filter.correctTime     = jobproperties.TileRecFlags.correctTime()
                theTileRawChannelBuilderOpt2Filter.NoiseFilterTools= NoiseFilterTools
                theTileRawChannelBuilderOpt2Filter.BestPhase       = False # no point to use best phase with interations
                theTileRawChannelBuilderOpt2Filter.OF2 = True
                theTileRawChannelBuilderOpt2Filter.PedestalMode = 1
                theTileRawChannelBuilderOpt2Filter.MaxIterations = 5
                theTileRawChannelBuilderOpt2Filter.Minus1Iteration = True
                theTileRawChannelBuilderOpt2Filter.AmplitudeCorrection = False # don't need correction after iterations
                theTileRawChannelBuilderOpt2Filter.TimeCorrection    = False # don't need correction after iterations
                theTileRawChannelBuilderOpt2Filter.DSPContainer = TileRawChannelContainerDSP
      
                mlog.info(" adding now TileRawChannelBuilderOpt2Filter to the algorithm: %s", theTileRawChannelMaker.name())
      
                theTileRawChannelMaker.TileRawChannelBuilder += [theTileRawChannelBuilderOpt2Filter]
                self._TileRawChannelBuilderOpt2Filter = theTileRawChannelMaker.TileRawChannelBuilder['TileRawChannelBuilderOpt2Filter']
      
      
            if jobproperties.TileRecFlags.doTileOptATLAS():
                try:
                    from TileRecUtils.TileRecUtilsConf import TileRawChannelBuilderOpt2Filter
                    theTileRawChannelBuilderOptATLAS= TileRawChannelBuilderOpt2Filter("TileRawChannelBuilderOptATLAS")
                except Exception:
                    mlog.error("could not get handle to TileRawChannelBuilderOpt2Filter Quit")
                    traceback.print_exc()
                    return False
                
                # setup COOL to get OFCs
                if jobproperties.TileRecFlags.OfcFromCOOL():
                    theTileRawChannelBuilderOptATLAS.TileCondToolOfc = toolOfcCool
                else:
                    theTileRawChannelBuilderOptATLAS.TileCondToolOfc = toolOfcOnFly
                    
                #TileRawChannelBuilderOptATLAS Options:
                if globalflags.DataSource()=='data': # don't use the name which is used for reco data from DSP
                    if jobproperties.TileRecFlags.TileRawChannelContainer == "TileRawChannelCnt":
                        jobproperties.TileRecFlags.TileRawChannelContainer = "TileRawChannelFixed"
                    theTileRawChannelBuilderOptATLAS.TileRawChannelContainer = "TileRawChannelFixed"
                else:
                    jobproperties.TileRecFlags.TileRawChannelContainer = "TileRawChannelCnt"
                    theTileRawChannelBuilderOptATLAS.TileRawChannelContainer = "TileRawChannelCnt"
                theTileRawChannelBuilderOptATLAS.RunType         = jobproperties.TileRecFlags.TileRunType()
                theTileRawChannelBuilderOptATLAS.calibrateEnergy = jobproperties.TileRecFlags.calibrateEnergy()
                if jobproperties.TileRecFlags.BestPhaseFromCOOL(): # can't correct time and use best phase at the same time
                    theTileRawChannelBuilderOptATLAS.correctTime = False
                else:
                    theTileRawChannelBuilderOptATLAS.correctTime = jobproperties.TileRecFlags.correctTime()
                theTileRawChannelBuilderOptATLAS.BestPhase       = jobproperties.TileRecFlags.BestPhaseFromCOOL()
                theTileRawChannelBuilderOptATLAS.NoiseFilterTools= NoiseFilterTools
                theTileRawChannelBuilderOptATLAS.OF2 = True
                #theTileRawChannelBuilderOptATLAS.PedestalMode = 1 # not sure if we need this option here
                theTileRawChannelBuilderOptATLAS.MaxIterations = 1 # just one iteration
                theTileRawChannelBuilderOptATLAS.Minus1Iteration = False # assume that max sample is at t=0
                theTileRawChannelBuilderOptATLAS.AmplitudeCorrection = jobproperties.TileRecFlags.correctAmplitude()
                theTileRawChannelBuilderOptATLAS.TimeCorrection = jobproperties.TileRecFlags.correctTimeNI()
                theTileRawChannelBuilderOptATLAS.AmpMinForAmpCorrection = jobproperties.TileRecFlags.AmpMinForAmpCorrection()
                if jobproperties.TileRecFlags.TimeMaxForAmpCorrection() > jobproperties.TileRecFlags.TimeMinForAmpCorrection():
                    theTileRawChannelBuilderOptATLAS.TimeMinForAmpCorrection = jobproperties.TileRecFlags.TimeMinForAmpCorrection()
                    theTileRawChannelBuilderOptATLAS.TimeMaxForAmpCorrection = jobproperties.TileRecFlags.TimeMaxForAmpCorrection()

                theTileRawChannelBuilderOptATLAS.DSPContainer = TileRawChannelContainerDSP
                
                mlog.info(" adding now TileRawChannelBuilderOpt2Filter with name TileRawChannelBuilderOptATLAS to the algorithm: %s",
                          theTileRawChannelMaker.name())
                
                theTileRawChannelMaker.TileRawChannelBuilder += [theTileRawChannelBuilderOptATLAS]
                self._TileRawChannelBuilderOptATLAS = theTileRawChannelMaker.TileRawChannelBuilder['TileRawChannelBuilderOptATLAS']

            if jobproperties.TileRecFlags.doTileWiener():
                try:
                    from TileRecUtils.TileRecUtilsConf import TileRawChannelBuilderWienerFilter
                    theTileRawChannelBuilderWienerFilter= TileRawChannelBuilderWienerFilter()
                except Exception:
                    mlog.error("could not get handle to TileRawChannelBuilderWienerFilter Quit")
                    traceback.print_exc()
                    return False

                #TileRawChannelBuilderWienerFilter Options:
                jobproperties.TileRecFlags.TileRawChannelContainer           = "TileRawChannelWiener"
                theTileRawChannelBuilderWienerFilter.TileRawChannelContainer = "TileRawChannelWiener"
                theTileRawChannelBuilderWienerFilter.RunType                 = jobproperties.TileRecFlags.TileRunType()
                theTileRawChannelBuilderWienerFilter.calibrateEnergy         = jobproperties.TileRecFlags.calibrateEnergy()
                theTileRawChannelBuilderWienerFilter.correctTime             = jobproperties.TileRecFlags.correctTime()
                theTileRawChannelBuilderWienerFilter.NoiseFilterTools        = NoiseFilterTools
                theTileRawChannelBuilderWienerFilter.BestPhase               = False # no point to use best phase with interations
                theTileRawChannelBuilderWienerFilter.MC                      = globalflags.DataSource()!='data'
                theTileRawChannelBuilderWienerFilter.PedestalMode            = 1
                theTileRawChannelBuilderWienerFilter.MaxIterations           = 5
                theTileRawChannelBuilderWienerFilter.Minus1Iteration         = True
                theTileRawChannelBuilderWienerFilter.AmplitudeCorrection     = False # don't need correction after iterations
                theTileRawChannelBuilderWienerFilter.TimeCorrection          = False # don't need correction after iterations
                theTileRawChannelBuilderWienerFilter.DSPContainer            = TileRawChannelContainerDSP

                from LumiBlockComps.BunchCrossingCondAlgDefault import BunchCrossingCondAlgDefault
                BunchCrossingCondAlgDefault()


                ServiceMgr.TileInfoLoader.LoadWienerFilterWeights            = True

                mlog.info(" adding now TileRawChannelBuilderWienerFilter to the algorithm: %s", theTileRawChannelMaker.name())

                theTileRawChannelMaker.TileRawChannelBuilder += [theTileRawChannelBuilderWienerFilter]
                self._TileRawChannelBuilderWienerFilter = theTileRawChannelMaker.TileRawChannelBuilder['TileRawChannelBuilderWienerFilter']
            

            # now add algorithm to topSequence
            # this should always come at the end

            mlog.info(" now adding to topSequence")        

            if jobproperties.TileRecFlags.noiseFilter() == 2:
                # Instantiation of the C++ algorithm
                try:
                    from TileRecUtils.TileRecUtilsConf import TileRawCorrelatedNoise
                    theTileRawCorrelatedNoise=TileRawCorrelatedNoise("TileRCorreNoise")
                except Exception:
                    mlog.error("could not import TileRecUtils.TileRawCorrelatedNoise")
                    traceback.print_exc()
                    return False
                #theTileRawCorrelatedNoise.UseMeanFiles = False
                #theTileRawCorrelatedNoise.PMTOrder = True
                jobproperties.TileRecFlags.TileDigitsContainer = "NewDigitsContainer"
                topSequence += theTileRawCorrelatedNoise

            jobproperties.TileRecFlags.print_JobProperties('tree&value')

            theTileRawChannelMaker.TileDigitsContainer = jobproperties.TileRecFlags.TileDigitsContainer()
            topSequence += theTileRawChannelMaker

        else:
            mlog.info(" Disable all OF methods because readDigits flag set to False ")
            jobproperties.TileRecFlags.doTileFlat = False
            jobproperties.TileRecFlags.doTileFit = False
            jobproperties.TileRecFlags.doTileFitCool = False
            jobproperties.TileRecFlags.doTileOpt2 = False
            jobproperties.TileRecFlags.doTileOptATLAS = False
            jobproperties.TileRecFlags.doTileManyAmps = False
            jobproperties.TileRecFlags.doTileMF = False
            jobproperties.TileRecFlags.doTileOF1 = False
            jobproperties.TileRecFlags.doTileWiener = False
            jobproperties.TileRecFlags.OfcFromCOOL = False
            jobproperties.TileRecFlags.print_JobProperties('tree&value')

        return True

    def TileRChMaker(self):
        return self._TileRChMaker

    def TileRawChannelBuilderFitFilter(self):
        return self._TileRawChannelBuilderFitFilter

    def TileRawChannelBuilderFitFilterCool(self):
        return self._TileRawChannelBuilderFitFilterCool

    def TileRawChannelBuilderMF(self):
        return self._TileRawChannelBuilderMF

    def TileRawChannelBuilderOF1(self):
        return self._TileRawChannelBuilderOF1

    def TileRawChannelBuilderOpt2Filter(self):
        return self._TileRawChannelBuilderOpt2Filter

    def TileRawChannelBuilderOptATLAS(self):
        return self._TileRawChannelBuilderOptATLAS

    def TileRawChannelBuilderWienerFilter(self):
        return self._TileRawChannelBuilderWienerFilter

##    # would work only if one output object type
##    def outputKey(self):
##        return self._output[self._outputType]

##    def outputType(self):
##        return self._outputType


