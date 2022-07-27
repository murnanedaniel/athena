# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

from __future__ import division

# get the logger
from AthenaCommon.Logging import logging
digilog = logging.getLogger('Digi_tf')
digilog.info( '****************** STARTING DIGITIZATION *****************' )

include("SimuJobTransforms/CommonSkeletonJobOptions.py")

if hasattr(runArgs, "jobNumber"):
    if runArgs.jobNumber < 1:
        raise ValueError('jobNumber must be a postive integer. %s lies outside this range', str(runArgs.jobNumber))

from AthenaCommon.GlobalFlags import globalflags
if hasattr(runArgs,"geometryVersion"):
    # strip _VALIDATION
    digilog.info("stripping _VALIDATION")
    if runArgs.geometryVersion.endswith("_VALIDATION"):
        pos=runArgs.geometryVersion.find("_VALIDATION")
        globalflags.DetDescrVersion.set_Value_and_Lock( runArgs.geometryVersion[:pos] )
    else:
        globalflags.DetDescrVersion.set_Value_and_Lock( runArgs.geometryVersion )


#==============================================================
# Job definition parameters:
#==============================================================
from AthenaCommon.AthenaCommonFlags import athenaCommonFlags
#Jobs should stop if an include fails.
if hasattr(runArgs,"IgnoreConfigError"):
    athenaCommonFlags.AllowIgnoreConfigError=runArgs.IgnoreConfigError
else:
    athenaCommonFlags.AllowIgnoreConfigError=False

if hasattr(runArgs,"inputFile"): athenaCommonFlags.FilesInput.set_Value_and_Lock( runArgs.inputFile )
if hasattr(runArgs,"inputHITSFile"):
    globalflags.InputFormat.set_Value_and_Lock('pool')
    athenaCommonFlags.PoolHitsInput.set_Value_and_Lock( runArgs.inputHITSFile )
elif hasattr(runArgs,"inputHitsFile"):
    globalflags.InputFormat.set_Value_and_Lock('pool')
    athenaCommonFlags.PoolHitsInput.set_Value_and_Lock( runArgs.inputHitsFile )
elif not hasattr(runArgs,"testPileUpConfig"):
    raise RuntimeError("No inputHITSFile provided.")

# runNumber is MC channel number in reco
if hasattr(runArgs, 'runNumber'):
    # always set it in legacy config
    athenaCommonFlags.MCChannelNumber.set_Value(runArgs.runNumber)
    digilog.info('Got MC channel number %d from runNumber', athenaCommonFlags.MCChannelNumber())

digilog.info( '**** Transformation run arguments' )
digilog.info( str(runArgs) )

#==============================================================
# Job Configuration parameters:
#==============================================================
## Pre-exec
if hasattr(runArgs,"preExec"):
    digilog.info("transform pre-exec")
    for cmd in runArgs.preExec:
        digilog.info(cmd)
        exec(cmd)

## Pre-include
if hasattr(runArgs,"preInclude"):
    for fragment in runArgs.preInclude:
        include(fragment)

#--------------------------------------------------------------
# Override pile-up configuration on the command-line
#--------------------------------------------------------------
from Digitization.DigitizationFlags import digitizationFlags
PileUpConfigOverride=False
import math

## First check for depreacted command-line options
if hasattr(runArgs,"numberOfLowPtMinBias"):
    if not math.fabs(digitizationFlags.numberOfLowPtMinBias.get_Value()-runArgs.numberOfLowPtMinBias)<0.00000001: #FIXME comparing two floats
        digilog.info( "Changing digitizationFlags.numberOfLowPtMinBias from %s to %s", digitizationFlags.numberOfLowPtMinBias.get_Value(),runArgs.numberOfLowPtMinBias)
        digitizationFlags.numberOfLowPtMinBias=float(runArgs.numberOfLowPtMinBias)
        PileUpConfigOverride=True
if hasattr(runArgs,"numberOfHighPtMinBias"):
    if not math.fabs(digitizationFlags.numberOfHighPtMinBias.get_Value()-runArgs.numberOfHighPtMinBias)<0.00000001: #FIXME comparing two floats
        digilog.info( "Changing digitizationFlags.numberOfHighPtMinBias from %s to %s", digitizationFlags.numberOfHighPtMinBias.get_Value(),runArgs.numberOfHighPtMinBias)
        digitizationFlags.numberOfHighPtMinBias=float(runArgs.numberOfHighPtMinBias)
        PileUpConfigOverride=True
if hasattr(runArgs,"numberOfBeamHalo"):
    if not math.fabs(digitizationFlags.numberOfBeamHalo.get_Value()-runArgs.numberOfBeamHalo)<0.00000001: #FIXME comparing two floats
        digilog.info( "Changing digitizationFlags.numberOfBeamHalo from %s to %s", digitizationFlags.numberOfBeamHalo.get_Value(),runArgs.numberOfBeamHalo)
        digitizationFlags.numberOfBeamHalo=float(runArgs.numberOfBeamHalo)
        PileUpConfigOverride=True
if hasattr(runArgs,"numberOfBeamGas"):
    if not math.fabs(digitizationFlags.numberOfBeamGas.get_Value()-runArgs.numberOfBeamGas)<0.00000001: #FIXME comparing two floats
        digilog.info( "Changing digitizationFlags.numberOfBeamGas from %s to %s", digitizationFlags.numberOfBeamGas.get_Value(),runArgs.numberOfBeamGas)
        digitizationFlags.numberOfBeamGas=float(runArgs.numberOfBeamGas)
        PileUpConfigOverride=True
if hasattr(runArgs,"numberOfCavernBkg"):
    if not digitizationFlags.numberOfCavern.get_Value()==runArgs.numberOfCavernBkg:
        digilog.info( "Changing digitizationFlags.cavernEvents from %s to %s", digitizationFlags.numberOfCavern.get_Value(),runArgs.numberOfCavernBkg)
        digitizationFlags.numberOfCavern=runArgs.numberOfCavernBkg
        PileUpConfigOverride=True
if hasattr(runArgs,"bunchSpacing"):
    if  digitizationFlags.BeamIntensityPattern.statusOn:
        digilog.warning("Redefine bunch-structure with a fixed %s ns bunchSpacing. Overwriting the previous setting of %s", runArgs.bunchSpacing,digitizationFlags.BeamIntensityPattern.get_Value())
    else:
        digilog.info( "Setting up job to run with a fixed %s ns bunchSpacing.", runArgs.bunchSpacing)
    digitizationFlags.bunchSpacing = 25
    digitizationFlags.BeamIntensityPattern.createConstBunchSpacingPattern(int(runArgs.bunchSpacing)) #FIXME This runArg should probably inherit from argInt rather than argFloat
    digilog.info( "New bunch-structure = %s", digitizationFlags.BeamIntensityPattern.get_Value())
    from AthenaCommon.BeamFlags import jobproperties
    jobproperties.Beam.bunchSpacing = int(runArgs.bunchSpacing) #FIXME This runArg should probably inherit from argInt rather than argFloat
    PileUpConfigOverride=True
if hasattr(runArgs,"pileupInitialBunch"):
    if not (digitizationFlags.initialBunchCrossing.get_Value()==runArgs.pileupInitialBunch):
        digilog.info( "Changing digitizationFlags.initialBunchCrossing from %s to %s", digitizationFlags.initialBunchCrossing.get_Value(),runArgs.pileupInitialBunch)
        digitizationFlags.initialBunchCrossing=runArgs.pileupInitialBunch
        PileUpConfigOverride=True
if hasattr(runArgs,"pileupFinalBunch"):
    if not (digitizationFlags.finalBunchCrossing.get_Value()==runArgs.pileupFinalBunch):
        digilog.info( "Changing digitizationFlags.finalBunchCrossing from %s to %s", digitizationFlags.finalBunchCrossing.get_Value(),runArgs.pileupFinalBunch)
        digitizationFlags.finalBunchCrossing=runArgs.pileupFinalBunch
        PileUpConfigOverride=True
if hasattr(runArgs,"digiSteeringConf"):
    if not (digitizationFlags.digiSteeringConf.get_Value()==runArgs.digiSteeringConf+"PileUpToolsAlg"):
        digilog.info( "Changing digitizationFlags.digiSteeringConf from %s to %s", digitizationFlags.digiSteeringConf.get_Value(),runArgs.digiSteeringConf)
        digitizationFlags.digiSteeringConf=runArgs.digiSteeringConf+"PileUpToolsAlg"
        PileUpConfigOverride=True
if digitizationFlags.initialBunchCrossing > digitizationFlags.finalBunchCrossing:
    raise ValueError( "Initial bunch crossing should not be larger than the final one" )
if PileUpConfigOverride:
    digilog.info( "NB Some pile-up (re-)configuration was done on the command-line.")
del PileUpConfigOverride
#--------------------------------------------------------------
# Test Pile-up configuration
#--------------------------------------------------------------
from SimuJobTransforms.SimTransformUtils import pileUpCalc
if hasattr(runArgs,"testPileUpConfig"):
    nSignalEvts=1000
    if (athenaCommonFlags.EvtMax>0):
        nSignalEvts=float(athenaCommonFlags.EvtMax.get_Value())
    digilog.info( "#######################################################################" )
    digilog.info( "#######################################################################" )
    digilog.info( "##################  Testing Pile-Up Configuration  ####################" )
    digilog.info( "#######################################################################" )
    digilog.info( "#######################################################################" )
    digilog.info( "Calculating the required number of background events for a %s event job", nSignalEvts )
    refreshrate=1.0/150.0
    nbunches = 1 + digitizationFlags.finalBunchCrossing.get_Value() - digitizationFlags.initialBunchCrossing.get_Value()
    digilog.info( "#######################################################################\n")
    digilog.info( "-----------------------------------------------------------------------" )
    digilog.info( "Required number of background events (no. files assuming 2000 evts per file)" )
    digilog.info( "-----------------------------------------------------------------------" )
    if digitizationFlags.numberOfLowPtMinBias.get_Value()>0.0:
        refreshrate=1.0/150.0 #current default for this collection - FIXME this will be an under-estimate as in-time background is always refreshed
        if(digitizationFlags.doXingByXingPileUp()):
            refreshrate=1.0
        nBkg = pileUpCalc(nSignalEvts, refreshrate, digitizationFlags.numberOfLowPtMinBias.get_Value(),nbunches)
        digilog.info( "\t Low Pt Minimum bias events = %s\t(%s)", int(math.ceil(nBkg)), int(math.ceil(nBkg/5000.0)) )
    if digitizationFlags.numberOfHighPtMinBias.get_Value()>0.0:
        refreshrate=1.0 #current default for this collection
        if(digitizationFlags.doXingByXingPileUp()):
            refreshrate=1.0
        nBkg = pileUpCalc(nSignalEvts, refreshrate, digitizationFlags.numberOfHighPtMinBias.get_Value(),nbunches)
        digilog.info( "\t High Pt Minimum bias events = %s\t(%s)", int(math.ceil(nBkg)), int(math.ceil(nBkg/2000.0)) )
    if digitizationFlags.numberOfCavern.get_Value()>0.0:
        refreshrate=1.0/150.0 #current default for this collection
        if(digitizationFlags.doXingByXingPileUp()):
            refreshrate=1.0
        nBkg = pileUpCalc(nSignalEvts, refreshrate, digitizationFlags.numberOfCavern.get_Value(),nbunches)
        digilog.info( "\t               Cavern Background events = %s\t(%s)", int(math.ceil(nBkg)), int(math.ceil(nBkg/5000.0)) )
    if digitizationFlags.numberOfBeamHalo.get_Value()>0.0:
        refreshrate=1.0 #current default for this collection
        nBkg = pileUpCalc(nSignalEvts, refreshrate, digitizationFlags.numberOfBeamHalo.get_Value(),nbunches)
        digilog.info( "\t                       Beam Halo events = %s\t(%s)", int(math.ceil(nBkg)), int(math.ceil(nBkg/5000.0)) )
    if digitizationFlags.numberOfBeamGas.get_Value()>0.0:
        refreshrate=1.0 #current default for this collection
        nBkg = pileUpCalc(nSignalEvts, refreshrate, digitizationFlags.numberOfBeamGas.get_Value(),nbunches)
        digilog.info( "\t                        Beam Gas events = %s\t(%s)", int(math.ceil(nBkg)), int(math.ceil(nBkg/5000.0)) )
    digilog.info( "-----------------------------------------------------------------------" )
    digilog.info( "\n#######################################################################" )
    sys.exit()


#--------------------------------------------------------------
# Get the flags
#--------------------------------------------------------------
if hasattr(runArgs,"digiSeedOffset1"):
    digitizationFlags.rndmSeedOffset1=int(runArgs.digiSeedOffset1)
else:
    digilog.warning( 'digiSeedOffset1 not set' )
    digitizationFlags.rndmSeedOffset1=1

if hasattr(runArgs,"digiSeedOffset2"):
    digitizationFlags.rndmSeedOffset2=int(runArgs.digiSeedOffset2)
else:
    digilog.warning( 'digiSeedOffset2 not set' )
    digitizationFlags.rndmSeedOffset2=2

if hasattr(runArgs,"samplingFractionDbTag"): #FIXME change this to PhysicsList?
    digitizationFlags.physicsList=runArgs.samplingFractionDbTag

if hasattr(runArgs,"digiRndmSvc"):
    digitizationFlags.rndmSvc=runArgs.digiRndmSvc

if hasattr(runArgs,"conditionsTag"):
    if(runArgs.conditionsTag!='NONE'):
        digitizationFlags.IOVDbGlobalTag = runArgs.conditionsTag

if hasattr(runArgs,"PileUpPresampling"):
    digilog.info("Doing pile-up presampling")
    digitizationFlags.PileUpPresampling = runArgs.PileUpPresampling

#--------------------------------------------------------------
# Pileup configuration
#--------------------------------------------------------------
from SimuJobTransforms.SimTransformUtils import makeBkgInputCol,getInputColOffset
def HasInputFiles(runArgs, key):
    if hasattr(runArgs, key):
        cmd='runArgs.%s' % key
        if eval(cmd):
            return True
    return False


## Low Pt minbias set-up
if HasInputFiles(runArgs, "inputLowPtMinbiasHitsFile"):
    digitizationFlags.LowPtMinBiasInputCols = makeBkgInputCol(runArgs.inputLowPtMinbiasHitsFile,
                                                              digitizationFlags.numberOfLowPtMinBias.get_Value(), True, digilog)
if digitizationFlags.LowPtMinBiasInputCols.statusOn:
    digitizationFlags.doLowPtMinBias = True
else:
    digitizationFlags.doLowPtMinBias = False

## High Pt minbias set-up
if HasInputFiles(runArgs, "inputHighPtMinbiasHitsFile"):
    if(digitizationFlags.HighPtMinBiasInputColOffset.get_Value()<0):
        #Calculate a pseudo random offset into the collection from the jobNumber
        digitizationFlags.HighPtMinBiasInputColOffset = getInputColOffset(runArgs.inputHighPtMinbiasHitsFile, runArgs.jobNumber, digilog)
    digitizationFlags.HighPtMinBiasInputCols = makeBkgInputCol(runArgs.inputHighPtMinbiasHitsFile, 
                                                               digitizationFlags.numberOfHighPtMinBias.get_Value(), True, digilog,
                                                               digitizationFlags.HighPtMinBiasInputColOffset.get_Value())
if digitizationFlags.HighPtMinBiasInputCols.statusOn:
    digitizationFlags.doHighPtMinBias = True
else:
    digitizationFlags.doHighPtMinBias = False

## Cavern Background set-up
if HasInputFiles(runArgs, "inputCavernHitsFile"):
    digitizationFlags.cavernInputCols = makeBkgInputCol(runArgs.inputCavernHitsFile,
                                                        digitizationFlags.numberOfCavern.get_Value(), (not digitizationFlags.cavernIgnoresBeamInt.get_Value()), digilog)
if digitizationFlags.cavernInputCols.statusOn:
    digitizationFlags.doCavern = True
else:
    digitizationFlags.doCavern = False

## Beam Halo set-up
if HasInputFiles(runArgs, "inputBeamHaloHitsFile"):
    digitizationFlags.beamHaloInputCols = makeBkgInputCol(runArgs.inputBeamHaloHitsFile,
                                                          digitizationFlags.numberOfBeamHalo.get_Value(), True, digilog)
if digitizationFlags.beamHaloInputCols.statusOn:
    digitizationFlags.doBeamHalo = True
else:
    digitizationFlags.doBeamHalo = False

## Beam Gas set-up
if HasInputFiles(runArgs, "inputBeamGasHitsFile"):
    digitizationFlags.beamGasInputCols = makeBkgInputCol(runArgs.inputBeamGasHitsFile,
                                                         digitizationFlags.numberOfBeamGas.get_Value(), True, digilog)
if digitizationFlags.beamGasInputCols.statusOn:
    digitizationFlags.doBeamGas = True
else:
    digitizationFlags.doBeamGas = False

#--------------------------------------------------------------
# Other configuration: LVL1, turn off sub detectors, calo noise
#--------------------------------------------------------------

if hasattr(runArgs, "detectors"):
    if not 'DetFlags' in dir():
        from Digitization.DigitizationFlagsHelpers import setupDigitizationLegacyDetectorFlags
        DetFlags = setupDigitizationLegacyDetectorFlags(runArgs.detectors)

if hasattr(runArgs,"doAllNoise"):
    if runArgs.doAllNoise!="NONE":
        digilog.info('doAllNoise = %s: Overriding doInDetNoise, doCaloNoise and doMuonNoise', runArgs.doAllNoise)
        if runArgs.doAllNoise=="True":
            digitizationFlags.doInDetNoise=True
            digitizationFlags.doCaloNoise=True
            digitizationFlags.doMuonNoise=True
        else:
            digitizationFlags.doInDetNoise=False
            digitizationFlags.doCaloNoise=False
            digitizationFlags.doMuonNoise=False

if hasattr(runArgs,"LucidOn") or hasattr(runArgs,"ALFAOn") or hasattr(runArgs,"ZDCOn") or hasattr(runArgs,"AFPOn") or hasattr(runArgs,"FwdRegionOn") or hasattr(runArgs,"HGTDOn"):
    if not 'DetFlags' in dir():
        #if you configure one detflag, you're responsible for configuring them all!
        from AthenaCommon.DetFlags import DetFlags
        DetFlags.all_setOn()
        DetFlags.ALFA_setOff() #Default for now
        DetFlags.ZDC_setOff() #Default for now
        DetFlags.AFP_setOff() #Default for now
        DetFlags.FwdRegion_setOff() #Default for now
        checkHGTDOff = getattr(DetFlags, 'HGTD_setOff', None)
        if checkHGTDOff is not None:
            checkHGTDOff() #Default for now


    if hasattr(runArgs,"LucidOn"):
        if not runArgs.LucidOn:
            DetFlags.Lucid_setOff()
    if hasattr(runArgs,"ALFAOn"):
        if runArgs.ALFAOn:
            DetFlags.ALFA_setOn()
    if hasattr(runArgs,"ZDCOn"):
        if runArgs.ZDCOn:
            DetFlags.ZDC_setOn()
    if hasattr(runArgs, "AFPOn"):
        if runArgs.AFPOn:
            DetFlags.AFP_setOn()
    if hasattr(runArgs, "FwdRegionOn"):
        if runArgs.FwdRegionOn:
            DetFlags.FwdRegion_setOn()
    if hasattr(runArgs, "HGTDOn"):
        if runArgs.HGTDOn:
            checkHGTDOn = getattr(DetFlags, 'HGTD_setOn', None)
            if checkHGTDOn is not None:
                checkHGTDOn()
            else:
                digilog.warning('The HGTD DetFlag is not supported in this release')


def NoTriggerConfig(runArgs):
    if hasattr(runArgs,"triggerConfig"):
        if runArgs.triggerConfig!="NONE":
            return False
    return True

if NoTriggerConfig(runArgs):
    if not 'DetFlags' in dir():
        #if you configure one detflag, you're responsible for configuring them all!
        from AthenaCommon.DetFlags import DetFlags
        DetFlags.all_setOn()
        DetFlags.ALFA_setOff() #Default for now
        DetFlags.ZDC_setOff() #Default for now
    DetFlags.LVL1_setOff()
    DetFlags.Print()

## Output RDO File
if hasattr(runArgs,"outputRDOFile") or hasattr(runArgs,"outputRDO_FILTFile"):
    if hasattr(runArgs,"outputRDOFile"):
        if hasattr(runArgs,"outputRDO_FILTFile"):
            digilog.fatal("Both outputRDOFile and outputRDO_FILTFile specified - this configuration should not be used!")
            raise SystemError
        athenaCommonFlags.PoolRDOOutput.set_Value_and_Lock( runArgs.outputRDOFile )
    elif hasattr(runArgs,"outputRDO_FILTFile"):
        athenaCommonFlags.PoolRDOOutput.set_Value_and_Lock( runArgs.outputRDO_FILTFile )
    if hasattr(runArgs, "AddCaloDigi") and runArgs.AddCaloDigi:
        digilog.info("Will write out all LArDigitContainers and TileDigitsContainers to RDO file.")
        digitizationFlags.experimentalDigi+=["AddCaloDigi"]
else:
    digilog.info("no output file (outputRDOFile or outputRDO_FILTFile) specified - switching off output StreamRDO")
    if not 'DetFlags' in dir():
        #if you configure one detflag, you're responsible for configuring them all!
        from AthenaCommon.DetFlags import DetFlags
        DetFlags.all_setOn()
        DetFlags.ALFA_setOff() #Default for now
        DetFlags.ZDC_setOff() #Default for now
    DetFlags.writeRDOPool.all_setOff()

#--------------------------------------------------------------
# Go for it
#--------------------------------------------------------------
if hasattr(runArgs,"DataRunNumber"):
    if runArgs.DataRunNumber>0:
        digilog.info( 'Overriding run number to be: %s ', runArgs.DataRunNumber )
        digitizationFlags.dataRunNumber=runArgs.DataRunNumber

from AthenaCommon.AlgSequence import AlgSequence
topSeq = AlgSequence()

## Set Overall per-Algorithm time-limit on the AlgSequence
topSeq.TimeOut = 43200 * Units.s


if jobproperties.ConcurrencyFlags.NumThreads() == 0:
    try:
        timingOutput = "HITStoRDO_timings"
        if digitizationFlags.PileUpPresampling and 'LegacyOverlay' not in digitizationFlags.experimentalDigi():
            from OverlayCommonAlgs.OverlayFlags import overlayFlags
            timingOutput = overlayFlags.bkgPrefix() + timingOutput

        from RecAlgs.RecAlgsConf import TimingAlg
        topSeq += TimingAlg("DigiTimerBegin", TimingObjOutputName = timingOutput)
    except:
        digilog.warning('Could not add TimingAlg, no timing info will be written out.')
else:
    digilog.info("MT mode, not scheduling TimingAlg")

include ("Digitization/Digitization.py")

from PyUtils import AMITagHelper
AMITagHelper.SetAMITag(runArgs=runArgs)

#==========================================================
# Use ZLIB for compression of all temporary outputs
#==========================================================
from AthenaCommon.AppMgr import ServiceMgr as svcMgr
import AthenaPoolCnvSvc.AthenaPool  # noqa: F401
from AthenaPoolCnvSvc import PoolAttributeHelper as pah
Out = athenaCommonFlags.PoolRDOOutput()
if hasattr(runArgs, "outputRDOFile") and (runArgs.outputRDOFile.endswith('_000') or runArgs.outputRDOFile.startswith('tmp.')): # noqa: F821
    svcMgr.AthenaPoolCnvSvc.PoolAttributes += [ pah.setFileCompAlg( Out, 1 ) ]
    svcMgr.AthenaPoolCnvSvc.PoolAttributes += [ pah.setFileCompLvl( Out, 1 ) ]

# Set AutoFlush to 1 as per ATLASSIM-4274
# This helps with the overall HITtoRDO memory footprint
svcMgr.AthenaPoolCnvSvc.PoolAttributes += [ pah.setTreeAutoFlush( Out, "CollectionTree", 1 ) ]
svcMgr.AthenaPoolCnvSvc.PoolAttributes += [ pah.setTreeAutoFlush( Out, "POOLContainer", 1 ) ]
svcMgr.AthenaPoolCnvSvc.PoolAttributes += [ pah.setTreeAutoFlush( Out, "POOLContainerForm", 1 ) ]

## Post-include
if hasattr(runArgs,"postInclude"):
    for fragment in runArgs.postInclude:
        include(fragment)

## Post-exec
if hasattr(runArgs,"postExec"):
    digilog.info("transform post-exec")
    for cmd in runArgs.postExec:
        digilog.info(cmd)
        exec(cmd)
