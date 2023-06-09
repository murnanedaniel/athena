#still needs a lot of cleaning ... 

from RecExConfig.RecFlags import rec
from RecExConfig.RecAlgsFlags import recAlgs
from AthenaCommon.AthenaCommonFlags import athenaCommonFlags as acf

#if not acf.EvtMax.is_locked():
acf.EvtMax=15
if not ('OutputLevel' in dir()):
    rec.OutputLevel=INFO
#scan for RTT files (only if dsName and fileRange set)
include("TriggerTest/TrigScanFiles.py")
###############################

doTrigger=True
TriggerModernConfig=True
rec.doWriteAOD=True
rec.doWriteESD=False
rec.doWriteTAG=False
rec.doAOD=True
rec.doJetMissingETTag=True
recAlgs.doMissingET=False
rec.doMuonCombined=False
rec.doTau=False
rec.doBTagging=False
#rec.doESD.set_Value_and_Lock(False) 
rec.doESD=True
doTAG=False
rec.doCBNT=False 
#rec.doTruth=True

#-----------------------------------------------------------
include("RecExCond/RecExCommon_flags.py")
#-----------------------------------------------------------

# set up trigger monitoring                                                                                                                                                        
if not ('RunningRTT' in dir()):
    TriggerFlags.enableMonitoring = [ 'Validation', 'Time' , 'Log' ]
else:
    TriggerFlags.enableMonitoring = [ 'Validation', 'Time' ]

#------------ This is for ATN/RTT tests only ---------
#TriggerFlags.triggerMenuSetup = 'default'
TriggerFlags.readHLTconfigFromXML=False
TriggerFlags.readLVL1configFromXML=False
#if  ('menu' in dir()):
   # TiriggerFlags.triggerMenuSetup=menu 

#include("TrigEgammaValidation/TrigEgammaValidation_RTT_Chains.py") # uncomment to get used-defined trigger list

#egammatrigChainlist=electronChains() #Retrieve trigger menu to run from TrigEgammaValidation_RTT_Chains.py
#TriggerFlags.triggerMenuSetup="MC_pp_v5"
TriggerFlags.Slices_all_setOff()
TriggerFlags.EgammaSlice.setAll()
#TriggerFlags.EgammaSlice.signatures = egammatrigChainlist


TriggerFlags.doHLT=True
TriggerFlags.L1PrescaleSet  = '' 
TriggerFlags.HLTPrescaleSet = '' 
TriggerFlags.AODEDMSet="AODFULL"

#-------------end of flag for tests-------------------

#------------ This is a temporary fix ---------------
#from RecExConfig.RecConfFlags import recConfFlags
#recConfFlags.AllowIgnoreConfigError=False
#athenaCommonFlags.AllowIgnoreConfigError=False
#-------------end of temporary fix-------------------

#from ParticleBuilderOptions.AODFlags import AODFlags 
#AODFlags.FastSimulation=False 
# see comments in https://savannah.cern.ch/bugs/?83735
#AODFlags.MuonTrackSlimmer=False
#print AODFlags.Print()

#-----------------------------------------------------------
include("RecExCommon/RecExCommon_topOptions.py")
#-----------------------------------------------------------

# abort when there is an unchecked status code
StatusCodeSvc.AbortOnError=False

include("TriggerTest/TriggerTestCommon.py")
