###############################################################
#
# Job options file
#
## @file DoubleEventSelectorOverlayTest.py
##
## @brief For Athena POOL test: read an RDO and a HITS file, output a RDO file
##
## @author Miha Muskinja <miha.muskinja@cern.ch>
#
#==============================================================

# basic job configuration
import AthenaCommon.AtlasUnixStandardJob

# get a handle to the default top-level algorithm sequence
from AthenaCommon.AlgSequence import AlgSequence
topSequence = AlgSequence()

# get a handle to the ServiceManager
from AthenaCommon.AppMgr import ServiceMgr as svcMgr

# get a handle to the ApplicationManager
from AthenaCommon.AppMgr import theApp

#--------------------------------------------------------------
# Load POOL support for DoubleEventSelector
#--------------------------------------------------------------
import AthenaPoolCnvSvc.ReadAthenaPoolDouble

#--------------------------------------------------------------
# Set flags and load det descr
#--------------------------------------------------------------
from AthenaCommon.GlobalFlags  import globalflags
from RecExConfig.RecFlags      import rec
from OverlayCommonAlgs.OverlayFlags import overlayFlags

overlayFlags.isOverlayMT.set_Value_and_Lock(True)

# For general flags
rec.doAOD       = False
rec.doTrigger   = False
rec.doWriteTAG  = False
DetDescrVersion = "ATLAS-R2-2016-01-00-01"

#--------------------------------------------------------------
# Input options
#--------------------------------------------------------------
svcMgr.DoubleEventSelector.PrimaryInputCollections = [ "/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/OverlayMonitoringRTT/PileupPremixing/22.0/v1/RDO.merged-pileup-MT.100events.pool.root" ]
svcMgr.DoubleEventSelector.SecondaryaryInputCollections = [ "/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/OverlayMonitoringRTT/valid1.410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.simul.HITS.e4993_s3091/HITS.10504490._000765.pool.root.1" ]
svcMgr.DoubleEventSelector.OutputLevel = DEBUG

#--------------------------------------------------------------
# Remapping Service
#--------------------------------------------------------------
from SGComps import AddressRemappingSvc
AddressRemappingSvc.addInputRename("EventInfo","McEventInfo" ,"Sig_McEventInfo")
AddressRemappingSvc.addInputRename("McEventCollection","TruthEvent" ,"Sig_TruthEvent")
AddressRemappingSvc.addInputRename("RecoTimingObj","EVNTtoHITS_timings" ,"Sig_EVNTtoHITS_timings")
svcMgr.AddressRemappingSvc.OutputLevel = DEBUG

#--------------------------------------------------------------
# Event related parameters
#--------------------------------------------------------------
theApp.EvtMax = 10

#--------------------------------------------------------------
# Algorithms
#--------------------------------------------------------------
from AthenaCommon import CfgGetter
topSequence += CfgGetter.getAlgorithm("CopyMcEventCollection")
topSequence += CfgGetter.getAlgorithm("CopyTimings")

#--------------------------------------------------------------
# Athena EventLoop Manager
#--------------------------------------------------------------
AthenaEventLoopMgr = Service( "AthenaEventLoopMgr" )
AthenaEventLoopMgr.UseSecondaryEventNumber = True
AthenaEventLoopMgr.OutputLevel = INFO

#--------------------------------------------------------------
# DEBUG messaging
#--------------------------------------------------------------
svcMgr.ProxyProviderSvc.OutputLevel = DEBUG
svcMgr.AthenaPoolAddressProviderSvcPrimary.OutputLevel = DEBUG
svcMgr.AthenaPoolAddressProviderSvcSecondary.OutputLevel = DEBUG
svcMgr.DoubleEventSelector.OutputLevel = DEBUG

#--------------------------------------------------------------
# Output options
#--------------------------------------------------------------
from AthenaPoolCnvSvc.WriteAthenaPool import AthenaPoolOutputStream
# TODO: noTag=True to avoid warning, needs EventInfo overlay implemented
Stream1 = AthenaPoolOutputStream( "Stream1", asAlg=True, noTag=True )
Stream1.OutputLevel = INFO

Stream1.OutputFile  = "OutputRDO.root"
# List of DO's to write out
Stream1.ItemList =  []
Stream1.ItemList += ["McEventCollection#TruthEvent"]
Stream1.ItemList += ["RecoTimingObj#EVNTtoHITS_timings"]

#--------------------------------------------------------------
# Set output level threshold (2=DEBUG, 3=INFO, 4=WARNING, 5=ERROR, 6=FATAL )
#--------------------------------------------------------------
svcMgr.MessageSvc = Service( "MessageSvc" )
svcMgr.MessageSvc.OutputLevel = WARNING
svcMgr.MessageSvc.debugLimit  = 100000

# No stats printout
include( "AthenaPoolTest/NoStats_jobOptions.py" )

#==============================================================
#
# End of job options file
#
###############################################################
