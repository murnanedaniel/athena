# Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
#########################################################################
#
## @file EventInfoOverlayLegacyTest.py
##
## @brief xAOD::EventInfo overlay event loop test:
## read an RDO and a HITS file, output a RDO file
## Test legacy HITS with the old EventInfo
##
## @author Tadej Novak <tadej.novak@cern.ch>
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
# Set flags and load det descr
#--------------------------------------------------------------
from AthenaCommon.GlobalFlags  import globalflags
from RecExConfig.RecFlags      import rec
from OverlayCommonAlgs.OverlayFlags import overlayFlags

# Set that we are running MC+MC overlay in MT mode
globalflags.isOverlay.set_Value_and_Lock(True)
overlayFlags.isDataOverlay.set_Value_and_Lock(False)
overlayFlags.isOverlayMT.set_Value_and_Lock(True)

#--------------------------------------------------------------
# Load POOL support for DoubleEventSelector
#--------------------------------------------------------------
import AthenaPoolCnvSvc.ReadAthenaPoolDouble

# For general flags
rec.doAOD       = False
rec.doTrigger   = False
rec.doWriteTAG  = False
DetDescrVersion = "ATLAS-R2-2016-01-00-01"

# the correct tag should be specified
from IOVDbSvc.CondDB import conddb
conddb.setGlobalTag("OFLCOND-SDR-BS7T-04-00")

#--------------------------------------------------------------
# Input options
#--------------------------------------------------------------
svcMgr.EventSelector.InputCollections = [ "/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/OverlayTests/PresampledPileUp/22.0/Run2/large/mc20_13TeV.900149.PG_single_nu_Pt50.digit.RDO.e8307_s3482_s3136_d1715/RDO.26811908._031801.pool.root.1" ]
svcMgr.SecondaryEventSelector.InputCollections  = [ "/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/OverlayTests/Special/TestCase_xAODEventInfo.root" ]

#--------------------------------------------------------------
# Remapping Service
#--------------------------------------------------------------
from SGComps import AddressRemappingSvc
AddressRemappingSvc.addInputRename("xAOD::EventInfo", "EventInfo", "Sig_EventInfo")
AddressRemappingSvc.addInputRename("xAOD::EventAuxInfo", "EventInfoAux.", "Sig_EventInfoAux.")

#--------------------------------------------------------------
# Event related parameters
#--------------------------------------------------------------
theApp.EvtMax = 100

#--------------------------------------------------------------
# Algorithms
#--------------------------------------------------------------
# Beam spot conditions
from AthenaCommon.AlgSequence import AthSequencer
condSeq = AthSequencer("AthCondSeq")
from IOVDbSvc.CondDB import conddb
conddb.addFolderSplitOnline("INDET", "/Indet/Onl/Beampos", "/Indet/Beampos", className="AthenaAttributeList")
from BeamSpotConditions.BeamSpotConditionsConf import BeamSpotCondAlg
condSeq += BeamSpotCondAlg("BeamSpotCondAlg")

# Run the overlay
from AthenaCommon import CfgGetter
EventInfoOverlay = CfgGetter.getAlgorithm("EventInfoOverlay")
EventInfoOverlay.OutputLevel = DEBUG
topSequence += EventInfoOverlay

#--------------------------------------------------------------
# EventLoop
#--------------------------------------------------------------
from AthenaCommon.ConcurrencyFlags import jobproperties as jp
nThreads = jp.ConcurrencyFlags.NumThreads()
if nThreads > 0:
    from AthenaServices.AthenaServicesConf import AthenaHiveEventLoopMgr
    EventLoop = AthenaHiveEventLoopMgr()
else:
    from AthenaServices.AthenaServicesConf import AthenaEventLoopMgr
    EventLoop = AthenaEventLoopMgr()
EventLoop.RequireInputAttributeList = True
EventLoop.UseSecondaryEventNumber = True
svcMgr += EventLoop

#--------------------------------------------------------------
# Output options
#--------------------------------------------------------------
from AthenaPoolCnvSvc.WriteAthenaPool import AthenaPoolOutputStream
Stream = AthenaPoolOutputStream( "Stream" )
Stream.OutputFile  = "OutputRDO.root"
Stream.ItemList =  [ 'xAOD::EventInfo#EventInfo', 'xAOD::EventAuxInfo#EventInfoAux.' ]

#--------------------------------------------------------------
# Set output level threshold (2=DEBUG, 3=INFO, 4=WARNING, 5=ERROR, 6=FATAL )
#--------------------------------------------------------------
svcMgr.MessageSvc = Service( "MessageSvc" )
svcMgr.MessageSvc.OutputLevel = WARNING
svcMgr.MessageSvc.debugLimit  = 100000

# No stats printout
include( "AthenaPoolTest/NoStats_jobOptions.py" )
