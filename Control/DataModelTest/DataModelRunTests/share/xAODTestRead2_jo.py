#
# $Id$
#
# File: share/xAODTestRead2_jo.py
# Author: snyder@bnl.gov
# Date: May 2014
# Purpose: Test reading xAOD objects.
#          Read output of xAODTestRead_jo.py.
#

## basic job configuration (for generator)
import AthenaCommon.AtlasUnixStandardJob

## get a handle to the default top-level algorithm sequence
from AthenaCommon.AlgSequence import AlgSequence
topSequence = AlgSequence()

## get a handle to the ServiceManager
from AthenaCommon.AppMgr import ServiceMgr as svcMgr

## get a handle to the ApplicationManager
from AthenaCommon.AppMgr import theApp

#--------------------------------------------------------------
# Load POOL support
#--------------------------------------------------------------
import AthenaPoolCnvSvc.WriteAthenaPool
import AthenaPoolCnvSvc.ReadAthenaPool

#--------------------------------------------------------------
# Define input
#--------------------------------------------------------------
svcMgr.EventSelector.InputCollections        = [ "xaoddata2.root" ]

#--------------------------------------------------------------
# Event related parameters
#--------------------------------------------------------------
theApp.EvtMax = 20

#--------------------------------------------------------------
# Application:
#--------------------------------------------------------------

from DataModelTestDataRead.DataModelTestDataReadConf import \
     DMTest__xAODTestReadCVec, \
     DMTest__xAODTestRead
topSequence += DMTest__xAODTestReadCVec ("xAODTestReadCVec")
topSequence += DMTest__xAODTestRead ("xAODTestRead")
topSequence += DMTest__xAODTestReadCVec ("xAODTestReadCVec_copy",
                                         CVecKey = "copy_cvec")
topSequence += DMTest__xAODTestRead ("xAODTestRead_copy",
                                     ReadPrefix = "copy_")
topSequence += DMTest__xAODTestReadCVec ("xAODTestReadCVec_scopy",
                                         CVecKey = "scopy_cvec")
topSequence += DMTest__xAODTestRead ("xAODTestRead_scopy",
                                     ReadPrefix = "scopy_")

# Note: can't autoload these.
import ROOT
import cppyy
cppyy.loadDictionary("libDataModelTestDataCommonDict")
cppyy.loadDictionary("libDataModelTestDataReadDict")
ROOT.DMTest.B
ROOT.DMTest.setConverterLibrary ('libDataModelTestDataReadCnvPoolCnv.so')


#--------------------------------------------------------------
# Set output level threshold (2=DEBUG, 3=INFO, 4=WARNING, 5=ERROR, 6=FATAL )
#--------------------------------------------------------------
svcMgr.MessageSvc.OutputLevel = 3
svcMgr.MessageSvc.debugLimit  = 100000
svcMgr.ClassIDSvc.OutputLevel = 3

# No stats printout
ChronoStatSvc = Service( "ChronoStatSvc" )
ChronoStatSvc.ChronoPrintOutTable = FALSE
ChronoStatSvc.PrintUserTime       = FALSE
ChronoStatSvc.StatPrintOutTable   = FALSE

#svcMgr.ExceptionSvc.Catch = "None"

# Explicitly specify the output file catalog
# to avoid races when running tests in parallel.
PoolSvc = Service( "PoolSvc" )
PoolSvc.WriteCatalog = "file:xAODTestRead2_catalog.xml"
