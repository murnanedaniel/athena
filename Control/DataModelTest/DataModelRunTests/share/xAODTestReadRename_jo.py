#
# $Id$
#
# File: DataModelRunTests/share/xAODTestReadRename_jo.py
# Author: snyder@bnl.gov
# Date: Aug 2016
# Purpose: Test reading xAOD objects data with renaming on input.
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
import AthenaPoolCnvSvc.ReadAthenaPool

include ('DataModelRunTests/loadReadDicts.py')

#--------------------------------------------------------------
# Define input
#--------------------------------------------------------------
svcMgr.EventSelector.InputCollections        = [ "xaoddata.root" ]

from SGComps import AddressRemappingSvc
AddressRemappingSvc.addInputRename ('DMTest::CVec', 'cvec', 'cvec_renamed')
AddressRemappingSvc.addInputRename ('DMTest::CAuxContainer',
                                    'cvecAux.', 'cvec_renamedAux.')

AddressRemappingSvc.addInputRename ('DMTest::CVec', 'cvec.dInt1',
                                    'cvec_renamed.dInt1_renamed')
AddressRemappingSvc.addInputRename ('DMTest::C', 'cinfo.dInt1',
                                    'cinfo.dInt1_renamed')
AddressRemappingSvc.addInputRename ('DMTest::C', 'cinfo.dInt1Base',
                                    'cinfo.dInt1_renamedBase')

#--------------------------------------------------------------
# Event related parameters
#--------------------------------------------------------------
theApp.EvtMax = 20

#--------------------------------------------------------------
# Application:
#--------------------------------------------------------------

from DataModelTestDataCommon.DataModelTestDataCommonConf import \
     DMTest__xAODTestReadDecor
from DataModelTestDataRead.DataModelTestDataReadConf import \
     DMTest__xAODTestReadCVec


topSequence += DMTest__xAODTestReadCVec ('xAODTestReadCVec',
                                         CVecKey = 'cvec_renamed')
topSequence += DMTest__xAODTestReadDecor ('xAODTestReadDecor',
                                          CVecName = 'cvec_renamed',
                                          DecorName = 'dInt1_renamed')


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

# Avoid races when running tests in parallel.
if 'FILECATALOG' not in globals():
    FILECATALOG = 'xAODTestReadRename_catalog.xml'
include ('DataModelRunTests/setCatalog.py')
