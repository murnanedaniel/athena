#
# Job Options for reading LAr HV Pathologies from COOL using local SQLite database
#
from time import strptime,time
from calendar import timegm

if "GlobalTag" not in dir():
    GlobalTag = 'CONDBR2-BLKPA-2017-03'

if "inputsqlite" not in dir():
   dbConn="COOLOFL_LAR/CONDBR2"
else:   
   dbConn="sqlite://;schema="+inputsqlite+";dbname=CONDBR2"   

if "foldertag" not in dir():
   foldertag="LARHVPathologiesOflPathologies-RUN2-UPD1-00"

if "mode" not in dir():
   mode=0

if "date" not in dir():
    date="2015-05-29:12:00:00"

if "RunNumber" not in dir():
   RunNumber=999999

if "LumiBlock" not in dir():
   LumiBlock=1

if "TimeStamp" not in dir():
   try:
       ts=strptime(date+'/UTC','%Y-%m-%d:%H:%M:%S/%Z')
       TimeStamp=int(timegm(ts))*1000000000
   except ValueError:
       print("ERROR in time specification, use e.g. 2007-05-25:14:01:00")

## basic job configuration
import AthenaCommon.AtlasUnixGeneratorJob

## get a handle to the default top-level algorithm sequence
from AthenaCommon.AlgSequence import AlgSequence
topSequence = AlgSequence()

## get a handle to the ServiceManager
#from AthenaCommon.AppMgr import ServiceMgr as svcMgr

## get a handle to the ApplicationManager
#from AthenaCommon.AppMgr import theApp

from AthenaCommon.DetFlags import DetFlags
DetFlags.all_setOff()
DetFlags.LAr_setOn()
DetFlags.Calo_setOn()
DetFlags.digitize.all_setOff()

from AthenaCommon.GlobalFlags import globalflags
globalflags.DetGeo = 'atlas'
globalflags.DataSource = 'data'
globalflags.InputFormat = 'bytestream'
globalflags.DatabaseInstance="CONDBR2"

from AthenaCommon.GlobalFlags import jobproperties
jobproperties.Global.DetDescrVersion='ATLAS-R2-2015-04-00-00'

from AtlasGeoModel import SetGeometryVersion
from AtlasGeoModel import GeoModelInit
from AtlasGeoModel import SetupRecoGeometry

svcMgr.IOVDbSvc.GlobalTag = GlobalTag

include( "CaloDetMgrDetDescrCnv/CaloDetMgrDetDescrCnv_joboptions.py")
include( "CaloIdCnv/CaloIdCnv_joboptions.py" )
include( "LArDetDescr/LArDetDescr_joboptions.py" )
include( "LArConditionsCommon/LArIdMap_comm_jobOptions.py" )

conddb.addFolder("LAR_OFL","/LAR/IdentifierOfl/HVLineToElectrodeMap",className="AthenaAttributeList")
from AthenaCommon.AlgSequence import AthSequencer
condseq = AthSequencer("AthCondSeq")

from LArRecUtils.LArRecUtilsConf import LArHVIdMappingAlg
hvmapalg = LArHVIdMappingAlg(ReadKey="/LAR/IdentifierOfl/HVLineToElectrodeMap",WriteKey="LArHVIdMap")
condseq += hvmapalg
#--------------------------------------------------------------
# Access to IOVSvc, IOVDbSvc and CondDBMySQLCnvSvc
#--------------------------------------------------------------
import IOVDbSvc.IOVDb

from LArCondUtils.LArCondUtilsConf import LArHVPathologyDbAlg
LArHVPathologyDbAlg = LArHVPathologyDbAlg()
LArHVPathologyDbAlg.Folder = "/LAR/HVPathologiesOfl/Pathologies"
LArHVPathologyDbAlg.Mode=mode
if 'OutFile' in dir():
   LArHVPathologyDbAlg.OutFile=OutFile
topSequence += LArHVPathologyDbAlg

if 'inputsqlite' in dir():   
   svcMgr.IOVDbSvc.dbConnection  = "sqlite://;schema="+inputsqlite+";dbname=CONDBR2"
   svcMgr.IOVDbSvc.Folders += [ "/LAR/HVPathologiesOfl/Pathologies<tag>"+foldertag+"</tag>" ]
else:  
   conddb.addFolder("LAR_OFL","/LAR/HVPathologiesOfl/Pathologies<tag>"+foldertag+"</tag>")
    
from AthenaCommon.AlgSequence import AthSequencer
condseq = AthSequencer("AthCondSeq")

from LArRecUtils.LArRecUtilsConf import LArHVIdMappingAlg
hvmapalg = LArHVIdMappingAlg(ReadKey="/LAR/IdentifierOfl/HVLineToElectrodeMap",WriteKey="LArHVIdMap")
condseq += hvmapalg

svcMgr.MessageSvc.OutputLevel = 4
svcMgr.MessageSvc.debugLimit  = 99999999
svcMgr.MessageSvc.infoLimit   = 99999999

LArHVPathologyDbAlg.OutputLevel     = 3
svcMgr.IOVDbSvc.OutputLevel         = 3

#--------------------------------------------------------------
# The following turns off the ability to set EventSelector parameters
# below
#ApplicationMgr.EvtSel = "NONE"
#--------------------------------------------------------------
# Set the following when reading back to adjust the run/event numbers
svcMgr.EventSelector.RunNumber         = RunNumber
svcMgr.EventSelector.EventsPerRun      = 1
svcMgr.EventSelector.FirstEvent        = 1
svcMgr.EventSelector.EventsPerLB       = 1
svcMgr.EventSelector.FirstLB           = LumiBlock
svcMgr.EventSelector.InitialTimeStamp  = int(TimeStamp/1e9)
svcMgr.EventSelector.TimeStampInterval = 5
theApp.EvtMax                   = 1
