
from future import standard_library
standard_library.install_aliases()
import subprocess

###########################################################################
#
# <Marco.Delmastro@cern.ch>
#
# Example jobOptions to compute OFC cali in all 3 gains
#
# Last update: 09/12/2008 <Fabien.Tarrade@cern.ch>
#
###########################################################################

include("LArCalibProcessing/LArCalib_Flags.py")

###########################################################################
#                 Input selection (Autocorr, CaliWave)
###########################################################################

# RunNumber trigger IOV if reading from COOL is eneabled
if not 'RunNumber' in dir():
   RunNumber = '0000000'

if not 'SubDet' in dir():
   SubDet = "Barrel"  

if not 'Partition' in dir():
   if (SubDet=='EndCap'):
      Partition      = "EB-EC*"
   else:
      Partition      = "EB-EMB*"
      
## AutoCorr

if not 'ReadAutoCorrFromCOOL' in dir():
   ReadAutoCorrFromCOOL = True 

if not 'InputAutoCorrPoolDir' in dir():
   InputAutoCorrPoolDir = subprocess.getoutput("pwd")

if not 'InputAutoCorrPoolFileName' in dir():
   InputAutoCorrPoolFileName = "LArAutoCorr.pool.root"

## CaliWave

if not 'ReadCaliWaveFromCOOL' in dir():
   ReadCaliWaveFromCOOL = True

if not 'InputCaliWavePoolDir' in dir():
   InputCaliWavePoolDir = subprocess.getoutput("pwd")

if not 'InputCaliWavePoolFileName' in dir():
   InputCaliWavePoolFileName = "LArCaliWave.pool.root"


## Options   

if not 'StripsXtalkCorr' in dir():
   StripsXtalkCorr = True

if not 'ChannelSelection' in dir():
   # Read all
   ChannelSelection = " "
   ## Example to read only cool for Barrel C : Strips->Back
   #ChannelSelection = "<channelSelection>0,35:66</channelSelection>"

from string import *
def DBConnectionFile(sqlitefile):
   return "sqlite://;schema="+sqlitefile+";dbname=CONDBR2"

###########################################################################
#                             OFC properties
###########################################################################

if not 'Nsamples' in dir():
   Nsamples = 5

if not 'Nphases' in dir():
   Nphases = 50

if not 'Dphases' in dir():
   Dphases = 1
   
if not 'Ndelays' in dir():
   Ndelays = 24
	
if not 'ContainerKey' in dir():
   ContainerKey = "LArCaliWave"
   
if not 'OFCKey' in dir():
   OFCKey = "LArOFC"

if not 'ShapeKey' in dir():
   ShapeKey = "LArShape"
	
if not 'Normalize' in dir():
   Normalize = True

if not 'TimeShift' in dir() :
   TimeShift = False

if not 'TimeShiftByIndex' in dir() :
   TimeShiftByIndex = -1

if not 'FillShape' in dir():
   FillShape = False # Do not fill a LArShapeComplete object for calibration OFC!

if not 'DumpOFC' in dir():
   DumpOFC = False

###########################################################################
#                              OFC output
###########################################################################

if not 'GroupingType' in dir():
   GroupingType = "ExtendedSubDetector"
   
if not 'WriteNtuple' in dir():
   WriteNtuple = LArCalib_Flags.WriteNtuple

if not 'WritePoolFile' in dir():
   WritePoolFile = LArCalib_Flags.WritePoolFile

if not 'WriteIOV' in dir():
   WriteIOV = LArCalib_Flags.WriteIOV

if not 'IOVBegin' in dir():
   IOVBegin = int(RunNumber)

if not 'IOVEnd' in dir():
   IOVEnd = LArCalib_Flags.IOVEnd

if not 'DBConnectionCOOL' in dir():   
   DBConnectionCOOL = "oracle://ATLAS_COOLPROD;schema=ATLAS_COOLOFL_LAR;dbname=CONDBR2;"
   
if not 'OutputOFCRootFileDir' in dir():
   OutputOFCRootFileDir = subprocess.getoutput("pwd")
   
if not 'OutputOFCPoolFileDir' in dir():
   OutputOFCPoolFileDir = subprocess.getoutput("pwd")

if not 'OutputShapePoolFileDir' in dir():
   OutputShapePoolFileDir = subprocess.getoutput("pwd")

OFCFileTag = str(RunNumber)+"_"+Partition.replace("*","")

if (StripsXtalkCorr):
   OFCFileTag += "_StripsXtalkCorr"

OFCFileTag += "_"+str(Nsamples)+"samples"

if (Dphases>1):
   OFCFileTag += "_"+str(Dphases)+"Dphase"

if not 'OutputOFCRootFileName' in dir():
   OutputOFCRootFileName = "LArOFCCali_"+OFCFileTag + ".root"
   
if not 'OutputOFCPoolFileName' in dir():
   OutputOFCPoolFileName = "LArOFCCali_"+OFCFileTag + ".pool.root"

if not 'OutputShapePoolFileName' in dir():
   OutputShapePoolFileName = "LArShapeCali_"+OFCFileTag + ".pool.root"

if not 'LArCalibFolderOutputTag' in dir():
   LArCalibFolderOutputTag = LArCalib_Flags.tagSuffix  

if not 'OutputDB' in dir():
   OutputDB = LArCalib_Flags.OutputDB

if 'OutputSQLiteFile' in dir():
   OutputDB = DBConnectionFile(OutputSQLiteFile)

if not 'OutputObjectSpecOFC' in dir():
   if   ( ContainerKey == "LArCaliWave" ):
      if ( not StripsXtalkCorr ):
         OutputObjectSpecOFC   = "LArOFCComplete#"  +OFCKey  +"#"+ LArCalib_Flags.LArOFCCaliFolder	
         OutputObjectSpecTagOFC    = LArCalibFolderTag(LArCalib_Flags.LArOFCCaliFolder,LArCalibFolderOutputTag)
      else:
         OutputObjectSpecOFC   = "LArOFCComplete#"  +OFCKey  +"#"+ LArCalib_Flags.LArOFCCaliFolderXtlk
         OutputObjectSpecTagOFC    = LArCalibFolderTag(LArCalib_Flags.LArOFCCaliFolderXtlk,LArCalibFolderOutputTag)
			
   elif ( ContainerKey == "LArMasterWave" ):
      if ( not StripsXtalkCorr ):
         OutputObjectSpecOFC   = "LArOFCComplete#"+OFCKey+"#"+ LArCalib_Flags.LArOFCMasterWaveFolder
         OutputObjectSpecTagOFC    = LArCalibFolderTag(LArCalib_Flags.LArOFCMasterWaveFolder,LArCalibFolderOutputTag)
      else:
         OutputObjectSpecOFC   = "LArOFCComplete#"+OFCKey+"#"+ LArCalib_Flags.LArOFCMasterWaveFolderXtlk
         OutputObjectSpecTagOFC    = LArCalibFolderTag(LArCalib_Flags.LArOFCMasterWaveFolderXtlk,LArCalibFolderOutputTag)
			
if not 'OutputObjectSpecShape' in dir():
   if   ( ContainerKey == "LArCaliWave" ):
      if ( not StripsXtalkCorr ):
         OutputObjectSpecShape = "LArShapeComplete#"+ShapeKey+"#"+ LArCalib_Flags.LArShapeCaliWaveFolder		
         OutputObjectSpecTagShape  = LArCalibFolderTag(LArCalib_Flags.LArShapeCaliWaveFolder,LArCalibFolderOutputTag)
      else:
         OutputObjectSpecShape = "LArShapeComplete#"+ShapeKey+"#"+ LArCalib_Flags.LArShapeCaliWaveFolderXtlk
         OutputObjectSpecTagShape  = LArCalibFolderTag(LArCalib_Flags.LArShapeCaliWaveFolderXtlk,LArCalibFolderOutputTag)

   elif ( ContainerKey == "LArMasterWave" ):
         if ( not StripsXtalkCorr ):
            OutputObjectSpecShape = "LArShapeComplete#"+ShapeKey+"#"+LArCalib_Flags.LArShapeMasterWaveFolder
            OutputObjectSpecTagShape  = LArCalibFolderTag(LArCalib_Flags.LArShapeMasterWaveFolder,LArCalibFolderOutputTag)		
         else:
            OutputObjectSpecShape = "LArShapeComplete#"+ShapeKey+"#"+ LArCalib_Flags.LArShapeMasterWaveFolderXtlk
            OutputObjectSpecTagShape  = LArCalibFolderTag(LArCalib_Flags.LArShapeMasterWaveFolderXtlk,LArCalibFolderOutputTag)
   
if ( ReadAutoCorrFromCOOL ):      
   if 'InputAutoCorrSQLiteFile' in dir():
      InputDBConnectionAutoCorr = DBConnectionFile(InputAutoCorrSQLiteFile)
   else:
      InputDBConnectionAutoCorr = DBConnectionCOOL

if ( ReadCaliWaveFromCOOL ):      
   if 'InputCaliWaveSQLiteFile' in dir():
      InputDBConnectionCaliWave = DBConnectionFile(InputCaliWaveSQLiteFile)
   else:
      InputDBConnectionCaliWave = DBConnectionCOOL

## Bad Channel   
   
if not 'ReadBadChannelFromCOOL' in dir():
   ReadBadChannelFromCOOL = True   

if ( ReadBadChannelFromCOOL ):      
   if 'InputBadChannelSQLiteFile' in dir():
      InputDBConnectionBadChannel = DBConnectionFile(InputBadChannelSQLiteFile)
   else:
      #InputDBConnectionBadChannel = "oracle://ATLAS_COOLPROD;schema=ATLAS_COOLONL_LAR;dbname=CONDBR2;"
      InputDBConnectionBadChannel = "COOLOFL_LAR/CONDBR2"            

###########################################################################
# Print summary
###########################################################################

OFCLog = logging.getLogger( "OFCLog" )
OFCLog.info( " ======================================================== " )
OFCLog.info( " ***                 LAr OFC summary                  *** " )
OFCLog.info( " ======================================================== " )
OFCLog.info( " RunNumber                          = "+str(RunNumber) )
if ( ReadAutoCorrFromCOOL ):
   OFCLog.info( " InputDBConnectionAutoCorr          = "+InputDBConnectionAutoCorr )
else :
   OFCLog.info( " InputAutoCorrPoolFileName          = "+InputAutoCorrPoolFileName )
if 'AutoCorrLArCalibFolderTag' in dir() :   
   OFCLog.info( " AutoCorrLArCalibFolderTag          = "+AutoCorrLArCalibFolderTag )
if ( ReadCaliWaveFromCOOL ):
   OFCLog.info( " InputDBConnectionCaliWave          = "+InputDBConnectionCaliWave )
else :
   OFCLog.info( " InputCaliWavePoolFileName          = "+InputCaliWavePoolFileName )
if ( (ReadAutoCorrFromCOOL or ReadCaliWaveFromCOOL) and ChannelSelection != " " ):
   OFCLog.info( " ChannelSelection                   = "+ChannelSelection )
if 'CaliWaveLArCalibFolderTag' in dir() :   
   OFCLog.info( " CaliWaveLArCalibFolderTag          = "+CaliWaveLArCalibFolderTag )
OFCLog.info( " OutputOFCRootFullFileName          = "+OutputOFCRootFileDir+"/"+OutputOFCRootFileName )
OFCLog.info( " OutputOFCPoolFullFileName          = "+OutputOFCPoolFileDir+"/"+OutputOFCPoolFileName )
OFCLog.info( " OutputShapePoolFileName            = "+OutputShapePoolFileDir+"/"+OutputShapePoolFileName )
OFCLog.info( " OutputObjectSpecOFC                = "+str(OutputObjectSpecOFC) )
OFCLog.info( " OutputObjectSpecTagOFC             = "+str(OutputObjectSpecTagOFC) )
OFCLog.info( " OutputObjectSpecShape              = "+str(OutputObjectSpecShape) )
OFCLog.info( " OutputObjectSpecTagShape           = "+str(OutputObjectSpecTagShape) )
OFCLog.info( " IOVBegin                           = "+str(IOVBegin) )
OFCLog.info( " IOVEnd                             = "+str(IOVEnd) )
OFCLog.info( " LArCalibOutputDB                   = "+OutputDB )
OFCLog.info( " ======================================================== " )
###########################################################################

###########################################################################
#                                                                         #
#                       Global settings                                   #
#                                                                         #
###########################################################################

include( "AthenaCommon/Atlas_Gen.UnixStandardJob.py" )

#
# Provides ByteStreamInputSvc name of the data file to process in the offline context
#

## get a handle to the default top-level algorithm sequence
from AthenaCommon.AlgSequence import AlgSequence 
topSequence = AlgSequence()

## get a handle to the ApplicationManager, to the ServiceManager and to the ToolSvc
from AthenaCommon.AppMgr import (theApp, ServiceMgr as svcMgr,ToolSvc)

include("LArCalibProcessing/LArCalib_MinimalSetup.py")

###########################################################################
#                                                                         #
#    Input conditions data (AutoCorr, CaliWave) for COOL and/or POOL      #
#                                                                         #
###########################################################################

include("AthenaPoolCnvSvc/AthenaPool_jobOptions.py")
include("LArCondAthenaPool/LArCondAthenaPool_joboptions.py")

from IOVDbSvc.CondDB import conddb
PoolFileList     = []

include ("LArCalibProcessing/LArCalib_BadChanTool.py")

if not 'InputBadChannelSQLiteFile' in dir():
   OFCLog.info( "Read Bad Channels from Oracle DB")
else :   
   OFCLog.info( "Read Bad Channels from SQLite file") 

if 'BadChannelsLArCalibFolderTag' in dir() :
   BadChannelsTagSpec = LArCalibFolderTag (BadChannelsFolder,BadChannelsLArCalibFolderTag) 
   conddb.addFolder("",BadChannelsFolder+"<tag>"+BadChannelsTagSpec+"</tag>"+"<dbConnection>"+InputDBConnectionBadChannel+"</dbConnection>")
else :
   conddb.addFolder("",BadChannelsFolder+"<dbConnection>"+InputDBConnectionBadChannel+"</dbConnection>")

if 'MissingFEBsLArCalibFolderTag' in dir() :
   MissingFEBsTagSpec = LArCalibFolderTag (MissingFEBsFolder,MissingFEBsLArCalibFolderTag)   
   conddb.addFolder("",MissingFEBsFolder+"<tag>"+MissingFEBsTagSpec+"</tag>"+"<dbConnection>"+InputDBConnectionBadChannel+"</dbConnection>")
else :
   conddb.addFolder("",MissingFEBsFolder+"<dbConnection>"+InputDBConnectionBadChannel+"</dbConnection>")

## define the DB Gobal Tag :
svcMgr.IOVDbSvc.GlobalTag   = LArCalib_Flags.globalFlagDB
try:
   svcMgr.IOVDbSvc.DBInstance=""
except: 
   pass

from LArCalibProcessing.LArCalibCatalogs import larCalibCatalogs
svcMgr.PoolSvc.ReadCatalog += larCalibCatalogs


if ( ReadAutoCorrFromCOOL or ReadCaliWaveFromCOOL ):
   if not 'InputAutoCorrSQLiteFile' in dir():
      OFCLog.info( "Read AutoCorr from  Oracle DB" )
   else :
      OFCLog.info( "Read AutoCorr from SQLite file" )
      
   if not 'InputCaliWaveSQLiteFile'  in dir():
      OFCLog.info( "Read CaliWave from Oracle DB" )
   else :
      OFCLog.info( "Read CaliWave from SQLite file" )

if ( ReadAutoCorrFromCOOL ):
   AutoCorrFolder  = LArCalib_Flags.LArAutoCorrFolder
   if 'AutoCorrLArCalibFolderTag' in dir() :
      AutoCorrTagSpec = LArCalibFolderTag(AutoCorrFolder,AutoCorrLArCalibFolderTag)
      conddb.addFolder("",AutoCorrFolder+"<tag>"+AutoCorrTagSpec+"</tag>"+"<dbConnection>"+InputDBConnectionAutoCorr+"</dbConnection>"+ ChannelSelection)
   else :
      conddb.addFolder("",AutoCorrFolder+"<dbConnection>"+InputDBConnectionAutoCorr+"</dbConnection>"+ ChannelSelection)
else:
   if 'InputAutoCorrPoolFileName' in dir():
      OFCLog.info( "Read AutoCorr from POOL file" )
      PoolFileList += [ InputAutoCorrPoolDir+"/"+InputAutoCorrPoolFileName ]
   else:
      OFCLog.info( "No PoolFileList found! Please list the POOL files containing AutoCorrelation Matrix or read from COOL." )
      theApp.exit(-1)      
      
if ( ReadCaliWaveFromCOOL ):
   if ( not StripsXtalkCorr ):
      CaliWaveFolder = LArCalib_Flags.LArCaliWaveFolder
   else:
      CaliWaveFolder = LArCalib_Flags.LArCaliWaveFolderXtlk
   if 'CaliWaveLArCalibFolderTag' in dir() :  
      CaliWaveTagSpec = LArCalibFolderTag(CaliWaveFolder,CaliWaveLArCalibFolderTag)
      conddb.addFolder("",CaliWaveFolder+"<tag>"+CaliWaveTagSpec+"</tag>"+"<dbConnection>"+InputDBConnectionCaliWave+"</dbConnection>" + ChannelSelection)
   else :
      conddb.addFolder("",CaliWaveFolder+"<dbConnection>"+InputDBConnectionCaliWave+"</dbConnection>" + ChannelSelection)
   
else:
   if 'InputCaliWavePoolFileName' in dir():
      OFCLog.info( "Read CaliWave from POOL file" )
      PoolFileList += [ InputCaliWavePoolDir+"/"+InputCaliWavePoolFileName ]
   else:
      OFCLog.info( "No PoolFileList found! Please list the POOL files containing CaliWave or read from COOL." )
      theApp.exit(-1)

if ( len(PoolFileList)>0 ):   

   from AthenaCommon.ConfigurableDb import getConfigurable
   svcMgr += getConfigurable( "ProxyProviderSvc" )()
   svcMgr.ProxyProviderSvc.ProviderNames += [ "CondProxyProvider" ]
   
   svcMgr += getConfigurable( "CondProxyProvider" )()
   svcMgr.CondProxyProvider.InputCollections += PoolFileList
   
###########################################################################
#                            OFC computation
###########################################################################

from LArCalibUtils.LArCalibUtilsConf import LArOFCAlg
LArCaliOFCAlg = LArOFCAlg("LArCaliOFCAlg")
LArCaliOFCAlg.ReadCaliWave = True
LArCaliOFCAlg.KeyList   = [ ContainerKey ]
LArCaliOFCAlg.Nphase    = Nphases
LArCaliOFCAlg.Dphase    = Dphases
LArCaliOFCAlg.Ndelay    = Ndelays
LArCaliOFCAlg.Nsample   = Nsamples
LArCaliOFCAlg.Normalize = Normalize
LArCaliOFCAlg.TimeShift = TimeShift
LArCaliOFCAlg.TimeShiftByIndex = TimeShiftByIndex
LArCaliOFCAlg.Verify    = True
LArCaliOFCAlg.FillShape = FillShape
if ( DumpOFC ) :
	LArCaliOFCAlg.DumpOFCfile = "LArOFCCali.dat"
LArCaliOFCAlg.GroupingType = GroupingType
topSequence+=LArCaliOFCAlg

from LArCalibUtils.LArCalibUtilsConf import LArAutoCorrDecoderTool
theLArAutoCorrDecoderTool = LArAutoCorrDecoderTool()
ToolSvc += theLArAutoCorrDecoderTool

###########################################################################

if ( WriteNtuple ) :

   from LArCalibTools.LArCalibToolsConf import LArOFC2Ntuple
   LArOFC2Ntuple = LArOFC2Ntuple("LArOFC2Ntuple")
   LArOFC2Ntuple.ContainerKey = OFCKey 	   
   topSequence+=LArOFC2Ntuple

   if ( FillShape ):
      from LArCalibTools.LArCalibToolsConf import LArShape2Ntuple
      LArShape2Ntuple = LArShape2Ntuple("LArShape2Ntuple")
      LArShape2Ntuple.ContainerKey = ShapeKey 	   
      topSequence+=LArShape2Ntuple
   
   theApp.HistogramPersistency = "ROOT"
   from GaudiSvc.GaudiSvcConf import NTupleSvc
   if os.path.exists(OutputOFCRootFileDir+"/"+OutputOFCRootFileName): 
      os.remove(OutputOFCRootFileDir+"/"+OutputOFCRootFileName)  
   svcMgr += NTupleSvc()
   svcMgr.NTupleSvc.Output = [ "FILE1 DATAFILE='"+OutputOFCRootFileDir+"/"+OutputOFCRootFileName+"' OPT='NEW'" ]
   
if (  WritePoolFile ) :
    
   from RegistrationServices.OutputConditionsAlg import OutputConditionsAlg

   if os.path.exists(OutputOFCPoolFileDir+"/"+OutputOFCPoolFileName): 
      os.remove(OutputOFCPoolFileDir+"/"+OutputOFCPoolFileName)  
   OutputConditionsAlgOFC = OutputConditionsAlg("OutputConditionsAlgOFC",OutputOFCPoolFileDir+"/"+OutputOFCPoolFileName,
                                                [OutputObjectSpecOFC],[OutputObjectSpecTagOFC],WriteIOV)
   OutputConditionsAlgOFC.Run1     = IOVBegin
   if IOVEnd>0:
      OutputConditionsAlgOFC.Run2  = IOVEnd
   
   if ( FillShape ):
      
      if os.path.exists(OutputShapePoolFileDir+"/"+OutputShapePoolFileName): 
         os.remove(OutputShapePoolFileDir+"/"+OutputShapePoolFileName) 
      OutputConditionsAlgShape = OutputConditionsAlg("OutputConditionsAlgShape",OutputShapePoolFileDir+"/"+OutputShapePoolFileName,
                                                     [OutputObjectSpecShape],[OutputObjectSpecTagShape],WriteIOV)
      OutputConditionsAlgShape.Run1     = IOVBegin
      if IOVEnd>0:
         OutputConditionsAlgShape.Run2  = IOVEnd

   svcMgr.IOVDbSvc.dbConnection  = OutputDB
   
   from RegistrationServices.RegistrationServicesConf import IOVRegistrationSvc
   svcMgr += IOVRegistrationSvc()
   svcMgr.IOVRegistrationSvc.OutputLevel = DEBUG
   svcMgr.IOVRegistrationSvc.RecreateFolders = False
   
###########################################################################
#               Use EventSelector to select IOV                           #
###########################################################################	

from McEventSelector.McEventSelectorConf import McEventSelector
svcMgr += McEventSelector("EventSelector")
svcMgr.EventSelector.RunNumber	= int(RunNumber)
svcMgr.EventSelector.EventsPerRun      = 1
svcMgr.EventSelector.FirstEvent	= 1
svcMgr.EventSelector.InitialTimeStamp  = 0
svcMgr.EventSelector.TimeStampInterval = 1

##########################################################################
#          don't remove otherwise infinite loop                          #
##########################################################################

theApp.EvtMax = 1

###########################################################################

svcMgr.MessageSvc.OutputLevel  = WARNING
svcMgr.MessageSvc.defaultLimit = 10000
svcMgr.MessageSvc.Format       = "% F%20W%S%7W%R%T %0W%M"

svcMgr+=CfgMgr.AthenaEventLoopMgr(OutputLevel = VERBOSE)

from AthenaCommon.AppMgr import theAuditorSvc
from AthenaCommon.ConfigurableDb import getConfigurable
theAuditorSvc += getConfigurable("MemStatAuditor")(OutputLevel = DEBUG)
theAuditorSvc += getConfigurable("ChronoAuditor")()
theAuditorSvc += getConfigurable("NameAuditor")()

###########################################################################


