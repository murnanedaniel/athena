###################################################################
#
# Job options file for accessing LArRawConditions objects from COOL
#
#==================================================================

include.block ( "LArConditionsCommon/LArConditionsCommon_MC_jobOptions.py" )

from LArConditionsCommon.LArCondFlags import larCondFlags 

LArDBConnection = ""
LArDB = "LAR_OFL"

if larCondFlags.LArDBConnection.statusOn  :
  LArDBConnection = larCondFlags.LArDBConnection()
  LArDB = ""

from AthenaCommon.AlgSequence import AthSequencer
condSeq = AthSequencer("AthCondSeq")

larCondFlags.config_ElecCalibMC()

if svcMgr.MessageSvc.OutputLevel <= DEBUG :
  print larCondFlags

# POOL Converters
#include( "LArCondAthenaPool/LArCondAthenaPool_joboptions.py" )
#include ("LArRawConditions/LArRawConditionsDict_joboptions.py")


# Access to IOVSvc and IOVDbSvc
# Must list the folders to be used for reading
#--------------------------------------------------------------

from IOVDbSvc.CondDB import conddb

larCondDBFolders = ["/LAR/ElecCalibMC/Ramp",
                    "/LAR/ElecCalibMC/AutoCorr",
                    "/LAR/ElecCalibMC/DAC2uA",
                    "/LAR/ElecCalibMC/Pedestal",
                    "/LAR/ElecCalibMC/Noise",
                    "/LAR/ElecCalibMC/fSampl",
                    "/LAR/ElecCalibMC/uA2MeV",
                    "/LAR/ElecCalibMC/MinBias",
                    "/LAR/ElecCalibMC/MinBiasAverage"]

if larCondFlags.useMCShape():
    larCondDBFolders += ["/LAR/ElecCalibMC/Shape"]

include( "LArConditionsCommon/LArIdMap_MC_jobOptions.py" )

from LArBadChannelTool.LArBadChannelToolConf import LArBadChannelCondAlg,LArBadFebCondAlg
larCondDBFolders += [("CondAttrListCollection","/LAR/BadChannels/BadChannels")]
condSeq+=LArBadChannelCondAlg(ReadKey="/LAR/BadChannels/BadChannels")
larCondDBFolders += [("AthenaAttributeList","/LAR/BadChannels/MissingFEBs")]
condSeq+=LArBadFebCondAlg(ReadKey="/LAR/BadChannels/MissingFEBs")

## these may be conditional. 
if larCondFlags.hasMphys() :
  larCondDBFolders += ["/LAR/ElecCalibMC/MphysOverMcal"]

# HV Scale Corr
if larCondFlags.hasHVCorr() :
  larCondDBFolders += [ ('LArHVScaleCorrComplete', '/LAR/ElecCalibMC/HVScaleCorr') ]


## fill them all 
for i in larCondDBFolders :
  className = None
  if type(i) == type(()):
    className, i = i
  conddb.addFolder(LArDB,i+LArDBConnection, className=className)
  ## allow override
  larCondFlags.addTag(i,conddb)  

## apply hierarchical tag                    
larCondFlags.addTag('/LAR/ElecCalibMC',conddb)
