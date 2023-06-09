#--------------------------------------------------------------
# jobOptions for writing RDO to POOL
#--------------------------------------------------------------
# Load POOL support
include( "AthenaPoolCnvSvc/WriteAthenaPool_jobOptions.py" )

# check dictionary
include( "AthenaSealSvc/AthenaSealSvc_joboptions.py" )

#AthenaSealSvc.CheckDictionary = true;
# Define the output Db parameters (the default value are shown)
Stream1.OutputFile = "MuonPool.root"

# Converters:
include( "EventAthenaPool/EventAthenaPool_joboptions.py" )
include( "MuonEventAthenaPool/MuonEventAthenaPool_joboptions.py" )
#include( "MuonEventAthenaPool/MuonEventAthenaPoolDict_joboptions.py")

# list of output objects key
# EventInfo
Stream1 = Algorithm( "Stream1" )
Stream1.ItemList+=["EventInfo#*"]
# MDT
Stream1.ItemList+=["MdtCsmContainer#*"]
# RPC
Stream1.ItemList+=["RpcPadContainer#*"]
# TGC
Stream1.ItemList+=["TgcRdoContainer#*"]
# CSC
Stream1.ItemList+=["CscRawDataContainer#*"]
