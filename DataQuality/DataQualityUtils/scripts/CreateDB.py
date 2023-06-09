#!/usr/bin/env python

# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
# +-----------------------------------------------------+
# | Vassilis D. Kazazakis, CreateDB.py v1.0		|
# | Creates the appropriate folder structure in COOL db	|
# +-----------------------------------------------------+

import sys, os
import DataQualityUtils.DBInfo_Histo as DBInfo

from PyCool import cool, coral

# user & password information should be handled by each individual program
# if you have authentication file, this works

home = os.environ.get('HOME')
os.environ['CORAL_AUTH_PATH'] = home+"/private"

#userName = "ATLAS_COOLOFL_GLOBAL"
#password = "*********" # no need for a pwd 

dropExistingDB = False

connectString  = "sqlite://" +	DBInfo.getWriteServer()
connectString += ";schema=" + 	DBInfo.getSchema()
#connectString += ";user=" + 	userName
connectString += ";dbname=" + 	DBInfo.getDbName()
#connectString += ";password=" + password

if len(sys.argv) > 1 and sys.argv[1] == "createdb":
  dbSvc = cool.DatabaseSvcFactory.databaseService()
  db = dbSvc.createDatabase( connectString )

if len(sys.argv) > 1 and sys.argv[1] == "drop":
	dropExistingDB = True

folderNames = DBInfo.getFolders()

dbSvc = cool.DatabaseSvcFactory.databaseService()
print connectString
db = dbSvc.openDatabase( connectString, False) 

for folderName in folderNames:
  folderList = folderName.split("/")[1:]
  folder = ""
  for name in folderList[:-1]:
          folder = folder + "/" + name
          if db.existsFolderSet(folder):
                  print "Folderset", folder, "already exists."
                  continue
          print "Creating folderset", folder,
          db.createFolderSet(folder)
          print "... Done"

  if dropExistingDB:
          if db.existsFolder(folderName):
                  print "Droping", folderName,
                  folder = db.dropNode(folderName)
                  print "... Done"
  if db.existsFolder(folderName):
          folder = db.getFolder(folderName)
  else:
          spec = cool.RecordSpecification()
          #spec.extend("filename", cool.StorageType.String255)
          #spec.extend("guid", cool.StorageType.Int32)
          spec.extend("Code", cool.StorageType.Int32)
          spec.extend("deadFrac", cool.StorageType.Float)
          spec.extend("Thrust",cool.StorageType.Float)
          # heaven forgive me ...
          if 'SHIFT' in folderName:
            spec.extend("Comment",cool.StorageType.String255)
          print "Creating folder", folderName, 
          # Deprecated/dropped:  folder = db.createFolder(folderName, spec, "", cool.FolderVersioning.MULTI_VERSION, False)
          folderSpec=cool.FolderSpecification(cool.FolderVersioning.MULTI_VERSION, spec)
          folder = db.createFolder(folderName, folderSpec, "", False)
          print "... Done"

  print "Creating channels."
  channelDict = DBInfo.getChannelDict()
  for channelName in channelDict.keys():
          try:
                  folder.createChannel(channelDict[channelName], channelName)
                  print "Created channel:", channelName
          except:
                  print "Channel", channelName, "already exists."
