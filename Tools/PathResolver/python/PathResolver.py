# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration


#usage is:
#from PathResolver import PathResolver
#PathResolver.FindCalibFile("blah")

#if working in athena release, we will need to start the ApplicationMgr so that the appropriate MessageSvc is created
try:
   #check if Gaudi already started
   from GaudiPython.Bindings import _gaudi as GaudiAppMgr
   if GaudiAppMgr==None:
      print "PathResolver: Starting Gaudi for access to correct instance of MessageSvc"
      from GaudiPython import AppMgr
      AppMgr() #note this is not the same thing as 'theApp' in joboptions, but 'theApp' will use this
      #the line above starts the application mgr, and _gaudi will become set to that
except ImportError:
   pass

import cppyy
cppyy.loadDict('libPathResolver')
cppyy.loadDict('libPathResolverDict') #why do we have to load both the library and the dict :-( To be fixed by ROOT...
FindCalibFile = cppyy.gbl.PathResolverFindCalibFile
FindCalibDirectory = cppyy.gbl.PathResolverFindCalibDirectory
SetOutputLevel = cppyy.gbl.PathResolverSetOutputLevel

