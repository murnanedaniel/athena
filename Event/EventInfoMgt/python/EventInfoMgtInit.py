# Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration

## @file EventInfoMgtInit.py
## @brief Configurable for TagInfoMgr service initialization
## @author RD Schaffer <R.D.Schaffer@cern.ch>
###############################################################
#
# The jobOptions for TagInfoMgr
#
#==============================================================
#
# Required libs:

from __future__ import print_function

def _loadBasicEventInfoMgt():
    """Loads the basic services for EventInfoMgt"""

    import os
    import string
    class EventInfoMgtInit:

        def __init__(self):
            self.release = self.getRelease()

        #----------------------------------------------------
        #
        # Executing a shell command
        #
        def execute (self, cmd):
            #print ('> ' + cmd)
            r = os.popen(cmd)
            lines = []
            for line in r.readlines():
                #print (line)
                line = string.rstrip (line)
                lines.append (line)
            r.close()
            return lines

        # Get the project name and it version version, the top level
        # project, AtlasOffline, has the release number
        def getRelease (self):
            try:
                #print ("EventInfoMgtInit.getRelease: get project, version")
                project = os.environ ['AtlasProject']
                version = os.environ ['AtlasVersion']
                #print ("EventInfoMgtInit.getRelease: project, version",project, version)
                return project + '-' + version
            except Exception:
                # These variables can be missing during CI builds,
                # so don't complain if they're not there.
                #print ("EventInfoMgtInit getRelease: except caught")
                #print (sys.exc_info()[0])
                #print (sys.exc_info()[1])
                pass
            return "Unknown-Unknown"


    ## make sure we have been -at least- minimally correctly configured
    import AthenaCommon.AtlasUnixStandardJob # noqa: F401

    from AthenaCommon.Logging import logging
    from AthenaCommon.AppMgr  import ServiceMgr as svcMgr

    msg = logging.getLogger( 'loadBasicEventInfoMgt' )
    msg.debug( "Loading basic services for EventInfoMgt..." )

    #from EventInfoMgt.EventInfoMgtConf import TagInfoMgr
    from EventInfoMgt.EventInfoMgtConf import TagInfoMgr
    tagInfoMgr = TagInfoMgr()
    svcMgr += tagInfoMgr

    # Add in extra tag for the release number:
    evtMgt  = EventInfoMgtInit()
    release = evtMgt.release
    print ("EventInfoMgtInit: Got release version ",release)
    svcMgr.TagInfoMgr.ExtraTagValuePairs = {"AtlasRelease" : release }

    # Add TagInfoMgr to services activated by the framework 
    from AthenaCommon.AppMgr import theApp
    theApp.CreateSvc += [ tagInfoMgr.name() ]

    # include basic IOVDbSvc functionality for /TagInfo
    from IOVDbSvc import IOVDb   # noqa: F401

    msg.debug( "Loading basic services for EventInfoMgt... [DONE]" )
    return

## load basic services configuration at module import
_loadBasicEventInfoMgt()

## clean-up: avoid running multiple times this method
del _loadBasicEventInfoMgt

    
