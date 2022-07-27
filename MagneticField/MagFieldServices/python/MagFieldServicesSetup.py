# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

# JobOption fragment to set up the magnetic field services and algorithms 
# Valerio Ippolito - Harvard University

# inspired by https://svnweb.cern.ch/trac/atlasoff/browser/MuonSpectrometer/MuonCnv/MuonCnvExample/trunk/python/MuonCalibConfig.py

from AthenaCommon.Logging import logging
logging.getLogger().info("Importing %s", __name__)

from AthenaCommon.AthenaCommonFlags import athenaCommonFlags
from AthenaCommon import CfgMgr

#--------------------------------------------------------------

def AtlasFieldCacheCondAlg(name="AtlasFieldCacheCondAlg",**kwargs):
  if athenaCommonFlags.isOnline():
    kwargs.setdefault( "UseDCS", False )
    # currents for scaling wrt to map currents
    kwargs.setdefault( "UseSoleCurrent", 7730 )
    kwargs.setdefault( "UseToroCurrent", 20400 )
    kwargs.setdefault( "LockMapCurrents", True )
  else:
    kwargs.setdefault( "UseDCS", True )
  return CfgMgr.MagField__AtlasFieldCacheCondAlg(name,**kwargs)

def AtlasFieldMapCondAlg(name="AtlasFieldMapCondAlg",**kwargs):
  if athenaCommonFlags.isOnline():
    # online does not use DCS
    kwargs.setdefault( "UseMapsFromCOOL", False )

  return CfgMgr.MagField__AtlasFieldMapCondAlg(name,**kwargs)



def H8FieldSvc(name="H8FieldSvc",**kwargs):
  return CfgMgr.MagField__H8FieldSvc(name,**kwargs)

def GetFieldCacheCondAlg(name="AtlasFieldCacheCondAlg",**kwargs):
    return AtlasFieldCacheCondAlg(name, **kwargs)
    
def GetFieldMapCondAlg(name="AtlasFieldMapCondAlg",**kwargs):
    return AtlasFieldMapCondAlg(name, **kwargs)
