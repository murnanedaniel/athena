# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration

from TrigmuComb.TrigmuCombConf import *
from TrigmuComb.TrigmuCombMonitoring import *
from AthenaCommon.GlobalFlags import globalflags
from AthenaCommon.AppMgr import ServiceMgr
from TrigTimeMonitor.TrigTimeHistToolConfig import TrigTimeHistToolConfig
from MuonByteStream.MuonByteStreamFlags import muonByteStreamFlags
#from TrigMuonBackExtrapolator.TrigMuonBackExtrapolatorConfig import *
from AthenaCommon.AppMgr import ToolSvc
#from CombinedMuonRefit.CombinedMuonRefitConf import ExtrapolateMuonToIPTool

#ToolSvc += MuonBackExtrapolatorForAlignedDet()
#ToolSvc += MuonBackExtrapolatorForMisalignedDet()
#ToolSvc += MuonBackExtrapolatorForData()
#ToolSvc += ExtrapolateMuonToIPTool()

class TrigmuCombConfig (muComb):

    __slots__ = []
    
    def __new__( cls, *args, **kwargs ):
      newargs = ['%s_%s_%s' % (cls.getType(),args[0],args[1]) ] + list(args)
      if len(args) == 3:
        if args[2] == 'OLDextr':
          newargs = ['%s_%s_%s_%s' % (cls.getType(),args[0],args[1],args[2]) ] + list(args)
      return super( TrigmuCombConfig, cls ).__new__( cls, *newargs, **kwargs )

    def __init__( self, name, *args, **kwargs ):
      super( TrigmuCombConfig, self ).__init__( name )

      self.UseAthenaFieldService = True
      self.MuCombStrategy = 0

      self.UseBackExtrapolatorG4 = True
      if len(args) == 3:
        if args[2] == 'OLDextr':
          self.UseBackExtrapolatorG4 = False

      self.MinPtTRK = 0.

      if globalflags.DetDescrVersion().find('CSC')!=-1:
        self.IDSCANBarrelRes  = [0.02169,0.0004186]
        self.IDSCANEndcap1Res = [0.03054,0.000325]
        self.IDSCANEndcap2Res = [0.03557,0.0005383]
        self.IDSCANEndcap3Res = [0.04755,0.0007718]
        self.IDSCANEndcap4Res = [0.06035,0.0001145]

      self.IDalgo='InDetTrigTrackingxAODCnv_Muon_FTF'
          
      validation = TrigmuCombValidationMonitoring()
      online     = TrigmuCombOnlineMonitoring()
      cosmic     = TrigmuCombCosmicMonitoring()
      
      time = TrigTimeHistToolConfig("Time")
      time.TimerHistLimits = [0, 5]
      
      self.AthenaMonTools = [ validation, online, cosmic, time ]
