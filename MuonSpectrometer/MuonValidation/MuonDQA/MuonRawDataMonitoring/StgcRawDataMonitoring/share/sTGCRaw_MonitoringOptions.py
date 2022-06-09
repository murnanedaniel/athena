###### Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration #######

from AthenaCommon.AthenaCommonFlags import athenaCommonFlags
isTier0Flag = not athenaCommonFlags.isOnline()

if not 'MuonDQAFlags' in dir():
    print("MuonDQAFlags.py: MuonDQAFlags not yet imported - I import them now")
    from MuonDQAMonFlags.MuonDQAFlags import MuonDQAFlags as MuonDQAFlags

mmRawMonMan = AthenaMonManager(name="sTGCRawMonManager",
                                FileKey             = DQMonFlags.monManFileKey(),
                                Environment         = DQMonFlags.monManEnvironment(),
                                OutputLevel         = muonOutputLevel)

############# MMRawDataValAlg #############
from sTGCRawDataMonitoring.sTGCRawDataMonitoringConf import sTGCRawDataValAlg

#This it the main tool w/ detailed histograms
sTGCRawDataValAlg_main = sTGCRawDataValAlg(name='sTGCRawDataValAlg_main',
                                       Title='sTGC',
                                       OutputLevel = INFO
                                       )

sTGCRawMonMan.AthenaMonTools += [ sTGCRawDataValAlg_main ]


topSequence += sTGCRawMonMan
print sTGCRawMonMan
##################################################################
