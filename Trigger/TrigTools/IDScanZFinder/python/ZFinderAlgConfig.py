# Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration

from TrigEDMConfig.TriggerEDMRun3 import recordable
from IDScanZFinder.IDScanZFinderConf import TrigZFinderAlg
from IDScanZFinder.IDScanZFinderConf import TrigZFinder


MinBiasZFinderAlg = TrigZFinderAlg("TrigZFinderAlg", vertexKey=recordable("HLT_vtx_z"))
MinBiasZFinderAlg.ZFinderTools += [TrigZFinder("default")]
postfix=-1
def tool_name():
    global postfix
    postfix += 1
    return "ZFindersLike"+str(postfix)
# Default: TripletMode=0, MinZBinSize=0.2, PhiBinSize=0.20, NumberOfPeaks=1, UseOnlyPixels=False, MaxLayer=?
MinBiasZFinderAlg.ZFinderTools += [TrigZFinder(tool_name(), TripletMode=1, MinZBinSize=0.2, PhiBinSize=0.20, NumberOfPeaks=5, UseOnlyPixels=True, MaxLayer=3)]
MinBiasZFinderAlg.ZFinderTools += [TrigZFinder(tool_name(), TripletMode=1, MinZBinSize=0.2, PhiBinSize=0.30, NumberOfPeaks=5, UseOnlyPixels=True, MaxLayer=3)]
MinBiasZFinderAlg.ZFinderTools += [TrigZFinder(tool_name(), TripletMode=1, MinZBinSize=0.2, PhiBinSize=0.40, NumberOfPeaks=5, UseOnlyPixels=True, MaxLayer=3)]
MinBiasZFinderAlg.ZFinderTools += [TrigZFinder(tool_name(), TripletMode=1, MinZBinSize=0.2, PhiBinSize=0.50, NumberOfPeaks=5, UseOnlyPixels=True, MaxLayer=3)]

MinBiasZFinderAlg.ZFinderTools += [TrigZFinder(tool_name(), TripletMode=1, MinZBinSize=0.5, PhiBinSize=0.20, NumberOfPeaks=5, UseOnlyPixels=True, MaxLayer=3)]
MinBiasZFinderAlg.ZFinderTools += [TrigZFinder(tool_name(), TripletMode=1, MinZBinSize=1.5, PhiBinSize=0.20, NumberOfPeaks=5, UseOnlyPixels=True, MaxLayer=3)]
MinBiasZFinderAlg.ZFinderTools += [TrigZFinder(tool_name(), TripletMode=1, MinZBinSize=2.5, PhiBinSize=0.20, NumberOfPeaks=5, UseOnlyPixels=True, MaxLayer=3)]
MinBiasZFinderAlg.ZFinderTools += [TrigZFinder(tool_name(), TripletMode=1, MinZBinSize=3.5, PhiBinSize=0.20, NumberOfPeaks=5, UseOnlyPixels=True, MaxLayer=3)]

MinBiasZFinderAlg.ZFinderTools += [TrigZFinder(tool_name(), TripletMode=1, MinZBinSize=0.5, PhiBinSize=0.50, NumberOfPeaks=5, UseOnlyPixels=True, MaxLayer=3)]
MinBiasZFinderAlg.ZFinderTools += [TrigZFinder(tool_name(), TripletMode=1, MinZBinSize=1.5, PhiBinSize=0.50, NumberOfPeaks=5, UseOnlyPixels=True, MaxLayer=3)]
MinBiasZFinderAlg.ZFinderTools += [TrigZFinder(tool_name(), TripletMode=1, MinZBinSize=2.5, PhiBinSize=0.50, NumberOfPeaks=5, UseOnlyPixels=True, MaxLayer=3)]
MinBiasZFinderAlg.ZFinderTools += [TrigZFinder(tool_name(), TripletMode=1, MinZBinSize=3.5, PhiBinSize=0.50, NumberOfPeaks=5, UseOnlyPixels=True, MaxLayer=3)]

MinBiasZFinderAlg.ZFinderTools += [TrigZFinder(tool_name(), TripletMode=1, MinZBinSize=3.5, PhiBinSize=0.30, NumberOfPeaks=5, UseOnlyPixels=True, MaxLayer=3)]
MinBiasZFinderAlg.ZFinderTools += [TrigZFinder(tool_name(), TripletMode=1, MinZBinSize=3.5, PhiBinSize=0.40, NumberOfPeaks=5, UseOnlyPixels=True, MaxLayer=3)]


from AthenaMonitoringKernel.GenericMonitoringTool import GenericMonitoringTool
monTool = GenericMonitoringTool('MonTool')

monTool.defineHistogram( 'ZVertex', path='EXPERT', type='TH1F', title='Vertex Z distribution;z [mm];Entries',
                         xbins=400, xmin=-200, xmax=200 )
monTool.defineHistogram( 'ZVertexWeight', path='EXPERT', type='TH1F', title='Vertex Weight;Weight;Entries',
                         xbins=100, xmin=0.0, xmax=100 )
                         
MinBiasZFinderAlg.MonTool = monTool