# Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration

from AthenaConfiguration.ComponentFactory import CompFactory

def HLTResultMTMakerCfg(name="HLTResultMTMaker"):
   from AthenaMonitoringKernel.GenericMonitoringTool import GenericMonitoringTool

   m = CompFactory.HLTResultMTMaker(name)

   # ROBs/SubDets which are enabled but not necessarily part of the ROS-ROB map
   from libpyeformat_helper import SourceIdentifier,SubDetector
   subdets = [
      SubDetector.TDAQ_CALO_CLUSTER_PROC_ROI,
      SubDetector.TDAQ_CALO_JET_PROC_ROI,
      SubDetector.TDAQ_HLT,
      SubDetector.TDAQ_FTK,
      SubDetector.TDAQ_CALO_TOPO_PROC,
      SubDetector.TDAQ_CALO_DIGITAL_PROC,
      SubDetector.TDAQ_CALO_FEAT_EXTRACT_ROI,
   ]
   m.ExtraEnabledSubDets = []
   for subdetId in subdets:
      m.ExtraEnabledSubDets.append( int(subdetId) )

   def addROBs(dest,subdet,modules):
      for moduleId in modules:
         dest.append(SourceIdentifier(subdet,moduleId).code())

   m.ExtraEnabledROBs = []
   addROBs(m.ExtraEnabledROBs, SubDetector.TDAQ_CALO_CLUSTER_PROC_ROI, [0xa8, 0xa9, 0xaa, 0xab])
   addROBs(m.ExtraEnabledROBs, SubDetector.TDAQ_CALO_JET_PROC_ROI,     [0xac, 0xad])
   addROBs(m.ExtraEnabledROBs, SubDetector.TDAQ_MUON_CTP_INTERFACE,    [0x01])
   addROBs(m.ExtraEnabledROBs, SubDetector.TDAQ_CTP,                   [0x01])
   addROBs(m.ExtraEnabledROBs, SubDetector.TDAQ_CALO_TOPO_PROC,        [0x81, 0x91, 0x82, 0x92])

   # Configure HLT result monitoring histograms
   m.MonTool = GenericMonitoringTool('MonTool', HistPath='HLTFramework/'+name)
   m.MonTool.defineHistogram('TIME_makeResult', path='EXPERT', type='TH1F', title='makeResult() call time;Time [ms];Events',
                             xbins=200, xmin=0, xmax=50 )

   return m

def TriggerEDMSerialiserToolCfg(name="TriggerEDMSerialiserTool"):
   from AthenaCommon.Configurable import Configurable
   Configurable.configurableRun3Behavior += 1

   from TrigEDMConfig.DataScoutingInfo import getFullHLTResultID

   TriggerEDMSerialiserTool = CompFactory.TriggerEDMSerialiserTool  
   # Configuration helper methods
   def addCollection(self, typeNameAux, moduleIds):
      moduleIdsRepresentation = ','.join([f"{id}" for id in moduleIds])
      self.CollectionsToSerialize.append( f"{typeNameAux};{moduleIdsRepresentation}")

   def addCollectionToMainResult(self, typeNameAux):
      self.addCollection(typeNameAux,moduleIds=[getFullHLTResultID()])

   def addCollectionListToResults(self, typeNameAuxList, moduleIds):
      for typeNameAux in typeNameAuxList:
         self.addCollection(typeNameAux, moduleIds)

   def addCollectionListToMainResult(self, typeNameAuxList):
      self.addCollectionListToResults(typeNameAuxList,moduleIds=[getFullHLTResultID()])
   # Add the helper methods to the TriggerEDMSerialiserTool python class
   TriggerEDMSerialiserTool.addCollection = addCollection
   TriggerEDMSerialiserTool.addCollectionToMainResult = addCollectionToMainResult
   TriggerEDMSerialiserTool.addCollectionListToResults = addCollectionListToResults
   TriggerEDMSerialiserTool.addCollectionListToMainResult = addCollectionListToMainResult

   # Create and return a serialiser tool object
   serialiser = TriggerEDMSerialiserTool(name)

   from TrigEDMConfig.TriggerEDMRun3 import tpMap
   tpTool = CompFactory.TrigSerTPTool()
   tpTool.TPMap = tpMap()
   serialiser.TPTool = tpTool

   from TrigEDMConfig.DataScoutingInfo import TruncationThresholds as truncThresholds
   serialiser.TruncationThresholds = truncThresholds

   # Configure monitoring histograms
   from AthenaMonitoringKernel.GenericMonitoringTool import GenericMonitoringTool
   serialiser.MonTool = GenericMonitoringTool('MonTool', HistPath='HLTFramework/'+name)
   serialiser.MonTool.defineHistogram('Truncation_ModuleId', path='EXPERT', type='TH1F',
                                      title='Module IDs of truncated HLT results;Module ID;Num of truncated results',
                                      xbins=20, xmin=0, xmax=20)
   serialiser.MonTool.defineHistogram('Truncation_TotalSize', path='EXPERT', type='TH1F',
                                      title='Total size of truncated HLT result;Size [kB];Num of truncated results',
                                      xbins=200, xmin=0, xmax=5000)
   serialiser.MonTool.defineHistogram('Truncation_LargestName', path='EXPERT', type='TH1F',
                                      title='Name of the largest collection;;Num of truncated results',
                                      xbins=1, xmin=0, xmax=1)
   serialiser.MonTool.defineHistogram('Truncation_LargestSize', path='EXPERT', type='TH1F',
                                      title='Size of the largest collection;Size [kB];Num of truncated results',
                                      xbins=200, xmin=0, xmax=5000)

   Configurable.configurableRun3Behavior -= 1
   return serialiser

def StreamTagMakerToolCfg(name="StreamTagMakerTool"):

   stmaker = CompFactory.StreamTagMakerTool(name)
   # Extra configuration may come here

   return stmaker


def TriggerBitsMakerToolCfg(name="TriggerBitsMakerTool"):

   bitsmaker = CompFactory.TriggerBitsMakerTool(name)
   # Extra configuration may come here

   return bitsmaker
