# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
# file to simulate the HypoTool configuration of the signatures

from AthenaConfiguration.ComponentFactory import CompFactory

def TestHypoTool(name, prop, threshold_value):
    value  =  int(threshold_value)*1000
    UseThisLinkName="initialRoI"   
    HLTTest__TestHypoTool=CompFactory.getComp("HLTTest::TestHypoTool") 
    return HLTTest__TestHypoTool(name, Threshold=value, Property=prop, LinkName=UseThisLinkName)

def MuTestHypoTool(chainDict):
    name = chainDict['chainName']
    threshold = getThreshold(chainDict) 
    return TestHypoTool(name,prop="pt", threshold_value=threshold)

def ElTestHypoTool(chainDict):
    name = chainDict['chainName']
    threshold = getThreshold(chainDict) 
    return TestHypoTool(name,prop="et", threshold_value=threshold)

def GammTestHypoTool(chainDict):
    name = chainDict['chainName']
    threshold = getThreshold(chainDict) 
    return TestHypoTool(name,prop="et", threshold_value=threshold)


def MuTest2HypoTool(chainDict):
    name = chainDict['chainName']
    threshold = getThreshold(chainDict) 
    return TestHypoTool(name,prop="pt2", threshold_value=threshold)

def ElTest2HypoTool(chainDict):
    name = chainDict['chainName']
    threshold = getThreshold(chainDict) 
    return TestHypoTool(name,prop="et", threshold_value=threshold)


def getThreshold(chainDict):
    name = chainDict['chainParts'][0]['chainPartName']
    from TriggerMenuMT.HLT.Config.Utility.DictFromChainName import getChainThresholdFromName
    return getChainThresholdFromName( name.split("_"), "TestChain")



def dimuDrComboHypoTool(chainDict):
    from DecisionHandling.DecisionHandlingConf import DeltaRRoIComboHypoTool
    name = chainDict['chainName']
    tool= DeltaRRoIComboHypoTool(name)
    tool.DRcut=0.3
    return tool

