# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration

import copy

def getGroomedJetsConfig():
    
    dictsConf = []

    GroomedDictsTrimming = [
                 #    { 'Type' : 'Trimming', 'args' : { 'SmallR' : 0.3, 'PtFrac' : 0.01} },
                     { 'Type' : 'Trimming' , 'args' : {} }, # default ={ 'SmallR' : 0.3, 'PtFrac' : 0.05} },
                 #    { 'Type' : 'Trimming', 'args' : { 'SmallR' : 0.3, 'PtFrac' : 0.03} },

                 #    { 'Type' : 'Trimming', 'args' : { 'SmallR' : 0.2, 'PtFrac' : 0.01} },
                 #    { 'Type' : 'Trimming', 'args' : { 'SmallR' : 0.2, 'PtFrac' : 0.03} },
                 #    { 'Type' : 'Trimming', 'args' : { 'SmallR' : 0.2, 'PtFrac' : 0.05} }
                   ]

    GroomedDictsFiltering = [
                     { 'Type' : 'BDRSFiltering', 'args' : { 'minSplitR' : 0., 'massFraction' : 0.20 } },
                     { 'Type' : 'BDRSFiltering', 'args' : { 'minSplitR' : 0., 'massFraction' : 0.33 } },
                     { 'Type' : 'BDRSFiltering', 'args' : { 'minSplitR' : 0., 'massFraction' : 0.67 } }
                   ]

    GroomedDictsPruning = [
                     { 'Type' : 'Pruning', 'args' : { 'Name' : 'PrunedKt', 'RcutFactor' : 0.5, 'Zcut' : 0.05 } },
                     { 'Type' : 'Pruning', 'args' : { 'Name' : 'PrunedKt', 'RcutFactor' : 0.5, 'Zcut' : 0.10 } },
                     { 'Type' : 'Pruning', 'args' : { 'Name' : 'PrunedKt', 'RcutFactor' : 0.5, 'Zcut' : 0.15 } },
    
                     { 'Type' : 'Pruning', 'args' : { 'Name' : 'PrunedKt', 'RcutFactor' : 1., 'Zcut' : 0.05 } },
                     { 'Type' : 'Pruning', 'args' : { 'Name' : 'PrunedKt', 'RcutFactor' : 1., 'Zcut' : 0.10 } },
                     { 'Type' : 'Pruning', 'args' : { 'Name' : 'PrunedKt', 'RcutFactor' : 1., 'Zcut' : 0.15 } }
                   ]

    GroomedDictsKtSubJets = [ { 'Type' : 'KtSubJets' ,'args' : {} } ] # default = { 'NSubJets' : 3  }  }  ]


#### AntiKt R=0.4 ####

    ParentDict = {
                  'JetFinder'  :  'AntiKt',
                  'JetdR'      :   0.4,
                  'JetInput'   :   'Track',
                 }

    dictsConf += [ [ParentDict, copy.deepcopy(GroomedDictsTrimming)]  ]

##
    ParentDict = {
                  'JetFinder'  :  'AntiKt',
                  'JetdR'      :   0.4,
                  'JetInput'   :   'Truth',
                 }

    dictsConf += [ [ParentDict, copy.deepcopy(GroomedDictsTrimming)]  ]

##
    ParentDict = {
                  'JetFinder'  :  'AntiKt',
                  'JetdR'      :   0.4,
                  'JetInput'   :   'LCTopo',
                 }

    dictsConf += [ [ParentDict, copy.deepcopy(GroomedDictsTrimming)]  ]


#### AntiKt R=0.6 ####

    ParentDict = {
                  'JetFinder'  :  'AntiKt',
                  'JetdR'      :   0.6,
                  'JetInput'   :   'Track',
                 }

    dictsConf += [ [ParentDict, []] ] #copy.deepcopy(GroomedDictsTrimming)]  ]

##
    ParentDict = {
                  'JetFinder'  :  'AntiKt',
                  'JetdR'      :   0.6,
                  'JetInput'   :   'Truth',
                 }

    dictsConf += [ [ParentDict,  []] ] #copy.deepcopy(GroomedDictsTrimming)]  ]

##
    ParentDict = {
                  'JetFinder'  :  'AntiKt',
                  'JetdR'      :   0.6,
                  'JetInput'   :   'LCTopo',
                 }

    dictsConf += [ [ParentDict, []] ] #copy.deepcopy(GroomedDictsTrimming)]  ]



#### AntiKt R=1.0 ####

    ParentDict = {
                  'JetFinder'  :  'AntiKt',
                  'JetdR'      :   1.0,
                  'JetInput'   :   'Track',
                 }

    dictsConf += [ [ParentDict, copy.deepcopy(GroomedDictsTrimming+GroomedDictsKtSubJets)]  ]  

##
    ParentDict = {
                  'JetFinder'  :  'AntiKt',
                  'JetdR'      :   1.0,
                  'JetInput'   :   'Truth',
                 }

    dictsConf += [ [ParentDict, copy.deepcopy(GroomedDictsTrimming+GroomedDictsKtSubJets)]  ]

##
    ParentDict = {
                  'JetFinder'  :  'AntiKt',
                  'JetdR'      :   1.0,
                  'JetInput'   :   'LCTopo',
                 }

    dictsConf += [ [ParentDict, copy.deepcopy(GroomedDictsTrimming+GroomedDictsKtSubJets)]  ]


#### C/A R=1.2####
    
    ParentDict = {
                  'JetFinder'  :  'CamKt',
                  'JetdR'      :   1.2,
                  'JetInput'   :   'Truth',
                 }

    #dictsConf += [ [ParentDict, copy.deepcopy(GroomedDictsTrimming+GroomedDictsKtSubJets)]  ]
    dictsConf += [ [ParentDict,  []] ] #copy.deepcopy(GroomedDictsTrimming)]  ]

##
    ParentDict = {
                  'JetFinder'  :  'CamKt',
                  'JetdR'      :   1.2,
                  'JetInput'   :   'Track',
                 }

    #dictsConf += [ [ParentDict, copy.deepcopy(GroomedDictsTrimming+GroomedDictsKtSubJets)]  ]
    dictsConf += [ [ParentDict,  []] ] #copy.deepcopy(GroomedDictsTrimming)]  ]
##
    ParentDict = {
                  'JetFinder'  :  'CamKt',
                  'JetdR'      :   1.2,
                  'JetInput'   :   'LCTopo',
                 }

    #dictsConf += [ [ParentDict, copy.deepcopy(GroomedDictsTrimming+GroomedDictsKtSubJets)]  ]
    dictsConf += [ [ParentDict,  []] ] #copy.deepcopy(GroomedDictsTrimming)]  ]  

#### AntiKt R=1.2 ####
    '''
    ParentDict = {
                  'JetFinder'  :  'AntiKt',
                  'JetdR'      :   1.2,
                  'JetInput'   :   'Track',
                 }

    dictsConf += [ [ParentDict, copy.deepcopy(GroomedDictsTrimming+GroomedDictsKtSubJets)]  ]

##
    ParentDict = {
                  'JetFinder'  :  'AntiKt',
                  'JetdR'      :   1.2,
                  'JetInput'   :   'Truth',
                 }

    dictsConf += [ [ParentDict, copy.deepcopy(GroomedDictsTrimming+GroomedDictsKtSubJets)]  ]

##
    ParentDict = {
                  'JetFinder'  :  'AntiKt',
                  'JetdR'      :   1.2,
                  'JetInput'   :   'LCTopo',
                 }

    dictsConf += [ [ParentDict, copy.deepcopy(GroomedDictsTrimming+GroomedDictsKtSubJets)]  ]


#### C/A R=1.0####

    ParentDict = {
                  'JetFinder'  :  'CamKt',
                  'JetdR'      :   1.0,
                  'JetInput'   :   'Truth',
                 }

    dictsConf += [ [ParentDict, copy.deepcopy(GroomedDictsTrimming+GroomedDictsKtSubJets)]  ]

##
    ParentDict = {
                  'JetFinder'  :  'CamKt',
                  'JetdR'      :   1.0,
                  'JetInput'   :   'Track',
                 }

    dictsConf += [ [ParentDict, copy.deepcopy(GroomedDictsTrimming+GroomedDictsKtSubJets)]  ]

##
    ParentDict = {
                  'JetFinder'  :  'CamKt',
                  'JetdR'      :   1.0,
                  'JetInput'   :   'LCTopo',
                 }

    dictsConf += [ [ParentDict, copy.deepcopy(GroomedDictsTrimming+GroomedDictsKtSubJets)]  ]
    '''

    return dictsConf
