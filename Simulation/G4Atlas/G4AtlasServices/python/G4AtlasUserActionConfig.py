# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration

from G4AtlasApps.SimFlags import simFlags
from AthenaCommon.BeamFlags import jobproperties
from AthenaCommon import CfgGetter,CfgMgr,Logging

# actions to be run at begin of run
def defaultBoRActions():

    defaultUA=[]

    if not simFlags.ISFRun:
        defaultUA+=['G4SimTimer']

    if jobproperties.Beam.beamType() == 'cosmics' and not simFlags.CavernBG.statusOn:
        defaultUA+=['CosmicPerigeeAction']
        if simFlags.ISFRun:
            CfgGetter.getPublicTool('CosmicPerigeeAction').AllowMods=True

    if simFlags.StoppedParticleFile.statusOn:
        defaultUA+=['StoppedParticleAction']

    return defaultUA

# actions to be run at end of run
def defaultEoRActions():

    defaultUA=[]

    if not simFlags.ISFRun:
        defaultUA+=['G4SimTimer']

    defaultUA+=['G4TrackCounter']

    if simFlags.StoppedParticleFile.statusOn:
        defaultUA+=['StoppedParticleAction']
        defaultUA+=['G4CosmicFilter']

    if jobproperties.Beam.beamType() == 'cosmics':
        defaultUA+=['G4CosmicFilter']

    return defaultUA

# begin of event
def defaultBoEActions():

    defaultUA=[]

    if not simFlags.ISFRun:
        defaultUA+=['G4SimTimer']

    defaultUA+=['G4TrackCounter']

    if jobproperties.Beam.beamType() == 'cosmics' and not simFlags.CavernBG.statusOn:
        defaultUA+=['CosmicPerigeeAction']
        if simFlags.ISFRun:
            CfgGetter.getPublicTool('CosmicPerigeeAction').AllowMods=True

    if simFlags.StoppedParticleFile.statusOn:
        defaultUA+=['StoppedParticleAction']


    return defaultUA

# end of event
def defaultEoEActions():

    defaultUA=[]

    if not simFlags.ISFRun:
        defaultUA+=['G4SimTimer']

    if jobproperties.Beam.beamType() == 'cosmics' and not simFlags.CavernBG.statusOn:
        defaultUA+=['CosmicPerigeeAction']
        if simFlags.ISFRun:
            CfgGetter.getPublicTool('CosmicPerigeeAction').AllowMods=True
    
    if simFlags.CavernBG.statusOn and simFlags.CavernBG.get_Value() == 'Read':

        defaultUA+=['HitWrapper']

    if simFlags.StoppedParticleFile.statusOn:
        defaultUA+=['StoppedParticleAction']
        defaultUA+=['G4CosmicFilter']
        if simFlags.ISFRun:
             CfgGetter.getPublicTool('G4CosmicFilter').VolumeName='StoppingPositions'
    
    if jobproperties.Beam.beamType() == 'cosmics':
        defaultUA+=['G4CosmicFilter']

    return defaultUA

# stepping
def defaultSteActions():

    defaultUA=[]

    if jobproperties.Beam.beamType() == 'cosmics' and not simFlags.CavernBG.statusOn:
        defaultUA+=['CosmicPerigeeAction']
        if simFlags.ISFRun:
            CfgGetter.getPublicTool('CosmicPerigeeAction').AllowMods=True

    if simFlags.PhysicsList == 'QGSP_BERT_HP':
        defaultUA+=['PhotonKiller']

    return defaultUA

# PreUserTracking
def defaultBoTActions():

    defaultUA=[]

    defaultUA+=['G4TrackCounter']

    return defaultUA

# PostUserTracking
def defaultEoTActions():

    defaultUA=[]

    return defaultUA

def defaultStaClaActions():
    return []

def defaultStaPrepareActions():
    return []

def defaultStaNewStageActions():
    return []

class G4AtlasUAStore:

    def __init__(self):
        from AthenaCommon.CfgGetter import getService
        self.systemActions=[]
        self.log=Logging.logging.getLogger('G4UAStore')
        self.theSvc=getService('UserActionSvc')
        # map roles to attributes
        self.roleMap={'BeginOfRun': 'BeginOfRunActions',
                      'EndOfRun': 'EndOfRunActions',
                      'BeginOfEvent': 'BeginOfEventActions',
                      'EndOfEvent': 'EndOfEventActions',
                      'BeginOfTracking': 'PreTrackingActions',
                      'EndOfTracking': 'PostTrackingActions',
                      'Step': 'SteppingActions',
                      'Classification': 'StackingActionsClassification',
                      'PrepareNewEvent': 'StackingActionsPrepareNewEvent',
                      'NewStage':  'StackingActionsNewStage'}

        # set default actions. separate method, in case we want to move it somewhere else
        self.setDefaults()

    def checkAppState(self):
        from AthenaCommon.AppMgr import theApp, AthAppMgr
        if theApp.state() > AthAppMgr.State.CONFIGURED:
            self.log.fatal('Attempt to add a User Action when athena is already initialized. Check you configuration')
            return False
        return True

    # method to add a system action, i.e. an action that is prepended to the list
    # forceoverride bypasses the check on pre-existing system actions
    def addSystemAction(self,actionTool,roles,forceoverride=False):
        if self.checkAppState():
            #print "G4AtlasUAStore: adding system action",actionTool," with roles ",roles

            for role in roles:
                if not role in self.roleMap.keys():
                    self.log.fatal('Unkown role '+role+' requested for action ')
                else:
                    if not actionTool in getattr(self.theSvc, self.roleMap[role]):
                        # allowing only one system action per role
                        if (not role in self.systemActions) or forceoverride :
                            setattr(self.theSvc,self.roleMap[role],[actionTool]+getattr(self.theSvc,self.roleMap[role]))
                            self.systemActions+=[role]
                        else:
                            self.log.fatal('Requesting a system action for '+role+', which already has one. Check your configuration.')

                
    def addAction(self,actionTool,roles):
        if self.checkAppState():
            #print "G4AtlasUAStore: adding action",actionTool," with roles ",roles

            for role in roles:
                if not role in self.roleMap.keys():
                    self.log.fatal('Unkown role '+role+' requested for action ')
                else:
                    if not actionTool in getattr(self.theSvc, self.roleMap[role]):                        
                        setattr(self.theSvc,self.roleMap[role],getattr(self.theSvc,self.roleMap[role])+[actionTool])



    def setDefaults(self):
        
        self.theSvc.BeginOfRunActions=defaultBoRActions()
        self.theSvc.EndOfRunActions=defaultEoRActions()
        self.theSvc.BeginOfEventActions=defaultBoEActions()
        self.theSvc.EndOfEventActions=defaultEoEActions()
        self.theSvc.SteppingActions=defaultSteActions()
        self.theSvc.PreTrackingActions=defaultBoTActions()
        self.theSvc.PostTrackingActions=defaultEoTActions()
        self.theSvc.StackingActionsClassification=defaultStaClaActions()
        self.theSvc.StackingActionsPrepareNewEvent=defaultStaPrepareActions()
        self.theSvc.StackingActionsNewStage=defaultStaNewStageActions()

        # ADS FIXME: is this LAr action needed also for ISF?
        self.addSystemAction('LArHitsEventAction',['BeginOfRun','EndOfRun','BeginOfEvent','EndOfEvent'])
        self.addSystemAction('AthenaStackingAction',['Classification'])

        # add them using internal functions so we keep track of the existence of these system actions
        if not simFlags.ISFRun:
            self.addSystemAction('AthenaTrackingAction',['BeginOfTracking','EndOfTracking'])
            self.addSystemAction('MCTruthSteppingAction',['Step'])

UAStore=G4AtlasUAStore()

