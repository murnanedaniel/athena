
## Enable the looper killer, with options for overlay
#  Don't kill the event, but flag it as bad
try:
    if (hasattr(simFlags, 'UseV2UserActions') and simFlags.UseV2UserActions()):
        # this configures the MT LooperKiller
        from G4UserActions import G4UserActionsConfig
        try:
            G4UserActionsConfig.addLooperKillerTool() #FIXME should be LooperKillerEventOverlay
        except AttributeError:
            atlasG4log.warning("Could not add the MT-version of the LooperKiller")
    else:
        # this configures the non-MT looperKiller
        try:
            from G4AtlasServices.G4AtlasUserActionConfig import UAStore
        except ImportError:
            from G4AtlasServices.UserActionStore import UAStore
        from AthenaCommon.CfgGetter import getPublicToolClone
        # create a special instance of the LooperKiller, with specific configuration
        lkAction = getPublicToolClone("LooperKillerEventOverlay", "LooperKiller", PrintSteps=10, MaxSteps=1000000, VerboseLevel=1, AbortEvent=0, SetError=1)
        UAStore.addAction(lkAction,['Step'])
except:
    # Pre UserAction Migration
    def use_looperkiller():
        from G4AtlasApps import PyG4Atlas, AtlasG4Eng
        lkAction = PyG4Atlas.UserAction('G4UserActions', 'LooperKiller', ['BeginOfRun', 'EndOfRun', 'BeginOfEvent', 'EndOfEvent', 'Step'])
        #lkAction.set_Properties( {"PrintSteps":"0", "MaxSteps":"10", "VerboseLevel":"0", "AbortEvent":"0", "SetError":"1"} ) #to kill most everything
        lkAction.set_Properties( {"PrintSteps":"10", "MaxSteps":"1000000", "VerboseLevel":"1", "AbortEvent":"0", "SetError":"1"} )
        AtlasG4Eng.G4Eng.menu_UserActions.add_UserAction(lkAction)

    simFlags.InitFunctions.add_function("postInit", use_looperkiller)
