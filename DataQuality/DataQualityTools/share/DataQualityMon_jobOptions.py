#**************************************************************
# jopOptions file for Combined Monitoring in Athena
#**************************************************************

#Make m_trigDecTool available:
TrigDecisionTool= monTrigDecTool if DQMonFlags.useTrigger() else "",

# Global monitoring checks to make sure all triggers are firing. The following triggers
# are monitored. Triggers are listed here:
# https://twiki.cern.ch/twiki/bin/viewauth/Atlas/ExpressStream#Physics_pp_v2_menu_collisions
listOfTriggers = ['EF_g20_loose', 'EF_tauNoCut', 'EF_mu15', 'EF_2mu4_Upsimumu',
    'EF_2mu10_loose', 'EF_tauNoCut_L1TAU50', 'EF_2e12_medium', 'EF_mu4_L1J10_matched',
    'EF_tau16_IDTrkNoCut', 'EF_2mu4_Jpsimumu_IDTrkNoCut', 'EF_e20_medium_IDTrkNoCut',
    'EF_L1J10_firstempty_NoAlg', 'EF_L1J30_firstempty_NoAlg', 'EF_rd0_filled_NoAlg',
    'EF_rd0_empty_NoAlg', 'EF_tauNoCut_cosmic', 'EF_j240_a4tc_EFFS', 'EF_fj30_a4tc_EFFS',
    'EF_tau50_IDTrkNoCut', 'EF_xe20_noMu', 'EF_mbMbts_1_eff', 'EF_2e5_tight', 'EF_2mu10',
    'EF_b10_IDTrkNoCut', 'EF_mu15_mu10_EFFS', 'EF_j30_a4tc_EFFS', 'EF_mu20_IDTrkNoCut',
    'EF_InDetMon_FS', 'EF_2mu13_Zmumu_IDTrkNoCut', 'EF_mu20_muCombTag_NoEF']

# Import AlgSequence
from AthenaCommon.AlgSequence import AlgSequence
topSequence = AlgSequence()

from AthenaCommon.AthenaCommonFlags import athenaCommonFlags
isOnline=False
isOffline=True
if athenaCommonFlags.isOnline==True:
    isOnline=True
    isOffline=False


from RecExConfig.RecFlags import rec
if rec.doHeavyIon():
    JetCollectionKey='AntiKt4HIJets'
else:
    JetCollectionKey='AntiKt4EMTopoJets'

# add isolation variables for IsolationSelection
if DQMonFlags.monManEnvironment() == 'tier0ESD' or DQMonFlags.monManEnvironment() == 'tier0':
    if not hasattr(topSequence,"IsolationBuilderNonprompt_All_MaxWeight1000"):
        from IsolationAlgs.IsoUpdatedTrackCones import GetUpdatedIsoTrackCones
        from TrackVertexAssociationTool.getTTVAToolForReco import getTTVAToolForReco
        WP='Nonprompt_All_MaxWeight'
        topSequence += GetUpdatedIsoTrackCones(TTVATool=getTTVAToolForReco(WP,WorkingPoint=WP))

from AthenaCommon.JobProperties import jobproperties
if not 'InDetKeys' in dir():
    from InDetRecExample.InDetKeys import InDetKeys
    
IDTrkContNames = [ InDetKeys.Tracks() ]

isCosmics=False
isBeam=True

if jobproperties.Beam.beamType()=='cosmics':
    isCosmics=True
    isBeam=False

if DQMonFlags.doDataFlowMon():
    from DataQualityTools.DQTDataFlowMonAlg import DQTDataFlowMonAlgConfigOld
    topSequence += DQTDataFlowMonAlgConfigOld(DQMonFlags)

if DQMonFlags.doGlobalMon():
    if DQMonFlags.monManEnvironment != 'tier0ESD':
        # Import Det Synch tool
        if DQMonFlags.monManEnvironment in ('tier0Raw', 'tier0') and globalflags.DataSource.get_Value() != 'geant4':
            from DataQualityTools.DQTDetSynchMonAlg import DQTDetSynchMonAlgConfigOld
            topSequence += DQTDetSynchMonAlgConfigOld(DQMonFlags)

        # Background Monitoring
        if DQMonFlags.useTrigger():
            from DataQualityTools.DQTBackgroundMon import DQTBackgroundMonAlgConfig
            topSequence += DQTBackgroundMonAlgConfig(DQMonFlags,isOld=True)

        # Default values
        MinSCTHits=5
        MinPtCut=4000

        #For now, to increase statistics in cosmics data taking
        if athenaCommonFlags.isOnline==True:
            MinSCTHits=0
            MinPtCut=500

        if not rec.doMuon:
            try:
                svcMgr.ByteStreamAddressProviderSvc.TypeNames.remove("RpcPadContainer/RPCPAD")
            except:
                printfunc ('RPCPAD cannot be removed')

    if isBeam==True and (DQMonFlags.monManEnvironment != 'tier0Raw') and rec.doInDet() and DQMonFlags.useTrigger():

        topSequence += AthenaMonManager( "GlobalMonPhysicsManager" )
        ManagedAthenaGlobalPhysMon = topSequence.GlobalMonPhysicsManager
        ManagedAthenaGlobalPhysMon.FileKey             = DQMonFlags.monManFileKey()
        ManagedAthenaGlobalPhysMon.ManualDataTypeSetup = DQMonFlags.monManManualDataTypeSetup()
        ManagedAthenaGlobalPhysMon.DataType            = DQMonFlags.monManDataType()
        ManagedAthenaGlobalPhysMon.Environment         = DQMonFlags.monManEnvironment()

        from DataQualityTools.DQTGlobalWZFinderAlg import DQTGlobalWZFinderAlgConfig
        topSequence += DQTGlobalWZFinderAlgConfig(DQMonFlags)

        from DataQualityTools.DQTLumiMonAlg import DQTLumiMonAlgConfig
        topSequence += DQTLumiMonAlgConfig(DQMonFlags, isOld=True)
