#
#  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
#

'''@file tauMonitorAlgorithm.py
@date 2018-01-11
@brief tau python configuration for the Run III AthenaMonitoring package
'''

from math import pi

def tauMonitoringConfig(inputFlags,**kwargs):
    '''Function to configures some algorithms in the monitoring system.'''


    if not (inputFlags.Common.isOnline == 'online' or inputFlags.Input.isMC):
        kwargs.setdefault('useReadyFilterTool', True)
    else:
        kwargs.setdefault('useReadyFilterTool', False)

    from AthenaConfiguration.ComponentFactory import CompFactory
    from AthenaMonitoring.BadLBFilterToolConfig import LArBadLBFilterToolCfg
    from AthenaMonitoring.AtlasReadyFilterConfig import AtlasReadyFilterCfg

    ### STEP 1 ###

    # The following class will make a sequence, configure algorithms, and link
    # them to GenericMonitoringTools
    from AthenaMonitoring import AthMonitorCfgHelper
    cfgHelper = AthMonitorCfgHelper(inputFlags, monName='tauMonitorAlgFamily')
    cfg = cfgHelper.result()

    if kwargs['useReadyFilterTool'] and 'ReadyFilterTool' not in kwargs:
        readyFilterTool = cfg.popToolsAndMerge(AtlasReadyFilterCfg(inputFlags))
        kwargs['ReadyFilterTool'] = readyFilterTool


    ### STEP 2 ###
    # Adding an algorithm to the helper. Here, we will use the example 
    # algorithm in the AthenaMonitoring package. Just pass the type to the 
    # helper. Then, the helper will instantiate an instance and set up the 
    # base class configuration following the inputFlags. The returned object 
    # is the algorithm.


    tauMonitorAlgorithm = CompFactory.tauMonitorAlgorithm
    tauMonAlgBA = cfgHelper.addAlgorithm( tauMonitorAlgorithm, name='tauMonAlgBA', addFilterTools = [LArBadLBFilterToolCfg(inputFlags)])
    tauMonAlgCR = cfgHelper.addAlgorithm( tauMonitorAlgorithm, name='tauMonAlgCR', addFilterTools = [LArBadLBFilterToolCfg(inputFlags)])
    tauMonAlgEC = cfgHelper.addAlgorithm( tauMonitorAlgorithm, name='tauMonAlgEC', addFilterTools = [LArBadLBFilterToolCfg(inputFlags)])
    tauMonAlgGlobal = cfgHelper.addAlgorithm( tauMonitorAlgorithm, name='tauMonAlgGlobal', addFilterTools = [LArBadLBFilterToolCfg(inputFlags)])
    tauMonAlgTauTrig1 = cfgHelper.addAlgorithm( tauMonitorAlgorithm, name='tauMonAlgTauTrig1', addFilterTools = [LArBadLBFilterToolCfg(inputFlags)])
    tauMonAlgTauTrig2 = cfgHelper.addAlgorithm( tauMonitorAlgorithm, name='tauMonAlgTauTrig2', addFilterTools = [LArBadLBFilterToolCfg(inputFlags)])
    tauMonAlgTauTrig3 = cfgHelper.addAlgorithm( tauMonitorAlgorithm, name='tauMonAlgTauTrig3', addFilterTools = [LArBadLBFilterToolCfg(inputFlags)])
    tauMonAlgTauTrig4 = cfgHelper.addAlgorithm( tauMonitorAlgorithm, name='tauMonAlgTauTrig4', addFilterTools = [LArBadLBFilterToolCfg(inputFlags)])
    tauMonAlgTauTrig5 = cfgHelper.addAlgorithm( tauMonitorAlgorithm, name='tauMonAlgTauTrig5', addFilterTools = [LArBadLBFilterToolCfg(inputFlags)])
    tauMonAlgTauTrig6 = cfgHelper.addAlgorithm( tauMonitorAlgorithm, name='tauMonAlgTauTrig6', addFilterTools = [LArBadLBFilterToolCfg(inputFlags)])
    tauMonAlgTauTrig7 = cfgHelper.addAlgorithm( tauMonitorAlgorithm, name='tauMonAlgTauTrig7', addFilterTools = [LArBadLBFilterToolCfg(inputFlags)])
    tauMonAlgEleTrig = cfgHelper.addAlgorithm( tauMonitorAlgorithm, name='tauMonAlgEleTrig', addFilterTools = [LArBadLBFilterToolCfg(inputFlags)])
    tauMonAlgJetTrig = cfgHelper.addAlgorithm( tauMonitorAlgorithm, name='tauMonAlgJetTrig', addFilterTools = [LArBadLBFilterToolCfg(inputFlags)])




    ### STEP 3 ###
    # Edit properties of a algorithm
    # some generic property
    # exampleMonAlg.RandomHist = True
    # to enable a trigger filter, for example:
    #exampleMonAlg.TriggerChain = 'HLT_mu26_ivarmedium'

    tauMonAlgTauTrig1.TriggerChain = "HLT_tau25_mediumRNN_tracktwoMVABDT_L1TAU12IM"
    tauMonAlgTauTrig2.TriggerChain = "HLT_tau160_mediumRNN_tracktwoMVABDT_L1TAU100"
    tauMonAlgTauTrig3.TriggerChain = "HLT_tau25_mediumRNN_tracktwoMVABDT_L1eTAU20M"
    tauMonAlgTauTrig4.TriggerChain = "HLT_tau160_mediumRNN_tracktwoMVABDT_L1eTAU140"
    tauMonAlgTauTrig5.TriggerChain = "HLT_tau80_mediumRNN_tracktwoMVA_tau35_mediumRNN_tracktwoMVA_03dRAB30_L1TAU60_DR-TAU20ITAU12I"
    tauMonAlgTauTrig6.TriggerChain = "HLT_tau35_mediumRNN_tracktwoMVABDT_tau25_mediumRNN_tracktwoMVABDT_03dRAB30_L1DR-TAU20ITAU12I-J25"
    tauMonAlgTauTrig7.TriggerChain = "HLT_tau200_mediumRNN_tracktwoMVABDT_L1TAU100"

    tauMonAlgEleTrig.TriggerChain = "HLT_e[2-9][0-9]_.*"
    tauMonAlgJetTrig.TriggerChain = "HLT_j[2-9][0-9]_.*"




    tauMonAlgBA.etaMin = -1.
    tauMonAlgBA.etaMax = 1.3
    tauMonAlgCR.etaMin = 1.3
    tauMonAlgCR.etaMax = 1.7
    tauMonAlgEC.etaMin = 1.7
    tauMonAlgEC.etaMax = 3.5
    tauMonAlgGlobal.etaMin = -100
    tauMonAlgGlobal.etaMax = 100
    tauMonAlgTauTrig1.etaMin = -100
    tauMonAlgTauTrig1.etaMax = 100
    tauMonAlgTauTrig2.etaMin = -100
    tauMonAlgTauTrig2.etaMax = 100
    tauMonAlgTauTrig3.etaMin = -100
    tauMonAlgTauTrig3.etaMax = 100
    tauMonAlgTauTrig4.etaMin = -100
    tauMonAlgTauTrig4.etaMax = 100
    tauMonAlgTauTrig5.etaMin = -100
    tauMonAlgTauTrig5.etaMax = 100
    tauMonAlgTauTrig6.etaMin = -100
    tauMonAlgTauTrig6.etaMax = 100
    tauMonAlgTauTrig7.etaMin = -100
    tauMonAlgTauTrig7.etaMax = 100

    tauMonAlgEleTrig.etaMin = -100
    tauMonAlgEleTrig.etaMax = 100
    tauMonAlgJetTrig.etaMin = -100
    tauMonAlgJetTrig.etaMax = 100

    tauMonAlgBA.kinGroupName = 'tauMonKinGroupBA'
    tauMonAlgCR.kinGroupName = 'tauMonKinGroupCR'
    tauMonAlgEC.kinGroupName = 'tauMonKinGroupEC'
    tauMonAlgGlobal.kinGroupName = 'tauMonKinGroupGlobal'
    tauMonAlgTauTrig1.kinGroupName = 'tauMonKinGroupTauTrig1'
    tauMonAlgTauTrig2.kinGroupName = 'tauMonKinGroupTauTrig2'
    tauMonAlgTauTrig3.kinGroupName = 'tauMonKinGroupTauTrig3'
    tauMonAlgTauTrig4.kinGroupName = 'tauMonKinGroupTauTrig4'
    tauMonAlgTauTrig5.kinGroupName = 'tauMonKinGroupTauTrig5'
    tauMonAlgTauTrig6.kinGroupName = 'tauMonKinGroupTauTrig6'
    tauMonAlgTauTrig7.kinGroupName = 'tauMonKinGroupTauTrig7'

    tauMonAlgEleTrig.kinGroupName = 'tauMonKinGroupEleTrig'
    tauMonAlgJetTrig.kinGroupName = 'tauMonKinGroupJetTrig'

    ### STEP 4 ###
    # Add some tools. N.B. Do not use your own trigger decion tool. Use the
    # standard one that is included with AthMonitorAlgorithm.

    # Add a generic monitoring tool (a "group" in old language). The returned 
    # object here is the standard GenericMonitoringTool.
    myKinGroupBA = cfgHelper.addGroup(alg=tauMonAlgBA, name='tauMonKinGroupBA', topPath='Tau/TauB/' )
    myKinGroupCR = cfgHelper.addGroup(alg=tauMonAlgCR, name='tauMonKinGroupCR', topPath='Tau/TauCR/' )
    myKinGroupEC = cfgHelper.addGroup(alg=tauMonAlgEC, name='tauMonKinGroupEC', topPath='Tau/TauE/' )
    myKinGroupGlobal = cfgHelper.addGroup(alg=tauMonAlgGlobal, name='tauMonKinGroupGlobal', topPath='Tau/' )
    myKinGroupTauTrig1 = cfgHelper.addGroup(alg=tauMonAlgTauTrig1, name='tauMonKinGroupTauTrig1', topPath='Tau/Trigger/tauTrigger1/' )
    myKinGroupTauTrig2 = cfgHelper.addGroup(alg=tauMonAlgTauTrig2, name='tauMonKinGroupTauTrig2', topPath='Tau/Trigger/tauTrigger2/' )
    myKinGroupTauTrig3 = cfgHelper.addGroup(alg=tauMonAlgTauTrig3, name='tauMonKinGroupTauTrig3', topPath='Tau/Trigger/tauTrigger3/' )
    myKinGroupTauTrig4 = cfgHelper.addGroup(alg=tauMonAlgTauTrig4, name='tauMonKinGroupTauTrig4', topPath='Tau/Trigger/tauTrigger4/' )
    myKinGroupTauTrig5 = cfgHelper.addGroup(alg=tauMonAlgTauTrig5, name='tauMonKinGroupTauTrig5', topPath='Tau/Trigger/tauTrigger5/' )
    myKinGroupTauTrig6 = cfgHelper.addGroup(alg=tauMonAlgTauTrig6, name='tauMonKinGroupTauTrig6', topPath='Tau/Trigger/tauTrigger6/' )
    myKinGroupTauTrig7 = cfgHelper.addGroup(alg=tauMonAlgTauTrig7, name='tauMonKinGroupTauTrig7', topPath='Tau/Trigger/tauTrigger7/' )
    myKinGroupEleTrig = cfgHelper.addGroup(alg=tauMonAlgEleTrig, name='tauMonKinGroupEleTrig', topPath='Tau/Trigger/EleTrig/' )
    myKinGroupJetTrig = cfgHelper.addGroup(alg=tauMonAlgJetTrig, name='tauMonKinGroupJetTrig', topPath='Tau/Trigger/JetTrig/' )


    naming= {
            'BA': "Tau_TauB_",
            'CR': "Tau_TauCR_",
            'EC': "Tau_TauE_",
            'Global': "",
            'EleTrig': "emTriggered_",
            'JetTrig': "jetTriggered_",
            'TauTrig1': "tauTriggered1_",
            'TauTrig2': "tauTriggered2_",
            'TauTrig3': "tauTriggered3_",
            'TauTrig4': "tauTriggered4_",
            'TauTrig5': "tauTriggered5_",
            'TauTrig6': "tauTriggered6_",
            'TauTrig7': "tauTriggered7_",

            }

    def namer(variables,histogramName, folder, algorithm):
        if folder != "":
            #This wouldn't have been necessary except the calo folder in the naming convention:
            #This adds the Tau prefix to the histograms in the calo folder
            if postfix == "Global" and (folder =="Calo" or folder=="Track" or folder=="Track_allTracks" or folder=="Track_leadTrack"  or folder=="Cluster" ):
                return variables+";"+"Tau_"+naming[algorithm]+folder+"_"+histogramName
            else:
                return variables+";"+naming[algorithm]+folder+"_"+histogramName
        else:
            return variables+";"+naming[algorithm]+folder+histogramName

    # Add a GMT for the other example monitor algorithm
    # anotherGroup = cfgHelper.addGroup(anotherTauMonAlg,'tauMonitor')

    ### STEP 5 ###
    # Configure histograms
    for itup in [(myKinGroupBA,'BA'),
                 (myKinGroupCR,'CR'),
                 (myKinGroupEC,'EC'),
                 (myKinGroupGlobal,'Global'),
                 (myKinGroupEleTrig,'EleTrig'),
                 (myKinGroupJetTrig,'JetTrig'),
                 (myKinGroupTauTrig1, 'TauTrig1'),
                 (myKinGroupTauTrig2, 'TauTrig2'),
                 (myKinGroupTauTrig3, 'TauTrig3'),
                 (myKinGroupTauTrig4, 'TauTrig4'),
                 (myKinGroupTauTrig5, 'TauTrig5'),
                 (myKinGroupTauTrig6, 'TauTrig6'),
                 (myKinGroupTauTrig7, 'TauTrig7')]:

        (igroup, postfix) = itup

        PHIMIN = -pi
        PHIMAX = pi

        folder = ""

        if(postfix =="BA" or postfix =="CR" or postfix=="EC" or postfix == "Global" or postfix.startswith('TauTrig') or postfix=="EleTrig" or postfix =="JetTrig"):


            #potentialHigh Pt replacement
            igroup.defineHistogram(namer('tauEtaEt15,tauPhiEt15','tauPhiVsEta_et15','',postfix), type='TH2F', title='Eta Vs phi (Et>15GeV);Eta;Phi', 
               xbins=30,xmin=-2.55,xmax=2.55,ybins=32,ymin=PHIMIN,ymax=PHIMAX)

            igroup.defineHistogram(namer('nTauCandidates', 'nTauCandidates',"",postfix), title='Number of tau candidates;Number of Taus per Event', 
                                   xbins=61, xmin=-0.5, xmax=60.5,path=folder)


            igroup.defineHistogram(namer('tauCharge','tauCharge','',postfix), title='Charge of tau candidates;Charge;Number of Candidates',
                                   xbins=11, xmin=-5.5, xmax=5.5,path=folder)


            igroup.defineHistogram(namer('RNNJetScore','RNNJetScore','',postfix), title='RNN Jet Score',
                                   xbins=100, xmin=0, xmax=1,path=folder)

            igroup.defineHistogram(namer('RNNJetScoreSigTrans','RNNJetScoreSigTrans','',postfix), title='RNN Jet Score Sig Trans',
                                   xbins=48, xmin=0, xmax=1.1,path=folder)



            igroup.defineHistogram(namer('tauEt','tauEt','',postfix), title='Et of tau candidates;Transverse Energy (GeV);Number of Candidates',
                                   xbins=60, xmin=0., xmax=300.,path=folder)

            igroup.defineHistogram(namer('tauEtaEXP','tauEtaEXP','',postfix), title='Eta of tau candidates;Eta;Number of Candidates',
                                   xbins=51, xmin=-2.55, xmax=2.55,path=folder)


            
            igroup.defineHistogram(namer('tauEta','tauEta','',postfix), title='Eta of tau candidates;Eta;Number of Candidates',
                                   xbins=51, xmin=-2.55, xmax=2.55,path=folder)
            
            igroup.defineHistogram(namer('LB,tauEta','tauEtaVsLB','',postfix), type='TH2F', title='Tau Eta Vs Lumiblock;Lumiblock;eta', 
                                   xbins=1200,xmin=0,xmax=1200,ybins=51,ymin=-2.55,ymax=2.55,path=folder)
            
            igroup.defineHistogram(namer('tauPhi','tauPhi','',postfix), title='Phi of tau candidates;Phi;Number of Candidates',
                                   xbins=65, xmin=PHIMIN-0.098174/2., xmax=PHIMAX+0.098174/2.,path=folder)
            
            igroup.defineHistogram(namer('tauEta,tauPhi','tauPhiVsEta','',postfix), type='TH2F', title='Phi Vs Eta;Eta;Phi', 
                                   xbins=30,xmin=-2.55,xmax=2.55,ybins=32,ymin=PHIMIN,ymax=PHIMAX,path=folder)

            
            igroup.defineHistogram(namer('LB,tauPhi','tauPhiVsLB','',postfix), type='TH2F', title='Tau Phi vs Lumiblock;Lumiblock;Phi', 
                                   xbins=1200,xmin=0,xmax=1200,ybins=80,ymin=PHIMIN,ymax=PHIMAX,path=folder)

        if(postfix =="BA" or postfix =="CR" or postfix=="EC" or postfix.startswith('TauTrig')):

            igroup.defineHistogram(namer('tauPhiRNNLoose','phi','Identification_RNNLoose',postfix), title='Phi of tau candidates ( RNNLoose) ;Phi;Number of Candidates',
                xbins=65, xmin=PHIMIN-0.098174/2., xmax=PHIMAX+0.098174/2., path=folder+"Identification/RNNLoose" )

            igroup.defineHistogram(namer('tauEtaRNNLoose','eta','Identification_RNNLoose',postfix), title='Eta of tau candidates ( RNNLoose) ;Eta;Number of Candidates',
                xbins=51, xmin=-2.55, xmax=2.55, path=folder+"Identification/RNNLoose")

            igroup.defineHistogram(namer('tauEtRNNLoose','et','Identification_RNNLoose',postfix), title='Et of tau candidates;Transverse Energy (GeV);Number of Candidates',
                xbins=60, xmin=0., xmax=300.,path=folder+"Identification/RNNLoose")

            igroup.defineHistogram(namer('NumTracksRNNLoose','NumTracks','Identification_RNNLoose',postfix), title='Number Of Tracks for Tau Candidates (RNNLoose);Number Of Tracks;Number Of Candidates',
                xbins=21, xmin=-0.5, xmax=20.5,path=folder+"Identification/RNNLoose")

            igroup.defineHistogram(namer('tauPhiRNNMedium','phi','Identification_RNNMedium',postfix), title='Phi of tau candidates ( RNNMedium) ;Phi;Number of Candidates',
                xbins=65, xmin=PHIMIN-0.098174/2., xmax=PHIMAX+0.098174/2., path=folder+"Identification/RNNMedium" )

            igroup.defineHistogram(namer('tauEtaRNNMedium','eta','Identification_RNNMedium',postfix), title='Eta of tau candidates ( RNNMedium) ;Eta;Number of Candidates',
                xbins=51, xmin=-2.55, xmax=2.55, path=folder+"Identification/RNNMedium")

            igroup.defineHistogram(namer('tauEtRNNMedium','et','Identification_RNNMedium',postfix), title='Et of tau candidates;Transverse Energy (GeV);Number of Candidates',
                xbins=60, xmin=0., xmax=300.,path=folder+"Identification/RNNMedium")

            igroup.defineHistogram(namer('NumTracksRNNMedium','NumTracks','Identification_RNNMedium',postfix), title='Number Of Tracks for Tau Candidates (RNNMedium);Number Of Tracks;Number Of Candidates',
                xbins=21, xmin=-0.5, xmax=20.5,path=folder+"Identification/RNNMedium")

            igroup.defineHistogram(namer('tauPhiRNNTight','phi','Identification_RNNTight',postfix), title='Phi of tau candidates ( RNNTight) ;Phi;Number of Candidates',
                xbins=65, xmin=PHIMIN-0.098174/2., xmax=PHIMAX+0.098174/2., path=folder+"Identification/RNNTight" )

            igroup.defineHistogram(namer('tauEtaRNNTight','eta','Identification_RNNTight',postfix), title='Eta of tau candidates ( RNNTight) ;Eta;Number of Candidates',
                xbins=51, xmin=-2.55, xmax=2.55, path=folder+"Identification/RNNTight")

            igroup.defineHistogram(namer('tauEtRNNTight','et','Identification_RNNTight',postfix), title='Et of tau candidates;Transverse Energy (GeV);Number of Candidates',
                xbins=60, xmin=0., xmax=300.,path=folder+"Identification/RNNTight")

            igroup.defineHistogram(namer('NumTracksRNNTight','NumTracks','Identification_RNNTight',postfix), title='Number Of Tracks for Tau Candidates (RNNTight);Number Of Tracks;Number Of Candidates',
                xbins=21, xmin=-0.5, xmax=20.5,path=folder+"Identification/RNNTight")




            igroup.defineHistogram(namer('tauPhiEt15RNNLoose','phi','Identification_RNNLoose15GeV',postfix), title='Phi of tau candidates (Et>15, RNNLoose) ;Phi;Number of Candidates',
                xbins=65, xmin=PHIMIN-0.098174/2., xmax=PHIMAX+0.098174/2., path=folder+"Identification/RNNLoose15GeV" )

            igroup.defineHistogram(namer('tauEtaEt15RNNLoose','eta','Identification_RNNLoose15GeV',postfix), title='Eta of tau candidates (Et>15, RNNLoose) ;Eta;Number of Candidates',
                xbins=51, xmin=-2.55, xmax=2.55, path=folder+"Identification/RNNLoose15GeV")

            igroup.defineHistogram(namer('nClustersEt15RNNLoose','nCluster','Identification_RNNLoose15GeV',postfix), title='Number Of CaloTopoClusters (Et>15,RNNLoose);Number Of Clusters;Number Of Candidates',
                xbins=40, xmin=0., xmax=40.,path=folder+"Identification/RNNLoose15GeV" )

            igroup.defineHistogram(namer('NumTracksEt15RNNLoose','NumTracks','Identification_RNNLoose15GeV',postfix), title='Number Of Tracks for Tau Candidates (Et>15,RNNLoose);Number Of Tracks;Number Of Candidates',
                xbins=21, xmin=-0.5, xmax=20.5,path=folder+"Identification/RNNLoose15GeV")

            igroup.defineHistogram(namer('tauEtEt15RNNLoose','et','Identification_RNNLoose15GeV',postfix), title='Et of tau candidates;Transverse Energy (GeV);Number of Candidates',
                                   xbins=60, xmin=0., xmax=300.,path=folder+"Identification/RNNLoose15GeV")

            igroup.defineHistogram(namer('panModeEt15RNNLoose','panMode','Identification_RNNLoose15GeV',postfix), title='tau decay mode from panTau upon JetRNNSigMedium;mode',
                                   xbins=5, xmin=0., xmax=5., path=folder+"Identification/RNNLoose15GeV", xlabels=["1p0n","1p1n","1pXn","3p0n","3pXn"])

            igroup.defineHistogram(namer('tauEtaEt15RNNLoose,tauPhiEt15RNNLoose','tauPhiVsEta_et15_RNNLoose','',postfix), type='TH2F', title='Phi vs Eta (Et>15, RNNLoose) ;Eta;Phi', 
               xbins=30,xmin=-2.55,xmax=2.55,ybins=32,ymin=PHIMIN,ymax=PHIMAX)


            igroup.defineHistogram(namer('jetSeedEta','jetSeedEta','Calo',postfix), title='Calorimeter eta of tau candidates;Eta;Numbers of Candidates',path=folder+"Calo",
            xbins=50, xmin=-2.5, xmax=2.5 )

            igroup.defineHistogram(namer('jetSeedPhi','jetSeedPhi','Calo',postfix), title='Calorimeter phi of tau candidates;Phi;Numbers of Candidates',path=folder+"Calo",
            xbins=64, xmin=PHIMIN, xmax=PHIMAX )

            igroup.defineHistogram(namer('jetSeedPt','jetSeedPt','Calo',postfix), title='Calorimeter pt of tau candidates;Et;Numbers of Candidates',path=folder+"Calo",
            xbins=100, xmin=0., xmax=300 )
            
            igroup.defineHistogram(namer('centFrac,LB','CentFracVsLB','Calo',postfix), type='TH2F', title='Centrality Fraction vs Lumiblock;Centrality Fraction;Lumiblock', path=folder+"Calo", 
            xbins=51,xmin=0,xmax=1.02,ybins=1200,ymin=0.,ymax=1200.)

            igroup.defineHistogram('centFrac,LB', type='TH2F', title='Centrality Fraction vs Lumiblock;Centrality Fraction;Lumiblock', path=folder+"Calo", 
            xbins=51,xmin=0,xmax=1.02,ybins=1200,ymin=0.,ymax=1200.)

            igroup.defineHistogram('isolFrac,LB', type='TH2F', title='Isolation Fraction vs Lumiblock;Isolation Fraction;Lumiblock', path=folder+"Calo", 
            xbins=51,xmin=0,xmax=1.02,ybins=1200,ymin=0.,ymax=1200.)



            igroup.defineHistogram(namer('muonVeto','muonVeto','Identification',postfix), title='Muon Veto',path=folder+"Identification",
            xbins=2, xmin=-0.5, xmax=1.5, xlabels=["False","True"] )

            igroup.defineHistogram(namer('tauRNNLoose','tauRNNLoose','Identification',postfix), title='Identification Flag: tauRNNLoose',path=folder+"Identification",
            xbins=2, xmin=-0.5, xmax=1.5 , xlabels=["False","True"])

            igroup.defineHistogram(namer('tauRNNMedium','tauRNNMedium','Identification',postfix), title='Identification Flag: tauRNNMedium',path=folder+"Identification",
            xbins=2, xmin=-0.5, xmax=1.5 , xlabels=["False","True"])

            igroup.defineHistogram(namer('tauRNNTight','tauRNNTight','Identification',postfix), title='Identification Flag: tauRNNTight',path=folder+"Identification",
            xbins=2, xmin=-0.5, xmax=1.5, xlabels=["False","True"])


            #Substructure Histograms
            igroup.defineHistogram(namer('panModeSubstructure','panMode','SubStructure',postfix), title='tau decay mode from PanTau upon JetBDTSigMedium;mode',
            xbins=8, xmin=0., xmax=8., path=folder+"SubStructure")

            igroup.defineHistogram(namer('BDTScoreAsP0','BDTScoreAsP0','SubStructure',postfix), title='BDT score indentifying pion zero ; BDT score',path=folder+"SubStructure",
            xbins=15, xmin=-0, xmax=1.2 )

            igroup.defineHistogram(namer('EMFracTrk','EMFracTrk','SubStructure',postfix), title='Ratio of pt to shot electromagnetic energy for associated tracks; track pt ratio in EM',path=folder+"SubStructure",
            xbins=15, xmin=0, xmax=1.5 )

            igroup.defineHistogram(namer('EfracL2EffCluster','EfracL2EffCluster','SubStructure',postfix), title='Energy fraction of leading two effective clusters in shot; energy fraction',path=folder+"SubStructure",
            xbins=15, xmin=0, xmax=1.5 )

            igroup.defineHistogram(namer('EisoEffCluster','EisoEffCluster','SubStructure',postfix), title='Isolation Energy after correction in effective clusters ; isolation energy (GeV)',path=folder+"SubStructure",
            xbins=10, xmin=0, xmax=50. )

            igroup.defineHistogram(namer('InvMassEffClusters','InvMassEffClusters','SubStructure',postfix), title='Invariant mass of effective clusters in shot; invariant mass (GeV)',path=folder+"SubStructure",
            xbins=40, xmin=0, xmax=8. )

            igroup.defineHistogram(namer('nNeutPFO','nNeutPFO','SubStructure',postfix), title='_NumNeutPFO;Number of neutral ParticleFlow objects ; PFO number',path=folder+"SubStructure",
            xbins=20, xmin=0, xmax=20. )

            igroup.defineHistogram(namer('nShot','nShot','SubStructure',postfix), title='number of shots ; shot number ',path=folder+"SubStructure",
            xbins=20, xmin=0, xmax=20. )

            igroup.defineHistogram(namer('panEta','panEta','SubStructure',postfix), title='tau Eta from PanTau ; substructure Eta',path=folder+"SubStructure",
            xbins=16, xmin=-3.2, xmax=3.2 )

            igroup.defineHistogram(namer('panPhi','panPhi','SubStructure',postfix), title='tau Phi from PanTau ; substructure Phi',path=folder+"SubStructure",
            xbins=16, xmin=PHIMIN, xmax=PHIMAX )

            igroup.defineHistogram(namer('panPt','panPt','SubStructure',postfix), title='tau Pt from PanTau ; substructure pt (GeV)',path=folder+"SubStructure",
            xbins=20, xmin=0, xmax=200 )

            igroup.defineHistogram(namer('pt3','shots_pt3','SubStructure',postfix), title='weighted cell pt in 3x3 window in shots; pt3 (GeV) ',path=folder+"SubStructure",
            xbins=24, xmin=0, xmax=12 )

            igroup.defineHistogram(namer('PSSFrac','PSSFracEffCluster','SubStructure',postfix), title='Energy fraction for PreSampling and sampling layers in effective clusters in shot;Sampling Energy Fraction',path=folder+"SubStructure",
            xbins=10, xmin=0, xmax=1. )

            igroup.defineHistogram(namer('ptRatioEflowApprox','ptRatioApprox','SubStructure',postfix), title='Ratio of pt to shot total energy for associated tracks; track pt ratio',
                xbins=20, xmin=0., xmax=2., path=folder+"SubStructure" )





        
        if(postfix =="BA" or postfix =="CR" or postfix=="EC" or postfix.startswith('TauTrig') or postfix=="Global"):
            #Cluster Histograms
            igroup.defineHistogram(namer('clusterLogEt','logEt','Cluster',postfix), title='log(cluster E_{T});log(cluster E_{T});Entries',path=folder+"Cluster",
            xbins=50, xmin=1.5, xmax=6.5 )
            igroup.defineHistogram(namer('clusterEta','eta','Cluster',postfix), title='cluster eta;cluster eta;Entries',path=folder+"Cluster",
            xbins=50, xmin=-2.6, xmax=2.6 )
            igroup.defineHistogram(namer('clusterPhi','phi','Cluster',postfix), title='cluster phi;cluster phi;Entries',path=folder+"Cluster",
            xbins=32, xmin=-3.2, xmax=3.2 )
            igroup.defineHistogram(namer('clusterSecondR','secondR','Cluster',postfix), title='cluster SECOND_R;cluster SECOND_R;Entries',path=folder+"Cluster",
            xbins=50, xmin=-1.0, xmax=6.0 )
            igroup.defineHistogram(namer('clusterSecondLambda','secondLambda','Cluster',postfix), title='cluster SECOND_LAMBDA;cluster SECOND_LAMBDA;Entries',path=folder+"Cluster",
            xbins=50, xmin=2.0, xmax=6.5 )
            igroup.defineHistogram(namer('clusterCenterLambda','centerLambda','Cluster',postfix), title='cluster CENTER_LAMBDA;cluster CENTER_LAMBDA;Entries',path=folder+"Cluster",
            xbins=50, xmin=1.0, xmax=4.0 )
            igroup.defineHistogram(namer('clusterFirstEngDens','firstEngDens','Cluster',postfix), title='cluster FIRST_ENG_DENS;cluster FIRST_ENG_DENS;Entries',path=folder+"Cluster",
            xbins=50, xmin=-6.0, xmax=1.0 )
            igroup.defineHistogram(namer('clusterEMproba','EMproba','Cluster',postfix), title='cluster EM_PROBABILITY;cluster EM_PROBABILITY;Entries',path=folder+"Cluster",
            xbins=20, xmin=0.0, xmax=1.0 )

            igroup.defineHistogram(namer('clustersMeanCenterLambda','meanCenterLambda','Cluster',postfix), title='clustersMeanCenterLambda;clustersMeanCenterLambda;Entries',path=folder+"Cluster",
            xbins=40, xmin=0., xmax=2500.0 )
            igroup.defineHistogram(namer('clustersMeanFirstEngDens','meanFirstEngDens','Cluster',postfix), title='clustersMeanFirstEngDens;clustersMeanFirstEngDens;Entries',path=folder+"Cluster",
            xbins=40, xmin=-8.7, xmax=-5.5 )
            igroup.defineHistogram(namer('clustersMeanEMProbability','meanEMProbability','Cluster',postfix), title='clustersMeanEMProbability;clustersMeanEMProbability;Entries',path=folder+"Cluster",
            xbins=20, xmin=0.0, xmax=1.0 )
            igroup.defineHistogram(namer('clustersMeanSecondLambda','meanSecondLambda','Cluster',postfix), title='clustersMeanSecondLambda;clustersMeanSecondLambda;Entries',path=folder+"Cluster",
            xbins=20, xmin=0., xmax=600000 )
            igroup.defineHistogram(namer('clustersMeanPresamplerFrac','meanPresamplerFrac','Cluster',postfix), title='clustersMeanPresamplerFrac;clustersMeanPresamplerFrac;Entries',path=folder+"Cluster",
            xbins=20, xmin=0., xmax=0.2 )
            igroup.defineHistogram(namer('clustersPFOEngRelDiff','PFOEngRelDiff','Cluster',postfix), title='PFOEngRelDiff;PFOEngRelDiff;Entries',path=folder+"Cluster",
            xbins=20, xmin=-1.0, xmax=1.0 )
            #Calo Histograms
            igroup.defineHistogram(namer('EMRadius','EMRadius','Calo',postfix), title='Uncalibrated EM Radius;EM Radius;Number Of Candidates',path=folder+"Calo",
            xbins=50, xmin=0., xmax=1.2 )

            igroup.defineHistogram(namer('hadRadius','hadRadius','Calo',postfix), title='Hadronic Radius of tau candidates;Hadronic Radius; Number Of Candidates',path=folder+"Calo",
            xbins=50, xmin=0., xmax=1. )

            igroup.defineHistogram(namer('stripWidth2','stripWidth2','Calo',postfix), title='Strip Width of tau candidates;Strip Width;Number of Candidates',path=folder+"Calo",
            xbins=50, xmin=-0.1, xmax=0.12 )

            igroup.defineHistogram(namer('isolFrac','isolFrac','Calo',postfix), title='Isolation Fraction;Et Isolation Fraction;Number of Candidates',path=folder+"Calo",
            xbins=51, xmin=0.0, xmax=1.02 )

            igroup.defineHistogram(namer('nStrip','nStrip','Calo',postfix), title='Number of strip cells of tau candidates;Number of Strip Cells;Number of Candidates',path=folder+"Calo",
            xbins=56, xmin=-0.5, xmax=55.5 )

            igroup.defineHistogram(namer('etEMAtEMScale','etEMAtEMScale','Calo',postfix), title='EM energy at the EM scale;EM Et (GeV) );Number of Candidates',path=folder+"Calo",
            xbins=50, xmin=0., xmax=200. )

            igroup.defineHistogram(namer('etHadAtEMScale','etHadAtEMScale','Calo',postfix), title='Hadronic Energy at the EM Scale;Had Et (GeV) );Number of Candidates',path=folder+"Calo",
            xbins=50, xmin=0., xmax=200. )
            
            igroup.defineHistogram(namer('centFrac,LB','centFracVsLB','Calo',postfix), type='TH2F', title='Centrality Fraction vs Lumiblock;Centrality Fraction;Lumiblock', path=folder+"Calo", 
            xbins=51,xmin=0,xmax=1.02,ybins=1200,ymin=0.,ymax=1200.)

            igroup.defineHistogram(namer('isolFrac,LB','isolFracVsLB','Calo',postfix), type='TH2F', title='Isolation Fraction vs Lumiblock;Isolation Fraction;Lumiblock', path=folder+"Calo", 
            xbins=51,xmin=0,xmax=1.02,ybins=1200,ymin=0.,ymax=1200.)

            igroup.defineHistogram(namer('centFrac',"centralityFraction","Calo",postfix), title='Centrality Fraction;Centrality Fraction;Number Of Candidates',
                xbins=70, xmin=0., xmax=1.4, path=folder+"Calo")

            igroup.defineHistogram(namer('dRmax',"dRmax","Calo",postfix), title='Maximum delta R;delta R;Number of Candidates',
                xbins=42, xmin=0., xmax=0.21, path=folder+"Calo" )

            igroup.defineHistogram(namer('tauEta','eta','Calo',postfix), title='Eta of tau candidates;Eta;Number of Candidates',
                xbins=51, xmin=-2.55, xmax=2.55,path=folder+"Calo")

            igroup.defineHistogram(namer('tauPhi','phi','Calo',postfix), title='Phi of tau candidates;Phi;Number of Candidates',
                                   xbins=65, xmin=PHIMIN, xmax=PHIMAX,path=folder+"Calo")

            igroup.defineHistogram(namer('trkAvgDist','trkAvgDist','Calo',postfix), title='Average Track Distance from Calorimeter Seed;Distance (mm);Number of Candidates',
                xbins=20, xmin=0., xmax=1., path=folder+"Calo" )

            #Track Histograms!

            #track variables used for track classification and RNN tau ID
            igroup.defineHistogram(namer('trackLogSeedJetPt','logSeedJetPt','Track_allTracks',postfix), title='log(ptJetSeed);log(ptJetSeed);Entries',path=folder+"Track/allTracks",
            xbins=50, xmin=10.0, xmax=14.5 )
            igroup.defineHistogram(namer('trackLogPt','logPt','Track_allTracks',postfix), title='track log(p_{T});track log(p_{T});Entries',path=folder+"Track/allTracks",
            xbins=50, xmin=6.0, xmax=15.0 )
            igroup.defineHistogram(namer('trackEta','eta','Track_allTracks',postfix), title='track eta;track eta;Entries',path=folder+"Track/allTracks",
            xbins=50, xmin=-2.6, xmax=2.6 )
            igroup.defineHistogram(namer('trackd0TJVA','d0TJVA','Track_allTracks',postfix), title='track tanh(d0TJVA);track tanh(d0TJVA);Entries',path=folder+"Track/allTracks",
            xbins=50, xmin=-0.1, xmax=0.1 )
            igroup.defineHistogram(namer('trackZ0SinthetaTJVA','Z0SinthetaTJVA','Track_allTracks',postfix), title='track z0SinthetaTJVA;track z0SinthetaTJVA;Entries',path=folder+"Track/allTracks",
            xbins=50, xmin=-100.0, xmax=100.0 )
            igroup.defineHistogram(namer('trackD0SigTJVA','d0SigTJVA','Track_allTracks',postfix), title='track d0SigTJVA;track d0SigTJVA;Entries',path=folder+"Track/allTracks",
            xbins=50, xmin=-10.0, xmax=10.0 )
            igroup.defineHistogram(namer('trackZ0sinthetaSigTJVA','Z0sinthetaSigTJVA','Track_allTracks',postfix), title='track z0sinthetaSigTJVA;track z0sinthetaSigTJVA;Entries',path=folder+"Track/allTracks",
            xbins=50, xmin=-10.0, xmax=10.0 )
            igroup.defineHistogram(namer('trackCharge','charge','Track_allTracks',postfix), title='track charge;track charge;Entries',path=folder+"Track/allTracks",
            xbins=8, xmin=-4.0, xmax=4.0 )
            igroup.defineHistogram(namer('trackqOverP','qOverP','Track_allTracks',postfix), title='track qOverP;track qOverP;Entries',path=folder+"Track/allTracks",
            xbins=50, xmin=-2.0, xmax=2.0 )
            igroup.defineHistogram(namer('trackLogRConv','logRConv','Track_allTracks',postfix), title='track log(rConv);track log(rConv);Entries',path=folder+"Track/allTracks",
            xbins=50, xmin=0.0, xmax=7.0 )
            igroup.defineHistogram(namer('trackTanhRConvII','tanhRConvII','Track_allTracks',postfix), title='track tanh(rConvII/500);track tanh(rConvII/500);Entries',path=folder+"Track/allTracks",
            xbins=50, xmin=-0.4, xmax=0.4 )
            igroup.defineHistogram(namer('trackPtRatioSeedJet','ptRatioSeedJet','Track_allTracks',postfix), title='(track p_{T})/(ptJetSeed);(track p_{T})/(ptJetSeed);Entries',path=folder+"Track/allTracks",
            xbins=50, xmin=0.0, xmax=0.1 )
            igroup.defineHistogram(namer('trackdRJetSeedAxis','dRJetSeedAxis','Track_allTracks',postfix), title='#DeltaR(tau,track);#DeltaR(tau,track);Entries',path=folder+"Track/allTracks",
            xbins=50, xmin=0.0, xmax=0.5 )
            igroup.defineHistogram(namer('trackNInnermostPixHits','nInnermostPixHits','Track_allTracks',postfix), title='track nInnermostPixHits;track nInnermostPixHits;Entries',path=folder+"Track/allTracks",
            xbins=5, xmin=0, xmax=5.0 )
            igroup.defineHistogram(namer('trackNPixHits','nPixHits','Track_allTracks',postfix), title='track nPixHits;track nPixHits;Entries',path=folder+"Track/allTracks",
            xbins=10, xmin=0.0, xmax=10.0 )
            #nPixelSharedHits;track
            igroup.defineHistogram(namer('trackNSiHits','nSiHits','Track_allTracks',postfix), title='track nSiHits;track nSiHits;Entries',path=folder+"Track/allTracks",
            xbins=20, xmin=0.0, xmax=20.0 )
            #nSCTSharedHits
            #nTRTHits
            igroup.defineHistogram(namer('trackeProbabilityHT','eProbabilityHT','Track_allTracks',postfix), title='track eProbabilityHT;track eProbabilityHT;Entries',path=folder+"Track/allTracks",
            xbins=50, xmin=0.0, xmax=1.0 )
            igroup.defineHistogram(namer('trackeProbabilityNN','eProbabilityNN','Track_allTracks',postfix), title='track eProbabilityNN;track eProbabilityNN;Entries',path=folder+"Track/allTracks",
            xbins=50, xmin=0.0, xmax=1.0 )
            igroup.defineHistogram(namer('trackeProbabilityHTorNN','eProbabilityHTorNN','Track_allTracks',postfix), title='track eProbabilityHTorNN;track eProbabilityHTorNN;Entries',path=folder+"Track/allTracks",
            xbins=50, xmin=0.0, xmax=1.0 )

            igroup.defineHistogram(namer('trackIdScoreCharged','idScoreCharged','Track_allTracks',postfix), title='track RNN idScoreCharged;track idScoreCharged;Entries',path=folder+"Track/allTracks",
            xbins=50, xmin=0.0, xmax=1.0 )
            igroup.defineHistogram(namer('trackIdScoreIso','idScoreIso','Track_allTracks',postfix), title='track RNN idScoreIso;track idScoreIso;Entries',path=folder+"Track/allTracks",
            xbins=50, xmin=0.0, xmax=1.0 )
            igroup.defineHistogram(namer('trackIdScoreConv','idScoreConv','Track_allTracks',postfix), title='track RNN idScoreConv;track idScoreConv;Entries',path=folder+"Track/allTracks",
            xbins=50, xmin=0.0, xmax=1.0 )
            igroup.defineHistogram(namer('trackIdScoreFake','idScoreFake','Track_allTracks',postfix), title='track RNN idScoreFake;track idScoreFake;Entries',path=folder+"Track/allTracks",
            xbins=50, xmin=0.0, xmax=1.0 )

            #original Track Histograms
            igroup.defineHistogram(namer('d0','d0','Track_leadTrack',postfix), title='Track d0;Transverse Impact Parameter (mm)',path=folder+"Track/leadTrack",
            xbins=50, xmin=-5.0, xmax=5.0 )

            igroup.defineHistogram(namer('etaTrack','eta','Track_leadTrack',postfix), title='Track Eta;Eta',path=folder+"Track/leadTrack",
            xbins=51, xmin=-2.55, xmax=2.55 )

            igroup.defineHistogram(namer('leadTrkPt','leadTrkPt','Track_leadTrack',postfix), title='Pt of Leading track;Transverse Momentum (GeV)',path=folder+"Track/leadTrack",
            xbins=50, xmin=0., xmax=200. )

            igroup.defineHistogram(namer('nHighPtTaus','nHighPtTaus','Track',postfix), title='Number of High Pt tau candidates;Number of Taus per Event', 
                    xbins=15, xmin=-0.5, xmax=14.5,path=folder+"Track")

            igroup.defineHistogram(namer('numberOfTRTHighThresholdHits','numberOfTRTHighThresholdHits','Track_leadTrack',postfix), title='Number of TRT High Threshold Hits;Number of High Threshold TRT Hits',path=folder+"Track/leadTrack",
            xbins=15, xmin=-0.5, xmax=14.5 )

            igroup.defineHistogram(namer('numberOfTRTHighThresholdOutliers','numberOfTRTHighThresholdOutliers','Track_leadTrack',postfix), title="Number of TRT High Threshold Outliers;Number of TRT High Threshold Outliers",path=folder+"Track/leadTrack",
            xbins=26, xmin=-0.5, xmax=25.5 )

            igroup.defineHistogram(namer('numberOfTRTHits','numberOfTRTHits','Track_leadTrack',postfix), title='Number of TRT Low Threshold Hits;Number of Low Threshold TRT Hits',path=folder+"Track/leadTrack",
            xbins=101, xmin=-0.5, xmax=100.5 )

            igroup.defineHistogram(namer('numberOfTRTOutliers','numberOfTRTOutliers','Track_leadTrack',postfix), title="Number of TRT Low Threshold Outliers;Number of TRT Low Threshold Outliers",path=folder+"Track/leadTrack",
            xbins=31, xmin=-0.5, xmax=25.5)
            
            igroup.defineHistogram(namer('ptTrack','pT','Track_leadTrack',postfix), title='Track pT;Transverse Momentum (GeV)',path=folder+"Track/leadTrack",
            xbins=60, xmin=0, xmax=300 )

            igroup.defineHistogram(namer('phiTrack','phi','Track_leadTrack',postfix), title='Track Phi;Phi',path=folder+"Track/leadTrack",
            xbins=64, xmin=PHIMIN, xmax=PHIMAX )

            igroup.defineHistogram(namer('trkWidth2','trkWidth2','Track',postfix), title='Weighted Track Width;Momentum-Weighted Width of Track System',path=folder+"Track",
            xbins=50, xmin=0.0, xmax=0.1 )

            igroup.defineHistogram(namer('ipZ0SinThetaSigLeadTrk','ipZ0SinThetaSigLeadTrk','Track_leadTrack',postfix), title='Impact Parameter z0 Sine Theta Significance of Leading Track;Z0SinTheta Significance',path=folder+"Track/leadTrack",
            xbins=50, xmin=-10.0, xmax=10.0 )

            igroup.defineHistogram(namer('numberOfPixelHits','numberOfPixelHits','Track_leadTrack',postfix), title='Number of Pixel Hits;Number of Pixel Hits',path=folder+"Track/leadTrack",
            xbins=25, xmin=-0.5, xmax=25.5 )

            igroup.defineHistogram(namer('numberOfSCTHits','numberOfSCTHits','Track_leadTrack',postfix), title='Number of SCT Hits;Number of SCT Hits',path=folder+"Track/leadTrack",
            xbins=26, xmin=0, xmax=26 )

            igroup.defineHistogram(namer('etOverPtLeadTrack','etOverPtLeadTrack','Track_leadTrack',postfix), title='Et over Pt of lead track of tau candidates;Et/Pt;Number of Candidates',
            xbins=50, xmin=0., xmax=12., path=folder+"Track/leadTrack" )

            igroup.defineHistogram(namer('ipSigLeadTrk','ipSigLeadTrk','Track_leadTrack',postfix), title='Impact Parameter Significance of Leading Track;Transverse Impact Parameter Significance',
            xbins=100, xmin=-5., xmax=5., path=folder+"Track/leadTrack" )

            igroup.defineHistogram(namer('massTrkSys','massTrkSys','Track',postfix), title='Mass of the Track System;Invariant Mass (GeV)',
                xbins=30, xmin=-1., xmax=5., path=folder+"Track" )

            igroup.defineHistogram(namer('trFlightPathSig','trFlightPathSig','Track',postfix),title='Track Transverse Flight Path Significance;Transverse Flight Path Significance',
                xbins=125, xmin=-4.5, xmax=8., path=folder+"Track" )

            igroup.defineHistogram(namer('z0','z0','Track_leadTrack',postfix), title='Track z0;Longitudinal Impact Parameter (mm)',path=folder+"Track/leadTrack",
            xbins=50, xmin=-35.0, xmax=35.0 )




        if postfix == 'Global' or postfix.startswith('TauTrig') or postfix == "EleTrig" or postfix=="JetTrig":
            igroup.defineHistogram(namer('LB','nTauPerLB','',postfix), title='Total number of tau candidates per LB;Luminosity Block);Number of Candidates', 
                    xbins=1000, xmin=0., xmax=1000. ,path=folder)

            igroup.defineHistogram(namer('nHighPtTauCandidates','nHightPtTauCandidates','',postfix), title='Number of High Pt tau candidates;Number of Taus per Event', 
                    xbins=15, xmin=-0.5, xmax=14.5,path=folder)

            igroup.defineHistogram(namer('NumTracks','tauNumTracks','',postfix), title='Number Of Tracks for Tau Candidates;Number Of Tracks;Number Of Candidates',
                      xbins=21, xmin=-0.5, xmax=20.5,path=folder)

            igroup.defineHistogram(namer('nClusters','nCluster','',postfix), title='Number Of CaloTopoClusters;Number Of Clusters;Number Of Candidates',
                      xbins=40, xmin=0., xmax=40. ,path=folder)

            igroup.defineHistogram(namer('tauEta,tauEt','tauEtVsEta','',postfix), type='TH2F', title='Tau Et Vs Eta;#eta;Transverse Energy (Gev)',
                      xbins=40,xmin=-2.55,xmax=2.55 ,ybins=300,ymin=0,ymax=300,path=folder)

            igroup.defineHistogram(namer('tauPhi,tauEt','tauEtVsPhi','',postfix), type='TH2F', title='EtVsPhi;Phi;Transverse Energy (GeV);',
                      xbins=80,xmin=PHIMIN,xmax=PHIMAX,ybins=300,ymin=-0.,ymax=300.,path=folder)


    ### STEP 6 ###
    # Finalize. The return value should be a tuple of the ComponentAccumulator
    # and the sequence containing the created algorithms. If we haven't called
    # any configuration other than the AthMonitorCfgHelper here, then we can 
    # just return directly (and not create "result" above)

    accumulator = cfgHelper.result()
    cfg.merge(accumulator)
    return cfg
    
    # # Otherwise, merge with result object and return
    # acc = cfgHelper.result()
    # result.merge(acc)
    # return result

if __name__=='__main__':
    # Setup logs
    from AthenaCommon.Logging import log
    from AthenaCommon.Constants import INFO
    log.setLevel(INFO)

    # Set the Athena configuration flags
    from AthenaConfiguration.AllConfigFlags import ConfigFlags


    nightly = '/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/CommonInputs/'
    file = 'data16_13TeV.00311321.physics_Main.recon.AOD.r9264/AOD.11038520._000001.pool.root.1'
    ConfigFlags.Input.Files = [nightly+file]

    ConfigFlags.Input.isMC = False
    ConfigFlags.Output.HISTFileName = 'tauRegions.root'
    
    ConfigFlags.lock()

    # Initialize configuration object, add accumulator, merge, and run.
    from AthenaConfiguration.MainServicesConfig import MainServicesCfg 
    from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
    cfg = MainServicesCfg(ConfigFlags)
    cfg.merge(PoolReadCfg(ConfigFlags))

    exampleMonitorAcc = tauMonitoringConfig(ConfigFlags)   # calls the main function above
    cfg.merge(exampleMonitorAcc)

    # If you want to turn on more detailed messages ...
    exampleMonitorAcc.getEventAlgo('tauMonAlgBA').OutputLevel = 2 # DEBUG
    exampleMonitorAcc.getEventAlgo('tauMonAlgCR').OutputLevel = 2 # DEBUG
    exampleMonitorAcc.getEventAlgo('tauMonAlgEC').OutputLevel = 2 # DEBUG
    exampleMonitorAcc.getEventAlgo('tauMonAlgGlobal').OutputLevel = 2 # DEBUG
    exampleMonitorAcc.getEventAlgo('tauMonAlgTauTrig1').OutputLevel = 2 # DEBUG
    exampleMonitorAcc.getEventAlgo('tauMonAlgTauTrig2').OutputLevel = 2 # DEBUG
    exampleMonitorAcc.getEventAlgo('tauMonAlgTauTrig3').OutputLevel = 2 # DEBUG
    exampleMonitorAcc.getEventAlgo('tauMonAlgTauTrig4').OutputLevel = 2 # DEBUG
    exampleMonitorAcc.getEventAlgo('tauMonAlgTauTrig5').OutputLevel = 2 # DEBUG
    exampleMonitorAcc.getEventAlgo('tauMonAlgTauTrig6').OutputLevel = 2 # DEBUG
    exampleMonitorAcc.getEventAlgo('tauMonAlgTauTrig7').OutputLevel = 2 # DEBUG
    exampleMonitorAcc.getEventAlgo('tauMonAlgEleTrig').OutputLevel = 2 # DEBUG
    exampleMonitorAcc.getEventAlgo('tauMonAlgJetTrig').OutputLevel = 2 # DEBUG

    cfg.printConfig(withDetails=True) # set True for exhaustive info

    cfg.run() #use cfg.run(20) to only run on first 20 events
