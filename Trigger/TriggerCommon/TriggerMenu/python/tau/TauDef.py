# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration

""" Tau slice signatures """

__authors__  = "P.O. DeViveiros, FTK modifications by P. McNamara"
__version__ = "" 
__doc__="Implementation of Tau slice in new TM framework\n now with FTK chains"

from AthenaCommon.Logging import logging
logging.getLogger().info("Importing %s",__name__)
logTauDef = logging.getLogger("TriggerMenu.tau.TauDef")


from TriggerMenu.menu.HltConfig import L2EFChainDef, mergeRemovingOverlap
from AthenaCommon.SystemOfUnits import GeV

####################
# Helper functions #
####################

def L1InputTE(name):
    tmpte = name.split('_')[1]
    #tmpte = 'HA'+tmpte.strip('TAU')
    tmpte = tmpte.replace("TAU", "HA");
    logTauDef.debug('L1 input TE : %s'%tmpte)
    return tmpte

##############
# Main Class #
##############

class L2EFChain_tau(L2EFChainDef):

    def __init__(self, chainDict, theHypos):

        self.L2sequenceList   = []
        self.EFsequenceList   = []
        self.L2signatureList  = []
        self.EFsignatureList  = []
        self.TErenamingDict   = {}

        self.hypoProvider = theHypos

        self.chainPart = chainDict['chainParts']
        
        #if chainDict['stream'] == ['Combined']:
        #    self.chainL1Item = self.chainPart['L1item']
        #else:
        #    self.chainL1Item = chainDict['L1item']


        if self.chainPart['L1item']:
            self.chainL1Item = self.chainPart['L1item']
        else:
            self.chainL1Item = chainDict['L1item']            


        
        self.chainCounter = chainDict['chainCounter']       
        self.L2Name = 'L2_'+self.chainPart['chainPartName']
        self.EFName = 'EF_'+self.chainPart['chainPartName']
        self.chainName = chainDict['chainName']
        self.chainPartName = self.chainPart['chainPartName']

        self.L2InputTE = L1InputTE(self.chainL1Item)

        # Book-keeping for updating sequences
        self.currentItem = self.L2InputTE
        self.currentStep = 0

        selection    = self.chainPart['selection']
        preselection = self.chainPart['preselection']

        if 'r1' in selection or 'r1' in preselection:
            # Run-I set-up
            self.setup_tauChainRunOne()
        else:
            # Run-II set-up
            self.setup_tauChain()
       
        L2EFChainDef.__init__(self, self.chainName, self.L2Name, self.chainCounter, chainDict['L1item'], self.EFName, self.chainCounter, self.L2InputTE)

    def defineSequences(self):

        for sequence in self.L2sequenceList:
            self.addL2Sequence(*sequence)

        for sequence in self.EFsequenceList:
            self.addEFSequence(*sequence)
                
    def defineSignatures(self):
       
        for signature in self.L2signatureList:
            self.addL2Signature(*signature)

        for signature in self.EFsignatureList:
            self.addEFSignature(*signature)

    def defineTErenaming(self):
        self.TErenamingMap = self.TErenamingDict



############################### DEFINE GROUPS OF CHAINS HERE ##############################
#
    # create calorimeter preselection sequence without TrigTauRec and Hypo sequences
    def addCaloSequence(self,threshold,selection,preselection):
        # Run-II calo-based approach
        from TrigCaloRec.TrigCaloRecConfig import TrigCaloCellMaker_tau
        from TrigCaloRec.TrigCaloRecConfig import TrigCaloClusterMaker_topo
        from TrigTauHypo.TrigTauHypoConf import HLTTauCaloRoiUpdater

        cellmaker         = TrigCaloCellMaker_tau()
        clustermaker_topo = TrigCaloClusterMaker_topo()
        caloroiupdater    = HLTTauCaloRoiUpdater()

        # Run topoclustering
        self.EFsequenceList += [[[ self.currentItem ],
                                 [cellmaker, clustermaker_topo, caloroiupdater],
                                 self.continueChain('EF', 'clf0')]]
        
    #create the TrigTauRec Calorimeter only sequence    
    def addTrigTauRecCaloOnlySequence(self,threshold,selection,preselection):
        
        # Run TrigTauRec, calorimeter only (to get proper calibration, and cell-based vars)
        from TrigTauRec.TrigTauRecConfig import TrigTauRecMerged_TauCaloOnly
        caloRec = TrigTauRecMerged_TauCaloOnly()

        self.EFsequenceList += [[[ self.currentItem ],
                                 [caloRec],
                                 self.continueChain('EF', 'calorec')]]
        
    #create the Calorimeter hypo (selection) sequence    
    def addCaloHypoSequence(self,threshold,selection,preselection):    
        theHLTPre   = self.hypoProvider.GetHypo('L2', threshold, selection, 'calo', preselection)
        # Run the calo-based pre-selection
        self.EFsequenceList += [[[ self.currentItem ],
                                 [theHLTPre],
                                 self.continueChain('EF', 'calopre')]]
        
    
    
    #create the tracking sequence and tracking selection sequence    
    def addTrackPreselectionSequences(self,threshold,selection,preselection,idperf,trkprec):        
        # Run-II tracking-based approach
        theHLTTrackPre   = self.hypoProvider.GetHypo('L2', threshold, selection, 'id', preselection)

        # Get the necessary fexes
        from TrigInDetConf.TrigInDetSequence import TrigInDetSequence
        if preselection == 'FTK':
            from TrigInDetConf.TrigInDetFTKSequence import TrigInDetFTKSequence
            # use [:] so the list trkprec is modified by this function
            [trkfast, trkprec[:]] = TrigInDetFTKSequence("Tau","tau",sequenceFlavour=[""]).getSequence()
        elif preselection == 'FTKRefit':
            from TrigInDetConf.TrigInDetFTKSequence import TrigInDetFTKSequence
            [trkfast, trkprec[:]] = TrigInDetFTKSequence("Tau","tau",sequenceFlavour=["refit"]).getSequence()
        else:
            [trkfast, trkprec[:]] = TrigInDetSequence("Tau", "tau", "IDTrig").getSequence()
        # Use cosmic-specific tracking algorithm
        if selection == 'cosmic':
            [trkfast] = TrigInDetSequence("Cosmics", "cosmics", "IDTrig", "FTF").getSequence()

        # Run fast-tracking
        self.EFsequenceList += [[[ self.currentItem ],
                                 trkfast,
                                 self.continueChain('EF', 'trfast')]]
            
        # Run the track-based pre-selection
        # Only cut if we're not in idperf
        if not idperf:
            self.EFsequenceList += [[[ self.currentItem ],
                                     [theHLTTrackPre],
                                     self.continueChain('EF', 'trackpre')]]
            
    #create the vertexing selection sequence        
    def addVertexPreselectionSequence(self,threshold,selection,preselection,idperf):   
        
        #import FEX and xAOD Conversion algorithm   
        from TrigFTK_RecAlgs.TrigFTK_RecAlgs_Config import TrigFTK_VxPrimary_EF
        theTrigFTK_VxPrimary_EF = TrigFTK_VxPrimary_EF("TauFTKVertex", "Tau")
        theTrigFTK_VxPrimary_EF.vxContainerName = 'PrimVxFTK'
        theTrigFTK_VxPrimary_EF.getVertexContainer = False
        theTrigFTK_VxPrimary_EF.useRawTracks = False
        theTrigFTK_VxPrimary_EF.useRefittedTracks = False
        
        from InDetTrigParticleCreation.InDetTrigParticleCreationConf import InDet__TrigVertexxAODCnv 
        theInDet__TrigVertexxAODCnv = InDet__TrigVertexxAODCnv()
        theInDet__TrigVertexxAODCnv.InputVxContainerKey = 'PrimVxFTK'
        theInDet__TrigVertexxAODCnv.OutputVxContainerKey = 'PrimVtxFTK'

        vertexAlgorithms = [theTrigFTK_VxPrimary_EF, theInDet__TrigVertexxAODCnv]
        
        #selection is applied only if not an idperf chain
        if not idperf:
            from TrigTauHypo.TrigTauHypoConf import HLTVertexPreSelHypo
            theHLTVertexPreSelHypo = HLTVertexPreSelHypo()
            vertexAlgorithms.append(theHLTVertexPreSelHypo)
            
        self.EFsequenceList += [[[ self.currentItem ],
                                     vertexAlgorithms,
                                     self.continueChain('EF', 'vertpre')]]
        
    #create the TrigTauRec preselection sequence       
    def addTrigTauRecTauPreselectionSequence(self,threshold,selection,preselection,idperf):              
        # Run TrigTauRec to store pre-selected taus
        from TrigTauRec.TrigTauRecConfig import TrigTauRecMerged_TauPreselection
        recPreselection = TrigTauRecMerged_TauPreselection()

        self.EFsequenceList += [[[ self.currentItem ],
                                 [recPreselection],
                                 self.continueChain('EF', 'storepre')]]
    
    #create the TrigTauRec FTK preselection sequence    
    def addTrigTauRecTauFTKSequence(self,threshold,selection,preselection,idperf):            
        # Run TrigTauRec to store pre-selected taus
        from TrigTauRec.TrigTauRecConfig import TrigTauRecMerged_TauFTK
        recPreselection = TrigTauRecMerged_TauFTK()

        self.EFsequenceList += [[[ self.currentItem ],
                                 [recPreselection],
                                 self.continueChain('EF', 'storepre')]]
        
    #create the two step tracking sequences
    def addTwoStepTrackingSequence(self,threshold,selection,preselection,idperf,trkprec): 
        # Get the necessary fexes
        from TrigInDetConf.TrigInDetSequence import TrigInDetSequence
        # use [:] so the list trkprec is modified by this function
        [trkcore, trkiso, trkprec[:]] = TrigInDetSequence("Tau", "tau", "IDTrig", "2step").getSequence()

        # Get the HLTTrackTauHypo_rejectNoTracks
        from TrigTauHypo.TrigTauHypoBase import HLTTrackTauHypo_rejectNoTracks
        tauRejectEmpty = HLTTrackTauHypo_rejectNoTracks("TauRejectEmpty")

        # Here we load our new tau-specific RoI Updater
        from TrigTauHypo.TrigTauHypoConf import HLTTauTrackRoiUpdater
        tauRoiUpdater = HLTTauTrackRoiUpdater()
        # This will add up to a tolerance of 5 mm due to the extra 3mm tolerance from the FTF
        # tauRoiUpdater.z0HalfWidth = 2.0 # Temporarily widened to 10 mm
        tauRoiUpdater.z0HalfWidth = 7.0

        #ftracks = trkcore+[tauRoiUpdater]+trkiso
        if not idperf:
            ftracks = trkcore+[tauRejectEmpty, tauRoiUpdater]+trkiso
        else :
            ftracks = trkcore+[tauRoiUpdater]+trkiso

        # Run fast-tracking
        self.EFsequenceList += [[[ self.currentItem ],
                                 ftracks,
                                 self.continueChain('EF', 'trfasttwo')]]
        
    def addTwoStepTrackingSelectionSequence(self,threshold,selection,preselection,idperf): 
        theHLTTrackPre   = self.hypoProvider.GetHypo('L2', threshold, selection, 'id', preselection)
        
        # Run the track-based pre-selection
        # Only cut if we're not in idperf
        if not idperf:
            self.EFsequenceList += [[[ self.currentItem ],
                                     [theHLTTrackPre],
                                     self.continueChain('EF', 'trackpre')]]




    # Increment the step counter, set the proper adjusted TE name and return it
    # Also update the Signature list and TErenamingDict as it goes along
    def continueChain(self, level, identifier):
        lastLevel = self.currentItem.split('_')[0]
        if lastLevel != level:
            self.currentStep = 1

        self.currentItem = level+'_tau_step'+str(self.currentStep)

        if level=='L2':
            self.L2signatureList += [ [[self.currentItem]] ]
        elif level=='EF':
            self.EFsignatureList += [ [[self.currentItem]] ]

        self.TErenamingDict[self.currentItem] = mergeRemovingOverlap(level+'_', self.chainPartName+'_'+identifier)
            

        self.currentStep += 1
        return self.currentItem


    # Define the full tau chain
    # Self-configuring based on the dictionary parameters
    def setup_tauChain(self):

        threshold   = self.chainPart['threshold']
        calibration = self.chainPart['calib']
        recoAlg     = self.chainPart['recoAlg'] 
        selection   = self.chainPart['selection']
        preselection= self.chainPart['preselection']
        idperf      = "idperf" in self.chainPart['trkInfo']


        # Cleaner if-statements
        # Strategies which need calorimeter pre-selection
        needsCaloPre  = ['calo', 'ptonly', 'mvonly', 'caloonly',
                         'track', 'trackonly', 'tracktwo',
                         'trackcalo', 'tracktwocalo','tracktwo2015']
        # Strategies which need fast-track finding
        needsTrackTwoPre = ['tracktwo', 'tracktwoonly', 'tracktwocalo','tracktwo2015']
        needsTrackPre    = ['track', 'trackonly', 'trackcalo', 'FTK', 'FTKRefit']
        # Strategies which need Run-II final hypo
        needsRun2Hypo = ['calo', 'ptonly', 'mvonly', 'caloonly',
                         'trackonly', 'track', 'tracktwo', 'tracktwocalo', 'trackcalo', 'FTK', 'FTKRefit','tracktwo2015']
        fastTrackingUsed = needsTrackPre + needsTrackTwoPre
        
        #Set the default values
        from TrigInDetConf.TrigInDetSequence import TrigInDetSequence
        [trkcore, trkprec] = TrigInDetSequence("Tau", "tau", "IDTrig").getSequence()

        # Temporary hack to handle naming scheme
        if 'r1' in selection:
            preselection = 'r1'
            selection = selection.replace('r1', '')
    
        # Overrule the final EF selection
        if idperf:
            selection = 'perf'

        #Create FTK chain or other chains
        if preselection == 'FTK' or preselection == 'FTKRefit':
            self.addTrackPreselectionSequences(threshold, selection, preselection, idperf, trkprec)
            self.addVertexPreselectionSequence(threshold, selection, preselection, idperf)
            self.addCaloSequence(threshold, selection, preselection)
            self.addTrigTauRecTauFTKSequence(threshold,selection,preselection,idperf)
            self.addCaloHypoSequence(threshold,selection,preselection)
        else:
            # Calorimeter
            if preselection in needsCaloPre:
                self.addCaloSequence(threshold, selection, preselection)
                self.addTrigTauRecCaloOnlySequence(threshold,selection,preselection)
                self.addCaloHypoSequence(threshold,selection,preselection)
            # Two step fast-tracking
            if preselection in needsTrackTwoPre:
                self.addTwoStepTrackingSequence(threshold,selection,preselection,idperf, trkprec)
                if preselection != 'tracktwo':
                    self.addTwoStepTrackingSelectionSequence(threshold,selection,preselection,idperf)
                    self.addTrigTauRecTauPreselectionSequence(threshold,selection,preselection,idperf)
                else:
                    self.addTrigTauRecTauPreselectionSequence(threshold,selection,preselection,idperf)
                    self.addTwoStepTrackingSelectionSequence(threshold,selection,preselection,idperf)
            # One step fast-tracking
            if preselection in needsTrackPre:
                self.addTrackPreselectionSequences(threshold, selection, preselection, idperf, trkprec)   
                self.addTrigTauRecTauPreselectionSequence(threshold,selection,preselection,idperf)

        if preselection in needsRun2Hypo:
            # Only run tracking and tau-rec : no need for topoclustering
            if preselection == 'caloonly' or preselection == 'trackonly' or selection == 'cosmic':
                theEFHypo       = self.hypoProvider.GetHypo('EF', threshold, 'perf', '', 'r1')
            else: 
                theEFHypo       = self.hypoProvider.GetHypo('EF', threshold, selection, '', 'r1')


            # Change track selection if we're running on cosmics...
            if selection == 'cosmic':
                from TrigTauRec.TrigTauRecCosmicsConfig import TrigTauRecCosmics_Tau2012
                recmerged_2012    = TrigTauRecCosmics_Tau2012()
            else:
                from TrigTauRec.TrigTauRecConfig import TrigTauRecMerged_TauPrecision
                recmerged_2012    = TrigTauRecMerged_TauPrecision()

            efidinsideout = trkprec

            # Only run the fast-tracking if it wasn't run at pre-selection
            # Is two-step preselection good enough?
            if preselection not in fastTrackingUsed:
                efidinsideout = trkcore+trkprec


            # Precision tracking
            self.EFsequenceList += [[[ self.currentItem ],
                                     efidinsideout,
                                     self.continueChain('EF', 'tr')]]


            # TrigTauRec and Hypo (no BDT)
            if selection == 'dikaon' or selection == 'dikaontight':
                self.EFsequenceList += [[[ self.currentItem ],
                                         [recmerged_2012, theEFHypo],
                                         self.continueChain('EF', 'effinal')]]                
            else:
            # TrigTauRec, BDT and Hypo
                from TrigTauDiscriminant.TrigTauDiscriGetter import TrigTauDiscriGetter2015
                efmv              = TrigTauDiscriGetter2015()
                self.EFsequenceList += [[[ self.currentItem ],
                                         [recmerged_2012, efmv, theEFHypo],
                                         self.continueChain('EF', 'effinal')]]

    def setup_tauChainRunOne(self):
        
        threshold   = self.chainPart['threshold']
        calibration = self.chainPart['calib']
        recoAlg     = self.chainPart['recoAlg']
        selection   = self.chainPart['selection']
        preselection= self.chainPart['preselection']
        idperf      = "idperf" in self.chainPart['trkInfo']

        # Handle Run-II naming scheme
        if 'r1' in selection:
            preselection = 'r1'
            selection = selection.replace('r1', '')

        # Overrule the final EF selection
        if idperf:
            selection = 'perf'    

        if preselection == 'r1':
            # Try new hypo extraction method
            theL2CaloHypo   = self.hypoProvider.GetHypo('L2', threshold, selection, 'calo', 'r1')
            theL2IDHypo     = self.hypoProvider.GetHypo('L2', threshold, selection, 'id', 'r1')
            theL2FinalHypo  = self.hypoProvider.GetHypo('L2', threshold, selection, '', 'r1')
            # Get the necessary fexes
            from TrigT2CaloTau.TrigT2CaloTauConfig import T2CaloTau_Tau_Med
            from TrigL2SiTrackFinder.TrigL2SiTrackFinder_Config import TrigL2SiTrackFinder_TauB
            from TrigT2IDTau.T2IDTauConfig import T2IDTau_Tau_1GeV_dZ02_dR0103
            from TrigT2Tau.T2TauFinalConfig import T2TauFinal_Tau_dR03_1GeV_dZ02

            t2calo_2012 = T2CaloTau_Tau_Med()
            l2sitrkfinder_tauB = TrigL2SiTrackFinder_TauB()
            t2id_2012 = T2IDTau_Tau_1GeV_dZ02_dR0103()
            t2final_2012 = T2TauFinal_Tau_dR03_1GeV_dZ02()

            if idperf:
                from TrigL2SiTrackFinder.TrigL2SiTrackFinder_Config import TrigL2SiTrackFinder_TauA,TrigL2SiTrackFinder_TauB,TrigL2SiTrackFinder_TauC
                l2sitrkfinder_tauA = TrigL2SiTrackFinder_TauA()
                l2sitrkfinder_tauC = TrigL2SiTrackFinder_TauC()

            # L2 configuration
            self.L2sequenceList += [[[ self.currentItem ],
                                     [t2calo_2012, theL2CaloHypo],
                                     self.continueChain('L2', 'calo')]]

            if idperf:
                self.L2sequenceList += [[[ self.currentItem ],
                                         [l2sitrkfinder_tauA, l2sitrkfinder_tauB, l2sitrkfinder_tauC, t2id_2012, theL2IDHypo],
                                         self.continueChain('L2', 'id')]]
            else:
                self.L2sequenceList += [[[ self.currentItem ],
                                         [l2sitrkfinder_tauB, t2id_2012, theL2IDHypo],
                                         self.continueChain('L2', 'id')]]

            self.L2sequenceList += [[[ self.currentItem ],
                                     [t2final_2012, theL2FinalHypo],
                                     self.continueChain('L2', 'l2final')]]

        if preselection == 'r1':

            theEFHypo       = self.hypoProvider.GetHypo('EF', threshold, selection, '', 'r1')

            # Get the necessary fexes
            from TrigCaloRec.TrigCaloRecConfig import TrigCaloCellMaker_tau
            from TrigCaloRec.TrigCaloRecConfig import TrigCaloClusterMaker_topo
            from InDetTrigRecExample.EFInDetConfig import  TrigEFIDInsideOut_Tau
            from TrigTauRec.TrigTauRecConfig import TrigTauRecMerged_Tau2012
            from TrigTauDiscriminant.TrigTauDiscriGetter import TrigTauDiscriGetter

            cellmaker         = TrigCaloCellMaker_tau()
            clustermaker_topo = TrigCaloClusterMaker_topo()
            efidinsideout     = TrigEFIDInsideOut_Tau().getSequence()
            recmerged_2012    = TrigTauRecMerged_Tau2012()
            efmv              = TrigTauDiscriGetter()

            self.EFsequenceList += [[[ self.currentItem ],
                                     [cellmaker, clustermaker_topo],
                                     self.continueChain('EF', 'clf0')]]

            self.EFsequenceList += [[[ self.currentItem ],
                                     efidinsideout,
                                     self.continueChain('EF', 'tr')]]

            self.EFsequenceList += [[[ self.currentItem ],
                                     [recmerged_2012, efmv, theEFHypo],
                                     self.continueChain('EF', 'effinal')]]
            

    # Prototype for TwoStep configuration
    def setup_tauChainTwoStep(self):

        threshold   = self.chainPart['threshold']
        calibration = self.chainPart['calib']
        recoAlg     = self.chainPart['recoAlg'] 
        selection   = self.chainPart['selection']
        preselection= self.chainPart['preselection']
        idperf      = "idperf" in self.chainPart['trkInfo']
        
        # Overrule the final EF selection
        if idperf:
            selection = 'perf'

        if preselection == 'calo' or preselection == 'ptonly' or preselection == 'mvonly' or preselection == 'caloonly' or preselection == 'track' or preselection == 'trackonly':
            # Test 2015 approach
            logTauDef.info("Calo-based pre-selection configuration is not quite ready yet!")
            logTauDef.info("Very preliminary version!!")

            theHLTPre   = self.hypoProvider.GetHypo('L2', threshold, selection, 'calo', preselection)

            from TrigCaloRec.TrigCaloRecConfig import TrigCaloCellMaker_tau
            from TrigCaloRec.TrigCaloRecConfig import TrigCaloClusterMaker_topo

            cellmaker         = TrigCaloCellMaker_tau()
            clustermaker_topo = TrigCaloClusterMaker_topo()

            # Run topoclustering
            self.EFsequenceList += [[[ self.currentItem ],
                                     [cellmaker, clustermaker_topo],
                                     self.continueChain('EF', 'clf0')]]
                

        if preselection == 'track' or preselection == 'trackonly' or (preselection != 'r1' and idperf):
    
            theHLTTrackPre  = self.hypoProvider.GetHypo('L2', threshold, selection, 'id', preselection)

            # Get the necessary fexes
            from TrigInDetConf.TrigInDetSequence import TrigInDetSequence
            [trkfast, trkprec] = TrigInDetSequence("Tau", "tau", "IDTrig").getSequence()
            
            # Run fast-tracking
            self.EFsequenceList += [[[ self.currentItem ],
                                     trkfast,
                                     self.continueChain('EF', 'trfast')]]
            
        # Here we're running the TrigTauRec based on all the stuff that ran before.  Uh-oh, this is dangerous...
        # For now, assume fast-tracking is always run
        from TrigTauRec.TrigTauRecConfig import TrigTauRecMerged_TauPreselection
        
        recPreselection = TrigTauRecMerged_TauPreselection()

        self.EFsequenceList += [[[ self.currentItem ],
                                 recPreselection,
                                 self.continueChain('EF', 'prefinal')]]


        # Insert generic hypothesis
        
        
        
        if preselection == 'calo' or preselection == 'ptonly' or preselection == 'mvonly' or preselection == 'caloonly' or preselection == 'trackonly' or preselection == 'track':
            # Only run tracking and tau-rec : no need for topoclustering
            if preselection == 'caloonly' or preselection == 'trackonly' or selection == 'cosmic':
                theEFHypo       = self.hypoProvider.GetHypo('EF', threshold, 'perf', '', 'r1')
            else: 
                theEFHypo       = self.hypoProvider.GetHypo('EF', threshold, selection, '', 'r1')

            # Get the necessary fexes
            from TrigInDetConf.TrigInDetSequence import TrigInDetSequence
            [trkfast, trkprec] = TrigInDetSequence("Tau", "tau", "IDTrig").getSequence()

            from TrigTauDiscriminant.TrigTauDiscriGetter import TrigTauDiscriGetter2015

            # Change track selection if we're running on cosmics...
            if selection == 'cosmic':
                from TrigTauRec.TrigTauRecCosmicsConfig import TrigTauRecCosmics_Tau2012
                recmerged_2012    = TrigTauRecCosmics_Tau2012()
            else:
                from TrigTauRec.TrigTauRecConfig import TrigTauRecMerged_TauPrecision
                recmerged_2012    = TrigTauRecMerged_TauPrecision()

            # Only run the fast-tracking if it wasn't run at pre-selection
            if preselection != 'track' and preselection != 'trackonly' and not idperf:
                efidinsideout     = trkfast+trkprec
            else:
                efidinsideout     = trkprec
            
            efmv              = TrigTauDiscriGetter2015()


            self.EFsequenceList += [[[ self.currentItem ],
                                     efidinsideout,
                                     self.continueChain('EF', 'tr')]]

            self.EFsequenceList += [[[ self.currentItem ],
                                     [recmerged_2012, efmv, theEFHypo],
                                     self.continueChain('EF', 'effinal')]]
