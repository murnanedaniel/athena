# Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration

from PyUtils.Decorators import memoize
from six import add_metaclass

# Configure the scheduler
from AthenaCommon.AlgScheduler import AlgScheduler
AlgScheduler.ShowControlFlow( True )
AlgScheduler.ShowDataFlow( True )

from TriggerJobOpts.TriggerFlags import TriggerFlags

from TriggerMenuMT.HLTMenuConfig.Menu.TriggerConfigHLT  import TriggerConfigHLT
from TriggerMenuMT.HLTMenuConfig.Menu.HLTCFConfig import makeHLTTree
from TriggerMenuMT.HLTMenuConfig.Menu.DictFromChainName import dictFromChainName
from TriggerMenuMT.HLTMenuConfig.Menu.ChainDictTools import splitInterSignatureChainDict
from TriggerMenuMT.HLTMenuConfig.Menu.MenuPrescaleConfig import MenuPrescaleConfig, applyHLTPrescale
from TriggerMenuMT.HLTMenuConfig.Menu.ChainMerging import mergeChainDefs
from TriggerMenuMT.HLTMenuConfig.Menu.MenuAlignmentTools import MenuAlignment
from TriggerMenuMT.HLTMenuConfig.CommonSequences import EventBuildingSequenceSetup

from AthenaCommon.Logging import logging
log = logging.getLogger( __name__ )

_func_to_modify_signatures = None
        
class Singleton(type):
    _instances = {}
    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


# for now we make this a singleton because calling menu generation twice leads to problems
@add_metaclass(Singleton)
class GenerateMenuMT(object):

    @staticmethod
    def overwriteSignaturesWith(f):
        """
        == Function to generate menu for certain signatures only
        """
        log.info('In overwriteSignaturesWith ')
        global _func_to_modify_signatures
        if _func_to_modify_signatures is not None:
            log.warning('Updating the function to modify signatures from %s to %s',
                        _func_to_modify_signatures.__name__, f.__name__)
        _func_to_modify_signatures = f


    def __init__(self):
        self.chainsInMenu = []
        self.listOfErrorChainDefs = []
        self.selectChainsForTesting = []
        
        self.allChainsForAlignment = []
        self.chainDicts = []
        self.combinationsInMenu = []
        self.combinationsInMenu = []
        self.alignmentGroupsToAlign = set()
        self.configLengthDict = {}
        
        self.signaturesOverwritten = False
        self.L1Prescales = None
        self.HLTPrescales = None

        self.availableSignatures = []
        self.signaturesToGenerate = []
        self.allSignatures = ['Egamma', 'Muon', 'Jet', 'Bjet', 'Bphysics', 'MET', 'Tau',
                              'HeavyIon', 'Beamspot', 'Cosmic', 'EnhancedBias',
                              'Monitor', 'Calib', 'Streaming', 'Combined', 'MinBias', 'UnconventionalTracking', 'Test'] #, AFP
        self.calibCosmicMonSigs = ['Streaming','Monitor','Beamspot','Cosmic', 'Calib', 'EnhancedBias']
        
        # flags
        self.doEgammaChains         = True
        self.doJetChains            = True
        self.doBjetChains           = True
        self.doMuonChains           = True
        self.doBphysicsChains       = True
        self.doMETChains            = True
        self.doTauChains            = True
        self.doAFPChains            = True
        self.doMinBiasChains        = True
        self.doHeavyIonChains       = True
        self.doCosmicChains         = True
        self.doCalibChains          = True
        self.doStreamingChains      = True
        self.doMonitorChains        = True
        self.doBeamspotChains       = True
        self.doEnhancedBiasChains   = True
        self.doUnconventionalTrackingChains   = True
        self.doCombinedChains       = True
        self.doTestChains           = True


    def setTriggerConfigHLT(self):
        """
        == Setup of TriggerConfigHLT, menu and prescale names
        """
        (HLTPrescales) = self.setupMenu()
        log.debug("   and prescales : %s", HLTPrescales)


    def setupMenu(self):
        """
        == Setup of menu in terms of prescales and overwrite function
        """
        # go over the slices and put together big list of signatures requested
        (self.L1Prescales, self.HLTPrescales) = MenuPrescaleConfig(TriggerConfigHLT)
        global _func_to_modify_signatures
        if _func_to_modify_signatures is not None:
            log.info('setupMenu:  Modifying trigger signatures in TriggerFlags with %s',
                     _func_to_modify_signatures.__name__)
            _func_to_modify_signatures()
            self.signaturesOverwritten = True
        else:
            log.debug("No signature overwrite flags are set")

        return (self.HLTPrescales)


    def generateL1Topo(self):
        """
        == Generates the L1Topo menu
        """
        if not TriggerFlags.readL1TopoConfigFromXML() and not TriggerFlags.readMenuFromTriggerDb():
            log.info('Generating L1 topo configuration for %s', TriggerFlags.triggerMenuSetup())
            from TriggerMenuMT.LVL1MenuConfig.TriggerConfigL1Topo import TriggerConfigL1Topo
            self.trigConfL1Topo = TriggerConfigL1Topo( outputFile = TriggerFlags.outputL1TopoConfigFile() )
            # build the menu structure
            self.trigConfL1Topo.generateMenu()
            log.debug('Topo Menu has %i trigger lines', len(self.trigConfL1Topo.menu))
            # write xml file
            self.trigConfL1Topo.writeXML()
        elif TriggerFlags.readLVL1configFromXML():
            log.info("ReadingL1TopocofnigFromXML currently not implemented")
        else:
            log.info("Doing nothing with L1Topo menu configuration...")


    def generateLVL1(self):
        """
        == Generates the LVL1 menu
        """
        if not TriggerFlags.readLVL1configFromXML() and not TriggerFlags.readMenuFromTriggerDb():
            log.info('Generating L1 configuration for %s', TriggerFlags.triggerMenuSetup() )
            from TriggerMenuMT.LVL1MenuConfig.TriggerConfigLVL1 import TriggerConfigLVL1
            self.trigConfL1 = TriggerConfigLVL1( outputFile = TriggerFlags.outputLVL1configFile())
            # build the menu structure
            self.trigConfL1.generateMenu()
            log.debug('Menu has %i items', len(self.trigConfL1.menu.items) )
            # write xml file
            self.trigConfL1.writeXML()
        elif TriggerFlags.readLVL1configFromXML():
            log.info("ReadingLVL1cofnigFromXML currently not implemented")
        else:
            log.info("Doing nothing with L1 menu configuration...")
    

    def getChainDicts(self):
        chainCounter = 0
        all_chain_dicts = []
        for chain in self.chainsInMenu:
            log.debug("Now processing chain: %s ", chain) 
            chainDict = dictFromChainName(chain)

            chainCounter += 1
            chainDict['chainCounter'] = chainCounter

            #set default chain prescale
            chainDict['prescale'] = 1
            all_chain_dicts += [chainDict]
        self.chainDicts = all_chain_dicts
        
        return 

    def importSignaturesToGenerate(self):
        # check if all the signature files can be imported files can be imported
        # and then import them!
        log.debug("[getSignaturesInMenu] signaturesToGenerate: %s",  self.signaturesToGenerate)

        for sig in self.signaturesToGenerate:
            log.debug("[getSignaturesInMenu] sig: %s", sig)
            
            try:
                if eval('self.do' + sig + 'Chains'):
                    if sig == 'Egamma':
                        sigFolder = sig
                        subSigs = ['Electron', 'Photon']
                    elif sig in self.calibCosmicMonSigs:
                        sigFolder = 'CalibCosmicMon'
                         #only import the CalibCosmicMon signatures that we need, not all of them!
                        subSigs = [sig]
                    else:
                        sigFolder = sig
                        subSigs = [sig]
                    for ss in subSigs:
                        if sigFolder == 'Combined':
                            continue
                        else:
                            #import the includes into the global namespace. Only import the signature we need!
                            #this is equivalent having this line at the beginning of the file:
                            #import TriggerMenuMT.HLTMenuConfig.[sig].Generate[sig]ChainDefs as Generate[sig]ChainDefs
                            import_module = 'TriggerMenuMT.HLTMenuConfig.' + sigFolder +'.Generate' + ss + 'ChainDefs'
                            globals()['Generate'+ss+'ChainDefs'] = __import__(import_module,fromlist=['Generate'+ss+'ChainDefs'])
                            if ss not in self.availableSignatures:
                                self.availableSignatures.append(ss)

            except ImportError:
                log.exception('Problems when importing ChainDef generating code for %s', sig)
                import traceback
                traceback.print_exc()

        log.debug('Available signature(s) for chain generation: %s', self.availableSignatures)
        return

    def generateChains(self):

        all_chains = []
        combinations_in_menu = []
        alignmentGroups_to_align = set()
        length_of_configs = {}
        

        for chainDict in self.chainDicts:

            log.debug("Next: getting chain configuration for chain %s ", chainDict['chainName']) 
            chainConfig,lengthOfChainConfigs = self.__generateChainConfig(chainDict)
            
            all_chains += [(chainDict,chainConfig,lengthOfChainConfigs)]
            
            #update the alignment group length dictionary if we have a longer number of steps
            #or the signature isn't registered in the dictionary yet
            for config_length, config_grp in lengthOfChainConfigs:
                if config_grp in length_of_configs:
                    if config_length > length_of_configs[config_grp]:
                        length_of_configs[config_grp] = config_length
                else:
                    length_of_configs[config_grp] = config_length

            # find the chains that contain more than one alignment group, to keep track
            # of what combinations do we need to deal with.
            # we're using sets here so we don't end up with duplicates
            if len(set(chainDict['alignmentGroups'])) > 1:  
                combinations_in_menu += [list(set(chainDict['alignmentGroups']))]
                for align_group in list(set(chainDict['alignmentGroups'])):
                    alignmentGroups_to_align.update([align_group])
                    
        self.allChainsForAlignment = all_chains
        self.combinationsInMenu = combinations_in_menu
        self.alignmentGroupsToAlign = alignmentGroups_to_align
        self.configLengthDict = length_of_configs

        return 

    @memoize
    def generateAllChainConfigs(self):
        """
        == Obtains chain configs for all chains in menu
        """

        # get all chain names from menu
        log.debug ("getting chains from Menu")
        self.getChainsFromMenu()

        # decoding of the chain name
        self.getChainDicts()
        
        #import the necessary signatures
        self.importSignaturesToGenerate()
        
        self.generateChains()

        # align event building sequences
        EventBuildingSequenceSetup.alignEventBuildingSteps(self.allChainsForAlignment)

        #dict of signature: set it belongs to
        #e.g. {'Electron': ['Electron','Muon','Photon']}        
        menuAlignment = MenuAlignment(self.combinationsInMenu,
                                      self.alignmentGroupsToAlign,
                                      self.configLengthDict)
        menuAlignment.analyse_combinations()

        # alignmentGroups_to_align = menuAlignment.groupsToAlign
        # lengthOfChainConfigs = self.configLengthDict
        # combinationsInMenu = menuAlignment.combinationsInMenu
        # alignmentGroup_sets_to_align = menuAlignment.setsToAlign

        log.debug('Aligning the following signatures with sets: %s',sorted(menuAlignment.sets_to_align))
        log.debug('Length of each of the alignment groups: %s',self.configLengthDict)


        for chainDict,chainConfig,lengthOfChainConfigs in self.allChainsForAlignment:
        
              # start by ordering electron, photon, muon by having e+mu, g+mu, e+g chains
              # desired ordering: electron, photon, muon, tau, jet, met, b-jet
              
              # lengthOfChainConfigs is something like this: [(4, 'Photon'), (5, 'Muon')]
              # needs to match up with the maximum number of steps in a signature in the menu (length_of_configs)
              # start with electron! Only need to add post-steps for combined electron chains if the max length in a combined chain
              # is greater than the number of electron steps combined chain. Assume that the max length of an electron chain occurs 
              # in a combined chain.
              
              alignmentGroups = chainDict['alignmentGroups']
            
              #parallel-merged single-signature chains or single signature chains. Anything that needs no splitting!
              if len(set(alignmentGroups)) == 1: 
                  alignedChainConfig = menuAlignment.single_align(chainDict, chainConfig)
                  TriggerConfigHLT.registerChain( chainDict, alignedChainConfig )

              elif len(alignmentGroups) == 2:
                  alignedChainConfig = menuAlignment.multi_align(chainDict, chainConfig, lengthOfChainConfigs)
                  TriggerConfigHLT.registerChain( chainDict, alignedChainConfig )              

              else: 
                  log.error("Menu can't deal with combined chains with more than two alignmentGroups at the moment. oops...")
                  raise NotImplementedError("more than three alignment groups still needs implementing in ChainMerging.py, ATR-22206")

              if not TriggerConfigHLT.isChainRegistered(chainDict['chainName']):
                log.error("Chain %s has not been registered in the menu!", chainDict['chainName'])
                import pprint
                pp = pprint.PrettyPrinter(indent=4, depth=8)
                log.error('The chain dictionary is: %s', pp.pformat(chainDict))
                raise Exception("Please fix the menu or the chain.")
        
        return TriggerConfigHLT.configsList()

    @memoize
    def getChainsFromMenu(self):
        """
        == Returns the list of chain names that are in the menu
        """
        log.debug('Setting TriggerConfigHLT to get the right menu')

        self.setTriggerConfigHLT()

        log.debug('Creating one big list of of enabled signatures and chains')

        chains = []
        ## we can already use new set of flags
        #from AthenaConfiguration.AllConfigFlags import ConfigFlags
        #from TriggerMenuMT.HLTMenuConfig.Menu.LS2_v1_newJO import setupMenu as setupMenuFlags
        #setupMenuFlags( ConfigFlags )
        #ConfigFlags.lock()

        for sig in self.allSignatures:
            if eval('TriggerFlags.' + sig + 'Slice.signatures()') and eval('self.do' + sig + 'Chains'):
                log.debug("Adding %s chains to the list of chains to be configured", sig)
                chains+= eval('TriggerFlags.' + sig + 'Slice.signatures()')
                self.signaturesToGenerate.append(sig)
            elif not eval('TriggerFlags.' + sig + 'Slice.signatures()'):
                log.debug('Signature %s is not switched on (no chains in menu)', sig)
            elif not eval('self.do' + sig + 'Chains'):
                log.debug('Signature %s is not switched on (disabled by flag)', sig)

        log.info("The following signature(s) is (are) enabled: %s", self.signaturesToGenerate)

        if self.selectChainsForTesting:
            log.info("Narrowing down the list of chains with the selectChainsForTesting list")
            selectedChains = [ch for ch in chains if ch.name in self.selectChainsForTesting]
            if len(selectedChains) < len(self.selectChainsForTesting):
                selectedNames = [ch.name for ch in selectedChains]
                missingNames = [name for name in self.selectChainsForTesting if name not in selectedNames]
                log.warning("The following chains were specified in selectChainsForTesting but were not found in the menu: %s", str(missingNames))
            chains = selectedChains

        if len(chains) == 0:
            log.warning("There seem to be no chains in the menu - please check")
        else:
            log.debug("The following chains were found in the menu: \n %s", '\n '.join(map(str,(chains) )))
            
        self.chainsInMenu = chains 
                                

    def __generateChainConfig(self, mainChainDict):
        """
        # Assembles the chain configuration and returns a chain object with (name, L1see and list of ChainSteps)
        """

        # split the the chainDictionaries for each chain and print them in a pretty way
        chainDicts = splitInterSignatureChainDict(mainChainDict)

        if log.isEnabledFor(logging.DEBUG):
            import pprint
            pp = pprint.PrettyPrinter(indent=4, depth=8)
            log.debug('dictionary is: %s', pp.pformat(chainDicts))

        # Loop over all chainDicts and send them off to their respective assembly code
        listOfChainConfigs = []
        tmp_lengthOfChainConfigs = []

        for chainPartDict in chainDicts:
            chainPartConfig = None
            currentSig = chainPartDict['signature']
            currentAlignGroup = None
            if len(chainPartDict['chainParts']) == 1:
                currentAlignGroup = chainPartDict['chainParts'][0]['alignmentGroup']
            
            chainName = chainPartDict['chainName']
            log.debug('Checking chainDict for chain %s in signature %s, alignment group %s' , chainName, currentSig, currentAlignGroup)

            sigFolder = ''
            if currentSig == 'Electron' or currentSig == 'Photon':
                sigFolder = 'Egamma'
            elif currentSig in self.calibCosmicMonSigs:
                sigFolder = 'CalibCosmicMon'
            else:
                sigFolder = currentSig

            if currentSig in self.availableSignatures and currentSig != 'Combined':
                try:                    
                    functionToCall ='Generate'+currentSig+'ChainDefs.generateChainConfigs(chainPartDict)' 
                    log.debug("Trying to get chain config for %s in folder %s using %s", currentSig, sigFolder, functionToCall)
                    chainPartConfig = eval(functionToCall)
                except RuntimeError:
                    log.exception( 'Problems creating ChainDef for chain\n %s ', chainName)
                    continue
            else:
                log.error('Chain %s ignored - Signature not available', chainPartDict['chainName'])
            log.debug("Chain %s chain configs: %s",chainPartDict['chainName'],chainPartConfig)
            listOfChainConfigs.append(chainPartConfig)
            tmp_lengthOfChainConfigs.append((chainPartConfig.nSteps,chainPartConfig.alignmentGroups))

        # this will be a list of lists for inter-sig combined chains and a list with one 
        # multi-element list for intra-sig combined chains
        # here, we flatten it accordingly (works for both cases!)
        lengthOfChainConfigs = []
        for nSteps, aGrps in tmp_lengthOfChainConfigs:
            if len(nSteps) != len(aGrps):
                log.error("Chain part has %s steps and %s alignment groups - these don't match!",nSteps,aGrps)
            else:
                for a,b in zip(nSteps,aGrps):
                    lengthOfChainConfigs.append((a,b))         
            
        ## if log.isEnabledFor(logging.DEBUG):
        ##     import pprint
        ##     pp = pprint.PrettyPrinter(indent=4, depth=8)
        ##     log.debug('mainChainDict dictionary: %s', pp.pformat(mainChainDict))


        # This part is to deal with combined chains between different signatures
        if len(listOfChainConfigs) == 0:
            log.error('No Chain Configuration found ')
            return False

        elif len(listOfChainConfigs)>1:
                log.debug("Merging strategy from dictionary: %s", mainChainDict["mergingStrategy"])
                theChainConfig = mergeChainDefs(listOfChainConfigs, mainChainDict)

                # This needs to be added for topological chains - needs implementation
                #doTopo = self.CheckIntraSignatureTopo(chainDicts) and chainDict["topo"]
                #if doTopo:
                # theChainDef = TriggerMenu.combined.generateCombinedChainDefs._addTopoInfo(theChainDef,chainDicts,listOfChainDefs)

        else:
            theChainConfig = listOfChainConfigs[0]

        # Configure event building strategy
        eventBuildType = mainChainDict['eventBuildType']
        if eventBuildType:
            log.debug('Configuring event building sequence %s for chain %s', eventBuildType, mainChainDict['chainName'])
            EventBuildingSequenceSetup.addEventBuildingSequence(theChainConfig, eventBuildType, mainChainDict)

        log.debug('ChainConfigs  %s ', theChainConfig)
        return theChainConfig,lengthOfChainConfigs

    @memoize
    def generateMT(self):
        """
        == Main function of the class which generates L1, L1Topo and HLT menu
        """
        log.info('GenerateMenuMT.py:generateMT ')

        # --------------------------------------------------------------------
        # L1 menu generation
        # - from the code, from DB and from xmls (if we want to maintain this)
        # currently implementing the generation from configuration code
        # --------------------------------------------------------------------
        #generateL1Topo()
        #generateLVL1()


        # --------------------------------------------------------------------
        # HLT menu generation
        # --------------------------------------------------------------------
        finalListOfChainConfigs = self.generateAllChainConfigs()
        log.debug("Length of FinalListofChainConfigs %s", len(finalListOfChainConfigs))

        log.debug("finalListOfChainConfig %s", finalListOfChainConfigs)
        for cc in finalListOfChainConfigs:
            log.info('ChainName %s', cc.name)
            log.info('  ChainSteps %s', cc.steps)
            for step in cc.steps:
                log.info(step)

        makeHLTTree(newJO=False, triggerConfigHLT = TriggerConfigHLT)
        # the return values used for debugging, might be removed later

        # Having built the Menu add prescales for disabling items (e.g. MC production)
        applyHLTPrescale(TriggerConfigHLT, self.HLTPrescales, self.signaturesOverwritten)

        from TriggerMenuMT.HLTMenuConfig.Menu.HLTMenuJSON import generateJSON
        generateJSON()

        from TriggerMenuMT.HLTMenuConfig.Menu.HLTPrescaleJSON import generateJSON as generatePrescaleJSON
        generatePrescaleJSON()

        return finalListOfChainConfigs
