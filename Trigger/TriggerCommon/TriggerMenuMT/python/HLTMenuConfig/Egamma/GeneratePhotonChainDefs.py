# Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration

from ..Menu.ChainDictTools import splitChainDict
from ..Menu.ChainMerging import mergeChainDefs
from .PhotonChainConfiguration import PhotonChainConfiguration

import pprint
from AthenaCommon.Logging import logging
log = logging.getLogger( 'TriggerMenuMT.HLTMenuConfig.Egamma.generateChainConfigs' )
log.info("Importing %s",__name__)



def generateChainConfigs( chainDict ):
    log.debug('dictionary is: %s\n', pprint.pformat(chainDict))

    
    listOfChainDicts = splitChainDict(chainDict)
    listOfChainDefs = []

    for subChainDict in listOfChainDicts:
        log.debug('Assembling subChainsDict %s for chain %s', len(listOfChainDefs), subChainDict['chainName'] )  
        Photon = PhotonChainConfiguration(subChainDict).assembleChain() 

        listOfChainDefs += [Photon]
        

    if len(listOfChainDefs)>1:
        theChainDef = mergeChainDefs(listOfChainDefs, chainDict)
    else:
        theChainDef = listOfChainDefs[0]

    return theChainDef



    

    
