# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration

"""
Tools configurations for ISF_FastCaloSimServices
KG Tan, 04/12/2012
"""

from AthenaCommon import CfgMgr

#### FastCaloSimSvc
def getFastCaloSimSvc(name="ISF_FastCaloSimSvc", **kwargs):
    from ISF_FastCaloSimServices.ISF_FastCaloSimJobProperties import ISF_FastCaloSimFlags
    kwargs.setdefault("BatchProcessMcTruth"              , False                                             )
    kwargs.setdefault("SimulateUndefinedBarcodeParticles", False                                             )
    kwargs.setdefault("Identifier"                       , 'FastCaloSim'                                     )
    kwargs.setdefault("CaloCellsOutputName"              , ISF_FastCaloSimFlags.CaloCellsName()              )
    kwargs.setdefault("PunchThroughTool"                 , 'ISF_PunchThroughTool'             )
    kwargs.setdefault("DoPunchThroughSimulation"         , False                                             )
    kwargs.setdefault("ParticleBroker"                   , 'ISF_ParticleBrokerSvc'               )
    kwargs.setdefault("CaloCellMakerTools_setup"         , [ 'ISF_EmptyCellBuilderTool' ] )
    kwargs.setdefault("CaloCellMakerTools_simulate"      , [ 'ISF_FastShowerCellBuilderTool' ])
    kwargs.setdefault("CaloCellMakerTools_release"       , [ #'ISF_AddNoiseCellBuilderTool',
                                                             'ISF_CaloCellContainerFinalizerTool',
                                                             'ISF_FastHitConvertTool' ])
    kwargs.setdefault("Extrapolator"                     , 'ISF_NITimedExtrapolator')
    # register the FastCaloSim random number streams
    from G4AtlasApps.SimFlags import simFlags
    if not simFlags.RandomSeedList.checkForExistingSeed(ISF_FastCaloSimFlags.RandomStreamName()):
        simFlags.RandomSeedList.addSeed( ISF_FastCaloSimFlags.RandomStreamName(), 98346412, 12461240 )
    return CfgMgr.ISF__FastCaloSimSvc(name, **kwargs )

#### Pileup FastCaloSim
def getFastCaloSimPileupSvc(name="ISF_FastCaloSimPileupSvc", **kwargs):
    from ISF_FastCaloSimServices.ISF_FastCaloSimJobProperties import ISF_FastCaloSimFlags
    kwargs.setdefault("CaloCellsOutputName"              , ISF_FastCaloSimFlags.CaloCellsName()+'PileUp'     )
    kwargs.setdefault("CaloCellMakerTools_simulate"      , [ 'ISF_PileupFastShowerCellBuilderTool' ])
    return getFastCaloSimSvc(name, **kwargs)


#### Legacy FastCaloSim
def getLegacyAFIIFastCaloSimSvc(name="ISF_LegacyAFIIFastCaloSimSvc", **kwargs):
    kwargs.setdefault("BatchProcessMcTruth" , True )
    kwargs.setdefault("ParticleBroker"                   , 'ISF_AFIIParticleBrokerSvc' )
    kwargs.setdefault("CaloCellMakerTools_simulate"      , [ 'ISF_LegacyFastShowerCellBuilderTool' ] )
    return getFastCaloSimSvc(name, **kwargs)

def getFastHitConvAlgFastCaloSimSvc(name="ISF_FastHitConvAlgFastCaloSimSvc",**kwargs):
    kwargs.setdefault("CaloCellMakerTools_release", [
                                                           #'ISF_AddNoiseCellBuilderTool',
                                                            'ISF_CaloCellContainerFinalizerTool'
                                                    ] )
    # setup FastCaloSim hit converter and add it to the alg sequence:
    # -> creates HITS from reco cells
    from AthenaCommon.AlgSequence import AlgSequence
    topSequence=AlgSequence()
    from AthenaCommon.CfgGetter import getAlgorithm
    topSequence+=getAlgorithm('ISF_FastHitConvAlg')
    return getFastCaloSimSvc(name,**kwargs)

def getFastHitConvAlgLegacyAFIIFastCaloSimSvc(name="ISF_FastHitConvAlgLegacyAFIIFastCaloSimSvc", **kwargs):
    kwargs.setdefault("BatchProcessMcTruth" , True )
    return getFastHitConvAlgFastCaloSimSvc(name, **kwargs)

#### FastCaloSimV2
def getFastCaloSimSvcV2(name="ISF_FastCaloSimSvcV2", **kwargs):
    from ISF_FastCaloSimServices.ISF_FastCaloSimJobProperties import ISF_FastCaloSimFlags
    # register the FastCaloSim random number streams
    from G4AtlasApps.SimFlags import simFlags
    if not simFlags.RandomSeedList.checkForExistingSeed(ISF_FastCaloSimFlags.RandomStreamName()):
        simFlags.RandomSeedList.addSeed( ISF_FastCaloSimFlags.RandomStreamName(), 98346412, 12461240 )
    return CfgMgr.ISF__FastCaloSimSvcV2(name, **kwargs )
