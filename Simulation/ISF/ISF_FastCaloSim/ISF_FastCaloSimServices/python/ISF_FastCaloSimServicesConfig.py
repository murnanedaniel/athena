# Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration

"""
Tools configurations for ISF_FastCaloSimServices
KG Tan, 04/12/2012
"""

from AthenaCommon import CfgMgr

from FastChainPileup.FastChain_jobProperties import FastChain_Flags
#### FastCaloSimSvc
def getFastCaloSimSvc(name="ISF_FastCaloSimSvc", **kwargs):
    kwargs.setdefault("SimulatorTool",  'ISF_FastCaloTool')
    kwargs.setdefault("Identifier",     'FastCaloSim')
    return CfgMgr.ISF__LegacySimSvc(name, **kwargs )

#### Pileup FastCaloSim
def getFastCaloSimPileupSvc(name="ISF_FastCaloSimPileupSvc", **kwargs):
    kwargs.setdefault("SimulatorTool",  'ISF_FastCaloPileupTool')
    return getFastCaloSimSvc(name, **kwargs)

#### Legacy FastCaloSim
def getLegacyAFIIFastCaloSimSvc(name="ISF_LegacyAFIIFastCaloSimSvc", **kwargs):
    kwargs.setdefault("SimulatorTool",  'ISF_LegacyAFIIFastCaloTool')
    return getFastCaloSimSvc(name, **kwargs)

def getFastHitConvAlgFastCaloSimSvc(name="ISF_FastHitConvAlgFastCaloSimSvc",**kwargs):
    kwargs.setdefault("CaloCellMakerTools_release", [
                                                           #'ISF_AddNoiseCellBuilderTool',
                                                            'ISF_CaloCellContainerFCSFinalizerTool'
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

def getFastCaloSimPileupOTSvc(name="ISF_FastCaloSimPileupOTSvc", **kwargs):
    from ISF_FastCaloSimServices.ISF_FastCaloSimJobProperties import ISF_FastCaloSimFlags
    kwargs.setdefault("BatchProcessMcTruth"              , False                                             )
    kwargs.setdefault("SimulateUndefinedBarcodeParticles", False                                             )
    kwargs.setdefault("Identifier"                       , 'FastCaloSim'                                     )
    kwargs.setdefault("CaloCellsOutputName"              , ISF_FastCaloSimFlags.CaloCellsName()+'PileUp'     )
    #kwargs.setdefault("PUWeights"                        , FastChain_Flags.FastChainPUWeights()  )
    kwargs.setdefault("PUWeights_lar_bapre"              , FastChain_Flags.FastChainPUWeights_lar_bapre()  )
    kwargs.setdefault("PUWeights_lar_hec"                , FastChain_Flags.FastChainPUWeights_lar_hec()  )
    kwargs.setdefault("PUWeights_lar_em"                 , FastChain_Flags.FastChainPUWeights_lar_em()  )
    kwargs.setdefault("PUWeights_tile"                   , FastChain_Flags.FastChainPUWeights_tile()  )
    kwargs.setdefault("CaloCellMakerTools_setup"         , [ 'ISF_EmptyCellBuilderTool' ] )
    kwargs.setdefault("CaloCellMakerTools_simulate"      , [ 'ISF_FastShowerCellBuilderTool' ])
    kwargs.setdefault("CaloCellMakerTools_release"       , [ #'ISF_AddNoiseCellBuilderTool',
                                                             'ISF_CaloCellContainerFCSFinalizerTool',
                                                             'ISF_FastHitConvertTool' ])
    kwargs.setdefault("Extrapolator"                     , 'ISF_NITimedExtrapolator')
    # register the FastCaloSim random number streams
    from G4AtlasApps.SimFlags import simFlags
    if not simFlags.RandomSeedList.checkForExistingSeed(ISF_FastCaloSimFlags.RandomStreamName()):
        simFlags.RandomSeedList.addSeed( ISF_FastCaloSimFlags.RandomStreamName(), 98346412, 12461240 )
    kwargs.setdefault("ParticleTruthSvc"                 , simFlags.TruthStrategy.TruthServiceName() )
    return CfgMgr.ISF__FastCaloSimSvcPU(name, **kwargs )

#### FastCaloSimV2
def getFastCaloSimV2ParamSvc(name="ISF_FastCaloSimV2ParamSvc", **kwargs):
    from ISF_FastCaloSimServices.ISF_FastCaloSimJobProperties import ISF_FastCaloSimFlags
    kwargs.setdefault("ParamsInputFilename"              , ISF_FastCaloSimFlags.ParamsInputFilename())
    kwargs.setdefault("ParamsInputObject"                , 'SelPDGID')
    return CfgMgr.ISF__FastCaloSimV2ParamSvc(name, **kwargs )

def getFastCaloSimSvcV2(name="ISF_FastCaloSimSvcV2", **kwargs):
    kwargs.setdefault("SimulatorTool",  'ISF_FastCaloSimV2Tool')
    kwargs.setdefault("Identifier",     'FastCaloSim')
    return CfgMgr.ISF__LegacySimSvc(name, **kwargs )

#### DNNCaloSim
def getDNNCaloSimSvc(name="ISF_DNNCaloSimSvc", **kwargs):
    from ISF_FastCaloSimServices.ISF_FastCaloSimJobProperties import ISF_FastCaloSimFlags

    kwargs.setdefault("CaloCellsOutputName"              , ISF_FastCaloSimFlags.CaloCellsName()   )
    kwargs.setdefault("CaloCellMakerTools_setup"         , [ 'ISF_EmptyCellBuilderTool' ] )
    kwargs.setdefault("CaloCellMakerTools_release"       , [ 'ISF_CaloCellContainerFinalizerTool',
                                                           'ISF_FastHitConvertTool' ]) #DR needed ?
    kwargs.setdefault("ParamsInputFilename"              , ISF_FastCaloSimFlags.ParamsInputFilename())
    kwargs.setdefault("FastCaloSimCaloExtrapolation"     , 'FastCaloSimCaloExtrapolation')

    # register the FastCaloSim random number streams
    from G4AtlasApps.SimFlags import simFlags
    if not simFlags.RandomSeedList.checkForExistingSeed(ISF_FastCaloSimFlags.RandomStreamName()):
        simFlags.RandomSeedList.addSeed( ISF_FastCaloSimFlags.RandomStreamName(), 98346412, 12461240 )

    kwargs.setdefault("RandomStream"                     , ISF_FastCaloSimFlags.RandomStreamName())
    kwargs.setdefault("RandomSvc"                        , simFlags.RandomSvc.get_Value() )

    return CfgMgr.ISF__DNNCaloSimSvc(name, **kwargs )
