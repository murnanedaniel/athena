# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory
from AthenaConfiguration.Enums import BunchStructureSource


def BunchCrossingCondAlgCfg(configFlags):
    BunchCrossingCondAlg=CompFactory.BunchCrossingCondAlg
    from IOVDbSvc.IOVDbSvcConfig import addFolders

    result=ComponentAccumulator()

    run1=(configFlags.IOVDb.DatabaseInstance=='COMP200')
    cfgsvc = None
    folder = ''

    if configFlags.Beam.BunchStructureSource == BunchStructureSource.MC:
        folder = "/Digitization/Parameters"
        from AthenaConfiguration.Enums import ProductionStep
        if configFlags.Common.ProductionStep not in [ProductionStep.Digitization, ProductionStep.PileUpPresampling, ProductionStep.Overlay, ProductionStep.FastChain]:
            result.merge(addFolders(configFlags,folder,None,className="AthenaAttributeList",tag='HEAD'))
        else:
            # Here we are in a job which runs digitization, so the
            # /Digitization/Parameters metadata is not present in the
            # input file and will be created during the job
            pass
    elif configFlags.Beam.BunchStructureSource == BunchStructureSource.FILLPARAMS:
        folder = '/TDAQ/OLC/LHC/FILLPARAMS'
        result.merge(addFolders(configFlags,folder,'TDAQ',className = 'AthenaAttributeList',tag='HEAD'))
    elif configFlags.Beam.BunchStructureSource == BunchStructureSource.TrigConf:
        from TrigConfxAOD.TrigConfxAODConfig import getxAODConfigSvc
        cfgsvc = result.getPrimaryAndMerge(getxAODConfigSvc(configFlags))
    elif configFlags.Beam.BunchStructureSource == BunchStructureSource.Lumi:
        from .LuminosityCondAlgConfig import LuminosityCondAlgCfg
        result.merge(LuminosityCondAlgCfg(configFlags))

    alg = BunchCrossingCondAlg('BunchCrossingCondAlgDefault',
                               Run1=run1,
                               FillParamsFolderKey=folder,
                               Mode=configFlags.Beam.BunchStructureSource.value,
                               TrigConfigSvc=cfgsvc 
                               )

    result.addCondAlgo(alg)

    return result



if __name__=="__main__":
    from AthenaConfiguration.AllConfigFlags import ConfigFlags
    from AthenaConfiguration.MainServicesConfig import MainServicesCfg
    ConfigFlags.Input.Files = []
    ConfigFlags.Input.isMC=False
    ConfigFlags.IOVDb.DatabaseInstance="CONDBR2"
    ConfigFlags.IOVDb.GlobalTag="CONDBR2-BLKPA-2017-05"
    ConfigFlags.lock()

    result=MainServicesCfg(ConfigFlags)

    McEventSelector=CompFactory.McEventSelector
    McCnvSvc=CompFactory.McCnvSvc
    EvtPersistencySvc=CompFactory.EvtPersistencySvc

    #event & time-stamp from the q431 test input
    mcevtsel=McEventSelector(RunNumber=330470,
                             EventsPerRun=1,
                             FirstEvent=1183722158,
                             FirstLB=310,
                             EventsPerLB=1,
                             InitialTimeStamp=1500867637,
                             TimeStampInterval=1
                         )

    result.addService(mcevtsel)
    result.setAppProperty("EvtSel",mcevtsel.getFullJobOptName())

    mccnvsvc=McCnvSvc()
    result.addService(mccnvsvc)


    result.addService(EvtPersistencySvc("EventPersistencySvc",CnvServices=[mccnvsvc.getFullJobOptName(),]))

    result.merge(BunchCrossingCondAlgCfg(ConfigFlags))

    BunchCrossingCondTest=CompFactory.BunchCrossingCondTest
    result.addEventAlgo(BunchCrossingCondTest(FileName="BCData1.txt"))

    result.run(1)

    #f=open("test.pkl","wb")
    #result.store(f)
    #f.close()
