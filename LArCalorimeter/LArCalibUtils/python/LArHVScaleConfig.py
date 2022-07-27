# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory
from AthenaConfiguration.Enums import LHCPeriod


def LArHVScaleCfg(configFlags):
    result=ComponentAccumulator()

    from IOVDbSvc.IOVDbSvcConfig import addFolders
    LArHVCondAlg=CompFactory.LArHVCondAlg

    if configFlags.Input.isMC:
        result.merge(addFolders(configFlags,["/LAR/Identifier/HVLineToElectrodeMap<tag>LARHVLineToElectrodeMap-001</tag>"], "LAR_OFL", className="AthenaAttributeList"))

        LArHVIdMappingAlg=CompFactory.LArHVIdMappingAlg
        hvmapalg = LArHVIdMappingAlg(ReadKey="/LAR/Identifier/HVLineToElectrodeMap",WriteKey="LArHVIdMap")
        result.addCondAlgo(hvmapalg)

        result.addCondAlgo(LArHVCondAlg(doHV=False, doAffectedHV=False))

    elif not configFlags.Common.isOnline:
        result.merge(addFolders(configFlags,["/LAR/DCS/HV/BARREl/I16"], "DCS_OFL", className="CondAttrListCollection"))
        result.merge(addFolders(configFlags,["/LAR/DCS/HV/BARREL/I8"],  "DCS_OFL", className="CondAttrListCollection"))

        result.merge(addFolders(configFlags,["/LAR/IdentifierOfl/HVLineToElectrodeMap"], "LAR_OFL", className="AthenaAttributeList"))
        result.merge(addFolders(configFlags,["/LAR/HVPathologiesOfl/Pathologies"], "LAR_OFL", className="AthenaAttributeList"))
        if configFlags.GeoModel.Run is not LHCPeriod.Run1:
            result.merge(addFolders(configFlags,["/LAR/HVPathologiesOfl/Rvalues"], "LAR_OFL", className="AthenaAttributeList"))

        from LArBadChannelTool.LArBadChannelConfig import LArBadChannelCfg, LArBadFebCfg
        result.merge(LArBadChannelCfg(configFlags))
        result.merge(LArBadFebCfg(configFlags))

        LArHVIdMappingAlg=CompFactory.LArHVIdMappingAlg
        hvmapalg = LArHVIdMappingAlg(ReadKey="/LAR/IdentifierOfl/HVLineToElectrodeMap",WriteKey="LArHVIdMap")
        result.addCondAlgo(hvmapalg)

        LArHVPathologyDbCondAlg=CompFactory.LArHVPathologyDbCondAlg
        hvpath = LArHVPathologyDbCondAlg(PathologyFolder="/LAR/HVPathologiesOfl/Pathologies",
                                         HVMappingKey="LArHVIdMap",
                                         HVPAthologyKey="LArHVPathology")
        result.addCondAlgo(hvpath)

        from LArConfiguration.LArElecCalibDBConfig import LArElecCalibDbCfg
        result.merge(LArElecCalibDbCfg(configFlags,["HVScaleCorr",]))

        if configFlags.GeoModel.Run is not LHCPeriod.Run1:
           hvcond = LArHVCondAlg(HVPathologies="LArHVPathology")
        else:   
           hvcond = LArHVCondAlg(HVPathologies="LArHVPathology",doR=False)

        hvcond.UndoOnlineHVCorr=True
        hvcond.keyOutputCorr= "LArHVScaleCorrRecomputed"
        result.addCondAlgo(hvcond)

    return result

if __name__=="__main__":
    from AthenaConfiguration.AllConfigFlags import ConfigFlags

    nThreads=1
    ConfigFlags.Concurrency.NumThreads = nThreads
    if nThreads>0:
        ConfigFlags.Scheduler.ShowDataDeps = True
        ConfigFlags.Scheduler.ShowDataFlow = True
        ConfigFlags.Scheduler.ShowControlFlow = True
        ConfigFlags.Concurrency.NumConcurrentEvents = nThreads

    ConfigFlags.Input.Files = ["myESD-data.pool.root"]
    ConfigFlags.lock()

    from AthenaConfiguration.MainServicesConfig import MainServicesCfg
    from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
    cfg=MainServicesCfg(ConfigFlags)
    cfg.merge(PoolReadCfg(ConfigFlags))

    cfg.merge( LArHVScaleCfg(ConfigFlags) )

    cfg.run(10)
