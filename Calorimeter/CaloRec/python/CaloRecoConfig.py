
# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory
from AthenaConfiguration.Enums import Format


def CaloRecoCfg(configFlags, clustersname=None):
    result = ComponentAccumulator()
    if configFlags.Input.Format is Format.BS:
        #Data-case: Schedule ByteStream reading for LAr & Tile
        from LArByteStream.LArRawDataReadingConfig import LArRawDataReadingCfg
        result.merge(LArRawDataReadingCfg(configFlags))

        from ByteStreamCnvSvc.ByteStreamConfig import ByteStreamReadCfg

        result.merge(ByteStreamReadCfg(configFlags,type_names=['TileDigitsContainer/TileDigitsCnt',
                                                               'TileRawChannelContainer/TileRawChannelCnt',
                                                               'TileMuonReceiverContainer/TileMuRcvCnt']))
        result.getService("ByteStreamCnvSvc").ROD2ROBmap=["-1"]
        if configFlags.Output.doWriteESD:
            from TileRecAlgs.TileDigitsFilterConfig import TileDigitsFilterOutputCfg
            result.merge(TileDigitsFilterOutputCfg(configFlags))
        else: #Mostly for wrapping in RecExCommon
            from TileRecAlgs.TileDigitsFilterConfig import TileDigitsFilterCfg
            result.merge(TileDigitsFilterCfg(configFlags))

        from LArROD.LArRawChannelBuilderAlgConfig import LArRawChannelBuilderAlgCfg
        result.merge(LArRawChannelBuilderAlgCfg(configFlags))

        from TileRecUtils.TileRawChannelMakerConfig import TileRawChannelMakerCfg
        result.merge(TileRawChannelMakerCfg(configFlags))

    if not configFlags.Input.isMC and not configFlags.Common.isOnline:
        from LArCellRec.LArTimeVetoAlgConfig import LArTimeVetoAlgCfg
        result.merge(LArTimeVetoAlgCfg(configFlags))

    if not configFlags.Input.isMC and not configFlags.Overlay.DataOverlay:
        from LArROD.LArFebErrorSummaryMakerConfig import LArFebErrorSummaryMakerCfg
        result.merge(LArFebErrorSummaryMakerCfg(configFlags))


    #Configure cell-building
    from CaloRec.CaloCellMakerConfig import CaloCellMakerCfg
    result.merge(CaloCellMakerCfg(configFlags))

    #Configure topo-cluster builder
    from CaloRec.CaloTopoClusterConfig import CaloTopoClusterCfg
    result.merge(CaloTopoClusterCfg(configFlags, clustersname=clustersname))

    #Configure forward towers:
    from CaloRec.CaloFwdTopoTowerConfig import CaloFwdTopoTowerCfg
    result.merge(CaloFwdTopoTowerCfg(configFlags,CaloTopoClusterContainerKey="CaloCalTopoClusters"))

    #Configure NoisyROSummary
    from LArCellRec.LArNoisyROSummaryConfig import LArNoisyROSummaryCfg
    result.merge(LArNoisyROSummaryCfg(configFlags))

    #Configure TileLookForMuAlg
    from TileMuId.TileMuIdConfig import TileLookForMuAlgCfg
    result.merge(TileLookForMuAlgCfg(configFlags))

    if not configFlags.Input.isMC and not configFlags.Overlay.DataOverlay:
        #Configure LArDigitsThinner:
        from LArROD.LArDigitThinnerConfig import LArDigitThinnerCfg
        result.merge(LArDigitThinnerCfg(configFlags))

    #Configure MBTSTimeDiff
    #Clients are BackgroundWordFiller and (deprecated?) DQTBackgroundMonTool
    #Consider moving to BackgroundWordFiller config
    if configFlags.Detector.GeometryMBTS:
        from TileRecAlgs.MBTSTimeDiffEventInfoAlgConfig import MBTSTimeDiffEventInfoAlgCfg
        result.merge(MBTSTimeDiffEventInfoAlgCfg(configFlags))


    #Configure AOD Cell-Thinning based on samplings:
    from CaloRec.CaloThinCellsBySamplingAlgConfig import CaloThinCellsBySamplingAlgCfg
    result.merge(CaloThinCellsBySamplingAlgCfg(configFlags,'StreamAOD', ['TileGap3']))
    

    return result


def CaloRecoDebuggingCfg(configFlags):
    result = ComponentAccumulator()

    result.addEventAlgo(CompFactory.DumpLArRawChannels(LArRawChannelContainerName="LArRawChannels_FromDigits"))
    result.addEventAlgo(CompFactory.CaloCellDumper())

    ClusterDumper = CompFactory.ClusterDumper
    result.addEventAlgo(ClusterDumper("TopoDumper", ContainerName="CaloCalTopoClusters", FileName="TopoCluster.txt"))
    result.addEventAlgo(ClusterDumper("FwdTopoDumper", ContainerName="CaloCalFwdTopoTowers", FileName="FwdTopoCluster.txt"))

    return result


if __name__=="__main__":
    from AthenaConfiguration.AllConfigFlags import ConfigFlags
    from AthenaCommon.Logging import log
    from AthenaCommon.Constants import DEBUG,WARNING
    log.setLevel(DEBUG)

    ConfigFlags.Input.Files = ["/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/data17_13TeV.00330470.physics_Main.daq.RAW._lb0310._SFO-1._0001.data",]

    ConfigFlags.fillFromArgs()

    ConfigFlags.lock()

    from AthenaConfiguration.MainServicesConfig import MainServicesCfg
    acc = MainServicesCfg(ConfigFlags)

    acc.merge(CaloRecoCfg(ConfigFlags))


    CaloCellDumper=CompFactory.CaloCellDumper
    acc.addEventAlgo(CaloCellDumper())

    ClusterDumper=CompFactory.ClusterDumper
    acc.addEventAlgo(ClusterDumper("TopoDumper",ContainerName="CaloCalTopoClusters",FileName="TopoCluster.txt"))

    f=open("CaloRec.pkl","wb")
    acc.store(f)
    f.close()

    acc.run(10,OutputLevel=WARNING)
