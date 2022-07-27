# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory
from ByteStreamCnvSvc.ByteStreamConfig import ByteStreamReadCfg
from AthenaConfiguration.MainServicesConfig import MainServicesCfg
from LArConfiguration.LArConfigFlags import RawChannelSource 


def LArRawDataReadingCfg(configFlags, **kwargs):
    acc=ComponentAccumulator()
    from LArGeoAlgsNV.LArGMConfig import LArGMCfg
    acc.merge(LArGMCfg(configFlags))
    acc.merge(ByteStreamReadCfg(configFlags))

    if configFlags.Overlay.DataOverlay:
        kwargs.setdefault("LArDigitKey", configFlags.Overlay.BkgPrefix + "FREE")

    if configFlags.LAr.RawChannelSource is RawChannelSource.Calculated or configFlags.Overlay.DataOverlay:
        kwargs.setdefault("LArRawChannelKey", "")

    print('LArRawDataReadingCfg configFlags.LAr.RawChannelSource ',configFlags.LAr.RawChannelSource)

    acc.addEventAlgo(CompFactory.LArRawDataReadingAlg(**kwargs))
    return acc


if __name__=="__main__":

    from AthenaConfiguration.AllConfigFlags import ConfigFlags
    from AthenaCommon.Logging import log
    from AthenaCommon.Constants import DEBUG
    log.setLevel(DEBUG)

    from AthenaConfiguration.TestDefaults import defaultTestFiles
    ConfigFlags.LAr.doAlign=False
    ConfigFlags.Input.Files = defaultTestFiles.RAW
    ConfigFlags.lock()

    acc = MainServicesCfg( ConfigFlags )
    acc.merge(LArRawDataReadingCfg(ConfigFlags))
    
    DumpLArRawChannels=CompFactory.DumpLArRawChannels
    from LArCabling.LArCablingConfig import LArOnOffIdMappingCfg 
    acc.merge(LArOnOffIdMappingCfg(ConfigFlags))
    acc.addEventAlgo(DumpLArRawChannels(LArRawChannelContainerName="LArRawChannels",))

    acc.run(2,OutputLevel=DEBUG)

