# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

"""Define methods to construct configured TileHid2ReSrcIDCondAlg conditions algorithm"""

from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory

def TileHid2RESrcIDCondAlgCfg(flags, **kwargs):
    """Return component accumulator with configured TileHid2ReSrcIDCondAlg conditions algorithm"""

    forHLT = kwargs.get('ForHLT', False)
    hid2RESrcID = 'TileHid2RESrcIDHLT' if forHLT else 'TileHid2RESrcID'
    kwargs.setdefault('TileHid2RESrcID', hid2RESrcID)
    kwargs.setdefault('name', f'{hid2RESrcID}CondAlg')
    
    acc = ComponentAccumulator()

    from TileGeoModel.TileGMConfig import TileGMCfg
    acc.merge( TileGMCfg(flags) )

    TileHid2ReSrcIDCondAlg = CompFactory.TileHid2RESrcIDCondAlg
    acc.addCondAlgo( TileHid2ReSrcIDCondAlg(**kwargs) )

    return acc



if __name__ == "__main__":

    from AthenaConfiguration.AllConfigFlags import ConfigFlags
    from AthenaConfiguration.TestDefaults import defaultTestFiles
    from AthenaCommon.Logging import log
    from AthenaCommon.Constants import INFO
    
    # Test setup
    log.setLevel(INFO)

    ConfigFlags.Input.Files = defaultTestFiles.RAW
    ConfigFlags.lock()

    # Initialize configuration object, add accumulator, merge, and run.
    from AthenaConfiguration.MainServicesConfig import MainServicesCfg
    cfg = MainServicesCfg(ConfigFlags)

    from ByteStreamCnvSvc.ByteStreamConfig import ByteStreamReadCfg
    cfg.merge( ByteStreamReadCfg(ConfigFlags) )

    cfg.merge( TileHid2RESrcIDCondAlgCfg(ConfigFlags, ForHLT=True) )

    cfg.printConfig(withDetails = True, summariseProps = True)
    cfg.store( open('TileHid2ReSrcIDCondAlg.pkl','wb') )

    sc = cfg.run(3)

    import sys
    # Success should be 0
    sys.exit(not sc.isSuccess())
