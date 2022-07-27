#
#  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
#

'''
@file TileMonitoringConfig.py
@brief Python configuration of Tile Monitoring for the Run III
'''

def TileMonitoringCfg(flags):
    ''' Function to configure Tile Monitoring in the monitoring system for Run III.'''

    from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
    acc = ComponentAccumulator()

    from AthenaCommon.Logging import logging
    msg = logging.getLogger( 'TileMonitoringCfg' )

    environment = flags.DQ.Environment

    if environment in ('online', 'tier0', 'tier0Raw'):
        msg.info('Setup Tile Monitoring for RAW data due to environment: %s', environment)

        from TileMonitoring.TileDQFragMonitorAlgorithm import TileDQFragMonitoringConfig
        acc.merge( TileDQFragMonitoringConfig(flags) )

        from TileMonitoring.TileMBTSMonitorAlgorithm import TileMBTSMonitoringConfig
        acc.merge( TileMBTSMonitoringConfig(flags) )

        from TileMonitoring.TileDigiNoiseMonitorAlgorithm import TileDigiNoiseMonitoringConfig
        acc.merge( TileDigiNoiseMonitoringConfig(flags) )

    if environment in ('online', 'tier0', 'tier0ESD'):
        msg.info('Setup Tile Monitoring for ESD data due to environment: %s', environment)

        from TileMonitoring.TileCellMonitorAlgorithm import TileCellMonitoringConfig
        acc.merge( TileCellMonitoringConfig(flags) )

        from TileMonitoring.TileTowerMonitorAlgorithm import TileTowerMonitoringConfig
        acc.merge( TileTowerMonitoringConfig(flags) )

        from TileMonitoring.TileClusterMonitorAlgorithm import TileClusterMonitoringConfig
        acc.merge( TileClusterMonitoringConfig(flags) )

        from TileMonitoring.TileMuIdMonitorAlgorithm import TileMuIdMonitoringConfig
        acc.merge( TileMuIdMonitoringConfig(flags) )

        from TileMonitoring.TileJetMonitorAlgorithm import TileJetMonitoringConfig
        acc.merge( TileJetMonitoringConfig(flags) )

        if flags.IOVDb.DatabaseInstance == 'CONDBR2' and flags.DQ.triggerDataAvailable:
            from TileMonitoring.TileTMDBRawChannelMonitorAlgorithm import TileTMDBRawChannelMonitoringConfig
            acc.merge( TileTMDBRawChannelMonitoringConfig(flags, FillRawChannelHistograms = False, FillEfficiencyHistograms = True) )

        from AthenaConfiguration.Enums import BeamType
        if flags.Beam.Type in [BeamType.Cosmics, BeamType.SingleBeam]:
            from TileCosmicAlgs.TileMuonFitterConfig import TileMuonFitterCfg
            acc.merge(TileMuonFitterCfg(flags))

            from TileMonitoring.TileMuonFitMonitorAlgorithm import TileMuonFitMonitoringConfig
            acc.merge( TileMuonFitMonitoringConfig(flags) )

    return acc



if __name__=='__main__':

    from AthenaConfiguration.AllConfigFlags import ConfigFlags
    from AthenaCommon.Logging import log
    from AthenaCommon.Constants import INFO
    log.setLevel(INFO)

    ConfigFlags.Input.Files = ['/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/q431/22.0/v1/myESD.pool.root']
    ConfigFlags.Output.HISTFileName = 'TileMonitoringOutput.root'
    ConfigFlags.DQ.enableLumiAccess = False
    ConfigFlags.DQ.useTrigger = False
    ConfigFlags.DQ.Environment = 'tier0'
    ConfigFlags.Exec.MaxEvents = 3
    ConfigFlags.fillFromArgs()
    ConfigFlags.lock()

    # Initialize configuration object, add accumulator, merge, and run.
    from AthenaConfiguration.MainServicesConfig import MainServicesCfg
    from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
    acc = MainServicesCfg(ConfigFlags)
    acc.merge(PoolReadCfg(ConfigFlags))

    acc.merge( TileMonitoringCfg(ConfigFlags) )

    acc.printConfig(withDetails = True, summariseProps = True)
    ConfigFlags.dump()
    acc.store(open("TileMonitoring.pkl","wb"))

    sc = acc.run()
    import sys
    # Success should be 0
    sys.exit(not sc.isSuccess())
