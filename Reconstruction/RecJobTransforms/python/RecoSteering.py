# Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
from AthenaConfiguration.MainServicesConfig import MainServicesCfg
from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaCommon.Logging import logging
log = logging.getLogger("RecoSteering.py")

def RecoSteering(flags):
    """
    Generates configuration of the reconstructions

    This driver configures all reconstruction steps unconditionally.
    The selftest available below can be used for simple jobs, yet full functionality is achieved with tansforms that set many flags.
    """
    acc = ComponentAccumulator()
    # setup input    
    if flags.Input.Format == 'BS':
        from ByteStreamCnvSvc.ByteStreamConfig import ByteStreamReadCfg
        acc.merge(ByteStreamReadCfg(flags))
        log.info("---------- Configured BS reading")
    else:
        from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
        acc.merge(PoolReadCfg(flags))
        log.info("---------- Configured POOL reading")

    # calorimeter    
    from CaloRec.CaloRecoConfig import CaloRecoCfg
    acc.merge(CaloRecoCfg(flags))
    log.info("---------- Configured calorimeter reconstruction")

    # ID    
    from InDetConfig.TrackRecoConfig import TrackRecoCfg
    acc.merge(TrackRecoCfg(flags))
    log.info("---------- Configured tracking")

    # muons
    from MuonConfig.MuonReconstructionConfig import MuonReconstructionCfg
    acc.merge(MuonReconstructionCfg(flags))
    log.info("---------- Configured muon tracking")

    # physics objects
    # egamma
    # jets
    # btagging
    # pflow

    #setup output
    if flags.Output.doWriteESD:
        log.info("ESD ItemList: %s", acc.getEventAlgo("OutputStreamESD").ItemList)
        log.info("---------- Configured ESD writing")

    if flags.Output.doWriteAOD:
        log.info("ESD ItemList: %s", acc.getEventAlgo("OutputStreamAOD").ItemList)
        log.info("---------- Configured AOD writing")

    return acc

import sys
def _run(input):
    import datetime
    startStamp = datetime.datetime.now()
    from AthenaCommon.Configurable import Configurable
    Configurable.configurableRun3Behavior=True
    from AthenaConfiguration.AllConfigFlags import ConfigFlags as flags
    flags.Exec.MaxEvents = 5
    flags.Concurrency.NumThreads = 1
    flags.Input.Files = []
    #TODO these flags should be defaulted in the divier function above, 
    #TODO    but then we ought to have option to set them from command line should the parser be passed there too?
    flags.Detector.RecoIBL=True
    flags.Detector.RecoPixel=True
    flags.Detector.RecoSCT=True
    flags.Detector.RecoTRT=True

    flags.Detector.GeometryBCM=True
    flags.Detector.GeometryDBM=True
    flags.Detector.GeometryPixel=True
    flags.Detector.GeometrySCT=True
    flags.Detector.GeometryTRT=True

    flags.Detector.GeometryTile=True
    flags.Detector.GeometryLAr=True

    flags.Calo.TopoCluster.doTopoClusterLocalCalib=False
    flags.Output.ESDFileName="outputESD.pool.root"
    flags.Output.AODFileName="outputAOD.pool.root"
    parser = flags.getArgumentParser()
    args = flags.fillFromArgs(parser=parser)

    # test inputs
    if flags.Input.Files == []:
        from AthenaConfiguration.TestDefaults import defaultTestFiles
        if input == "RAW":
            flags.Input.Files = defaultTestFiles.RAW
        if input == "ESD":
            flags.Input.Files = defaultTestFiles.ESD

    flags.lock()

    acc = MainServicesCfg(flags)
    acc.merge(RecoSteering(flags), sequenceName="AthAlgSeq")
    confStamp = datetime.datetime.now()
    log.info("configured in %d seconds", (confStamp-startStamp).seconds )
    flags.dump()
    acc.printConfig(withDetails=True)

    if args.configOnly:
        with open(args.configOnly, "wb") as confFile:
            acc.store(confFile)
            log.info("configOnly option specified. Saved in: %s ... exiting now.", args.configOnly )
    
            sys.exit(0)
    # running        
    statusCode = acc.run()
    endStamp = datetime.datetime.now()
    log.info("total time spent in %d seconds (running %s seconds) ", (endStamp-startStamp).seconds, (endStamp-confStamp).seconds )
    return statusCode


if __name__ == "__main__":
    statusCode = None
    if "--RAW" in sys.argv:
        del sys.argv[sys.argv.index("--RAW")]
        statusCode = _run(input="RAW")

    if "--ESD" in sys.argv:    
        del sys.argv[sys.argv.index("--ESD")]
        statusCode = _run(input="ESD")

#TODO enable digest production once able to read the AOD
#        if statusCode.isSuccess():
#            import subprocess
#            subprocess.run("xAODDigest.py outputAOD.pool.root digetst.txt")

    sys.exit(not statusCode.isSuccess())
