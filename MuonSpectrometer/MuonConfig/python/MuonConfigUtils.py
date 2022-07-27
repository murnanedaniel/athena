# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

# This file is just for shared functions etc used by this package.
from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory

def SetupMuonStandaloneArguments():
    from argparse import ArgumentParser
    
    parser = ArgumentParser()
    parser.add_argument("-t", "--threads", dest="threads", type=int,
                        help="number of threads", default=1)
                        
    parser.add_argument("-o", "--output", dest="output", default='newESD.pool.root',
                        help="write ESD to FILE", metavar="FILE")
                        
    parser.add_argument("--run", help="Run directly from the python. If false, just stop once the pickle is written.",
                        action="store_true")
                        
    parser.add_argument("--forceclone", help="Override default cloneability of algorithms to force them to run in parallel",
                        action="store_true")
    parser.add_argument("-d","--debug", default=None, help="attach debugger (gdb) before run, <stage>: conf, init, exec, fini")
    parser.add_argument("--input", "-i", help="Input file to run the config", nargs="+",
                                        default= ['/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/RecExRecoTest/ESD.16747874._000011_100events.pool.root'])
    args = parser.parse_args()
    
    return args
    
def SetupMuonStandaloneConfigFlags(args):
    from AthenaConfiguration.AllConfigFlags import ConfigFlags
    # Keeping this commented out so we can easily switch to the default for testing against that.
    # from AthenaConfiguration.TestDefaults import defaultTestFiles
    # ConfigFlags.Input.Files = defaultTestFiles.ESD
    ConfigFlags.Input.Files = args.input
    
    ConfigFlags.Concurrency.NumThreads=args.threads
    ConfigFlags.Concurrency.NumConcurrentEvents=args.threads # Might change this later, but good enough for the moment.

    ConfigFlags.Detector.GeometryMDT   = True 
    ConfigFlags.Detector.GeometryTGC   = True
    ConfigFlags.Detector.GeometryCSC   = True     
    ConfigFlags.Detector.GeometryRPC   = True
    # TODO: disable these for now, to be determined if needed
    ConfigFlags.Detector.GeometryCalo  = False
    ConfigFlags.Detector.GeometryID    = False

    # FIXME This is temporary. I think it can be removed with some other refactoring
    ConfigFlags.Muon.makePRDs          = False

    ConfigFlags.Output.ESDFileName=args.output
    
    ConfigFlags.Input.isMC = True
    ConfigFlags.lock()
    ConfigFlags.dump()
    return ConfigFlags
    
def SetupMuonStandaloneCA(args,ConfigFlags):
    # When running from a pickled file, athena inserts some services automatically. So only use this if running now.
    if args.run:
        from AthenaConfiguration.MainServicesConfig import MainServicesCfg
        cfg = MainServicesCfg(ConfigFlags)
        msgService = cfg.getService('MessageSvc')
        msgService.Format = "S:%s E:%e % F%128W%S%7W%R%T  %0W%M"
    else:
        cfg=ComponentAccumulator()

    from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
    cfg.merge(PoolReadCfg(ConfigFlags))

    if ConfigFlags.Input.isMC:
        from xAODTruthCnv.xAODTruthCnvConfigNew import GEN_AOD2xAODCfg
        cfg.merge(GEN_AOD2xAODCfg(ConfigFlags))

    return cfg
    
def SetupMuonStandaloneOutput(cfg, ConfigFlags, itemsToRecord):
    # Set up output
    from OutputStreamAthenaPool.OutputStreamConfig import OutputStreamCfg

    cfg.merge( OutputStreamCfg( ConfigFlags, 'ESD', ItemList=itemsToRecord) )
    outstream = cfg.getEventAlgo("OutputStreamESD")
    outstream.ForceRead = True

    # Fix for ATLASRECTS-5151
    Trk__EventCnvSuperTool = CompFactory.Trk.EventCnvSuperTool
    cnvTool = Trk__EventCnvSuperTool(name = 'EventCnvSuperTool')
    cnvTool.MuonCnvTool.FixTGCs = True 
    cfg.addPublicTool(cnvTool)
