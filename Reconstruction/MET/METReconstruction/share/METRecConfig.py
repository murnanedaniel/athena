from AthenaCommon import Loggingfrom AthenaConfiguration.ComponentAccumulator import ComponentAccumulatorif __name__=="__main__":    # Setting needed for the ComponentAccumulator to do its thing    from AthenaCommon.Configurable import Configurable    Configurable.configurableRun3Behavior=True        # Set message levels    from AthenaCommon import Constants    msgLvl = "INFO"    from AthenaCommon.Logging import log    log.setLevel(msgLvl)        # Config flags steer the job at various levels    from AthenaConfiguration.AllConfigFlags import ConfigFlags    ConfigFlags.Input.isMC  = True    ConfigFlags.Input.Files = ["/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/ASG/mc16_13TeV.410501.PowhegPythia8EvtGen_A14_ttbar_hdamp258p75_nonallhad.merge.AOD.e5458_s3126_r9364_r9315/AOD.11182705._000001.pool.root.1"]    # Flags relating to multithreaded execution    nthreads=0    ConfigFlags.Concurrency.NumThreads =nthreads    if nthreads>0:    	ConfigFlags.Concurrency.NumThreads = 1    	ConfigFlags.Concurrency.NumConcurrentEvents = 1	    if ConfigFlags.Beam.Type == 'cosmics' or ConfigFlags.Beam.Type == 'singlebeam':# used to have " or not rec.doInDet()" on the end        ConfigFlags.MET.UseTracks=False        ConfigFlags.MET.DoPFlow=False        print "METReconstruction_jobOptions: detected cosmics/single-beam configuration -- switch off track-based MET reco"     ConfigFlags.lock()    # Get a ComponentAccumulator setting up the fundamental Athena job    from AthenaConfiguration.MainServicesConfig import MainServicesThreadedCfg     cfg=MainServicesThreadedCfg(ConfigFlags)     # Add the components for reading in pool files    from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg    cfg.merge(PoolReadCfg(ConfigFlags))    print "CHECKPOINT 1"        # Get Jet Inputs    from JetRecConfig.StandardJetDefs import EMTopoOrigin, LCTopoOrigin, CHSPFlow    from JetRecConfig import JetRecConfig    cfg1 = JetRecConfig.JetInputCfg( [EMTopoOrigin], ConfigFlags)    cfg1.printConfig()    cfg.merge( cfg1 )    cfg2 = JetRecConfig.JetInputCfg( [LCTopoOrigin], ConfigFlags)    cfg2.printConfig()    cfg.merge( cfg2 )    cfg3 = JetRecConfig.JetInputCfg( [CHSPFlow], ConfigFlags)    cfg3.printConfig()    cfg.merge( cfg3 )    print "CHECKPOINT 2"    from METReconstruction.METCfg_Track import METTrack_Cfg    cfg4=METTrack_Cfg(ConfigFlags)    cfg4.printConfig()    cfg.merge(cfg4)        # Start by just trying to add in MET Reconstruction based on METReconstruction_jobOptions.py    #from METReconstruction.METRecoFlags import metFlags    # from AthenaCommon.BeamFlags import jobproperties NO LONGER ALLOWED    # from RecExConfig.RecFlags import rec NO LONGER ALLOWED    #NEED TO CHANGE THIS TO DEPEND ON ConfigFlags.Beam.Type => for now ignore    # TJ: Best to start each of these from scratch as a new CA module,   # Can e.g. make files called METCaloConfig.py that just put the   # old alg into a CA.      # Rather than have N reco tools that get thrown into one alg later,   # have each CA generate its own METRecoAlg and add this to the sequence.       """    import METReconstruction.METConfig_Calo    import METReconstruction.METConfig_Track    if rec.doTruth():        import METReconstruction.METConfig_Truth        from METReconstruction.METRecoConfig import getMETRecoAlg    print "PICKING UP CHANGES"    metAlg = getMETRecoAlg('METReconstruction')    """   # Probably want to define one CA each for EMTopo, LCTopo and PFlow,   # then have a higher level one that merges in all three,   # then the top-level (i.e. this) can just pull in the all-associators CA    """    components_metAlg = ComponentAccumulator()    from AthenaCommon.AlgSequence import AthSequencer    components_metAlg.addSequence( AthSequencer('METReconstruction') ) #technically don't need a new sequence name for it    components_metAlg.addEventAlgo(metAlg,'METReconstruction')    cfg.merge(components_metAlg)    # Set up default configurations    import METReconstruction.METConfig_Associator    from METReconstruction.METAssocConfig import getMETAssocAlg    # Get the configuration directly from METRecoFlags    # Can also provide a dict of configurations or list of RecoTools or both    assocAlg = getMETAssocAlg('METAssociation')    components_assocAlg = ComponentAccumulator()    components_assocAlg.addSequence(AthSequencer('METAssociation') )    components_assocAlg.addEventAlgo(assocAlg,'METAssociation')    cfg.merge(components_assocAlg)    from METUtilities.METMakerConfig import getMETMakerAlg    for key,conf in metFlags.METAssocConfigs().iteritems():        if not conf.doTruth:            makerAlg = getMETMakerAlg(conf.suffix)            components_makerAlg=ComponentAccumulator()            components_makerAlg.addSequence(AthSequencer(conf.suffix) )            components_makerAlg.addEventAlgo(makerAlg,conf.suffix)	            cfg.merge(components_makerAlg)    """    print "Running final component accumulator"    cfg.printConfig()    cfg.run(maxEvents=10)