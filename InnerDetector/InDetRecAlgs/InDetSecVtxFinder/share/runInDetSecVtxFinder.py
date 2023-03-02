# Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
# Define method to construct configures Sec Vtx Finder alg
# attempted by N Ribaric (@LancasterUNI) neza.ribaric@cern.ch

if __name__ == "__main__":
    
    from AthenaCommon.Configurable import Configurable
    Configurable.configurableRun3Behavior = 1
    import AthenaCommon.Constants as Lvl

    # import the flags and set them
    from AthenaConfiguration.AllConfigFlags import ConfigFlags

    ConfigFlags.Exec.MaxEvents = 50

    # use one of the predefined files
    from AthenaConfiguration.TestDefaults import defaultTestFiles
    ConfigFlags.Input.Files = defaultTestFiles.AOD_RUN3_MC
    ConfigFlags.Input.isMC=True

    # lock the flags
    ConfigFlags.lock()

    # create basic infrastructure
    from AthenaConfiguration.MainServicesConfig import MainServicesCfg
    acc = MainServicesCfg(ConfigFlags)

    from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
    acc.merge(PoolReadCfg(ConfigFlags)) #!!!!! THIS IS IMPORTANT

    # add the algorithm to the configuration
    from InDetSecVtxFinder.InDetSecVtxFinderConfig import InDetSecVtxFinderAlgCfg
    acc.merge(InDetSecVtxFinderAlgCfg(ConfigFlags,
              name = "InDetSecVtxFinder",
              FinderTool = "AMVF",
              useTrackParticles = True,
              inputTrackParticles = "InDetTrackParticles",
              outputSecondaryVertices = "RecoSecVtx",
              doVertexMerging = False,
              OutputLevel = Lvl.INFO))

    # Contents
    from DerivationFrameworkCore.SlimmingHelper import SlimmingHelper
    TRUTH0SlimmingHelper = SlimmingHelper("TRUTH0SlimmingHelper", NamesAndTypes = ConfigFlags.Input.TypedCollections, ConfigFlags = ConfigFlags)
    TRUTH0SlimmingHelper.AppendToDictionary = {'EventInfo':'xAOD::EventInfo','EventInfoAux':'xAOD:EventAuxInfo',
                                               'TruthEvents':'xAOD::TruthEventContainer','TruthEventsAux':'xAOD::TruthEventAuxContainer',
                                               'TruthVertices':'xAOD::TruthVertexContainer','TruthVerticesAux':'xAOD::TruthVertexAuxContainer',
                                               'TruthParticles':'xAOD::TruthParticleContainer','TruthParticlesAux':'xAOD::TruthParticleAuxContainer'} 

    TRUTH0SlimmingHelper.AllVariables = [ 'EventInfo',
                                          'TruthEvents', 
                                          'TruthVertices',
                                          'TruthParticles']

    # Metadata
    TRUTH0MetaDataItems = [ "xAOD::TruthMetaDataContainer#TruthMetaData", "xAOD::TruthMetaDataAuxContainer#TruthMetaDataAux." ]

    # Create output stream 
    from OutputStreamAthenaPool.OutputStreamConfig import OutputStreamCfg
    TRUTH0ItemList = TRUTH0SlimmingHelper.GetItemList()
    TRUTH0ItemList+=["xAOD::VertexContainer#RecoSecVtx","xAOD::VertexContainer#RecoSecVtxAux."]
    acc.merge(OutputStreamCfg(ConfigFlags, "SecondayVertexOutput", ItemList=TRUTH0ItemList))

    # debug printout
    acc.printConfig(withDetails=True, summariseProps=True)

    # run the job
    status = acc.run()

    # report the execution status (0 ok, else error)
    import sys
    sys.exit(not status.isSuccess())
