# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory
from AthenaConfiguration.Enums import BeamType


if __name__=="__main__":
    # Setting needed for the ComponentAccumulator to do its thing
    from AthenaCommon.Configurable import Configurable
    Configurable.configurableRun3Behavior=True

    # Set message levels
    msgLvl = "WARNING"
    from AthenaCommon.Logging import log
    log.setLevel(msgLvl)

    # Config flags steer the job at various levels
    from AthenaConfiguration.AllConfigFlags import ConfigFlags
    ConfigFlags.Input.isMC  = True
    ConfigFlags.Input.Files = ["/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/ASG/mc16_13TeV.410501.PowhegPythia8EvtGen_A14_ttbar_hdamp258p75_nonallhad.merge.AOD.e5458_s3126_r9364_r9315/AOD.11182705._000001.pool.root.1"]

    # Flags relating to multithreaded execution
    nthreads = 0
    ConfigFlags.Concurrency.NumThreads = nthreads
    if nthreads > 0:
        ConfigFlags.Concurrency.NumThreads = 1
        ConfigFlags.Concurrency.NumConcurrentEvents = 1
    ConfigFlags.MET.UseTracks = True
    ConfigFlags.MET.DoPFlow = True
    if ConfigFlags.Beam.Type in [BeamType.Cosmics, BeamType.SingleBeam]: # used to have " or not rec.doInDet()" on the end
        ConfigFlags.MET.UseTracks = False
        ConfigFlags.MET.DoPFlow = False
        print("METReconstruction_jobOptions: detected cosmics/single-beam configuration -- switch off track-based MET reco")

    ConfigFlags.lock()

    # Get a ComponentAccumulator setting up the fundamental Athena job
    from AthenaConfiguration.MainServicesConfig import MainServicesCfg
    cfg=MainServicesCfg(ConfigFlags)

    # Add the components for reading in pool files
    from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
    cfg.merge(PoolReadCfg(ConfigFlags))
    StoreGateSvc=CompFactory.StoreGateSvc
    cfg.addService(StoreGateSvc("DetectorStore"))

    #Setup up general geometry
    modelConfig=ComponentAccumulator()
    from AtlasGeoModel.GeoModelConfig import GeoModelCfg
    modelConfig=GeoModelCfg(ConfigFlags)
    cfg.merge(modelConfig)

    from MagFieldServices.MagFieldServicesConfig import MagneticFieldSvcCfg
    cfg.merge(MagneticFieldSvcCfg(ConfigFlags))

    from TrkConfig.AtlasTrackingGeometrySvcConfig import TrackingGeometrySvcCfg
    cfg.merge(TrackingGeometrySvcCfg(ConfigFlags))

    from MuonConfig.MuonGeometryConfig import MuonGeoModelCfg
    cfg.merge(MuonGeoModelCfg(ConfigFlags))

    # Nowadays the jet calibration tool requires the EventInfo
    # to be decorated with lumi info, which is not in Run 2 AODs
    from LumiBlockComps.LuminosityCondAlgConfig import LuminosityCondAlgCfg
    cfg.merge(LuminosityCondAlgCfg(ConfigFlags))

    from AthenaConfiguration.ComponentFactory import CompFactory
    muWriter = CompFactory.LumiBlockMuWriter("LumiBlockMuWriter",LumiDataKey="LuminosityCondData")
    cfg.addEventAlgo(muWriter,"AthAlgSeq")

    # Get Jet Inputs
    from JetRecConfig.StandardJetConstits import stdConstitDic as cst
    from JetRecConfig import JetRecConfig
    for jetdef in [cst.EMTopoOrigin,cst.LCTopoOrigin,cst.EMPFlow]:
        cfg.merge(JetRecConfig.JetInputCfg([jetdef], ConfigFlags))

    # Need to rename the collections in the xAOD in order to avoid conflicts
    from SGComps.AddressRemappingConfig import InputRenameCfg
    cfg.merge(InputRenameCfg('xAOD::MissingETContainer','MET_Track','MET_Track_Old'))
    cfg.merge(InputRenameCfg('xAOD::MissingETAuxContainer','MET_TrackAux.','MET_Track_OldAux.'))
    cfg.merge(InputRenameCfg('xAOD::MissingETContainer','MET_EMTopo','MET_EMTopo_Old'))
    cfg.merge(InputRenameCfg('xAOD::MissingETAuxContainer','MET_EMTopoAux.','MET_EMTopo_OldAux.'))
    cfg.merge(InputRenameCfg('xAOD::MissingETContainer','MET_LocHadTopo','MET_LocHadTopo_Old'))
    cfg.merge(InputRenameCfg('xAOD::MissingETAuxContainer','MET_LocHadTopoAux.','MET_LocHadTopo_OldAux.'))

    from METReconstruction.METTrack_Cfg import METTrack_Cfg
    cfg.merge(METTrack_Cfg(ConfigFlags))

    from METReconstruction.METCalo_Cfg import METCalo_Cfg
    cfg.merge(METCalo_Cfg(ConfigFlags))

    if ConfigFlags.Input.isMC:
        from METReconstruction.METTruth_Cfg import METTruth_Cfg
        cfg.merge(METTruth_Cfg(ConfigFlags))

    from METReconstruction.METAssociatorCfg import METAssociatorCfg
    cfg.merge(METAssociatorCfg(ConfigFlags, "AntiKt4EMTopo"))
    cfg.merge(METAssociatorCfg(ConfigFlags, "AntiKt4LCTopo"))
    if ConfigFlags.MET.DoPFlow:
        cfg.merge(METAssociatorCfg(ConfigFlags, "AntiKt4EMPFlow"))

    outputlist = ["EventInfo#*"]
    outputlist+=["xAOD::MissingETContainer#"+"MET_Track","xAOD::MissingETAuxContainer#"+"MET_Track"+"Aux."]
    outputlist+=["xAOD::MissingETContainer#"+"MET_Track_Old","xAOD::MissingETAuxContainer#"+"MET_Track_Old"+"Aux."]
    outputlist+=["xAOD::MissingETContainer#"+"MET_EMTopo","xAOD::MissingETAuxContainer#"+"MET_EMTopo"+"Aux."]
    outputlist+=["xAOD::MissingETContainer#"+"MET_EMTopo_Old","xAOD::MissingETAuxContainer#"+"MET_EMTopo_Old"+"Aux."]
    outputlist+=["xAOD::MissingETContainer#"+"MET_AntiKt4EMPFlow","xAOD::MissingETAuxContainer#"+"MET_AntiKt4EMPFlow"+"Aux."]
    from OutputStreamAthenaPool.OutputStreamConfig import OutputStreamCfg
    cfg.merge(OutputStreamCfg(ConfigFlags,"xAOD",ItemList=outputlist))

    # Optionally, print the contents of the store every event
    cfg.getService("StoreGateSvc").Dump = False

    #cfg.printConfig()
    cfg.run(maxEvents=20)
