# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

if __name__=="__main__":

    from AthenaConfiguration.AllConfigFlags import ConfigFlags as cfgFlags

    cfgFlags.Concurrency.NumThreads=8
    cfgFlags.Exec.MaxEvents=100
    cfgFlags.Input.isMC=True
    cfgFlags.Input.Files= ["/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/PFlowTests/mc16_13TeV/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.recon.ESD.e6337_e5984_s3170_r12674/ESD.25732025._000034.pool.root.1"]
    cfgFlags.Output.AODFileName="output_AOD.root"
    cfgFlags.Output.doWriteAOD=True
    cfgFlags.lock()

    from AthenaConfiguration.MainServicesConfig import MainServicesCfg
    cfg=MainServicesCfg(cfgFlags)

    from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
    cfg.merge(PoolReadCfg(cfgFlags))

    from AtlasGeoModel.GeoModelConfig import GeoModelCfg
    cfg.merge(GeoModelCfg(cfgFlags))

    from LArGeoAlgsNV.LArGMConfig import LArGMCfg
    cfg.merge(LArGMCfg(cfgFlags))

    from TileGeoModel.TileGMConfig import TileGMCfg
    cfg.merge(TileGMCfg(cfgFlags))

    from PFlowUtils.PFlowCalibHitDecoratorCfg import PFlowCalibHitDecoratorCfg
    cfg.merge(PFlowCalibHitDecoratorCfg())
    cfg.getEventAlgo("PFlowCalibPFODecoratorAlgorithm").PFOWriteDecorHandleKey_NLeadingTruthParticles="JetETMissNeutralFlowElements.calpfo_NLeadingTruthParticleBarcodeEnergyPairs"

    cfg.run()
