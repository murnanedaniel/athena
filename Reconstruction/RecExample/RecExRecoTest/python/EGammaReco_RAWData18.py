# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

if __name__=="__main__":
    from AthenaConfiguration.AllConfigFlags import ConfigFlags
    
    ConfigFlags.Input.Files = ['/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/RecExRecoTest/data18_13TeV/data18_13TeV.00357750.physics_Main.daq.RAW/data18_13TeV.00357750.physics_Main.daq.RAW._lb0083._SFO-1._0001.data']
    from egammaConfig.egammaOnlyFromRawFlags import egammaOnlyFromRaw
    egammaOnlyFromRaw(ConfigFlags)
    ConfigFlags.lock()

    from RecJobTransforms.RecoSteering import RecoSteering
    acc = RecoSteering(ConfigFlags)

    with open("config.pkl", "wb") as file:
      acc.store(file)

    acc.run(100)
