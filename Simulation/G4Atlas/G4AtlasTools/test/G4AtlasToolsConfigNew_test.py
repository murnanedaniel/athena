#!/usr/bin/env python
"""Run tests on G4AtlasToolsConfigNew

Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
"""

if __name__ == '__main__':
  from AthenaConfiguration.MainServicesConfig import MainServicesCfg

  # Set up logging
  from AthenaCommon.Logging import log
  from AthenaCommon.Constants import DEBUG
  log.setLevel(DEBUG)


  #import config flags
  from AthenaConfiguration.AllConfigFlags import ConfigFlags
  from AthenaConfiguration.Enums import ProductionStep
  ConfigFlags.Common.ProductionStep = ProductionStep.Simulation
  ConfigFlags.Sim.ISFRun = True

  #Provide input
  from AthenaConfiguration.TestDefaults import defaultTestFiles
  inputDir = defaultTestFiles.d
  ConfigFlags.Input.Files = defaultTestFiles.EVNT

  from AthenaConfiguration.Enums import LHCPeriod
  ConfigFlags.GeoModel.Run = LHCPeriod.Run2

  detectors =['BCM', 'Pixel', 'SCT', 'TRT', 'LAr']
  # Setup detector flags
  from AthenaConfiguration.DetectorConfigFlags import setupDetectorsFromList
  setupDetectorsFromList(ConfigFlags, detectors, toggle_geometry=True)

  #ConfigFlags.GeoModel.AtlasVersion = "tb_Tile2000_2003"
  #ConfigFlags.GeoModel.AtlasVersion = "ctbh8"
  ConfigFlags.GeoModel.AtlasVersion = 'ATLAS-R2-2015-03-01-00'
  # Finalize
  ConfigFlags.lock()


  ## Initialize a new component accumulator
  cfg = MainServicesCfg(ConfigFlags)


  from G4AtlasTools.G4AtlasToolsConfigNew import SensitiveDetectorMasterToolCfg
  acc  = SensitiveDetectorMasterToolCfg(ConfigFlags)
  tool = cfg.popToolsAndMerge(acc)
  cfg.setPrivateTools(tool)

  cfg.printConfig(withDetails=True, summariseProps = True)
  ConfigFlags.dump()


  #cfg not being used so complains ...fine now!
  f=open("test.pkl","wb")
  cfg.store(f) #sets wasmerged = true
  f.close()


  print ("-----------------finished----------------------")
