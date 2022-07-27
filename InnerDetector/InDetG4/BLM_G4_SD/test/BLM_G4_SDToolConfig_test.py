#!/usr/bin/env python
"""Run tests on BLM_G4_SD configuration

Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
"""


if __name__ == '__main__':

  # Set up logging and config behaviour
  from AthenaCommon.Logging import log
  from AthenaCommon.Constants import DEBUG
  log.setLevel(DEBUG)


  #import config flags
  from AthenaConfiguration.AllConfigFlags import ConfigFlags
  ConfigFlags.Sim.ISFRun = True

  #Provide input
  from AthenaConfiguration.TestDefaults import defaultTestFiles
  inputDir = defaultTestFiles.d
  ConfigFlags.Input.Files = defaultTestFiles.EVNT

  # Finalize
  ConfigFlags.lock()

  # Setup the tool
  from BLM_G4_SD.BLM_G4_SDToolConfig import BLMSensorSDCfg
  cfg = BLMSensorSDCfg(ConfigFlags)
  cfg.printConfig(withDetails=True, summariseProps = True)
  ConfigFlags.dump()

  f=open("test.pkl","wb")
  cfg.store(f)
  f.close()

  print(cfg._privateTools)
  print("-----------------finished----------------------")
