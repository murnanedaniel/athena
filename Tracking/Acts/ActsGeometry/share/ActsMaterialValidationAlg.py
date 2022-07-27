"""
This job options file will run an example extrapolation using the
Acts tracking geometry and the Acts extrapolation toolchain.

Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
"""

# start from scratch with component accumulator

from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory

from ActsGeometry.ActsGeometryConfig import ActsExtrapolationToolCfg
from ActsGeometry.ActsGeometryConfig import ActsAlignmentCondAlgCfg

def ActsExtrapolationAlgCfg(configFlags, name = "ActsExtrapolationAlg", **kwargs):
  result = ComponentAccumulator()

  if "ExtrapolationTool" not in kwargs:
    extrapTool = ActsExtrapolationToolCfg(configFlags)
    kwargs["ExtrapolationTool"] = extrapTool.getPrimary()
    result.merge(extrapTool)

  ActsExtrapolationAlg = CompFactory.ActsExtrapolationAlg
  alg = ActsExtrapolationAlg(name, **kwargs)
  result.addEventAlgo(alg)

  return result

if "__main__" == __name__:
  from AthenaCommon.Logging import log
  from AthenaCommon.Constants import INFO
  from AthenaConfiguration.AllConfigFlags import ConfigFlags
  from AthenaConfiguration.MainServicesConfig import MainServicesCfg
  from ActsGeometry.ActsGeometryConfig import ActsMaterialTrackWriterSvcCfg

  ## Just enable ID for the moment.
  ConfigFlags.Input.isMC             = True
  ConfigFlags.GeoModel.AtlasVersion  = "ATLAS-R2-2016-01-00-01"
  ConfigFlags.IOVDb.GlobalTag        = "OFLCOND-SIM-00-00-00"
  ConfigFlags.Detector.GeometryBpipe = True
  ConfigFlags.Detector.GeometryID    = True
  ConfigFlags.Detector.GeometryPixel = True
  ConfigFlags.Detector.GeometrySCT   = True
  ConfigFlags.Detector.GeometryCalo  = True
  ConfigFlags.Detector.GeometryMuon  = False
  ConfigFlags.Detector.GeometryTRT   = True
  ConfigFlags.Acts.TrackingGeometry.MaterialSource = "material-maps.json"
  # ConfigFlags.Acts.TrackingGeometry.MaterialSource = "/eos/project-a/acts/public/MaterialMaps/ATLAS/material-maps.json"

  ConfigFlags.Concurrency.NumThreads = 10
  ConfigFlags.Concurrency.NumConcurrentEvents = 10

  ConfigFlags.lock()
  ConfigFlags.dump()

  cfg = MainServicesCfg(ConfigFlags)

  from BeamPipeGeoModel.BeamPipeGMConfig import BeamPipeGeometryCfg
  cfg.merge(BeamPipeGeometryCfg(ConfigFlags))

  alignCondAlgCfg = ActsAlignmentCondAlgCfg(ConfigFlags)

  cfg.merge(alignCondAlgCfg)

  cfg.merge(ActsMaterialTrackWriterSvcCfg(ConfigFlags,
                                          "ActsMaterialTrackWriterSvc",
                                          "MaterialTracks_mapped.root"))

  print('DEF WRITER : ')
  extrapol = ActsExtrapolationToolCfg(ConfigFlags,
                                      InteractionMultiScatering = True,
                                      InteractionEloss = True,
                                      InteractionRecord = True)
  cfg.merge(extrapol)
  
  alg = ActsExtrapolationAlgCfg(ConfigFlags,
                                OutputLevel=INFO,
                                NParticlesPerEvent=int(1e4),
                                EtaRange=[-2.5, 2.5],
                                PtRange=[20, 100],
                                WriteMaterialTracks = True,
                                ExtrapolationTool=extrapol.getPrimary())

  cfg.merge(alg)

  tgSvc = cfg.getService("ActsTrackingGeometrySvc")

  # Service will have removed TRT and Calo
  # We want them enabled for testing
  tgSvc.BuildSubDetectors += [
    "TRT",
    "Calo"
  ]
  # needed to construct the calo geometry in ACTS
  tgSvc.CaloVolumeBuilder = CompFactory.ActsCaloTrackingVolumeBuilder()


  cfg.printConfig()

  log.info("CONFIG DONE")

  cfg.run(100)

