#!/usr/bin/env python
"""Run ACTS geometry construction for ITk

Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
"""
from AthenaConfiguration.AllConfigFlags import initConfigFlags
flags = initConfigFlags()

from AthenaConfiguration.Enums import ProductionStep
flags.Common.ProductionStep = ProductionStep.Simulation
flags.GeoModel.AtlasVersion = "ATLAS-P2-RUN4-01-00-00"
flags.IOVDb.GlobalTag = "OFLCOND-MC16-SDR-14"
flags.GeoModel.Align.Dynamic = False
flags.Input.Files = ['/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/SimCoreTests/valid1.410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.evgen.EVNT.e4993.EVNT.08166201._000012.pool.root.1']

flags.Detector.GeometryITkPixel = True
flags.Detector.GeometryITkStrip = True
flags.Detector.GeometryBpipe = True
flags.Detector.GeometryCalo = False

flags.Concurrency.NumThreads = 64
flags.Concurrency.NumConcurrentEvents = 64

flags.Exec.MaxEvents = 200

flags.lock()
flags.dump()

from AthenaConfiguration.MainServicesConfig import MainServicesCfg
acc = MainServicesCfg( flags )

from ActsGeometry.ActsGeometryConfig import ActsExtrapolationAlgCfg, ActsTrackingGeometrySvcCfg

from AthenaCommon.Constants import INFO
tgSvc = ActsTrackingGeometrySvcCfg(flags,
                                   OutputLevel=INFO,
                                   ObjDebugOutput=True)
acc.merge(tgSvc)

alg = ActsExtrapolationAlgCfg(flags,
                              OutputLevel=INFO,
                              NParticlesPerEvent = int(100),
                              WritePropStep = True,
                              EtaRange = [-5, 5],
                              PtRange = [20, 100])

acc.merge(alg)
acc.printConfig()
acc.run()
