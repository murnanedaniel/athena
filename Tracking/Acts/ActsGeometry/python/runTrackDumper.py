#!/usr/bin/env python

from AthenaConfiguration.ComponentFactory import CompFactory
from AthenaConfiguration.Enums import LHCPeriod
from AthenaConfiguration.AllConfigFlags import initConfigFlags
flags = initConfigFlags()

flags.Input.Files = [
    "/eos/atlas/atlaslocalgroupdisk/dq2/rucio/mc15_14TeV/99/7f/AOD.30640765._000024.pool.root.1"
]

flags.lock()

from AthenaConfiguration.MainServicesConfig import MainServicesCfg
acc = MainServicesCfg(flags)
from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
acc.merge(PoolReadCfg(flags))

#  from InDetPhysValMonitoring.InDetPhysValMonitoringConfig import InDetPhysValMonitoringCfg
#  acc.merge(InDetPhysValMonitoringCfg(flags))


TrackDumperAlg = CompFactory.TrackDumperAlg
alg = TrackDumperAlg("TrackDumperAlg")
acc.addEventAlgo(alg)

acc.printConfig(withDetails=True)


#  # Execute and finish
sc = acc.run(maxEvents=10)

#  # Success should be 0
import sys
sys.exit(not sc.isSuccess())
