#!/bin/env python
# Set the Athena configuration flags
from AthenaConfiguration.AllConfigFlags import ConfigFlags
ConfigFlags.fillFromArgs()
ConfigFlags.lock()


from AthenaConfiguration.MainServicesConfig import MainServicesCfg
from AthenaConfiguration.ComponentFactory import CompFactory
from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
cfg = MainServicesCfg(ConfigFlags)
cfg.merge(PoolReadCfg(ConfigFlags))

from  AthenaMonitoring.TriggerInterface import TrigDecisionToolCfg

cfg.merge(TrigDecisionToolCfg(ConfigFlags))
checker = CompFactory.TrigDecisionChecker()
checker.WriteEventDecision = True
checker.WriteOutFilename = "trigger.counts.log"

cfg.addEventAlgo(checker)
cfg.run()
