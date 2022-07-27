# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

##############################################################

from AthenaConfiguration.AllConfigFlags import ConfigFlags
ConfigFlags.Scheduler.ShowDataDeps = True
ConfigFlags.Scheduler.ShowDataFlow = True
ConfigFlags.Scheduler.ShowControlFlow = True

ConfigFlags.lock()
ConfigFlags.dump()

##############################################################  

from AthenaConfiguration.MainServicesConfig import MainServicesCfg
from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
acc = MainServicesCfg( ConfigFlags )
acc.merge( PoolReadCfg( ConfigFlags ) )
acc.getService("MessageSvc").Format = "% F%80W%S%7W%R%T %0W%M"
acc.getService("MessageSvc").defaultLimit = 2147483647

from FlavourTaggingTests.FlavourTaggingConfiguration_HistoService import getHistoServiceCfg
acc.addService( getHistoServiceCfg( ConfigFlags.Output.HISTFileName ) )

from FlavourTaggingTests.FlavourTaggingConfiguration_PhysicsVariablePlots import getPhysicsVariablePlotsCfg
acc.merge( getPhysicsVariablePlotsCfg() )

acc.printConfig(withDetails = True, summariseProps = True)
acc.store( open("FTAG_PhysicsVariablePlots.pkl","wb") )
sc = acc.run()

import sys
sys.exit(not sc.isSuccess())

##############################################################  

