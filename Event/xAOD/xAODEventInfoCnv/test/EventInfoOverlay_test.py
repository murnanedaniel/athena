#!/usr/bin/env python
"""Run tests for EventInfo overlay

Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
"""
import sys

from AthenaConfiguration.AllConfigFlags import ConfigFlags
from AthenaConfiguration.MainServicesConfig import MainServicesCfg
from AthenaConfiguration.TestDefaults import defaultTestFiles
from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
from OverlayConfiguration.OverlayTestHelpers import \
    CommonTestArgumentParser, postprocessAndLockFlags, printAndRun
from xAODEventInfoCnv.xAODEventInfoCnvConfig import EventInfoOverlayCfg

# Argument parsing
parser = CommonTestArgumentParser("EventInfoOverlay_test.py")
args = parser.parse_args()

# Configure
ConfigFlags.Input.Files = defaultTestFiles.RDO_BKG_RUN2
ConfigFlags.Input.SecondaryFiles = ["/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/OverlayTests/Special/TestCase_xAODEventInfo.root"]
ConfigFlags.IOVDb.GlobalTag = "OFLCOND-MC16-SDR-16"
ConfigFlags.Overlay.DataOverlay = False
ConfigFlags.Output.RDOFileName = "myRDO.pool.root"
ConfigFlags.Output.RDO_SGNLFileName = "myRDO_SGNL.pool.root"

postprocessAndLockFlags(ConfigFlags, args)

# Function tests
accAlg = EventInfoOverlayCfg(ConfigFlags)
# reset to prevent errors on deletion
accAlg.__init__()

# Construct our accumulator to run
acc = MainServicesCfg(ConfigFlags)
acc.merge(PoolReadCfg(ConfigFlags))

# Print and run
sys.exit(printAndRun(acc, ConfigFlags, args))
