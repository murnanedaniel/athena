#!/usr/bin/env python
"""Run tests for overlay metadata

Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
"""
import sys

from AthenaConfiguration.AllConfigFlags import ConfigFlags
from AthenaConfiguration.MainServicesConfig import MainServicesCfg
from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
from xAODEventInfoCnv.xAODEventInfoCnvConfig import EventInfoOverlayCfg

from OverlayConfiguration.OverlayMetadata import overlayMetadataCheck, overlayMetadataWrite
from OverlayConfiguration.OverlayTestHelpers import defaultTestFlags, postprocessAndLockFlags, printAndRun, CommonTestArgumentParser

# Argument parsing
parser = CommonTestArgumentParser("OverlayMetadataConfig_test.py")
args = parser.parse_args()

# Configure
defaultTestFlags(ConfigFlags, args)
overlayMetadataCheck(ConfigFlags)
postprocessAndLockFlags(ConfigFlags, args)
ConfigFlags.initAll()
ConfigFlags.dump()
# Construct our accumulator to run
acc = MainServicesCfg(ConfigFlags)
acc.merge(PoolReadCfg(ConfigFlags))
acc.merge(overlayMetadataWrite(ConfigFlags))

# Add event info overlay for minimal output
acc.merge(EventInfoOverlayCfg(ConfigFlags))

# Print and run
sys.exit(printAndRun(acc, ConfigFlags, args))
