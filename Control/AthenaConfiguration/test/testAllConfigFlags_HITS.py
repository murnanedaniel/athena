#!/usr/bin/env python
# Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
from AthenaConfiguration.AllConfigFlags import initConfigFlags
from AthenaConfiguration.TestDefaults import defaultTestFiles
flags = initConfigFlags()
flags.Input.Files = defaultTestFiles.HITS_RUN2

# Load Detector flags
flags._loadDynaFlags("Detector")

# Load Sim flags
flags._loadDynaFlags("Sim")

# Load GeoModel flags
flags._loadDynaFlags("GeoModel")

# Init and print
flags.initAll()
flags.lock()
flags.dump()
