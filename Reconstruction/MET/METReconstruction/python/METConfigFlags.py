# Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaborationfrom AthenaConfiguration.AthConfigFlags import AthConfigFlagsdef createMETConfigFlags():    metConfigFlags=AthConfigFlags()    metConfigFlags.addFlag("MET.UseTracks",True)     metConfigFlags.addFlag("MET.DoPFlow",True)     return metConfigFlags