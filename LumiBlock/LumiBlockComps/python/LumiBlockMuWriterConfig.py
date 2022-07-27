# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
#
# File: LumiBlockComps/python/LumiBlockMuWriterConfig.py
# Created: May 2020, sss
# Purpose: Configure LumiBlockMuWriter.
#


from AthenaConfiguration.ComponentFactory import CompFactory
from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.Enums import BeamType


def LumiBlockMuWriterCfg (configFlags, name = 'LumiBlockMuWriter', seqName="AthAlgSeq"):
    result = ComponentAccumulator(seqName)

    if configFlags.Beam.Type is BeamType.Cosmics or configFlags.Input.isMC:
        condkey = ''
    else:
        from LumiBlockComps.LuminosityCondAlgConfig import LuminosityCondAlgCfg
        result.merge (LuminosityCondAlgCfg (configFlags))
        condkey = result.getCondAlgo ('LuminosityCondAlg').LuminosityOutputKey

    LumiBlockMuWriter = CompFactory.LumiBlockMuWriter # LumiBlockComps
    alg = LumiBlockMuWriter (name, LumiDataKey = condkey)
    #In the HLT we want to run LumiBlockMuWriter as a normal EventAlgo, but in a pre-event sequence (HLTBeginSeq)
    if  configFlags.Trigger.doHLT:
         result.addEventAlgo(alg)
    #For offline and particularly serial athena, add LumiBlockMuWriter to AthCondSeq to ensure it runs first (ATR-24721)
    else:
         result.addCondAlgo(alg)
    return result


if __name__ == "__main__":
    from AthenaConfiguration.AllConfigFlags import ConfigFlags
    from AthenaConfiguration.TestDefaults import defaultTestFiles
    ConfigFlags.Input.Files = []

    print ('--- collisions')
    flags1 = ConfigFlags.clone()
    flags1.Input.Files = defaultTestFiles.RAW
    flags1.lock()
    acc1 = LumiBlockMuWriterCfg (flags1)
    acc1.printCondAlgs (summariseProps=True)
    acc1.wasMerged()

    print ('--- cosmics')
    flags2 = ConfigFlags.clone()
    flags2.Beam.Type = BeamType.Cosmics
    flags2.lock()
    acc2 = LumiBlockMuWriterCfg (flags2)
    acc2.printCondAlgs (summariseProps=True)
    acc2.wasMerged()
