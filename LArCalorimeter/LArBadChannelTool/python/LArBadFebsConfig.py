# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory
from AthenaConfiguration.Enums import LHCPeriod
from IOVDbSvc.IOVDbSvcConfig import addFoldersSplitOnline

LArBadFebCondAlg = CompFactory.LArBadFebCondAlg


def LArKnownBadFebCfg(configFlags, tag=None):
    result=ComponentAccumulator()

    if configFlags.GeoModel.Run is LHCPeriod.Run1:
        rekey=""
    else:
        rekey="/LAR/BadChannels/KnownBADFEBs"
        result.merge(addFoldersSplitOnline(configFlags,"LAR","/LAR/BadChannels/KnownBADFEBs",
                                           f"/LAR/BadChannelsOfl/KnownBADFEBs<key>{rekey}</key>",tag=tag,
                                           className="AthenaAttributeList"))

    result.addCondAlgo(LArBadFebCondAlg("LArKnownBadFebAlg",ReadKey=rekey,WriteKey="LArKnownBadFEBs"))
    return result

def LArKnownMNBFebCfg(configFlags, tag=None):
    result=ComponentAccumulator()

    if configFlags.GeoModel.Run is LHCPeriod.Run1:
        rekey=""
    else:
        rekey="/LAR/BadChannels/KnownMNBFEBs"
        result.merge(addFoldersSplitOnline(configFlags,"LAR","/LAR/BadChannels/KnownMNBFEBs",
                                           f"/LAR/BadChannelsOfl/KnownMNBFEBs<key>{rekey}</key>",tag=tag,
                                           className="AthenaAttributeList"))

    result.addCondAlgo(LArBadFebCondAlg("LArKnownMNBFebAlg",ReadKey=rekey,WriteKey="LArKnownMNBFEBs"))
    return result



if __name__=="__main__":

    from AthenaConfiguration.AllConfigFlags import ConfigFlags
    from AthenaCommon.Logging import log
    from AthenaCommon.Constants import DEBUG
    log.setLevel(DEBUG)

    ConfigFlags.Input.isMC = False
    ConfigFlags.Input.Files = ["/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/Tier0ChainTests/data17_13TeV.00330470.physics_Main.daq.RAW._lb0310._SFO-1._0001.data"]
    ConfigFlags.lock()

    cfg=ComponentAccumulator()

    cfg.merge(LArKnownBadFebCfg(ConfigFlags))
    cfg.merge(LArKnownMNBFebCfg(ConfigFlags))
    f=open("LArBadFebCondAlgos.pkl","wb")
    cfg.store(f)
    f.close()
