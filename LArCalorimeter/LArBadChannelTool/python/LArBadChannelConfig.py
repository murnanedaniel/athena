# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory
from LArCabling.LArCablingConfig import LArOnOffIdMappingCfg
from IOVDbSvc.IOVDbSvcConfig import addFolders, addFoldersSplitOnline

def LArBadChannelCfg(configFlags, tag=None):

    

    result=LArOnOffIdMappingCfg(configFlags)
    rekey="/LAR/BadChannels/BadChannels"

    if configFlags.Overlay.DataOverlay:
        # TODO: move this in a better location
        result.merge(addFolders(configFlags, "/LAR/BadChannels/BadChannels", "LAR_OFL", className="CondAttrListCollection", tag="LARBadChannelsBadChannels-HECAQ3Missing", db="OFLP200"))
    elif configFlags.Input.isMC:
        result.merge(addFolders(configFlags,"/LAR/BadChannels/BadChannels","LAR_OFL",tag=tag,
                                className="CondAttrListCollection"))
    else:
        result.merge(addFoldersSplitOnline(configFlags,"LAR","/LAR/BadChannels/BadChannels",
                                        f"/LAR/BadChannelsOfl/BadChannels<key>{rekey}</key>",tag=tag,
                                        className="CondAttrListCollection"))  
    theLArBadChannelCondAlgo=CompFactory.LArBadChannelCondAlg(ReadKey=rekey)
    result.addCondAlgo(theLArBadChannelCondAlgo)
    return result


def LArBadFebCfg(configFlags, tag=None):
    result=ComponentAccumulator()
    rekey="/LAR/BadChannels/MissingFEBs"

    if configFlags.Overlay.DataOverlay:
        # TODO: move this in a better location
        result.merge(addFolders(configFlags, "/LAR/BadChannels/MissingFEBs", "LAR_OFL", className="AthenaAttributeList", tag="LArBadChannelsMissingFEBs-IOVDEP-04", db="OFLP200"))
    elif configFlags.Input.isMC:
        result.merge(addFolders(configFlags,"/LAR/BadChannels/MissingFEBs","LAR_OFL",tag=tag,
                                className="AthenaAttributeList"))
    else:
        result.merge(addFoldersSplitOnline(configFlags,"LAR","/LAR/BadChannels/MissingFEBs",
                                           f"/LAR/BadChannelsOfl/MissingFEBs<key>{rekey}</key>",tag=tag,
                                           className="AthenaAttributeList"))  
    result.addCondAlgo(CompFactory.LArBadFebCondAlg(ReadKey=rekey))
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
    
    cfg.merge(LArBadChannelCfg(ConfigFlags))
    cfg.merge(LArBadFebCfg(ConfigFlags))
    f=open("LArBCCondAlgos.pkl","wb")
    cfg.store(f)
    f.close()
    
