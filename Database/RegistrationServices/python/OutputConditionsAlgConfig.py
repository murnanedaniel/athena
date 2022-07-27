# Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
from AthenaConfiguration.ComponentFactory import CompFactory

def OutputConditionsAlgCfg(flags, name="OutputConditionsAlg",outputFile='condobjs.root', **kwargs):

    from AthenaPoolCnvSvc.PoolWriteConfig import PoolWriteCfg

    result=PoolWriteCfg(flags)

    kwargs.setdefault("WriteIOV",True)
    oca=CompFactory.OutputConditionsAlg(name,**kwargs)

    
    # create outputStream tool with given filename and pass to myOCA
    condstream=CompFactory.AthenaOutputStreamTool(name+"Tool")

    #To be fixed: OutputConditionsAlgo works with a string-name of a public tool
    oca.StreamName=name+"Tool"
    condstream.OutputFile=outputFile
    condstream.PoolContainerPrefix="ConditionsContainer"
    condstream.TopLevelContainerName = "<type>"
    condstream.SubLevelBranchName = "<key>"
    result.addPublicTool(condstream)

    result.addEventAlgo(oca)
    
    return result
