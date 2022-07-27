# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory
from AthenaConfiguration.Enums import ProductionStep
from AthenaCommon import Logging

def GeoModelCfg(configFlags):
    version=configFlags.GeoModel.AtlasVersion

    from AthenaCommon.AppMgr import release_metadata
    rel_metadata = release_metadata()
    relversion = rel_metadata['release'].split('.')
    if len(relversion) < 3:
        relversion = rel_metadata['base release'].split('.')

    result=ComponentAccumulator()

    #Get DetDescrCnvSvc (for identifier dictionaries (identifier helpers)
    from DetDescrCnvSvc.DetDescrCnvSvcConfig import DetDescrCnvSvcCfg
    result.merge(DetDescrCnvSvcCfg(configFlags))

    #TagInfoMgr used by GeoModelSvc but no ServiceHandle. Relies on string-name
    from EventInfoMgt.TagInfoMgrConfig import TagInfoMgrCfg
    result.merge(TagInfoMgrCfg(configFlags))

    gms=CompFactory.GeoModelSvc(AtlasVersion=version,
                                SupportedGeometry = int(relversion[0]))
    if configFlags.Common.ProductionStep == ProductionStep.Simulation:
        ## Protects GeoModelSvc in the simulation from the AlignCallbacks
        gms.AlignCallbacks = False
    result.addService(gms, primary=True, create=True)

    return result


if __name__ == "__main__":
    from AthenaConfiguration.AllConfigFlags import ConfigFlags

    ConfigFlags.Input.Files = []

    acc = GeoModelCfg( ConfigFlags )
    acc.store( open( "test.pkl", "wb" ) )
    Logging.log.info("All OK")
