# Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
#
# File: CoolLumiUtilities/python/OnlineLumiCalibrationCondAlgConfig.py
# Created: May 2019, sss
# Purpose: Configure OnlineLumiCalibrationCondAlg.
#


from __future__ import print_function


from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from IOVDbSvc.IOVDbSvcConfig import addFolders
from AthenaCommon.Logging import logging


def OnlineLumiCalibrationCondAlgCfg (configFlags):
    log = logging.getLogger ('OnlineLumiCalibrationCondAlgCfg')
    name = 'OnlineLumiCalibrationCondAlg'
    result = ComponentAccumulator()

    if configFlags.Input.isMC:
        return result

    # For both runs 1 and 2
    folder = '/TDAQ/OLC/CALIBRATIONS'
    result.merge (addFolders (configFlags, folder, 'TDAQ',
                              className='CondAttrListCollection'))

    from CoolLumiUtilities.CoolLumiUtilitiesConf import \
         OnlineLumiCalibrationCondAlg
    alg = OnlineLumiCalibrationCondAlg (name,
                                        CalibrationFolderInputKey = folder,
                                        LumiCalibOutputKey = 'OnlineLumiCalibrationCondData')

    result.addCondAlgo (alg)
    return result


if __name__ == "__main__":
    from AthenaCommon.Configurable import Configurable
    Configurable.configurableRun3Behavior=1
    from AthenaConfiguration.AllConfigFlags import ConfigFlags
    from AthenaConfiguration.TestDefaults import defaultTestFiles

    print ('--- data')
    flags1 = ConfigFlags.clone()
    flags1.Input.Files = defaultTestFiles.RAW
    flags1.lock()
    acc1 = OnlineLumiCalibrationCondAlgCfg (flags1)
    acc1.printCondAlgs(summariseProps=True)
    print ('IOVDbSvc:', acc1.getService('IOVDbSvc').Folders)
    acc1.wasMerged()

    print ('--- mc')
    flags2 = ConfigFlags.clone()
    flags2.Input.Files = defaultTestFiles.ESD
    flags2.lock()
    acc2 = OnlineLumiCalibrationCondAlgCfg (flags2)
    acc2.printCondAlgs(summariseProps=True)
    print ('IOVDbSvc:', acc2.getServices())
    acc2.wasMerged()
