# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

### This module contains functions which may need to peek at the input file metadata

from AthenaKernel.EventIdOverrideConfig import getMinMaxRunNumbers
## Get the logger
from AthenaCommon.Logging import logging
simMDlog = logging.getLogger('Sim_Metadata')

def fillAtlasMetadata(ConfigFlags, dbFiller):

    #add all ConfigFlags to the metadata
    #todo - only add certain ones?
    #in future this should be a ConfigFlags method...?
    for flag in sorted(ConfigFlags._flagdict): #only sim
        if "Sim" in flag:
            key = flag.split(".")[-1] #use final part of flag as the key
            value = eval("ConfigFlags."+flag)
            if not isinstance(value, str):
                value = str(value)
            dbFiller.addSimParam(key, value)
            simMDlog.info('SimulationMetaData: setting "%s" to be %s', key, value)

    dbFiller.addSimParam('G4Version', ConfigFlags.Sim.G4Version)
    dbFiller.addSimParam('RunType', 'atlas')
    dbFiller.addSimParam('beamType', ConfigFlags.Beam.Type.value)
    dbFiller.addSimParam('SimLayout', ConfigFlags.GeoModel.AtlasVersion)
    dbFiller.addSimParam('MagneticField', 'AtlasFieldSvc') # TODO hard-coded for now for consistency with old-style configuration.

    #---------
    ## Simulated detector flags: add each enabled detector to the simulatedDetectors list
    ConfigFlags._loadDynaFlags("Detector") # Ensure that Detector flags have been loaded
    from AthenaConfiguration.DetectorConfigFlags import allDetectors
    simDets = []
    for det in allDetectors:
        if det in ['Bpipe', 'Cavern']: # skip regions without sensitive detectors
            continue
        attrname = f'Detector.Geometry{det}'
        if ConfigFlags.hasFlag(attrname):
            testValue = ConfigFlags(attrname)
            if testValue:
                simDets.append(det)
        else:
            simMDlog.info("No flag called '%s' found in ConfigFlags", attrname)
    simDets.append('Truth') # FIXME Currently there is no way to switch off processing Truth containers

    simMDlog.info("Setting 'SimulatedDetectors' = %r", simDets)
    dbFiller.addSimParam('SimulatedDetectors', repr(simDets))

    ## Hard-coded simulation hit file magic number (for major changes)
    dbFiller.addSimParam('hitFileMagicNumber', '0') ##FIXME Remove this?

    if ConfigFlags.Sim.ISFRun:
        dbFiller.addSimParam('Simulator', ConfigFlags.Sim.ISF.Simulator.value)
        dbFiller.addSimParam('SimulationFlavour', ConfigFlags.Sim.ISF.Simulator.value.replace('MT', '')) # used by egamma
    else:
        # TODO hard-code for now, but set flag properly later
        dbFiller.addSimParam('Simulator', 'AtlasG4')
        dbFiller.addSimParam('SimulationFlavour', 'AtlasG4')


def writeSimulationParametersMetadata(ConfigFlags):
    from IOVDbMetaDataTools import ParameterDbFiller
    dbFiller = ParameterDbFiller.ParameterDbFiller()
    myRunNumber, myEndRunNumber = getMinMaxRunNumbers(ConfigFlags)
    simMDlog.debug('ParameterDbFiller BeginRun = %s', str(myRunNumber) )
    dbFiller.setBeginRun(myRunNumber)
    simMDlog.debug('ParameterDbFiller EndRun   = %s', str(myEndRunNumber) )
    dbFiller.setEndRun(myEndRunNumber)

    fillAtlasMetadata(ConfigFlags, dbFiller)

    #-------------------------------------------------
    # Make the MetaData Db
    #-------------------------------------------------
    dbFiller.genSimDb()

    from IOVDbSvc.IOVDbSvcConfig import IOVDbSvcCfg
    cfg = IOVDbSvcCfg(ConfigFlags)
    folder = "/Simulation/Parameters"
    dbConnection = "sqlite://;schema=SimParams.db;dbname=SIMPARAM"

    cfg.getService("IOVDbSvc").Folders += [ folder + "<dbConnection>" + dbConnection + "</dbConnection>" ]
    cfg.getService("IOVDbSvc").FoldersToMetaData += [ folder ]
    #cfg.getService("IOVSvc").partialPreLoadData = True #FIXME IOVSvc missing??
    return cfg
