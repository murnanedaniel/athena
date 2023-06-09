#! /usr/bin/env python

# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration

"""
Run ISF simulation on generated events and produce a HITS file.
"""

import os.path
import sys
import time

import logging

# Setup core logging here
from PyJobTransforms.trfLogger import msg
msg.info('logging set in %s' % sys.argv[0])

from PyJobTransforms.transform import transform
from PyJobTransforms.trfExe import athenaExecutor
from PyJobTransforms.trfArgs import addAthenaArguments, addDetectorArguments
from PyJobTransforms.trfDecorators import stdTrfExceptionHandler, sigUsrStackTrace
from SimuJobTransforms.simTrfArgs import addForwardDetTrfArgs, addCosmicsTrfArgs, addCommonSimTrfArgs, addCommonSimDigTrfArgs, addTrackRecordArgs, addSim_tfArgs


import PyJobTransforms.trfArgClasses as trfArgClasses

# Prodsys hack...
ListOfDefaultPositionalKeys=['--AFPOn', '--ALFAOn', '--CosmicFilterVolume', '--CosmicFilterVolume2', '--CosmicPtSlice', '--DBRelease', '--DataRunNumber', '--FwdRegionOn', '--LucidOn', '--ZDCOn', '--amiConfig', '--amiMetadataTag', '--asetup', '--athena', '--athenaopts', '--beamType', '--checkEventCount', '--command', '--conditionsTag', '--enableLooperKiller', '--env', '--eventAcceptanceEfficiency', '--execOnly', '--firstEvent', '--geometryVersion', '--ignoreErrors', '--ignoreFiles', '--ignorePatterns', '--imf', '--inputEVNTFile', '--inputEVNT_CAVERNFile', '--inputEVNT_COSMICSFile', '--jobNumber', '--maxEvents', '--outputEVNT_CAVERNTRFile', '--outputEVNT_COSMICSTRFile', '--outputHITSFile', '--physicsList', '--postExec', '--postInclude', '--preExec', '--preInclude', '--randomSeed', '--reportName', '--reportType', '--runNumber', '--showGraph', '--showPath', '--showSteps', '--simulator', '--skipEvents', '--skipFileValidation', '--skipInputFileValidation', '--skipOutputFileValidation', '--tcmalloc', '--useISF']

@stdTrfExceptionHandler
@sigUsrStackTrace
def main():

    msg.info('This is %s' % sys.argv[0])

    trf = getTransform()
    trf.parseCmdLineArgs(sys.argv[1:])
    trf.execute()
    trf.generateReport()

    msg.info("%s stopped at %s, trf exit code %d" % (sys.argv[0], time.asctime(), trf.exitCode))
    sys.exit(trf.exitCode)

def getTransform():
    executorSet = set()
    from SimuJobTransforms.SimTransformUtils import addSimulationSubstep, addSimulationArguments
    addSimulationSubstep(executorSet)
    trf = transform(executor = executorSet, description = 'ATLAS Simulation transform. Inputs must be EVNT else a single particle Generator job options must be specified. Outputs must be HITS or TrackRecords.')
    addAthenaArguments(trf.parser)
    addDetectorArguments(trf.parser)
    addSimulationArguments(trf.parser)
    return trf

## FIXME - not sure what the equivalent of the method below is in the new framework?

##     def doPreRunActions(self):
##         JobTransform.doPreRunActions(self)
##         if hasattr(self,'_maxEventsStrategy'):
##             self._maxEventsStrategy = 'ABORT'
##         else:
##             print "WARNING EVGENtoHITJobTransform has no attribute \'_maxEventsStrategy\'."


if __name__ == '__main__':
    main()
