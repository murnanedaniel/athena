#!/usr/bin/env python

# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration

## FTK Simulation Transform
#  Specialist version to do sim x 4 subregions and merge in one job
# @version $Id: TrigFTKSM4Un_tf.py 643664 2015-02-02 19:48:36Z jahreda $ 

import sys
import time

import logging

# Setup core logging here
from PyJobTransforms.trfLogger import msg
msg.info('logging set in %s' % sys.argv[0])

from PyJobTransforms.trfExitCodes import trfExit
from PyJobTransforms.transform import transform
from PyJobTransforms.trfExe import athenaExecutor
from PyJobTransforms.trfArgs import addAthenaArguments
from PyJobTransforms.trfDecorators import stdTrfExceptionHandler, sigUsrStackTrace

import PyJobTransforms.trfExceptions as trfExceptions
import PyJobTransforms.trfArgClasses as trfArgClasses

subregions = 4

from TrigFTKSim.FTKSimOptions import *

ListOfDefaultPositionalKeys=['--AMI', '--CachePath', '--CachedBank', '--DBBankLevel', '--DoRoadFile', '--FTKDoGrid', '--FTKForceAllInput', '--FTKSetupTag', '--FTKUnmergedInputPath', '--HWNDiff', '--HitWarrior', '--IBLMode', '--MakeCache', '--NBanks', '--NSubRegions', '--PixelClusteringMode', '--RoadFilesDir', '--SSFAllowExtraMiss', '--SSFMultiConnection', '--SSFNConnections', '--SSFTRDefn', '--SSFTRMaxEta', '--SSFTRMinEta', '--SaveRoads', '--SctClustering', '--SecondStageFit', '--SetAMSize', '--TRACKFITTER_MODE', '--TSPMinCoverage', '--TSPSimulationLevel', '--UseTSPBank', '--badmap_path', '--badmap_path_for_hit', '--bankpatterns', '--bankregion', '--execOnly', '--fit711constantspath', '--fitconstantspath', '--ignoreErrors', '--inputNTUP_FTKIPFile', '--inputNTUP_FTKTMP_0File', '--inputNTUP_FTKTMP_1File', '--inputNTUP_FTKTMP_2File', '--inputNTUP_FTKTMP_3File', '--inputTXT_FTKIPFile', '--loadHWConf_path', '--maxEvents', '--omitFileValidation', '--outputNTUP_FTKTMPFile', '--outputNTUP_FTKTMP_0File', '--outputNTUP_FTKTMP_1File', '--outputNTUP_FTKTMP_2File', '--outputNTUP_FTKTMP_3File', '--patternbank0path', '--patternbank1path', '--patternbank2path', '--patternbank3path', '--pmap_path', '--pmapcomplete_path', '--pmapunused_path', '--reportName', '--rmap_path', '--sectorpath', '--showGraph', '--showPath', '--showSteps', '--ssmap_path', '--ssmaptsp_path', '--ssmapunused_path', '--uploadtoami', '--validation', '--SaveTruthTree', '--postExec', '--postInclude', '--preExec', '--preInclude', '--firstEvent', '--Save1stStageTrks']


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


## Get the base transform with all arguments added
def getTransform():

    executorSet = set()
    for subregion in range(subregions):
        executorSet.add(athenaExecutor(name = 'FTKFullSimulationBank{0}'.format(subregion), 
                                       skeletonFile = 'TrigFTKSim/skeleton.FTKStandaloneSim.py',
                                       inData = ['NTUP_FTKIP','TXT_FTKIP'], outData = ['NTUP_FTKTMP_{0}'.format(subregion)],
                                       extraRunargs = {'banksubregion': [subregion]},
                                       # Need to ensure that the correct subregion is used
                                       runtimeRunargs = {'patternbankpath': 'runArgs.patternbank{0}path'.format(subregion),
                                                         'outputNTUP_FTKTMPFile': 'runArgs.outputNTUP_FTKTMP_{0}File'.format(subregion)}))
    executorSet.add(athenaExecutor(name = 'FTKSimulationMerge', 
                                   skeletonFile = 'TrigFTKSim/skeleton.FTKStandaloneMerge.py',
                                   inData = [tuple([ 'NTUP_FTKTMP_{0}'.format(subregion) for subregion in range(subregions) ])+('NTUP_FTKIP',)],
                                   outData = ['NTUP_FTKTMP'],
                                   extraRunargs = {'inputNTUP_FTKTMPFile': [ 'tmp.NTUP_FTKTMP_{0}'.format(subregion) for subregion in range(subregions)]},
                                   runtimeRunargs = {'MergeRegion': 'runArgs.bankregion[0]',
                                                     'FirstRegion': 'runArgs.bankregion[0]',
													 'TruthTrackTreeName' : "'truthtracks'",
													 'EvtInfoTreeName' : "'evtinfo'",
                                                     'SaveTruthTree' : '1'
													 },
                                   ))
    trf = transform(executor = executorSet, description = 'FTK Subregion simulate x {0} and merge.'.format(subregions))
 
    addFTKSimulationArgs(trf.parser)
    addTrigFTKSimOptions(trf.parser,nsubregions=subregions)
    addTrigFTKSimMergeOptions(trf.parser);
    addTrigFTKSimTFOptions(trf.parser)
    addTrigFTKSimRFOptions(trf.parser)
    return trf


def addFTKSimulationArgs(parser):
    # Add a specific FTK argument group
    parser.defineArgGroup('TrigFTKSim', 'Fast tracker simulation options')

    parser.add_argument('--bankregion', type=trfArgClasses.argFactory(trfArgClasses.argIntList, runarg=True),
                            help='Bank region number', group='TrigFTKSim', nargs='+')
    parser.add_argument('--sectorpath', type=trfArgClasses.argFactory(trfArgClasses.argList, runarg=True),
                            help='sectors path file for all the subregions', group='TrigFTKSim', nargs='+')

    # Add named parameters for each subregion
    for subregion in range(subregions):
        parser.add_argument('--patternbank{0}path'.format(subregion), type=trfArgClasses.argFactory(trfArgClasses.argList, runarg=True), 
                            help='Pattern bank path file, subregion {0}'.format(subregion), group='TrigFTKSim', nargs='+')

    
    parser.defineArgGroup('TrigFTKMerge', 'Fast tracker simulation merge options')

    # File handling
    parser.add_argument('--inputNTUP_FTKIPFile', 
                        type=trfArgClasses.argFactory(trfArgClasses.argNTUPFile, runarg=True, io='input', type='ntup_ftkiptmp', treeNames='ftkhits'),
                        help='FTK RDO file in ROOT format'.format(subregion), group='TrigFTKMerge', nargs='+')
    parser.add_argument('--inputTXT_FTKIPFile', 
                        type=trfArgClasses.argFactory(trfArgClasses.argFTKIPFile, runarg=True, io='input', type='txt_ftkip'), 
                        help='Wrapper files (in .dat.bz2 format)', group='TrigFTKSim', nargs='+')
    parser.add_argument('--outputNTUP_FTKTMPFile', '--outputNTUP_FTKFile', 
                        type=trfArgClasses.argFactory(trfArgClasses.argNTUPFile, runarg=True, io='output', type='ntup_ftk', treeNames='ftkdata'),
                        help='Subregion merged FTK file'.format(subregion), group='TrigFTKMerge',nargs='+')
    
    # The following for testing only
    for subregion in range(subregions):
        parser.add_argument('--inputNTUP_FTKTMP_{0}File'.format(subregion), 
                        type=trfArgClasses.argFactory(trfArgClasses.argNTUPFile, runarg=True, io='input', type='ntup_ftkiptmp', treeNames='ftkdata'),
                        help='FTK NTUP file from subregion {0} (for testing only)'.format(subregion), group='TrigFTKSim', nargs='+')
        parser.add_argument('--outputNTUP_FTKTMP_{0}File'.format(subregion), 
                        type=trfArgClasses.argFactory(trfArgClasses.argNTUPFile, runarg=True, io='output', type='ntup_ftkiptmp', treeNames='ftkdata'),
                        help='FTK NTUP file from subregion {0} (for testing only)'.format(subregion), group='TrigFTKSim', nargs='+')
      
        
    

if __name__ == '__main__':
    main()
