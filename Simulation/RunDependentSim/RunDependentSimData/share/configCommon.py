# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

include.block('RunDependentSimData/configCommon.py')
from Digitization.DigitizationFlags import jobproperties
if not 'logging' in dir(): import logging
digilog = logging.getLogger('Digi_trf')

import math
if 'runArgs' in dir():
    if hasattr(runArgs,"jobNumber") and hasattr(runArgs,"maxEvents"):
        trfJobNumber = runArgs.jobNumber
        trfMaxEvents = runArgs.maxEvents
        trfTotalEvents = runArgs.maxEvents
        trfSkipEvents = runArgs.skipEvents if hasattr(runArgs, "skipEvents") else 0

        # do executor step filtering
        if hasattr(runArgs, "totalExecutorSteps") and runArgs.totalExecutorSteps > 1:
            JobMaker = list(filter(lambda lb: 'step' not in lb or lb['step'] == runArgs.executorStep, JobMaker))
            if runArgs.totalExecutorSteps != len(runArgs.executorEventCounts):
                raise ValueError("Mismatch between total executor steps and event fractions size!")
            trfMaxEvents = runArgs.executorEventCounts[runArgs.executorStep]
            trfSkipEvents = runArgs.executorEventSkips[runArgs.executorStep]
            DoNotCorrectMaxEvents = True

        if runArgs.maxEvents==-1:
            raise SystemExit("maxEvents = %d is not supported! Please set this to the number of events per file times the number of files per job."%(runArgs.maxEvents,))
        if not 'DoNotCorrectMaxEvents' in dir():
            corrMaxEvents = math.ceil(float(trfMaxEvents)/100.0)*100.0 # round up to nearest 100 events..
        else:
            if not (hasattr(runArgs, "totalExecutorSteps") and runArgs.totalExecutorSteps > 1):
                digilog.warning('Using the actual number of HITS input events for this job -- not for production use!')
            corrMaxEvents = trfMaxEvents
    else: 
        raise SystemExit("Please provide jobNumber and maxEvents to runArgs.") 
else:
    #this is a test job not a trf job
    trfJobNumber=1
    trfMaxEvents=10
    trfTotalEvents=10
    corrMaxEvents=float(trfMaxEvents)
    
#We may need to repeat this run for long production jobs.
#NB: unlike vanilla variable-mu jobs, it's possible to waste
#  up to trfMaxEvents-1 events per complete run in prodsys if
#  the number of events specified by this run is not evenly
#  divisible by trfMaxEvents.
runMaxEvents=sum(lb['evts'] for lb in JobMaker)
digilog.info('There are %d events in this run.' % runMaxEvents)
jobsPerRun=int(math.ceil(float(runMaxEvents)/corrMaxEvents))
digilog.info('Assuming there are usually %d events per job. (Based on %d events in this job.)', corrMaxEvents, trfMaxEvents)
digilog.info('There must be %d jobs per run.' % jobsPerRun)

# Override event numbers with sequential ones if requested
sequentialEventNumbers = True if 'SequentialEventNumbers' in dir() and SequentialEventNumbers else False
if sequentialEventNumbers:
    digilog.info('All event numbers will be sequential.')

# Random mu sampling
randomMuSampling = True if 'RandomMuSampling' in dir() and RandomMuSampling else False
if randomMuSampling:
    digilog.info('Mu values will be sampled randomly from the set profile.')
    #Load needed tools 
    from RunDependentSimComps.RunDependentMCTaskIterator import getRandomlySampledRunLumiInfoFragment
    fragment=getRandomlySampledRunLumiInfoFragment(jobnumber=(trfJobNumber-1),task=JobMaker,maxEvents=trfMaxEvents,totalEvents=trfTotalEvents,skipEvents=trfSkipEvents,sequentialEventNumbers=sequentialEventNumbers)
else:
    #Load needed tools 
    from RunDependentSimComps.RunDependentMCTaskIterator import getRunLumiInfoFragment
    fragment=getRunLumiInfoFragment(jobnumber=(trfJobNumber-1),task=JobMaker,maxEvents=trfMaxEvents,totalEvents=trfTotalEvents,skipEvents=trfSkipEvents,sequentialEventNumbers=sequentialEventNumbers)

from RunDependentSimComps.RunLumiConfigTools import condenseRunLumiInfoFragment
digilog.info( 'Writing RunDMC trigger configuration fragment to file.  listOfRunsEvents = %s' %
              condenseRunLumiInfoFragment(fragment,"RunDMCTriggerRunsInfo.py"))

jobproperties.Digitization.RunAndLumiOverrideList=fragment
