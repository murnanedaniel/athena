# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration

"""
class Minder. Top of the Minder Heirarchy.

- Takes care of the state machine. This is a mixture of
a status object and a dispatcher.

The state machine  has states: queued, running, success, error,
it nows how to proceed from one state to the next, and species the actions to be
performed when the state changes. The way jobs are submited, and the names of
log files are examples of items that change with mode.

This class also has the following methods, which are common to all Minders,

poll
fromQueued
fromRunning
fromSuccess
fromError
setStatus
setStatusList
forceCloseDown
setDone

Many of these are overidden in derived classes.
"""

from MinderBase                     import MinderBase
from formatCollection               import formatCollection
from Defs                           import RTTdefs
from ModuleLoader                   import ModuleLoader
from exc2string2                    import exc2string2
try:
    from RTTTestRunner                  import RTTTestRunner  
except Exception, e:
    print str(e)
from UserLogger                     import UserLogger
from Tools                          import unique, timed_cmd
from Tools2                         import toRelPath, popenCmd
from BigFileIO                      import createBigFileIO
from RTTdict                        import RTTdict
from RTTSException                  import RTTCodingError
from DDDescriptor                   import DDDescriptor
from PPPaths                        import PPPaths

import logging
logger = logging.getLogger('rtt')

import os, shutil, fnmatch, shelve, time

# -------------------------------------------------------------------------
#from Tracer        import Tracer
# uncomment to monitor all method calls (shows args + return values)
#__metaclass__ = Tracer


from MethodTimer import MethodTimer



class Minder(MinderBase):
    class DiskSpaceCalculator:
        def __init__(self, runpath, respath, logger):
            self.runpath = runpath
            self.respath = respath
            self.logger  = logger # minder logger
            self.run = []
            self.res = []

        def du(self, path):
            try:
                ans = timed_cmd('du -hs %s' % self.runpath)
                return ans[0].split()[0].strip()
            except:
                self.logger.error('Unable to du %s' % path)
                self.logger.error('Answer was: %s' % ans)
                self.logger.error(exc2string2())
                return '-1'
            
        def __call__(self, phase):
            self.run.append((self.du(self.runpath), time.time(), phase))
            self.res.append((self.du(self.respath), time.time(), phase))

    def __init__(self, argBag, jDescriptor):
        """ One function of the minder is to bring together information
        from the job descriptor (gathered from the jobs config file),
        the job groups file (gathered from the job groups config file),
        and the Paths object (gathered from the RTT config file)

        For files to be transferred from the run directory to the results
        area, we use patterns (for file name matching). Lists of patterns
        are made at in this method, with a list per use (keep files). These lists will be used to identify
        filenames at the appropriat e moment, and converted to absolute
        path names. Source and destination will be kept in the dictionaries
        self.keepFiles
        
        """



        MinderBase.__init__(self, argBag.logDir,
                            argBag.package,
                            jDescriptor.identifiedName,
                            jDescriptor.jobSerialNumber,
                            argBag.elementCreator,
                            argBag.textNodeCreator)

        self.shareArea           = argBag.shareArea
        self.deleteNonKeepFiles  = argBag.deleteNonKeepFiles
        self.confFile            = argBag.confFile
        self.resultsBasePath     = argBag.resultsBasePath
        self.isNightly           = argBag.isNightly
        self.workDirs            = argBag.workDirs
        self.site                = argBag.site

        self.RTTLibDir           = argBag.RTTLibDir         # used for post processing
        self.RTTSrcDir           = argBag.RTTSrcDir         # used for post processing
        self.dCubeCfgFile        = argBag.dCubeCfgFile      # used for post processing
        self.installArea         = argBag.installArea       # used for post processing
        self.containerPackage    = argBag.containerPackage  # used for post processing
        self.runType             = argBag.runType           # used for post processing
        self.localRTTRun         = argBag.localRTTRun       # used for post processing
        self.distArea            = argBag.distArea          # used for post processing
        self.cmtPath             = argBag.cmtPath           # used for post processing
        self.cmtConfig           = argBag.cmtConfig         # used for post processing
        self.cmtLinesCmds        = argBag.cmtLinesCmds      # used for post processing
        self.topProject          = argBag.topProject        # used for post processing
        self.branch              = argBag.branch            # used for post processing
        self.release             = argBag.release           # used for post processing
        self.package             = argBag.package           # used for post processing
        self.archivers           = [a.duplicate() for a in argBag.archivers]
        
        self.runPath             = jDescriptor.runPath
        self.name                = jDescriptor.name
        self.jobDisplayName      = jDescriptor.jobDisplayName
        self.jobGroup            = jDescriptor.jobGroup
        self.jobDocString        = jDescriptor.jobDocString    
        self.jobDocURL           = jDescriptor.jobDocURL    
        self.rttPilotJob         = jDescriptor.rttPilotJob
        self.rttATNJob           = jDescriptor.rttATNJob
        self.log                 = jDescriptor.log
        self.elog                = jDescriptor.elog
        self.errorMessages       = jDescriptor.errorMessages
        self.resPath             = jDescriptor.resPath
        self.descDataForXMLNode  = jDescriptor.dataForXMLNode()
        self.keepFiles           = jDescriptor.keepFiles
        self.datasets            = jDescriptor.datasets
        self.missingDatasets     = jDescriptor.missingDatasets
        self.trendId             = jDescriptor.trendId
        self.batchWallTimeLimit  = jDescriptor.suggestedBatchWallTime

        # Combine some stuff from kit and job
        self.keepFilePatterns    = jDescriptor.keepFilePatterns + argBag.jobGroupKit['keepFilePatterns']
        self.auxFilePatterns     = jDescriptor.auxFilePatterns + argBag.jobGroupKit['auxFilePatterns']
        self.actionDescriptors   = jDescriptor.actions + argBag.jobGroupKit['actions']
        self.testDescriptors     = jDescriptor.tests + argBag.jobGroupKit['tests']
        
        self.testDBPath          = os.path.join(self.runPath, 'RTTtests.db')

        # This list holds those files that should never be copied or deleted from the run dir by this minder
        # Currently it will hold chain files (absolute paths).
        self.neverCopyAndNeverDelete = []

        # Filled in by the archivers
        self.keepFilesToVeto = []
        
        # Calculate disk space usage by the job. Stores history of usage by the job.
        self.diskspace = Minder.DiskSpaceCalculator(self.runPath,self.resPath, self.logger)

        # These jobs are NOT counted in  Launcher stats. Otherwise overide.
        self.jobWeight = 0 

        # bring together information on files to manipulate.
        # Keep and aux files are handled as patterns - the actual file
        # name resolution is deferred to later
        # Reference files references are complete file names due to
        # fact that they will go into keys before retirval - and patterns wont work
        

        self.finishedTests        = [] # a list of execcuted test objects.
        self.checkResults         = [] # results of tests which have been shelved.
        self.processingResult     = 'unavailable'
        self.postProcessingResult = 'unavailable'
        self.chainLength          = 0
        self.chainSuccesses       = 0
        self.chainFileCopier      = None # Might be reset by derived class
        
        # values for the commands are  worked out in derived class
        self.postProcessCommand   = argBag.postProcessCommand
        self.submitCommand        = argBag.submitCommand
        
        # declare dictionaries with source and destination paths
        # for files to be transferred from the run directories to the resulsts
        # directory.

        # declare lists of files - no  wild cards, as the elements of the list
        # are used in data base keys. Also,  no absolute paths for the same reason
        # (would break portability)
        
        
        # maximum no of times to retry a job if it enters the error state.
        self.errorStateMaxRetry = 3
        self.errorStateCurRetry = 0
          

        self.exitStatus  = 'Unknown'
        self.cpuTime     = 'Unknown'
        self.cpuTime2000 = 'Unknown'
        self.mem         = 'Unknown'
        self.vmem        = 'Unknown'
        self.wallTime    = 'Unknown'
        self.batchStatus = 'Unknown'
        

        self.scripts             = {}
        self.timedOut            = False
        # so complete file names go here here as well.

        self.handleFileNames()
        self.collectActionsAndTests()

        self.collectScripts(argBag.scriptWriter) # after collecting tests and actions
        # self.setupJobGroupDirs() # after making scripts

   
        jobGroupDirectoryMaker        = argBag.jobGroupDirectoryMaker
        jobGroupDirectoryMaker.logger = self.logger

        self.makeJobGroupDirs(jobGroupDirectoryMaker)
        self.setupJobGroupDirs(jobGroupDirectoryMaker)

        self.shelveTests()

    
    def makeJobGroupDirs(self, jobGroupDirectoryMaker):
        try:
            jobGroupDirectoryMaker.makeJobGroupDirs(self.runPath, self.resPath)
        except:
            msg  = 'Exception making directories by %s traceback:\n%s' % (self.identifiedName, exc2string2())
            self.logger.error(msg)
            raise RTTCodingError(msg)
        

    def unlink(self):
        'break circular references to allow garbage collection'

        self.stateEngine.unlink()
        del self.stateEngine
        self.xmlConverter.unlink()
        del self.xmlConverter
        MinderBase.unlink(self)

        
    def calcPostProcessingResult(self):
        if not (self.finishedTests or self.checks):
            self.postProcessingResult = 'no tests'
            return

        for t in self.finishedTests:
            if t.result != 0:
                self.postProcessingResult = 'error'
                return

        if 'error' in self.checkResults:
            self.postProcessingResult = 'error'
            return
        
        self.postProcessingResult = 'success'

    def collectScripts(self, scriptWriter):
        # post processing scripts are handled here, and so
        # apply to Watcher and Worker minders.
        # NOTE 26/01/09
        # self.postProcessingCommand is minder specific.
        # it is set in LinuxInteractive and LSFBatch, not yet set in WatcherMinder
        
        if not (self.tests or self.actions): return

        script = scriptWriter.postProcessorScript()

        if not self.postProcessCommand:
            msg = '%s using self.postProcessing command without it having been set to something useful' % self.__class__.__name__
            self.logger.error(msg)
            raise RTTCodingError(msg)
        
        dict = {'postProcessorScript': {'script':     script,
                                        'scriptName': 'postProcessorScript.sh',
                                        'cmd':         self.postProcessCommand
                                        }
                }
        self.scripts.update(dict)



# ------------------------------------------------------------------------

    def shelveTests(self):
        # associate string identifiers with each test and action
        # t.position is the user given order index
        npos = 0

        testIds = []
        for runOrder, instanceList in self.tests.items():
            for testInstance in instanceList:
                item = (testInstance, '%s_%d' % (testInstance.__class__.__name__,npos), runOrder)
                testIds.append(item)
                npos +=1

        actionIds = []
        for runOrder, instanceList in self.actions.items():
            for actionInstance in instanceList:
                item = (actionInstance, '%s_%d' % (actionInstance.__class__.__name__,npos), runOrder)
                actionIds.append(item)
                npos +=1
                
        if npos == 0: return

        # provide a database to store the tests
        db = shelve.open(self.testDBPath, 'n')
        db.close()

        # wrap the tests to allow transportation to the computing node
        wrappers = []
        
        isAction = False
        wrappers.extend([RTTTestRunner(self.testDBPath, t, isAction) for t in testIds])
        isAction = True
        wrappers.extend([RTTTestRunner(self.testDBPath, t, isAction) for t in actionIds])

        # Write the tests to a shelve db
        for w in wrappers:
            try:
                w.autoShelve()
            except:
                m = 'Error shelving %s\nTraceback: %s' % (w.testId, exc2string2())
                self.logger.error(m)
        
    # -------------------------------------------------------------------------
    # 09/03/21 PS This method is never called 
    # def hasChecks(self):
    #    answer = False

    #    if hasattr(self,'checks'):
    #        answer = len(self.checks) > 0 or answer

    #    if hasattr(self,'tests'):  
    #        answer = len(self.tests.keys()) > 0 or answer

    #    return answer

# ------------------------------------------------------------------------

    def collectActionsAndTests(self):
        """
        Test Descriptor is a list with four elements:
        position, modulename, classname, dictionary of arguments
        """


        # set up a "pruned descritpor". In 01/09, there was a large refactoring of the RTT
        # and the role of paths and descriptors was reduced. However user suppied test scripts may
        # be expecting these variables, so they are rebuilt here.
        # class DDDescriptor:pass
        # class Paths: pass
        def makePrunedDescriptor():
            desc = DDDescriptor()
            possibleAttrs = ['log', 'name', 'trendId', 'identifiedName', 'jobSerialNumber', 'jobDisplayName',
            'hashString', 'runPath', 'resPath', 'jobGroup', 'jobDocString', 'jobDocURL', 'package']
            [setattr(desc, p, getattr(self, p)) for p in possibleAttrs if hasattr(self, p)]


            paths = PPPaths()
            possibleAttrs =  ['installArea', 'containerPackage', 'cmtPath', 'shareArea',
                              'release', 'cmtHomeDir', 'localRTTRun', 'workDirs', 'resultsDirs', 'runType',
                              'build', 'originalBranch', 'branch', 'topProject', 'otherProject', 'cmtConfig',
                              'dCubeCfgFile', 'isNightly', 'RTTSrcDir', 'RTTLibDir', 'pathSegment', 'package',
                              'packageTag', 'fullPackageName', 'dCubeRefBaseDir', 'confFile', 'userHomeCmtDir',
                              'testProjectCmtDir', 'distArea', 'projectCMTDirs']
            [setattr(paths, p, getattr(self, p)) for p in possibleAttrs if hasattr(self, p)]

            desc.paths = paths
            return desc

        def process(descs, tgt):
            for descriptor in descs:
                # Create the pruned descriptor
                prunedDesc  = makePrunedDescriptor()
                # prunedDesc  = self.jDescriptor.prune()
                # prunedPaths = self.jDescriptor.paths.prune()
                # setattr(prunedDesc, 'paths', prunedPaths)

                # Add params to the dict
                descriptor.addParameter(('JobDescriptor', prunedDesc))
                descriptor.addParameter(('cmtLinesCmds', self.cmtLinesCmds[:]))
                instance = self.getClassFromDescriptor(descriptor)
                if instance != None:
                    tgt.setdefault(descriptor.position, []).append(instance)

        process(self.testDescriptors, self.tests)
        process(self.actionDescriptors, self.actions)
          
#-----------------------------------------------------------------------

    def getClassFromDescriptor(self, taDescriptor): # descriptor is a test or action descriptor
        self.logger.debug("Entering getClassFromDescriptor: %s" % taDescriptor.argDict)
        def runUserScriptInAtlasEnv():
            return taDescriptor.runInAtlasEnv!=""
        def createRTTdict(rtt_dict_keys, paramDict):
            self.logger.debug('Inside createRTTdict, paramDict: %s, rtt_dict_keys: %s' % (paramDict,
                                                                                          rtt_dict_keys))
            rtt_dict = {}
            for key in rtt_dict_keys:
                rtt_dict[key] = paramDict[key]
                del paramDict[key]
            user_dict = paramDict
            return RTTdict(rtt_dict, user_dict)
        
        # allow RTT library tools to be picked up from the local code base
        # rather than the installed shared area if so requested in the
        # RTT configuration file
        sharePathToUse = {'RttLibraryTools' : self.RTTLibDir}
        instance = None
        self.logger.debug("share path to use is: %s" % sharePathToUse)

        # RTTLibTools has been split into modules with one class per
        # module, with the module name = the class name. Needed this
        # for shelve to unshelve instances.

        moduleName = taDescriptor.moduleName
        sharePathIs = sharePathToUse.get(moduleName, self.shareArea)
        m = 'dict: %s, key: %s, value: %s' % (sharePathToUse, moduleName, sharePathIs)
        self.logger.debug(m)
        
        # now set module name = testname for shelving constraint
        moduleName = taDescriptor.testName
        className  = taDescriptor.testName
        self.logger.debug("modulename set = classname, i.e. %s" % moduleName)
        paramDict  = taDescriptor.argDict

        self.logger.debug("paramdict: %s" % paramDict)
        
        rtt_dict_keys = ['JobDescriptor']

        if runUserScriptInAtlasEnv():
            sharePathIs = sharePathToUse.get('RttLibraryTools')
            paramDict['userScript'] = moduleName
            paramDict['userScriptLoc'] = sharePathToUse.get(taDescriptor.moduleName,
                                                            self.shareArea)
            moduleName = className  = 'AtlasEnvSetup'
            
            rtt_dict_keys.extend(['userScript', 'userScriptLoc'])
          
        paramDict = createRTTdict(rtt_dict_keys, paramDict)
        
        # pickle the dictionary to the run path: actually gets dumped in JobGroupDirectoryMaker.setupRunDir
        self.rttArgDictPickle = {'what': paramDict, 'where': os.path.join(self.runPath,'rtt.argdict.cpickle')}
        
        # default is release share path if key not in dict
        # this happens when the user supplies own module, which is
        # picked up from the release area.

        self.logger.debug('module path ' + sharePathIs)
        self.logger.debug('moduleName  ' + moduleName)
        self.logger.debug('className   ' + className)
        self.logger.debug('paramDict   ' + str(paramDict))
        self.logger.debug('logger      ' + str(self.logger))

        try:            
            mL        = ModuleLoader(moduleName, sharePathIs, self.logger)
            instance  = mL.getClassInstance(className, paramDict)
            self.logger.debug('Created a test instance of class %s' % (instance.__class__.__name__))
        except Exception, e:
            self.logger.error('Failed to create a test instance')
            self.logger.error(exc2string2())
            self.logger.error(str(e))
            
        return instance
    
#------------------------------------------------------------------------

    def handleFileNames(self):
        """ method to collect file names, and resolve patterns where
        for those patterns for which it is possible at init."""
   
        # expand the wild cards - but do not create the full directory path
        # as the work sub directories have yet to be created.
        if not os.path.exists(self.shareArea):
            m = 'Cannot set self.auxfiles due to non-existent share directory: %s' % self.shareArea
            self.logger.fatal(m)
            raise RTTCodingError(m)

        # resolve auxFile patterns to file names
        auxFiles = []
        for pattern in self.auxFilePatterns:
            base, fnpattern = os.path.split(pattern)
            srcDir = os.path.normpath(os.path.join(self.shareArea, base))
            filesInShare = os.listdir(srcDir)
            auxFiles.extend([os.path.join(base,file) for file in filesInShare if fnmatch.fnmatch(file, fnpattern)])

        self.auxFiles = unique(auxFiles)

    # ------------------------------------------------------------------------
    # 09/03/31 This method is never called from RTT code base
    # def reportUserScriptError(self,message):
    #    filePath  = os.path.join(self.runPath,'Python_Script_Output.log')
    #    h = open(filePath,'a+')
    #    h.write(message)
    #    h.close()

    #-----------------------------------------------
    
    def readTestResults(self):
        #make UserLogger to transfer messages
        #from UserLogger import UserLogger
        ulogger = UserLogger().makeLogger(self.runPath, self.identifiedName)
        self.logger.debug('made logger for tests and actions :'+str(logger))

        db = shelve.open(self.testDBPath)

        dbValues = []
        try:
            dbValues.extend(db.values())
            self.finishedTests = db.values()
        except Exception, e:
            m  = 'Exception thrown trying to read values from %s\n' % str(self.testDBPath)
            m += str(e)            
            self.logger.error(m)
            # now close DB to prevent too many open files
            db.close()
        else:
            db.close()
        
        self.finishedActions = [t for t in self.finishedTests if t.isAction]
        [self.finishedTests.remove(t) for t in self.finishedActions]
        
        # for t in db.values():
        for t in dbValues:
            if t.error:
                self.logger.error('Error running test %s' % t.testId)
            else:
                self.logger.debug('No error running test %s' % t.testId)


            # collect any messages generated by the test
            t.transferLog(ulogger)

            delimeter  = "\n***********************************\n"
            delimeter += "           NEW ACTION/TEST         \n"
            delimeter += "***********************************\n"
            ulogger.debug(delimeter)

        self.calcPostProcessingResult()
       
    #--------------------------------------------------
    
    #def runMoniActions(self):
    #jobName  = self.jDescriptor.identifiedName
    #nActions = len(self.monActions)
    #msg      = 'Running %d monActions for job %s' % (nActions, jobName)
    #self.logger.debug(msg)
    #
    #for action in self.monActions:
    #className = str(action.__class__.__name__)
    #try:
    #self.logger.debug("Running monAction " + className)
    #dataDict=action.run()
    #self.monActionsData.append(dataDict)
    #
    # except Exception, e:
    # msg = "Could not run monitoring action " + className
    # self.logger.error(msg)
    # self.logger.error(str(e))
    # self.logger.error(exc2string2())
    # msg  = '-----------------------------------------\n'
    # msg += 'MoniAction: ' + className + ' could not be run!\n'
    # msg += exc2string2() + '\n'
    # msg += '-----------------------------------------\n\n'
    # self.reportUserScriptError(msg)
    # self.logger.debug("Running next available monAction")
    #       
    #
    # self.logger.debug('Minder moniData :'+str(self.monActionsData))
    # self.logger.debug("Finished running monActions for %s" % jobName)
      
    #--------------------------------------------------

    def fullResultsFileName(self, filename):
        return os.path.join(self.resPath, filename)

    #-----------------------------------------------
    
    def makeKeepFileEntry(self, file, infoString='',displayColor='', md5=''):
        "helper method for registerWildKeepFiles()"
        
        src = os.path.join(self.runPath, file)

        dest = {'keepFileString': self.fullResultsFileName(file),
                'infoString'    : infoString,
                'displayColor'  : displayColor,
                'md5sum'        : md5}
        
        self.keepFiles[src]=dest
        return dest['keepFileString']
    #-----------------------------------------------
    
    def registerWildKeepFiles(self):
        """
        Common implementation task.

        Obtain the wild card patterns for the current job group from
        a JobGroupKit.
        
        Collect all the wild carded keep files.
        These are files that are present in the run directory at when
        fromRunning is called, and which match a pattern in the pattern
        list.

        Give these files their full paths and add them to the keepFile
        dictionary.
        """

        for card in self.keepFilePatterns:

            keepString   = card['keepFileString'] # may contain subdirs
            infoString   = card['infoString']
            displayColor = card['displayColor']

            keepStringTokens = keepString.strip().split('/')

            if len(keepStringTokens) == 1: # current dir keepfile pattern
                wildFiles = fnmatch.filter(os.listdir(self.runPath), keepString)
                [self.makeKeepFileEntry(file, infoString, displayColor) for file in wildFiles]

            elif len(keepStringTokens) > 1: # subdirs
                matches = ['']
                for pathPart in keepStringTokens:
                    newMatches = []
                    for match in matches:
                        conts = os.listdir(os.path.join(self.runPath, match))
                        newMatches.extend([os.path.join(match, f) for f in fnmatch.filter(conts, pathPart)])

                    matches = newMatches
                [self.makeKeepFileEntry(m, infoString, displayColor) for m in matches]

        # now manually add the package configuration file to keep files
        self.makeKeepFileEntry(os.path.basename(self.confFile),
                               "Package XML test configuration file",
                               "#cc3333")


    #-----------------------------------------------

    def cleanSpace(self):
        self.removeNonKeepFiles()

        # If runpath != respath we need to delete keepfiles in the run dir
        # else they'll stick around.
        if self.runPath != self.resPath:
            allfiles = [os.path.join(self.runPath, o) for o in os.listdir(self.runPath)]
            toDelete = [a for a in allfiles if a not in self.neverCopyAndNeverDelete]
            self.logger.debug('Runpath!=Respath, deleting %d keepfiles from run dir' % len(toDelete))
            self.logger.debug('Files to delete:\n%s' % str(toDelete))
            for a in toDelete:
                os.system('rm -rf %s' % a) 

    # --------------------------------------------------------------------

    def getBigFiles(self):
        """Looks in the run/results dir for big files. Returns a list of them."""
        toMove = {self.resPath:[], self.runPath:[]}

        for tgt in toMove.keys():
            for root, dirs, files in os.walk(tgt):
                files = [os.path.join(root, f) for f in files]
                def isBig(f): # only files, no links
                    return os.path.isfile(f) and not os.path.islink(f) and os.path.getsize(f)>=self.bigFilesSize
                toMove[tgt].extend([f for f in files if isBig(f)])
        return toMove
    
    # --------------------------------------------------------------------

    def moveBigFile(self, src, dst):
        try:
            fileSize = os.path.getsize(src)
            shutil.move(src, dst)
            m = 'Successfully moved big file %s [%s bytes] to %s' % (os.path.basename(src),
                                                                     str(fileSize),
                                                                     dst)
            self.logger.info(m)
        except:
            m  = 'Unable to move big file %s to %s\n' % (os.path.basename(src), dst)
            m += exc2string2()
            m += str(os.listdir(os.path.dirname(src)))
            self.logger.error(m)
            
    # --------------------------------------------------------------------

    def makeReplacementKeepFile(self, bigFileSrc, bigFileDump):
        """If we move a big file from run path to separate volume,
        and that file just happened to be a keep file, then we need
        to replace it. Create a soft link in resPath to separate volume
        and add an expiry date."""

        linkName = os.path.basename(bigFileSrc)
        linkTgt  = os.path.join(toRelPath(os.path.dirname(bigFileSrc), bigFileDump), linkName)
        
        currDir = os.getcwd()
        os.system('cd %s;ln -s %s %s;cd %s' % (os.path.dirname(bigFileSrc), linkTgt, linkName, currDir))

        # add an md5sum to the keepfile to guard against mishap if the file pointed
        # at was overwritten (the web link would work, but you'd get the
        # overwriter, not the expected overwritee, and you would never know...
        # If md5sum different between now calculated and future web page one,
        # then grey-out the link.
        md5sum = popenCmd('md5sum %s' % bigFileSrc)
        self.logger.debug('Minder md5sum calculated is: %s' % str(md5sum))
        md5sum = md5sum[0].split()[0].strip()
        srcRelToRunPath = bigFileSrc.split(self.runPath)[1]
        if srcRelToRunPath.startswith('/'): srcRelToRunPath = srcRelToRunPath[1:]
        self.makeKeepFileEntry(srcRelToRunPath, md5=md5sum)
        
    # --------------------------------------------------------------------
    
    def moveBigFiles(self):
        """Move all files bigger than some value onto a separate volume."""
        if not self.bigFilesArea:
            self.logger.info('Moving of big files to a separate volume has not been requested.')
            return

        self.logger.info('Moving of big files to a separate volume is requested. Scanning...')
        
        if not os.path.exists(self.bigFilesArea):
            m = 'Cannot shift big files onto inexistent volume: %s' % self.bigFilesArea
            self.logger.error(m)
            return
        
        bigFiles = self.getBigFiles()

        if not [val for val in bigFiles.values() if val]:
            self.logger.info('No big files were found, returning.')
            return
        
        placeToDump = createBigFileIO(self.site, self.bigFilesArea, self.workDirs, self.isNightly).getJobDumpLocation(self)
        if not placeToDump:
            m = 'Unable to retrieve location of big files volume. Not moving big files.'
            self.logger.warning(m)
            return

        # We have files to move, let's move them
        for bigFileBaseDir, bigFiles in bigFiles.items():
            for bigFile in bigFiles:
                src = bigFile # file
                dst = placeToDump # directory
                self.moveBigFile(src, dst)
                # If big file origin is results path, replace with a soft link
                # to separate big file volume.
                if bigFileBaseDir == self.resPath:
                    self.makeReplacementKeepFile(bigFile, placeToDump)
            
    # --------------------------------------------------------------------                            

    def removeNonKeepFiles(self):
        """Delete all files in the run dir not explicitly kept."""

        if not self.deleteNonKeepFiles:
            self.logger.info('Deletion of non keep files (to save space) has not been requested.')
            return

        self.logger.info('Deletion of non-keep files is requested. Scanning...')
        
        toDelete = []
        for root, dirs, files in os.walk(self.runPath):
            toDelete.extend([os.path.join(root, f) for f in files if os.path.join(root, f) not in self.keepFiles.keys()])
            # os.walk ignores links, explicitly get them
            toDelete.extend([os.path.join(root, f) for f in os.listdir(root) if os.path.islink(os.path.join(root,f)) and os.path.join(root,f) not in self.keepFiles.keys()])

        for thing in os.listdir(self.runPath):
            self.logger.debug('removeNonKeepFiles::before delete: %s' % str(thing))
            
        for thing in toDelete:
            if thing in self.neverCopyAndNeverDelete:
                self.logger.debug('removeNonKeepFiles::Not deleting neverCopyAndNeverDelete file: \n%s' % thing)
                continue
            try:
                os.remove(thing)
                self.logger.debug('Deleted: %s' % thing)
            except:
                message  = 'Unable to delete non-keepfile: %s\n' % thing
                message += exc2string2()
                self.logger.error(message)

    # --------------------------------------------------------------------
    
    # copy files to be kept (log histos, ntuples...) to results dir
    def copyKeepFiles(self):
        
        # find if any of the wild keep cards have a match
        idName = str(self.identifiedName)
        try:
            self.registerWildKeepFiles()
        except Exception, e:
            msg = "Exception registering keep files for job: " + idName
            self.logger.error(msg)
            self.logger.error(str(e))
            self.logger.error(exc2string2())

        # Pop files in self.keepFiles dict that are in the self.keepFilesToVeto list
        popped = [self.keepFiles.pop(k) for k in self.keepFiles.keys() if k in self.neverCopyAndNeverDelete]
        self.logger.debug('%d neverCopyNeverDelete files were popped from the keepfile dict.' % len(popped))
        self.logger.debug('%s' % str(popped))

        popped = [self.keepFiles.pop(k) for k in self.keepFiles.keys() if k in self.keepFilesToVeto]
        self.logger.debug('%d keepFilesToVeto files were popped from the keepfile dict.' % len(popped))
        self.logger.debug('%s' % str(popped))

        msg  = 'About to transfer %d keep files' % len(self.keepFiles.keys())
        self.logger.debug(msg)

        for file in self.keepFiles.keys():
            if not os.path.exists(file):
                msg = "%s does not exist! Skipping transfer." % str(file)
                self.logger.error(msg)
                continue
                
            dictValue = self.keepFiles[file]
            desFile = dictValue['keepFileString'] # this is a filename

            fileName = os.path.basename(str(file))
            srcDir   = os.path.dirname(str(file))
            tgtDir   = os.path.dirname(str(desFile))

            if srcDir == tgtDir:
                msg  = 'Keepfile src and tgt dirs are same.\n'
                msg += '%s\n' % srcDir
                msg += '===> Not transferring file %s' % fileName
                self.logger.debug(msg)
                continue

            self.logger.debug("Moving file %s from %s to %s" % (fileName, srcDir, tgtDir))
            if not os.path.exists(os.path.dirname(desFile)):
                self.logger.debug('Destination dir does not exist. Making: %s' % os.path.dirname(desFile))
                try:
                    os.makedirs(os.path.dirname(desFile))
                except:
                    self.logger.error('Unable to create %s so as to copy keep file to it. Not copying it.' % os.path.dirname(desFile))
                    self.logger.error(exc2string2())
                    continue

            try:                    
                shutil.copy(file, os.path.dirname(desFile))
            except:
                # disk space problems?
                message = '***** COPYING OF KEEP FILES PROBLEM! *****\n'
                message += 'Unable to copy src file:\n'
                message += '   ' + file + '\n'
                message += 'to destination file:\n'
                message += '   ' + desFile + '\n'
                message += 'Doing command df on ' + desFile + ' yields the answer:\n'
                message += str(os.popen('df ' + desFile).read())
                message += '---------------------------------'                    
                self.logger.error(message)

        self.logger.debug("Finished copying keepfiles")

    # --------------------------------------------------------------------
        
    def runChecks(self):
        
        theJob = str(self.identifiedName)
        msg = 'Len self.checks is: %s for %s ' % (str(len(self.checks)),theJob)
        self.logger.debug(msg)
        
        for check in self.checks:

            check.setLogger(self.logger) # log to the minders logfile

            if not callable(check):
                msg = "Uncallable checker: %s for job %s" % (str(check),
                                                             theJob)
                self.logger.error(msg)
                self.checkResults.append('error')
                self.logger.error('Proceding to the next check')
                continue
            
            self.logger.debug('Job: %s calling check: %s'%(theJob, str(check)))
            status = ''
            try:
                rc = check(self) # int return code from check
                status = RTTdefs.status(rc) # convert to string
            except Exception, e:
                msg = 'Exception raised while executing %s' % str(check)
                self.logger.error(msg)
                self.logger.error(exc2string2())
                self.logger.error(str(e))
                self.checkResults.append('error')
                self.logger.error('Proceding to the next check')
                continue
            
            msg = 'Job: %s  check: %s return status %s' % (theJob,
                                                           str(check),
                                                           status)
            self.logger.debug(msg)
            
            self.checkResults.append(status)
        
        self.logger.debug(msg)

    # --------------------------------------------------------------------       
    """
    success if the job terminates, subclasses will override with a
    more precise meaning
    """
    def isSuccess(self): return self.isDone()
        
    # --------------------------------------------------------------------

    def setupJobGroupDirs(self, jobGroupDirectoryMaker):
        try:
            jobGroupDirectoryMaker.setupRunDir(self)
        except:
            msg  = 'Error  setting up run directory by %s: TraceBack:\n%s' % (
                self.identifiedName,
                exc2string2()
                )
            self.logger.warning(msg)
            raise RTTCodingError(msg)


    def dataForMonitoring(self):
        """
        return a dictionay of values for monitoring. repeats data in DOM
        document, but is faster to access.
        """
        dict = MinderBase.dataForMonitoring(self)
        
        dict['nTests']        = len(self.tests.keys())
        dict['done']          =   self.isDone()
        dict['nTestsSuccess'] = len([s for s in self.finishedTests if s.result == 0])
        dict['nTestsFailure'] = len([s for s in self.finishedTests if s.result != 0])
        dict['nRetries']      = self.errorStateCurRetry
        dict['ppFailure']     = (self.postProcessingResult == 'error')
        dict['ppSuccess']     = (self.postProcessingResult == 'success')

        return dict
    def __str__(self):
        s = '----------- Minder ---------------\n'
        s += ' done:                 '+str(self.isDone())+'\n'
        s += ' weight:               '+str(self.jobWeight)+'\n'
        s += ' runPath:              '+str(self.runPath)+'\n'
        s += ' keepFilePatterns:     '+formatCollection(self.keepFilePatterns)+'\n'
        s += ' auxFilePatterns:      '+formatCollection(self.auxFilePatterns)+'\n'
        s += ' keepFiles:            '+formatCollection(self.keepFiles)+'\n'
        s += ' auxFiles:             '+formatCollection(self.auxFiles)+'\n'
        s += ' actions:              '+formatCollection(self.actions)+'\n'
        s += ' tests:                '+formatCollection(self.tests)+'\n'
        s += ' scripts:              '+str(self.scripts.keys())+'\n'
        s += ' descriptor:\n'

        return s

    # --------------------------------------------------------------------
    
    def dump(self):
        self.logger.debug('|-------------------------------------------|')
        self.logger.debug('|                                           |')
        self.logger.debug('|            Minder  dump                   |')
        self.logger.debug('|                                           |')
        self.logger.debug('|-------------------------------------------|')
        self.logger.debug('\n'+self.__str__())
        
    def status(self):
        return self.stateEngine.state.state

    def printMethodTimes(self):
        
        m = '------ Minder  %s ---------\n%s' % (
            self.identifiedName,
            self.formatMethodTimes())


        #m += '\n-----State engine ---------  \n%s' % (
        #    self.stateEngine.formatMethodTimes())

        m += '\n-----XML Converter ---------  \n%s' % (
            self.xmlConverter.formatMethodTimes())
        
        self.logger.info(m)

# --------------------------------------------------------------------------
            
    def doPostProcessing(self):
        if self.processingResult == 'error': return False
        if 'postProcessorScript' in self.scripts.keys(): return True
        return False
            

    def setDone(self):
        self.cleanSpace()
        MinderBase.setDone(self)
