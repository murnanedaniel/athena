#!/bin/env python
#file RunDepTaskMaker.py
# Uses LumiCalc COOL query tools to get a dictionary of lumiblock configurations for use in the the EventIdModifierSvc.

import sys,getopt,os,itertools
from LumiBlockComps.LumiCalculator import lumiResult,RLBRange
from RunDependentSimComps.LumiResultsGetter import coolLumiResultsGetter

from AthenaCommon.Logging import logging
log = logging.getLogger("RunDepTaskMaker")
log.setLevel(logging.DEBUG)

#tool: format of dictionary entry
formatLine="'run':{run}, 'lb':{lb}, 'starttstamp':{starttstamp}, 'dt':{dt:.3f}, 'evts':_evts({evts:.0f}), 'mu':{mu:.3f}, 'force_new':False"

#tool: generator getting rid of bad elements (less than 1 event)
def getGoodLB(mylistIter):
    filteredout=0
    for x in mylistIter:
        if x.get('evts',0) >= 1:
            yield x
        else:
            filteredout += 1
    log.warning("Skipped %i lumiblocks because they would have contained less than one event.", filteredout)

#tool: generate the configuration preInclude file for digi trf.
def dressSkeleton( skeletonName, outName, dumpingIterator ):
    """ parse the skeleton.py file and create the file outName.
    outName is a copy of skeletonName with the space between #<xml comment escape> lines taken from dumpingIterator.
    Only good LB (contains nonzero 'evts' field) are dumped
    """
    from AthenaCommon.Include import FindFile,optionsPath
    sname = FindFile(skeletonName, optionsPath, os.R_OK )
    if not sname: raise IOError("Input skeleton file %s cannot be read." % skeletonName)
    import time
    s = open(sname,"r")
    outfile = open(outName,'w')
    outfile.write("#"*20)
    outfile.write("\n## File %s: autogenerated configuration file from command\n" % outName)
    outfile.write("##" + ' '.join(sys.argv) + "\n")
    outfile.write("## Created on %s\n" % time.asctime())
    outfile.write("#"*20+"\n")
    
    incomment = False
    dumpedAlready = False
    for line in s:
        if (line.find("#<!--") == 0): incomment = True
        if (line.find("#-->") == 0): incomment = False
        if not incomment:
            outfile.write(line)                        
        elif not dumpedAlready:
            for j in getGoodLB(dumpingIterator):
                outfile.write("   {"+formatLine.format(**j)+"},\n")
                pass
            dumpedAlready = True
            pass
        pass
    pass

class RunDepTaskDriver:
    def __init__(self):
        # defaults for using lumicalculator
        self.lumidbname='COOLOFL_TRIGGER/COMP200'
        self.statusdbname='COOLOFL_GLOBAL/COMP200'
        self.lumifoldername='/TRIGGER/OFLLUMI/LBLESTOFL'
        self.lumimethod='ATLAS_PREFERRED'
        self.generateMu=None      #for Beate's scenarios
        self.externaldict={}      #for Beate's scenarios
        self.tag='OflLumi-7TeV-002'
        self.outputfilename='configLumi.py'
        self.dbname='COOLONL_TRIGGER/COMP200' #can't change
        self.readoracle=False     #User can't change this
        self.runmin=0             #User can't change this
        self.runmax=0             #User can't change this
        self.lumimin=0            #User can't change this
        self.readxml=True         #User can't change this
        self.lumimax=(1 << 32)-2  #User can't change this
        self.skeletonName="RunDependentSimData/OverrideRunLBLumiDigitConfig.py" #User can't change this. 16.6.x.y
        self.trigger='L1_MBTS_1'  #Now needs changing
        self.useprescale=False    #User can't change this yet
        self.loglevel=0           
        self.statusreq=""         #User can't change this yet
        self.statustag=""         #User can't change this yet
        self.mcDatasetSize=10000  #New default! (assumes we override the mu distribution.)
        self.longpattern=False    #No need to cover lumiblocks from entire run, probably
        # process command line options and switches
        try:
            longopts=["database=","lumifolder=","lumimethod=","tag=","externalDict=", "nMC=","outfile=","trigger=", "debug","verbose","help", "longpattern"]
            opts,args=getopt.getopt(sys.argv[1:],'',longopts)
        except getopt.GetoptError,e:
            print e
            self.usage()
            sys.exit(1)
        self.filelist=args
        self.procopts(opts)
        self.execute()

    def usage(self):
        print ""
        print "usage: RunDepTaskMaker.py <options> [goodruns.xml]"
        print "Options are:"
        print "--database=<COOL databaseID>   : set COOL DB (%s) for luminosity " % self.lumidbname
        print "--tag=<COOL database tag>      : set COOL luminosity DB tag (%s)" % self.tag
        print "--trigger=<L1 trig>            : set high rate trigger (%s)" % self.trigger
        print "--lumifolder=<COOL foldername> : set COOL DB folder containing luminosity estimates(%s)" % self.lumifoldername
        print "--lumimethod=<method name>     : set name for lookup of lumi DB channel (%s)" % self.lumimethod
        print "--externalDict=<{mu:fraction}>: set python dictionary of ints/crossing (%s) -- lumimethod must be EXTERNAL" % self.externaldict
        print "--outfile=<filename.py>        : set output joboptions filename (%s)" %         self.outputfilename
        print "--nMC=<number of MC events>    : target number of events in MCProd task (%i) (or, repeat period for EXTERNAL lumimethod)" % self.mcDatasetSize
        print "--verbose : produce some output (IOV-level)"       
        print "--debug : produce maximum output (LB-level)"
        print "--longpattern : Try to use all of the LB in the goodrunlist (default is %s)" % self.longpattern
        print "--help : display this help text"
        

    def procopts(self,opts):
        "Process the command line parameters"
        for o,a in opts:
            if (o=='--database'): self.lumidbname=str(a)
            if (o=='--lumifolder'): self.lumifoldername=str(a)
            if (o=='--lumimethod'): self.lumimethod=str(a)
            if (o=='--tag'): self.tag=str(a)
            if (o=='--trigger'): self.trigger=str(a)
            if (o=='--externalDict'):
                self.externaldict=eval(str(a))
                self.generateMu = itertools.cycle(self.externaldict.items())
            if (o=="--outfile"): self.outputfilename=str(a)
            if (o=="--debug"): self.loglevel=2
            if (o=="--longpattern"): self.longpattern=True
            if (o=="--verbose"): self.loglevel=1
            if (o=='--nMC'): self.mcDatasetSize=int(a)
            if (o=='--help'):
                self.usage()
                sys.exit(0)
                pass
        return


    def trim(self, somelist):
        if self.lumimethod == 'EXTERNAL' : 
            s = len(somelist)
            mult = len(self.externaldict)
            if self.longpattern:
                print "Pattern has", mult, "entries. Truncating lumiblock list of length", s, "to", int(s/mult) * mult            
                del somelist[int(s/mult) * mult:]
            elif (s >= mult):
                print "Pattern has", mult, "entries. Truncating lumiblock list of length", s, "to", mult
                del somelist[mult:]
            pass
        pass
        
    def execute(self):
        
        from AthenaCommon.Include import FindFile,optionsPath
        sname = FindFile(self.skeletonName, optionsPath, os.R_OK )
        if not sname: raise IOError("Input skeleton file %s cannot be read." % self.skeletonName)

        # setup the tool
        log.info('Getting LumiResults using databases %s, %s:' % (self.dbname, self.lumidbname))

        lumitool=coolLumiResultsGetter(self.dbname,statusdbconn=self.statusdbname,coolfolderlumi=self.lumifoldername,useprescale=self.useprescale,
                                       lumidb=self.lumidbname, lumimethod=self.lumimethod,
                                       coolfoldertag=self.tag,
                                       loglevel=self.loglevel)
        rlblist=[]
        allresults=[]
        for file in self.filelist:
            rlblist    += lumitool.rlbFromFile(file, self.readxml)
            theresult   = lumitool.calcFromList(self.trigger,rlblist) #concat lists
            if (self.lumimethod == 'EXTERNAL') :                      #muss up the list for external profile
                lt=0
                l1=0
                newresult = []
                for t in theresult:
                    configItems = self.generateMu.next()             #config dict is {mu:fraction} so items are (mu, fraction)
                    newresult.append( (t[0],t[1],lumiResult(configItems[1],l1,0,0,lt,1,0),configItems[0],t[4]) )
                    pass
                allresults += newresult             
            else: allresults += theresult
            pass
        self.trim(allresults)                                                  #ensure integer number of pattern repeats
        if (self.lumimethod == 'EXTERNAL') : totl = 1.0                        #external profile puts mcDatasetSize * f events in each lb
        else: totl = sum([ a[2].intL for a in allresults if (a[2].intL > 0) ]) #filter nan to calculate total lumi
        lumi_averaged_mu = sum([ a[2].intL*(a[3]) for a in allresults if (a[2].intL > 0) and (a[3] > 0) ]) #filter nan to calculate avg. pileup  
        for a in allresults:
            if not (a[2].intL >= 0) : log.warning("Luminosity (%f) isn't a number: lumiblock %i[%i] skipped!!", a[2].intL, a[0], a[1])
            if not (a[3] >= 0) : log.warning("Mean number of collisions (%f) isn't a number: lumiblock %i[%i] skipped!!", a[3], a[0], a[1])
            pass        
        if totl == 0: log.warning("No luminosity in selected range.")        
        else:
            lumi_averaged_mu /= totl
            if self.lumimethod != 'EXTERNAL':
                log.info("Preparing a RunDMC task configuration object for %f nb^(-1) of integrated luminosity.", totl*1.0e-3)
            log.info("The average pileup in this task is %f", lumi_averaged_mu)
            f = self.mcDatasetSize/totl            
            JobMaker = [ {'run':a[0],
                          'lb':a[1],
                          'starttstamp':int(round(a[4] * 1e-09)),
                          'dt': (a[2]).livetime,
                          'evts': (a[2]).intL * f, 
                          'mu':a[3] }
                         for a in allresults ]
            l = len(JobMaker)
            log.info( "There are %i lumiblocks in the task.",l)
            if (l > 10):
                log.info( "Dumping first and last 5 lumiblocks:")
                for j in JobMaker[:5]+JobMaker[-5:]: print " ",j
            else:
                log.info( "Dumping allconfigured lumiblocks:")
                for j in JobMaker: print " ",j
            #now dump to outfile    
            allLB = iter(JobMaker)            
            dressSkeleton(self.skeletonName,self.outputfilename,allLB)

# main code
if __name__=='__main__':
    RunDepTaskDriver()


