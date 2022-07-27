####
# Setup VP1 jobOptions for running at P1
# --> do not rely on RecExCommon options (doVP1),
#     but setup things ourselves
#####

from AthenaCommon.AppMgr import ServiceMgr as svcMgr
from AthenaServices.AthenaServicesConf import OutputStreamSequencerSvc

outputStreamSequencerSvc = OutputStreamSequencerSvc()
outputStreamSequencerSvc.SequenceIncidentName = "EndEvent"
svcMgr += outputStreamSequencerSvc

### Add the algorithm producing VP1 events
from VP1AlgsEventProd.VP1AlgsEventProdConf import VP1EventProd
VP1EventProducer = VP1EventProd(InputPoolFile = StreamESD.OutputFile)
## =================== Added 09/03/15 by sjiggins ================= 
printfunc ("<<<<<<< VP1 Output File >>>>>>>")
printfunc ("OutputFile: %s" % StreamESD.OutputFile)
## ================================================================

#Write out files in the directory given by the stream name
VP1EventProducer.DestinationDirectory = "%s/.Unknown/" % OutputDirectory

#Set number of files large so deleting is doen by prune script
#Disable VP1 pruning by lshi on 19/May/2022
#The pruning is performed centrally by python/EventUtils.py
VP1EventProducer.MaxNumberOfFiles = -1

#Set the output level
if not 'VP1MsgLvl' in dir():
  VP1MsgLvl=WARNING
VP1EventProducer.OutputLevel=VP1MsgLvl

### Finally add this event producer to the main sequencer
from AthenaCommon.AlgSequence import AlgSequence
topSequence = AlgSequence()
topSequence += VP1EventProducer

### Finally print setup in debug mode
if VP1MsgLvl <= DEBUG:
  printfunc ("\n\n\t VP1 setup\n",VP1EventProducer,"\n\n")
