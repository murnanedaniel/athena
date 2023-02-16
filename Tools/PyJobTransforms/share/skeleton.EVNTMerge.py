#############################
## basic jobO configuration
include("PATJobTransforms/CommonSkeletonJobOptions.py")
## load pool support
import AthenaPoolCnvSvc.ReadAthenaPool
import AthenaPoolCnvSvc.WriteAthenaPool

## input
ServiceMgr.EventSelector.InputCollections = runArgs.inputEVNTFile

## output stream
from AthenaPoolCnvSvc.WriteAthenaPool import AthenaPoolOutputStream

outStream = AthenaPoolOutputStream("StreamEVGEN", runArgs.outputEVNT_MRGFile, noTag=True)

## copy everything from the input file
## must force reading of all input objects
outStream.TakeItemsFromInput = True


#==============================================================
# Job Configuration parameters:
#==============================================================
## Pre-exec
if hasattr(runArgs, "preExec"):
    for cmd in runArgs.preExec:
        exec(cmd)
        
## Pre-include
if hasattr(runArgs, "preInclude"):
    for fragment in runArgs.preInclude:
        include(fragment)

# Avoid command line preInclude for Event Service
if hasattr(runArgs, "eventService") and runArgs.eventService:
    import AthenaMP.EventService
        
## Post-include
if hasattr(runArgs, "postInclude"):
    for fragment in runArgs.postInclude:
        include(fragment)

## Post-exec
if hasattr(runArgs, "postExec"):
    for cmd in runArgs.postExec:
        exec(cmd)

#==============================================================
## Special configuration
from AthenaConfiguration.AutoConfigFlags import GetFileMD
hepmc_version = GetFileMD(ServiceMgr.EventSelector.InputCollections).get("hepmc_version", None)
if hepmc_version is not None:
    log = logging.getLogger('EVNTMerge')
    if hepmc_version == "2":
        hepmc_version = "HepMC2"
    elif hepmc_version == "3":
        hepmc_version = "HepMC3"
    log.info('Input file was produced with %s', hepmc_version)
    ServiceMgr.TagInfoMgr.ExtraTagValuePairs.update({"hepmc_version": hepmc_version})              
########## EOF ###############




