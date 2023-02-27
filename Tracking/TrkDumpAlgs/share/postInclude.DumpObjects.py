#----------------------------
# Algorithm
#----------------------------
from TrkDumpAlgs.TrkDumpAlgsConf import DumpObjects
dumpObjects = DumpObjects(name = "DumpObjects")
dumpObjects.NtupleFileName = '/DumpObjects/'
dumpObjects.NtupleTreeName = 'GNN4ITk'
dumpObjects.csvFile = True
dumpObjects.rootFile = True
topSequence += dumpObjects

#----------------------------
# Histogram and Tree Service
#----------------------------
def rootfile() :
    from AthenaCommon.AppMgr import ServiceMgr
    if not hasattr(ServiceMgr, 'THistSvc'):
        from GaudiSvc.GaudiSvcConf import THistSvc
        ServiceMgr += THistSvc()
    print("/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////")
    print ("Argument List:", str(sys.argv))
    print (sys.argv[3])
    print (dumpObjects.rootFile)
    print("/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////")
    if dumpObjects.rootFile :
        ServiceMgr.THistSvc.Output += ["DumpObjects DATAFILE='Dump_GNN4Itk.root' OPT='RECREATE'"]

