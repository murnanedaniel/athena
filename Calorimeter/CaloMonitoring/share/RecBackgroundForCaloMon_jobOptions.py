from RecExConfig.RecFlags import rec
from AthenaCommon.DetFlags import DetFlags
if rec.doInDet() and rec.doMuon() and rec.doCalo() and \
   DetFlags.detdescr.Muon_on() and DetFlags.detdescr.Calo_on() and DetFlags.detdescr.ID_on() :
   from RecBackgroundAlgs.RecBackgroundAlgsConf import BeamBackgroundFiller
   print "create BeamBackgroundFiller"
   BeamBackgroundFiller=BeamBackgroundFiller()
   print "set cutThetaCsc = 0"
   BeamBackgroundFiller.cutThetaCsc = 0
   print "set cutThetaMdtI = 90"
   BeamBackgroundFiller.cutThetaMdtI = 90
   print "set muonSegmentContainerKey = ConvertedMBoySegments "
   BeamBackgroundFiller.muonSegmentContainerKey = "ConvertedMBoySegments"
   print "add to top sequence"
   topSequence+=BeamBackgroundFiller

