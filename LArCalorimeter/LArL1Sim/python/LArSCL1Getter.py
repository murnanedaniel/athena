# Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration

# LArSCL1 creation from LArHits with LArSCL1Maker algorithm

from AthenaCommon.Logging import logging
from RecExConfig.Configured import Configured
from AthenaCommon.AlgSequence import AlgSequence
topSequence = AlgSequence()
import traceback

from IOVDbSvc.CondDB import conddb

def addLArFlatFolder (db, obj, calg, folder_base='/LAR/ElecCalibFlat/',qual=''):
    from AthenaCommon.AlgSequence import AthSequencer
    condSequence = AthSequencer("AthCondSeq")

    folder = folder_base + obj
    if not conddb.folderRequested(folder):
      conddb.addFolder(db, folder + qual,
                     className = 'CondAttrListCollection')
      condSequence += calg (ReadKey=folder, WriteKey='LAr'+obj+'SC')
    return

class LArSCL1Getter ( Configured )  :

        
    def configure(self):
        mlog = logging.getLogger( 'LArSCL1Getter::configure:' )
        mlog.info ('entering')        

        # get handle to upstream object
        try:
            from LArL1Sim.LArSCL1Getter import LArSCL1Getter
            theLArSCL1Getter=LArSCL1Getter()
        except Exception:
            mlog.error("could not get handle to LArSCL1Getter Quit")
            traceback.print_exc()
            return False

        if not theLArSCL1Getter.usable():
            if not self.ignoreConfigError():
                mlog.error("LArSCL1Getter unusable. Quit.")
                return False
            else:
                mlog.error("LArSCL1Getter unusable. Continue nevertheless")
                
        # Instantiation of the C++ algorithm
        try:        
            from LArL1Sim.LArL1SimConf import LArSCL1Maker                
        except Exception:
            mlog.error("could not import LArL1Sim.LArSCL1Maker")
            traceback.print_exc()
            return False

        from LArCabling.LArCablingAccess import LArOnOffIdMappingSC
        from LArRecUtils.LArRecUtilsConf import LArFlatConditionsAlg_LArNoiseSC_ as LArNoiseSCCondAlg
        from LArRecUtils.LArRecUtilsConf import LArFlatConditionsAlg_LArPedestalSC_ as LArPedestalSCFlatCondAlg
        from LArRecUtils.LArRecUtilsConf import LArFlatConditionsAlg_LArShapeSC_ as LArShapeSCCondAlg
        from LArRecUtils.LArRecUtilsConf import LArFlatConditionsAlg_LArfSamplSC_ as LArfSamplSCCondAlg
        

        LArOnOffIdMappingSC()
        addLArFlatFolder ('LAR_OFL', 'Shape', LArShapeSCCondAlg,'/LAR/ElecCalibMCSC/')
        addLArFlatFolder ('LAR_OFL', 'Pedestal', LArPedestalSCFlatCondAlg,'/LAR/ElecCalibMCSC/')
        addLArFlatFolder ('LAR_OFL', 'Noise', LArNoiseSCCondAlg,'/LAR/ElecCalibMCSC/')
        addLArFlatFolder ('LAR_OFL', 'fSampl', LArfSamplSCCondAlg,'/LAR/ElecCalibMCSC/')
        from LArRecUtils.LArAutoCorrNoiseSCCondAlgDefault import LArAutoCorrNoiseSCCondAlgDefault
        LArAutoCorrNoiseSCCondAlgDefault()
        from LArRecUtils.LArADC2MeVSCCondAlgDefault import LArADC2MeVSCCondAlgDefault
        LArADC2MeVSCCondAlgDefault()
        theLArSCL1Maker=LArSCL1Maker()
        from LArROD.LArRODFlags import larRODFlags
        theLArSCL1Maker.NSamples = larRODFlags.nSamples() + 2  # For consistency with LArAutoCorrNoiseSC - see ATLASSIM-5483
        from Digitization.DigitizationFlags import digitizationFlags
        if digitizationFlags.PileUpPresampling and 'LegacyOverlay' not in digitizationFlags.experimentalDigi():
            from OverlayCommonAlgs.OverlayFlags import overlayFlags
            theLArSCL1Maker.SCL1ContainerName = overlayFlags.bkgPrefix() + "LArDigitSCL2"
        else:
            theLArSCL1Maker.SCL1ContainerName = "LArDigitSCL2"

        self._LArSCL1Maker = theLArSCL1Maker

        from AthenaCommon.AlgSequence import AlgSequence
        topSequence = AlgSequence()

        # check if LArdigitization is run before. If yes, uses hit map from detector store produces from lardigitization
        from AthenaCommon.DetFlags import DetFlags
        if DetFlags.digitize.LAr_on():
            mlog.info("Using hit map from LArHitEMapMaker algoritm")
        else:
            mlog.info("digitmaker1 not found in topSequence, using own map in LArSCL1Maker")
            return False
        

        mlog.info(" now adding to topSequence")        
        topSequence += theLArSCL1Maker

        
        return True

    def LArSCL1Maker(self):
        return self._LArSCL1Maker
   
    def outputKey(cls):
       return cls._output[cls._outputType]

    def outputType(cls):
       return cls._outputType

    def outputTypeKey(self):
       return str(self.outputType()+"#"+self.outputKey())


