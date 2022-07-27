# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

from AthenaCommon.SystemOfUnits import TeV
from AthenaConfiguration.AthConfigFlags import AthConfigFlags
from AthenaConfiguration.AutoConfigFlags import GetFileMD, getInitialTimeStampsFromRunNumbers, getRunToTimestampDict, getSpecialConfigurationMetadata
from AthenaConfiguration.Enums import BeamType, Format, ProductionStep, BunchStructureSource
from PyUtils.moduleExists import moduleExists


def _addFlagsCategory (acf, name, generator, modName = None):
    if moduleExists (modName):
        return acf.addFlagsCategory (name, generator)
    return None


def _createCfgFlags():

    acf=AthConfigFlags()

    #Flags steering the job execution:
    from AthenaCommon.Constants import INFO
    acf.addFlag('Exec.OutputLevel',INFO) #Global Output Level
    acf.addFlag('Exec.MaxEvents',-1)
    acf.addFlag('Exec.SkipEvents',0)
    acf.addFlag('Exec.DebugStage','')

    acf.addFlag('ExecutorSplitting.TotalSteps', 0)
    acf.addFlag('ExecutorSplitting.Step', -1)
    acf.addFlag('ExecutorSplitting.TotalEvents', -1)

    #Flags describing the input data
    acf.addFlag('Input.Files', ["_ATHENA_GENERIC_INPUTFILE_NAME_",]) # former global.InputFiles
    acf.addFlag('Input.SecondaryFiles', []) # secondary input files for DoubleEventSelector
    acf.addFlag('Input.isMC', lambda prevFlags : "IS_SIMULATION" in GetFileMD(prevFlags.Input.Files).get("eventTypes", [])) # former global.isMC
    acf.addFlag('Input.OverrideRunNumber', False )
    acf.addFlag("Input.ConditionsRunNumber", -1) # Override the HITS file Run Number with one from a data run (TODO merge with Input.RunNumber)
    acf.addFlag('Input.RunNumber', lambda prevFlags : list(GetFileMD(prevFlags.Input.Files).get("runNumbers", []))) # former global.RunNumber
    acf.addFlag('Input.MCChannelNumber', lambda prevFlags : GetFileMD(prevFlags.Input.Files).get("mc_channel_number", 0))
    acf.addFlag('Input.LumiBlockNumber', lambda prevFlags : list(GetFileMD(prevFlags.Input.Files).get("lumiBlockNumbers", []))) # former global.RunNumber
    acf.addFlag('Input.TimeStamp', lambda prevFlags : getInitialTimeStampsFromRunNumbers(prevFlags.Input.RunNumber) if prevFlags.Input.OverrideRunNumber else [])
    # Configure EvtIdModifierSvc with a list of dictionaries of the form:
    # {'run': 152166, 'lb': 202, 'starttstamp': 1269948352889940910, 'evts': 1, 'mu': 0.005}
    acf.addFlag("Input.RunAndLumiOverrideList", [])
    # Job number
    acf.addFlag("Input.JobNumber", 1)
    acf.addFlag('Input.FailOnUnknownCollections', False)

    acf.addFlag('Input.ProjectName', lambda prevFlags : GetFileMD(prevFlags.Input.Files).get("project_name", "data17_13TeV")) # former global.ProjectName
    acf.addFlag('Input.TriggerStream', lambda prevFlags : GetFileMD(prevFlags.Input.Files).get("stream", "") if prevFlags.Input.Format == Format.BS
                                                          else GetFileMD(prevFlags.Input.Files).get("triggerStreamOfFile", "")) # former global.TriggerStream
    acf.addFlag('Input.Format', lambda prevFlags : Format.BS if GetFileMD(prevFlags.Input.Files).get("file_type", "BS") == "BS" else Format.POOL, enum=Format) # former global.InputFormat
    acf.addFlag('Input.ProcessingTags', lambda prevFlags : GetFileMD(prevFlags.Input.Files).get("processingTags", []) ) # list of names of streams written to this file
    acf.addFlag('Input.SpecialConfiguration', lambda prevFlags : getSpecialConfigurationMetadata(prevFlags.Input.Files, prevFlags.Input.SecondaryFiles))  # special Configuration options read from input file metadata

    def _inputCollections(inputFile):
        rawCollections = [type_key[1] for type_key in GetFileMD(inputFile).get("itemList", [])]
        collections = [col for col in rawCollections if not col.endswith('Aux.')]
        return collections

    def _typedInputCollections(inputFile):
        collections = ['%s#%s' % type_key for type_key in GetFileMD(inputFile).get("itemList", [])]
        return collections

    acf.addFlag('Input.Collections', lambda prevFlags : _inputCollections(prevFlags.Input.Files) )
    acf.addFlag('Input.SecondaryCollections', lambda prevFlags : _inputCollections(prevFlags.Input.SecondaryFiles) )
    acf.addFlag('Input.TypedCollections', lambda prevFlags : _typedInputCollections(prevFlags.Input.Files) )
    acf.addFlag('Input.SecondaryTypedCollections', lambda prevFlags : _typedInputCollections(prevFlags.Input.SecondaryFiles) )

    acf.addFlag('Concurrency.NumProcs', 0)
    acf.addFlag('Concurrency.NumThreads', 0 )
    acf.addFlag('Concurrency.NumConcurrentEvents', lambda prevFlags : prevFlags.Concurrency.NumThreads)
    acf.addFlag('Concurrency.DebugWorkers', False )

    acf.addFlag('Scheduler.CheckDependencies', True)
    acf.addFlag('Scheduler.ShowDataDeps', True)
    acf.addFlag('Scheduler.ShowDataFlow', True)
    acf.addFlag('Scheduler.ShowControlFlow', True)
    acf.addFlag('Scheduler.EnableVerboseViews', True)
    acf.addFlag('Scheduler.AutoLoadUnmetDependencies', True)

    acf.addFlag('MP.WorkerTopDir', 'athenaMP_workers')
    acf.addFlag('MP.OutputReportFile', 'AthenaMPOutputs')
    acf.addFlag('MP.Strategy', 'SharedQueue')
    acf.addFlag('MP.CollectSubprocessLogs', False)
    acf.addFlag('MP.PollingInterval', 100)
    acf.addFlag('MP.EventsBeforeFork', 0)
    acf.addFlag('MP.EventRangeChannel', 'EventService_EventRanges')
    acf.addFlag('MP.EvtRangeScattererCaching', False)
    acf.addFlag('MP.MemSamplingInterval', 0)
    """ Size of event chunks in the shared queue
        if chunk_size==-1, chunk size is set to auto_flush for files compressed with LZMA
        if chunk_size==-2, chunk size is set to auto_flush for files compressed with LZMA or ZLIB
        if chunk_size==-3, chunk size is set to auto_flush for files compressed with LZMA, ZLIB, or LZ4
        if chunk_size<=-4, chunk size is set to auto_flush
    """
    acf.addFlag('MP.ChunkSize', -1)
    acf.addFlag('MP.ReadEventOrders', False)
    acf.addFlag('MP.EventOrdersFile', 'athenamp_eventorders.txt')
    acf.addFlag('MP.UseSharedReader', False)
    acf.addFlag('MP.UseSharedWriter', False)
    acf.addFlag('MP.UseParallelCompression', True)

    acf.addFlag('Common.MsgSourceLength',50) #Length of the source-field in the format str of MessageSvc
    acf.addFlag('Common.isOnline', False ) #  Job runs in an online environment (access only to resources available at P1) # former global.isOnline
    acf.addFlag('Common.useOnlineLumi', lambda prevFlags : prevFlags.Common.isOnline ) #  Use online version of luminosity. ??? Should just use isOnline?
    acf.addFlag('Common.isOverlay', lambda prevFlags: (prevFlags.Common.ProductionStep == ProductionStep.Overlay or
                                                       (prevFlags.Common.ProductionStep == ProductionStep.FastChain and
                                                        prevFlags.Overlay.FastChain)))  # Enable Overlay
    acf.addFlag('Common.doExpressProcessing', False)
    acf.addFlag('Common.ProductionStep', ProductionStep.Default, enum=ProductionStep)

    def _checkProject():
        import os
        if "AthSimulation_DIR" in os.environ:
            return "AthSimulation"
        if "AthGeneration_DIR" in os.environ:
            return "AthGeneration"
        if "AthAnalysis_DIR" in os.environ:
            return "AthAnalysis"
        if "AthDerivation_DIR" in os.environ:
            return "AthDerivation"
        #TODO expand this method.
        return "Athena"
    acf.addFlag('Common.Project', _checkProject())

    # replace global.Beam*
    acf.addFlag('Beam.BunchSpacing', 25) # former global.BunchSpacing
    acf.addFlag('Beam.Type', lambda prevFlags : BeamType(GetFileMD(prevFlags.Input.Files).get('beam_type', 'collisions')), enum=BeamType)# former global.BeamType
    acf.addFlag("Beam.NumberOfCollisions", lambda prevFlags : 2. if prevFlags.Beam.Type is BeamType.Collisions else 0.) # former global.NumberOfCollisions
    acf.addFlag('Beam.Energy', lambda prevFlags : GetFileMD(prevFlags.Input.Files).get('beam_energy',7*TeV)) # former global.BeamEnergy
    acf.addFlag('Beam.estimatedLuminosity', lambda prevFlags : ( 1E33*(prevFlags.Beam.NumberOfCollisions)/2.3 ) *\
        (25./prevFlags.Beam.BunchSpacing)) # former flobal.estimatedLuminosity
    acf.addFlag('Beam.BunchStructureSource', lambda prevFlags: BunchStructureSource.MC if prevFlags.Input.isMC else BunchStructureSource.TrigConf)



    acf.addFlag('Output.EVNTFileName', '')
    acf.addFlag('Output.EVNT_TRFileName', '')
    acf.addFlag('Output.HITSFileName', '')
    acf.addFlag('Output.RDOFileName',  '')
    acf.addFlag('Output.RDO_SGNLFileName', '')
    acf.addFlag('Output.ESDFileName',  '')
    acf.addFlag('Output.AODFileName',  '')
    acf.addFlag('Output.HISTFileName', '')


    acf.addFlag('Output.doWriteRDO', lambda prevFlags: bool(prevFlags.Output.RDOFileName)) # write out RDO file
    acf.addFlag('Output.doWriteRDO_SGNL', lambda prevFlags: bool(prevFlags.Output.RDO_SGNLFileName)) # write out RDO_SGNL file
    acf.addFlag('Output.doWriteESD', lambda prevFlags: bool(prevFlags.Output.ESDFileName)) # write out ESD file
    acf.addFlag('Output.doWriteAOD', lambda prevFlags: bool(prevFlags.Output.AODFileName)) # write out AOD file
    acf.addFlag('Output.doWriteBS',  False) # write out RDO ByteStream file
    acf.addFlag('Output.doWriteDAOD',  False) # write out at least one DAOD file

    # Might move this elsewhere in the future.
    # Some flags from https://gitlab.cern.ch/atlas/athena/blob/master/Tracking/TrkDetDescr/TrkDetDescrSvc/python/TrkDetDescrJobProperties.py
    # (many, e.g. those that set properties of one tool are not needed)
    acf.addFlag('TrackingGeometry.MagneticFileMode', 6)
    acf.addFlag('TrackingGeometry.MaterialSource', 'COOL') # Can be COOL, Input or None

#Detector Flags:
    def __detector():
        from AthenaConfiguration.DetectorConfigFlags import createDetectorConfigFlags
        return createDetectorConfigFlags()
    acf.addFlagsCategory( "Detector", __detector )

#Simulation Flags:
    def __simulation():
        from SimulationConfig.SimConfigFlags import createSimConfigFlags
        return createSimConfigFlags()
    _addFlagsCategory (acf, "Sim", __simulation, 'SimulationConfig' )

#Digitization Flags:
    def __digitization():
        from Digitization.DigitizationConfigFlags import createDigitizationCfgFlags
        return createDigitizationCfgFlags()
    _addFlagsCategory(acf, "Digitization", __digitization, 'Digitization' )

#Overlay Flags:
    def __overlay():
        from OverlayConfiguration.OverlayConfigFlags import createOverlayConfigFlags
        return createOverlayConfigFlags()
    _addFlagsCategory(acf, "Overlay", __overlay, 'OverlayConfiguration' )

#Geo Model Flags:
    def __geomodel():
        from AthenaConfiguration.GeoModelConfigFlags import createGeoModelConfigFlags
        return createGeoModelConfigFlags()
    acf.addFlagsCategory( "GeoModel", __geomodel )

#Reco Flags:
    def __reco():
        from RecJobTransforms.RecoConfigFlags import createRecoConfigFlags
        return createRecoConfigFlags()
    _addFlagsCategory(acf, "Reco", __reco, 'RecJobTransforms')

#IOVDbSvc Flags:
    from IOVDbSvc.IOVDbAutoCfgFlags import getLastGlobalTag, getDatabaseInstanceDefault
    acf.addFlag("IOVDb.GlobalTag", getLastGlobalTag) # Retrieve last global tag used from metadata
    acf.addFlag("IOVDb.DatabaseInstance",getDatabaseInstanceDefault)
    # Run dependent simulation
    # map from runNumber to timestamp; migrated from RunDMCFlags.py
    acf.addFlag("IOVDb.RunToTimestampDict", lambda prevFlags: getRunToTimestampDict())
    acf.addFlag("IOVDb.DBConnection", lambda prevFlags : "sqlite://;schema=mycool.db;dbname=" + prevFlags.IOVDb.DatabaseInstance)

#PoolSvc Flags:
    acf.addFlag("PoolSvc.MaxFilesOpen", lambda prevFlags : 2 if prevFlags.MP.UseSharedReader else 0)


    def __bfield():
        from MagFieldConfig.BFieldConfigFlags import createBFieldConfigFlags
        return createBFieldConfigFlags()
    _addFlagsCategory(acf, "BField", __bfield, 'MagFieldConfig')

    def __lar():
        from LArConfiguration.LArConfigFlags import createLArConfigFlags
        return createLArConfigFlags()
    _addFlagsCategory(acf, "LAr", __lar, 'LArConfiguration' )

    def __tile():
        from TileConfiguration.TileConfigFlags import createTileConfigFlags
        return createTileConfigFlags()
    _addFlagsCategory(acf, 'Tile', __tile, 'TileConfiguration' )


    def __calo():
        from CaloRec.CaloConfigFlags import createCaloConfigFlags
        return createCaloConfigFlags()
    _addFlagsCategory(acf, 'Calo', __calo, 'CaloRec' )

#Random engine Flags:
    acf.addFlag("Random.Engine", "dSFMT") # Random service used in {"dSFMT", "Ranlux64", "Ranecu"}

    def __trigger():
        from TriggerJobOpts.TriggerConfigFlags import createTriggerFlags
        return createTriggerFlags(acf.Common.Project!='AthAnalysis')
    _addFlagsCategory(acf, "Trigger", __trigger, 'TriggerJobOpts' )

    def __indet():
        from InDetConfig.InDetConfigFlags import createInDetConfigFlags
        return createInDetConfigFlags()
    _addFlagsCategory(acf, "InDet", __indet, 'InDetConfig' )

    def __itk():
        from InDetConfig.ITkConfigFlags import createITkConfigFlags
        return createITkConfigFlags()
    _addFlagsCategory(acf, "ITk", __itk, 'InDetConfig' )

    def __acts():
        from ActsInterop.ActsConfigFlags import createActsConfigFlags
        return createActsConfigFlags()
    _addFlagsCategory(acf, "Acts", __acts, 'ActsInterop')

    def __hgtd():
        from HGTD_Config.HGTD_ConfigFlags import createHGTD_ConfigFlags
        return createHGTD_ConfigFlags()
    _addFlagsCategory(acf, "HGTD", __hgtd, 'HGTD_Config' )

    def __muon():
        from MuonConfig.MuonConfigFlags import createMuonConfigFlags
        return createMuonConfigFlags()
    _addFlagsCategory(acf, "Muon", __muon, 'MuonConfig' )

    def __muoncombined():
        from MuonCombinedConfig.MuonCombinedConfigFlags import createMuonCombinedConfigFlags
        return createMuonCombinedConfigFlags()
    _addFlagsCategory(acf, "MuonCombined", __muoncombined, 'MuonCombinedConfig' )

    def __egamma():
        from egammaConfig.egammaConfigFlags import createEgammaConfigFlags
        return createEgammaConfigFlags()
    _addFlagsCategory(acf, "Egamma", __egamma, 'egammaConfig' )

    def __met():
        from METReconstruction.METConfigFlags import createMETConfigFlags
        return createMETConfigFlags()
    _addFlagsCategory(acf,"MET",__met, 'METReconstruction')

    def __jet():
        from JetRecConfig.JetConfigFlags import createJetConfigFlags
        return createJetConfigFlags()
    _addFlagsCategory(acf,"Jet",__jet, 'JetRecConfig')

    def __tau():
        from tauRec.TauConfigFlags import createTauConfigFlags
        return createTauConfigFlags()
    _addFlagsCategory(acf, "Tau",__tau, 'tauRec')

    def __pflow():
        from eflowRec.PFConfigFlags import createPFConfigFlags
        return createPFConfigFlags()
    _addFlagsCategory(acf,"PF",__pflow, 'eflowRec')

    def __btagging():
        from JetTagConfig.BTaggingConfigFlags import createBTaggingConfigFlags
        return createBTaggingConfigFlags()
    _addFlagsCategory(acf,"BTagging",__btagging, 'JetTagConfig')

    def __hi():
        from HIRecConfig.HIRecConfigFlags import createHIRecConfigFlags
        return createHIRecConfigFlags()
    _addFlagsCategory(acf, "HeavyIon", __hi, "HIRecConfig")

    def __dq():
        from AthenaMonitoring.DQConfigFlags import createDQConfigFlags
        dqf = createDQConfigFlags()
        return dqf
    _addFlagsCategory(acf, "DQ", __dq, 'AthenaMonitoring' )

    def __perfmon():
        from PerfMonComps.PerfMonConfigFlags import createPerfMonConfigFlags
        return createPerfMonConfigFlags()
    _addFlagsCategory(acf, "PerfMon", __perfmon, 'PerfMonComps')

    return acf


ConfigFlags=_createCfgFlags()

if __name__=="__main__":
    import sys
    if len(sys.argv)>1:
        ConfigFlags.Input.Files = sys.argv[1:]
    else:
        ConfigFlags.Input.Files = [ "/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/CommonInputs/data16_13TeV.00311321.physics_Main.recon.AOD.r9264/AOD.11038520._000001.pool.root.1",]

    ConfigFlags.loadAllDynamicFlags()
    ConfigFlags.initAll()
    ConfigFlags.dump()
