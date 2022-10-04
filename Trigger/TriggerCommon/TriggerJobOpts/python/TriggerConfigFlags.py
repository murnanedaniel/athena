# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

import os

from AthenaConfiguration.AthConfigFlags import AthConfigFlags
from AthenaConfiguration.Enums import FlagEnum, Format
from AthenaCommon.SystemOfUnits import GeV
from AthenaCommon.Logging import logging

log=logging.getLogger('TriggerConfigFlags')

class ROBPrefetching(FlagEnum):
    # Enable mapping step InputMaker outputs as ROBPrefetchingAlg inputs
    StepRoI = 'StepRoI'
    # Enable mapping chains' first step to pre-HLT prefetching rules based on initial RoIs
    InitialRoI = 'InitialRoI'

def createTriggerFlags(doTriggerRecoFlags):
    flags = AthConfigFlags()

    # enables L1 simulation
    flags.addFlag('Trigger.doLVL1', lambda prevFlags: prevFlags.Input.isMC)

    # Run HLT selection algorithms
    flags.addFlag('Trigger.doHLT', False)

    # changes decoding of L1 so that allways all configured chains are enabled, testing mode
    flags.addFlag("Trigger.HLTSeeding.forceEnableAllChains", False)

    # Enable Run-3 LVL1 muon decoding
    flags.addFlag('Trigger.enableL1MuonPhase1', lambda prevFlags: prevFlags.Trigger.EDMVersion >= 3 or prevFlags.Detector.EnableMM or prevFlags.Detector.EnablesTGC)

    # Enable Run-3 LVL1 calo simulation and/or decoding
    # TODO: when Fex decoders work reliably, change it to: "lambda prevFlags: (prevFlags.Trigger.EDMVersion >= 3 or prevFlags.GeoModel.Run >= LHCPeriod.Run3) and not prevFlags.Trigger.doHLT"
    flags.addFlag('Trigger.enableL1CaloPhase1', False)

    # Enable L1Topo simulation to write inputs to txt
    flags.addFlag('Trigger.enableL1TopoDump', False)

    # Enable L1Topo simulation to BW run
    flags.addFlag('Trigger.enableL1TopoBWSimulation', False)

    # Enable Run-2 L1Calo simulation and/or decoding (possible even if enablePhase1 is True)
    flags.addFlag('Trigger.enableL1CaloLegacy', True)

    # Enable emulation tool for NSW-TGC coincidence
    flags.addFlag('Trigger.L1MuonSim.EmulateNSW', False)

    # Enable NSW MM trigger
    flags.addFlag('Trigger.L1MuonSim.doMMTrigger', False)

    # Enable NSW sTGC pad trigger
    flags.addFlag('Trigger.L1MuonSim.doPadTrigger', False)

    # Enable NSW sTGC strip trigger
    flags.addFlag('Trigger.L1MuonSim.doStripTrigger', False)

    # Enable Storing NSW Debug Ntuple
    flags.addFlag('Trigger.L1MuonSim.WriteNSWDebugNtuple', False)

    # Enable Storing MM branches in NSW Debug Ntuple
    flags.addFlag('Trigger.L1MuonSim.WriteMMBranches', False)

    # Enable Storing NSW Debug Ntuple
    flags.addFlag('Trigger.L1MuonSim.WritesTGCBranches', False)

    # Enable the veto mode of the NSW-TGC coincidence
    flags.addFlag('Trigger.L1MuonSim.NSWVetoMode', False)

    # Enable TGC-RPC BIS78 coincidence
    flags.addFlag('Trigger.L1MuonSim.doBIS78', False)

    # Offline CondDB tag for RPC/TGC coincidence window in rerunLVL1 on data
    flags.addFlag('Trigger.L1MuonSim.CondDBOffline', 'OFLCOND-MC16-SDR-RUN2-04')

    # Enable Inner Detector
    flags.addFlag('Trigger.doID', True)

    # Enable muon system
    flags.addFlag('Trigger.doMuon', True)

    # Enable calorimeters
    flags.addFlag('Trigger.doCalo', True)

    # Enable additional validation histograms
    flags.addFlag('Trigger.doValidationMonitoring', False)

    # Checks the validity of each Decision Object produced by a HypoAlg, including all of its
    # parents all the way back to the HLTSeeding. Potentially CPU expensive.
    # also enables per step decison printouts
    flags.addFlag('Trigger.doRuntimeNaviVal', False)

    # Select ROBPrefetching options, out of ROBPrefetching enum, empty list means no prefetching
    flags.addFlag('Trigger.ROBPrefetchingOptions', [ROBPrefetching.InitialRoI, ROBPrefetching.StepRoI])

    # if 1, Run1 decoding version is set; if 2, Run2; if 3, Run 3
    def EDMVersion(flags):
        """Determine Trigger EDM version based on the input file."""
        _log = logging.getLogger('TriggerConfigFlags.EDMVersion')

        default_version = -1     # intentionally invalid default value, ATR-22856

        if flags.Input.Format is Format.BS:
            _log.debug("Input format is ByteStream")

            if not any(flags.Input.Files) and flags.Common.isOnline:
                _log.info("Online reconstruction, no input file. Return default EDMVersion=%d", default_version)
                return default_version
            try:
                from TrigEDMConfig.Utils import getEDMVersionFromBS
            except ImportError:
                log.error("Failed to import TrigEDMConfig, analysing ByteStream files is not possible in this release!")
                raise


            version = getEDMVersionFromBS(flags.Input.Files[0])

            return version if version is not None else default_version

        else:
            # POOL files: decide based on HLT output type present in the file
            _log.debug("Input format is POOL -- determine from input file collections")
            collections = flags.Input.Collections
            if "HLTResult_EF" in collections:
                _log.info("Determined EDMVersion to be 1, because HLTResult_EF found in POOL file")
                return 1
            elif "TrigNavigation" in collections:
                _log.info("Determined EDMVersion to be 2, because TrigNavigation found in POOL file")
                return 2
            elif any("HLTNav_Summary" in s for s in collections):
                _log.info("Determined EDMVersion to be 3, because HLTNav_Summary.* found in POOL file")
                return 3
            elif not flags.Input.Collections:
                # Special case for empty input files (can happen in merge jobs on the grid)
                # The resulting version doesn't really matter as there's nothing to be done, but we want a valid configuration
                _log.warning("All input files seem to be empty, cannot determine EDM version. Guessing EDMVersion=3")
                return 3

        _log.warning("Could not determine EDM version from the input file. Return default EDMVersion=%d",
                     default_version)
        return default_version

    flags.addFlag('Trigger.EDMVersion', lambda prevFlags: EDMVersion(prevFlags))
    flags.addFlag('Trigger.doEDMVersionConversion', False)

    # Flag to control the scheduling of online Run 3 trigger navigation compactification into a single collection (uses slimming framework).
    flags.addFlag('Trigger.doOnlineNavigationCompactification', True)

    # Flag to control the scheduling of offline Run 3 trigger navigation slimming in RAWtoESD, RAWtoAOD, AODtoDAOD or RAWtoALL transforms.
    flags.addFlag('Trigger.doNavigationSlimming', True)

    # True if we have at least one input file, it is a POOL file, it has a metadata store, and the store has xAOD trigger configuration data
    # in either the run-2 or run-3 formats.
    def __trigConfMeta(flags):
        from AthenaConfiguration.AutoConfigFlags import GetFileMD
        md = GetFileMD(flags.Input.Files)
        return ("metadata_items" in md and any(('TriggerMenu' in key) for key in md["metadata_items"].keys()))

    # Flag to sense if trigger configuration POOL metadata is available on the job's input
    flags.addFlag('Trigger.InputContainsConfigMetadata', lambda prevFlags: __trigConfMeta(prevFlags))

    # Enables collection and export of detailed monitoring data of the HLT execution
    flags.addFlag('Trigger.CostMonitoring.doCostMonitoring', False)
    flags.addFlag('Trigger.CostMonitoring.chain', 'HLT_noalg_CostMonDS_L1All')
    flags.addFlag('Trigger.CostMonitoring.outputCollection', 'HLT_TrigCostContainer')
    flags.addFlag('Trigger.CostMonitoring.monitorAllEvents', False)
    flags.addFlag('Trigger.CostMonitoring.monitorROBs', True)

    # enable L1Muon ByteStream conversion / simulation
    flags.addFlag('Trigger.L1.doMuon', True)

    # enable ByteStream conversion / simulation of MUCTPI Topo TOBs
    flags.addFlag('Trigger.L1.doMuonTopoInputs', lambda prevFlags: prevFlags.Trigger.doLVL1)

    # enable L1Calo ByteStream conversion / simulation
    flags.addFlag('Trigger.L1.doCalo', True)

    # enable L1Calo Input ([ejg]Towers) ByteStream conversion / simulation
    flags.addFlag('Trigger.L1.doCaloInputs', lambda prevFlags: prevFlags.Trigger.enableL1CaloPhase1 and not prevFlags.Trigger.doHLT)

    # enable L1Topo ByteStream conversion / simulation
    flags.addFlag('Trigger.L1.doTopo', True)

    # enable CTP ByteStream conversion / simulation
    flags.addFlag('Trigger.L1.doCTP', True)

    # partition name used to determine online vs offline BS result writing
    flags.addFlag('Trigger.Online.partitionName', os.getenv('TDAQ_PARTITION') or '')

    # shortcut to check if job is running in a partition (i.e. partition name is not empty)
    flags.addFlag('Trigger.Online.isPartition', lambda prevFlags: len(prevFlags.Trigger.Online.partitionName)>0)

    # write BS output file
    flags.addFlag('Trigger.writeBS', False)

    # Write transient BS before executing HLT algorithms (for running on MC RDO with clients which require BS inputs)
    flags.addFlag('Trigger.doTransientByteStream', lambda prevFlags: True if prevFlags.Input.Format is Format.POOL and prevFlags.Trigger.doCalo else False)

    # list of EDM objects to be written to AOD
    flags.addFlag('Trigger.AODEDMSet', lambda flags: 'AODSLIM' if flags.Input.isMC else 'AODFULL')

    # list of objects to be written to ESD
    flags.addFlag('Trigger.ESDEDMSet', 'ESD')

    # to allow stroing extra EDM items via preExec
    flags.addFlag('Trigger.ExtraEDMList', [])

    # tag to be used for condutions used by HLT code
    flags.addFlag('Trigger.OnlineCondTag', 'CONDBR2-HLTP-2018-03')

    # geometry version used by HLT online
    flags.addFlag('Trigger.OnlineGeoTag', 'ATLAS-R2-2016-01-00-01')

    def __availableRecoMetadata(flags):
        systems = ['L1','HLT']
        # Online reco without input files
        if not any(flags.Input.Files) and flags.Common.isOnline:
            return systems
        # Makes no sense when running HLT
        elif flags.Trigger.doHLT:
            raise RuntimeError('Trigger.availableRecoMetadata is ill-defined if Trigger.doHLT==True')
        # RAW: check if keys are in COOL
        elif flags.Input.Format is Format.BS:
            from TrigConfigSvc.TriggerConfigAccess import getKeysFromCool
            keys = getKeysFromCool(flags.Input.RunNumber[0], lbNr = 1)  # currently only checking first file
            return ( (['L1'] if 'L1PSK' in keys else []) +
                     (['HLT'] if 'HLTPSK' in keys else []) )
        # POOL: metadata (do not distinguish L1/HLT yet, see discussions on GitLab commit f83ae2bc)
        else:
            return systems if flags.Trigger.InputContainsConfigMetadata else []

    # list of enabled trigger sub-systems in reconstruction: ['L1,'HLT']
    flags.addFlag('Trigger.availableRecoMetadata', lambda flags: __availableRecoMetadata(flags))

    # Offline flag, determins if the HLT trigger decision and the HLT result should be decoded in reconstruction in jobs with Reco.EnableTrigger
    flags.addFlag("Trigger.DecodeHLT", True)

    # the configuration source
    # see https://twiki.cern.ch/twiki/bin/view/Atlas/TriggerConfigFlag
    flags.addFlag('Trigger.triggerConfig', lambda flags: 'INFILE' if flags.Trigger.InputContainsConfigMetadata else 'FILE')

    # name of the trigger menu
    flags.addFlag('Trigger.triggerMenuSetup', 'MC_pp_run3_v1_BulkMCProd_prescale')

    # modify the slection of chains that are run (default run all), see more in GenerateMenuMT_newJO
    flags.addFlag('Trigger.triggerMenuModifier', ['all'])

    # debug output from control flow generation
    flags.addFlag('Trigger.generateMenuDiagnostics', False)

    # disable Consistent Prescale Sets, for testing only, useful when using selectChains (ATR-24744)
    flags.addFlag('Trigger.disableCPS', False)

    # Switch whether end-of-event sequence running extra algorithms for accepted events should be added
    flags.addFlag('Trigger.endOfEventProcessing.Enabled', True)

    # trigger reconstruction
    # Protection against import of packages not in the analysis release
    # Signature and other trigger reco flags should be handled here
    if doTriggerRecoFlags:
        flags.join( createTriggerRecoFlags() )

    return flags

def createTriggerRecoFlags():
    flags = AthConfigFlags()

    # enables the correction for pileup in cell energy calibration (should it be moved to some place where other calo flags are defined?)
    flags.addFlag('Trigger.calo.doOffsetCorrection', True )

    def __egamma():
        from TriggerMenuMT.HLT.Egamma.TrigEgammaConfigFlags import createTrigEgammaConfigFlags
        return createTrigEgammaConfigFlags()
    flags.addFlagsCategory('Trigger.egamma', __egamma )

    # muon offline reco flags varaint for trigger
    def __muonSA():
        from MuonConfig.MuonConfigFlags import createMuonConfigFlags
        muonflags = createMuonConfigFlags()
        muonflags.Muon.useTGCPriorNextBC=True
        muonflags.Muon.MuonTrigger=True
        muonflags.Muon.SAMuonTrigger=True
        return muonflags

    def __muon():
        from MuonConfig.MuonConfigFlags import createMuonConfigFlags
        muonflags = createMuonConfigFlags()
        muonflags.Muon.useTGCPriorNextBC=True
        muonflags.Muon.MuonTrigger=True
        return muonflags

    def __muonCombined():
        from MuonCombinedConfig.MuonCombinedConfigFlags import createMuonCombinedConfigFlags
        muonflags = createMuonCombinedConfigFlags()
        muonflags.MuonCombined.doCaloTrkMuId = False
        muonflags.MuonCombined.doSiAssocForwardMuons = False
        muonflags.MuonCombined.doStatisticalCombination = False
        muonflags.MuonCombined.doMuGirl = False
        return muonflags

    flags.addFlagsCategory('Trigger.Offline.SA', __muonSA, prefix=True)
    flags.addFlagsCategory('Trigger.Offline', __muon, prefix=True)
    flags.addFlagsCategory('Trigger.Offline.Combined', __muonCombined, prefix=True)

    def __tau():
        from TrigTauRec.TrigTauConfigFlags import createTrigTauConfigFlags
        return createTrigTauConfigFlags()
    flags.addFlagsCategory('Trigger.Offline.Tau', __tau )
    #TODO come back and use systematically the same

    def __idTrk():
        from TrigInDetConfig.TrigTrackingPassFlags import createTrigTrackingPassFlags
        return createTrigTrackingPassFlags()
    flags.addFlagsCategory( 'Trigger.InDetTracking', __idTrk )

    # NB: Longer term it may be worth moving these into a PF set of config flags, but right now the only ones that exist do not seem to be used in the HLT.
    # When we use component accumulators for this in the HLT maybe we should revisit this
    # PFO-muon removal option for the full-scan hadronic signatures.
    # Options are:
    #   "None": Do no PFO-muon removal
    #   "Calo": Use the calo-tagging tools from the muon slice
    #   "Iso" : Use the mainly isolation-based selections based on the MET associator package
    flags.addFlag("Trigger.FSHad.PFOMuonRemoval", "Calo")

    # the minimum pT threshold to use for the muon removal
    flags.addFlag("Trigger.FSHad.PFOMuonRemovalMinPt", 10 * GeV)

    # Switch on MC20 EOverP maps for the jet slice
    flags.addFlag("Trigger.Jet.doMC20_EOverP", True)

    # enable fast b-tagging for all fully calibrated HLT PFlow jets
    flags.addFlag("Trigger.Jet.fastbtagPFlow", True)

    # enables or disables the addition of VR track jet reconstruction sequence
    flags.addFlag("Trigger.Jet.doVRJets", False)

    # use online-derived calibration for HLT PFlow jets
    flags.addFlag("Trigger.Jet.useTriggerCalib", False)

    def __httFlags():
        """Additional function delays import"""
        from TrigHTTConfTools.HTTConfigFlags import createHTTConfigFlags
        return createHTTConfigFlags()
    flags.addFlagsCategory("Trigger.HTT", __httFlags, prefix=True )

    return flags

if __name__ == "__main__":
    import unittest

    class Tests(unittest.TestCase):
        def test_recoFlags(self):
            """Check if offline reco flags can be added to trigger"""
            from AthenaConfiguration.AllConfigFlags import ConfigFlags as flags
            flags.Trigger.Offline.Tau.doTauRec=False
            flags.Tau.doTauRec=True
            self.assertEqual(flags.Trigger.Offline.Tau.doTauRec, False, " dependent flag setting does not work")
            self.assertEqual(flags.Tau.doTauRec, True, " dependent flag setting does not work")

            newflags = flags.cloneAndReplace('Tau', 'Trigger.Offline.Tau')

            self.assertEqual(flags.Tau.doTauRec, True, " dependent flag setting does not work")
            self.assertEqual(newflags.Tau.doTauRec, False, " dependent flag setting does not work")
            newflags.dump()

    unittest.main()
