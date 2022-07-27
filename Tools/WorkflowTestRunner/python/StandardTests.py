# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
from typing import List

from .Checks import AODContentCheck, AODDigestCheck, FrozenTier0PolicyCheck
from .Inputs import input_EVNT, input_EVNT_AF3, input_HITS, \
    input_HITS_MC_overlay, input_RDO_BKG, \
    input_HITS_data_overlay, input_BS_SKIM, \
    input_HITS_minbias_low, input_HITS_minbias_high, input_HITS_neutrino
from .Test import TestSetup, WorkflowRun, WorkflowTest, WorkflowType


class QTest(WorkflowTest):
    """General workflow q-test."""

    def __init__(self, ID: str, run: WorkflowRun, type: WorkflowType, steps: List[str], setup: TestSetup, extra_args: str = "") -> None:
        if "maxEvents" not in extra_args:
            if type == WorkflowType.MCPileUpReco or run == WorkflowRun.Run4:
                extra_args += " --maxEvents 5"
            else:
                extra_args += " --maxEvents 20"

        if type == WorkflowType.MCPileUpReco:
            if "inputHITSFile" not in extra_args:
                extra_args += f" --inputHITSFile {input_HITS[run]}"
            if "inputRDO_BKGFile" not in extra_args:
                extra_args += " --inputRDO_BKGFile ../run_d*/myRDO.pool.root"

        threads = 1
        threads_argument = '--multithreaded'
        if setup.custom_threads is not None:
            threads = setup.custom_threads
            if threads <= 0:
                threads_argument = ''

        self.command = \
            (f"ATHENA_CORE_NUMBER={threads} Reco_tf.py {threads_argument} --AMIConfig {ID}"
             f" --imf False {extra_args}")

        self.output_checks = []
        # TODO: disable RDO comparison for now
        # if type == WorkflowType.MCReco:
        #     self.output_checks.append(FrozenTier0PolicyCheck(setup, "RDO", 10))
        self.output_checks.append(FrozenTier0PolicyCheck(setup, "ESD", 10))
        self.output_checks.append(FrozenTier0PolicyCheck(setup, "AOD", 30))

        self.digest_checks = []
        if "--CA" not in extra_args:
            self.digest_checks.append(AODContentCheck(setup))
            self.digest_checks.append(AODDigestCheck(setup))

        super().__init__(ID, run, type, steps, setup)


class SimulationTest(WorkflowTest):
    """Simulation workflow test."""

    def __init__(self, ID: str, run: WorkflowRun, type: WorkflowType, steps: List[str], setup: TestSetup, extra_args: str = "") -> None:
        if "maxEvents" not in extra_args:
            extra_args += " --maxEvents 20"

        if type == WorkflowType.AF3:
            input_file = input_EVNT_AF3[run]
        else:
            input_file = input_EVNT[run]

        self.command = \
            (f"Sim_tf.py --AMIConfig {ID}"
             f" --inputEVNTFile {input_file} --outputHITSFile myHITS.pool.root"
             f" --imf False {extra_args}")

            #  " --postExec 'ServiceMgr.AuditorSvc.FPEAuditor.NStacktracesOnFPE=500'"

        self.output_checks = [
            FrozenTier0PolicyCheck(setup, "HITS", 10)
        ]

        super().__init__(ID, run, type, steps, setup)


class OverlayTest(WorkflowTest):
    """MC overlay workflow test."""

    def __init__(self, ID: str, run: WorkflowRun, type: WorkflowType, steps: List[str], setup: TestSetup, extra_args: str = "") -> None:
        if "maxEvents" not in extra_args:
            extra_args += " --maxEvents 10"

        self.command = \
            (f"Overlay_tf.py --AMIConfig {ID}"
             f" --inputHITSFile {input_HITS_MC_overlay[run]} --inputRDO_BKGFile {input_RDO_BKG[run]} --outputRDOFile myRDO.pool.root"
             f" --imf False --athenaopts=\"--pmon=sdmonfp\" {extra_args}")

        # skip performance checks for now due to CA
        self.skip_performance_checks = True

        self.output_checks = [
            FrozenTier0PolicyCheck(setup, "RDO", 10)
        ]

        super().__init__(ID, run, type, steps, setup)


class DataOverlayTest(WorkflowTest):
    """Data overlay workflow test."""

    def __init__(self, ID: str, run: WorkflowRun, type: WorkflowType, steps: List[str], setup: TestSetup, extra_args: str = "") -> None:
        if "maxEvents" not in extra_args:
            extra_args += " --maxEvents 10"

        self.command = \
            (f"Overlay_tf.py --AMIConfig {ID}"
             f" --inputHITSFile {input_HITS_data_overlay[run]} --inputBS_SKIMFile {input_BS_SKIM[run]} --outputRDOFile myRDO.pool.root"
             " --triggerConfig 'Overlay=NONE'"  # disable trigger for now
             f" --imf False --athenaopts=\"--pmon=sdmonfp\" {extra_args}")

        self.output_checks = [
            FrozenTier0PolicyCheck(setup, "RDO", 10)
        ]

        super().__init__(ID, run, type, steps, setup)


class PileUpTest(WorkflowTest):
    """Digitization with pile-up workflow test."""

    def __init__(self, ID: str, run: WorkflowRun, type: WorkflowType, steps: List[str], setup: TestSetup, extra_args: str = "") -> None:
        if "maxEvents" not in extra_args:
            extra_args += " --maxEvents 5"

        self.command = \
            (f"Digi_tf.py --AMIConfig {ID} --jobNumber 1 --digiSeedOffset1 1 --digiSeedOffset2 1"
             f" --inputHITSFile {input_HITS_neutrino[run]} --inputHighPtMinbiasHitsFile {input_HITS_minbias_high[run]} --inputLowPtMinbiasHitsFile {input_HITS_minbias_low[run]} --outputRDOFile myRDO.pool.root"
             " --postExec 'FPEAuditor.NStacktracesOnFPE=500'"
             f" --imf False --athenaopts=\"--pmon=sdmonfp\" {extra_args}")

        self.output_checks = [
            FrozenTier0PolicyCheck(setup, "RDO", 5)
        ]

        super().__init__(ID, run, type, steps, setup)


class DerivationTest(WorkflowTest):
    """Derivations test."""

    def __init__(self, ID: str, run: WorkflowRun, type: WorkflowType, steps: List[str], setup: TestSetup, extra_args: str = "") -> None:
        if "maxEvents" not in extra_args:
            extra_args += " --maxEvents 10"

        threads = 2
        if setup.custom_threads is not None:
            threads = setup.custom_threads

        self.command = \
            (f"ATHENA_CORE_NUMBER={threads} Derivation_tf.py --AMIConfig {ID}"
             " --outputDAODFile myOutput.pool.root"
             " --formats PHYS"
             f" --imf False {extra_args}")

        # skip performance checks for now due to CA
        self.skip_performance_checks = True

        if threads == 0:
            # does not work with shared writer
            self.output_checks = [
                FrozenTier0PolicyCheck(setup, "DAOD_PHYS", 10)
            ]

        super().__init__(ID, run, type, steps, setup)
