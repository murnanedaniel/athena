# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration

import glob
import os
import shutil
import subprocess
from AthenaCommon import Logging
from ...decorators import timed
from ...utility import LHE, ProcessManager, SingleProcessThread

## Get handle to Athena logging
logger = Logging.logging.getLogger("PowhegControl")


@timed("MadSpin post-processing")
def MadSpin(process, powheg_LHE_output):
    """! Move output to correctly named file.

    @param process            External MadSpin process.
    @param powheg_LHE_output  Name of LHE file produced by PowhegBox.

    @author James Robinson <james.robinson@cern.ch>
    """
    process.expose()
    __construct_inputs(powheg_LHE_output, process)
    __run_executable(process.executable)
    __prepare_outputs(powheg_LHE_output)


def __construct_inputs(input_LHE_events, process):
    """! Construct MadSpin runcard.

    @param input_LHE_events  Input LHE file name.
    @param process           MadSpin process.

    @author James Robinson <james.robinson@cern.ch>
    """
    # Find insertion point for MadSpin header
    logger.info("Constructing MadSpin runcard header")
    pre_header  = "\n".join([elem for elem in [LHE.opening_tag(input_LHE_events), LHE.comment_block(input_LHE_events)] if elem])
    header_contents = LHE.header_block(input_LHE_events).replace("<header>", "").replace("</header>", "").strip("\n")
    post_header  = LHE.init_block(input_LHE_events)
    postamble  = LHE.postamble(input_LHE_events)

    # Write events to LHE file with MadSpin information
    with open("madspin_LHE_input.lhe", "wb") as f_madspin_LHE:
        f_madspin_LHE.write("{}\n".format(pre_header))
        f_madspin_LHE.write("<header>\n")
        f_madspin_LHE.write("<mgversion>\n")
        f_madspin_LHE.write("{}\n".format(os.environ["MADPATH"].split("/")[-1].split("v")[-1].split("_p")[0].replace("_", ".")))
        f_madspin_LHE.write("</mgversion>\n")
        f_madspin_LHE.write("<mg5proccard>\n")
        f_madspin_LHE.write("set group_subprocesses Auto\n")
        f_madspin_LHE.write("set ignore_six_quark_processes False\n")
        f_madspin_LHE.write("set loop_optimized_output True\n")
        f_madspin_LHE.write("set gauge unitary\n")
        f_madspin_LHE.write("set complex_mass_scheme False\n")
        f_madspin_LHE.write("import model {}\n".format(process.MadSpin_model))
        if int(process.MadSpin_nFlavours) == 4:
            f_madspin_LHE.write("define p = g u c d s u~ c~ d~ s~\n")
            f_madspin_LHE.write("define j = g u c d s u~ c~ d~ s~\n")
        elif int(process.MadSpin_nFlavours) == 5:
            f_madspin_LHE.write("define p = g u c d s b u~ c~ d~ s~ b~\n")
            f_madspin_LHE.write("define j = g u c d s b u~ c~ d~ s~ b~\n")
        else:
            raise ValueError("'MadSpin_nFlavours' must be set to 4 or 5")
        if process.MadSpin_taus_are_leptons:
            f_madspin_LHE.write("define l+ = e+ mu+ ta+\n")
            f_madspin_LHE.write("define l- = e- mu- ta-\n")
        else:
            f_madspin_LHE.write("define l+ = e+ mu+\n")
            f_madspin_LHE.write("define l- = e- mu-\n")
        f_madspin_LHE.write("{}\n".format(process.MadSpin_process))
        f_madspin_LHE.write("output tchan\n")
        f_madspin_LHE.write("</mg5proccard>\n")
        f_madspin_LHE.write("<mgruncard>\n")
        f_madspin_LHE.write("#0.01 = req_acc_FO ! needed for determining LO/NLO - see AGENE-1459\n")
        f_madspin_LHE.write("{} = nevents\n".format(LHE.event_counter(input_LHE_events)))
        f_madspin_LHE.write("1   =  lpp1     ! beam 1 type (0 = no PDF)\n")
        f_madspin_LHE.write("1   =  lpp2     ! beam 2 type (0 = no PDF)\n")
        f_madspin_LHE.write("{} =  ebeam1   ! beam 1 energy in GeV\n".format(process.beam_energy))
        f_madspin_LHE.write("{} =  ebeam2   ! beam 2 energy in GeV\n".format(process.beam_energy))
        f_madspin_LHE.write("{} =  bwcutoff\n".format(process.bwcutoff))
        f_madspin_LHE.write("</mgruncard>\n")
        f_madspin_LHE.write("<slha>\n")
        f_madspin_LHE.write("######################################################################\n")
        f_madspin_LHE.write("## PARAM_CARD AUTOMATICALY GENERATED BY MG5                       ####\n")
        f_madspin_LHE.write("######################################################################\n")
        f_madspin_LHE.write("###################################\n")
        f_madspin_LHE.write("## INFORMATION FOR LOOP\n")
        f_madspin_LHE.write("###################################\n")
        f_madspin_LHE.write("BLOCK LOOP #\n")
        f_madspin_LHE.write("      1 {:e} #   mu_r\n".format(process.mass_Z))
        f_madspin_LHE.write("###################################\n")
        f_madspin_LHE.write("## INFORMATION FOR MASS\n")
        f_madspin_LHE.write("###################################\n")
        f_madspin_LHE.write("BLOCK MASS #\n")
        f_madspin_LHE.write("      1 0.000000e+00 #   d : 0.0\n")
        f_madspin_LHE.write("      2 0.000000e+00 #   u : 0.0\n")
        f_madspin_LHE.write("      3 0.000000e+00 #   s : 0.0\n")
        f_madspin_LHE.write("      4 0.000000e+00 #   c : 0.0\n")
        f_madspin_LHE.write("      5 {:e} #   mb\n".format(process.mass_b))
        f_madspin_LHE.write("      6 {:e} #   mt\n".format(process.mass_t))
        f_madspin_LHE.write("      11 0.000000e+00 #   e- : 0.0\n")
        f_madspin_LHE.write("      12 0.000000e+00 #   ve : 0.0\n")
        f_madspin_LHE.write("      13 0.000000e+00 #   mu- : 0.0\n")
        f_madspin_LHE.write("      14 0.000000e+00 #   vm : 0.0\n")
        f_madspin_LHE.write("      15 {:e} #   mta\n".format(process.mass_tau))
        f_madspin_LHE.write("      16 0.000000e+00 #   vt : 0.0\n")
        f_madspin_LHE.write("      21 0.000000e+00 #   g : 0.0\n")
        f_madspin_LHE.write("      22 0.000000e+00 #   a : 0.0\n")
        f_madspin_LHE.write("      23 {:e} #   mz\n".format(process.mass_Z))
        f_madspin_LHE.write("      24 {:e} #   w+\n".format(process.mass_W))
        f_madspin_LHE.write("      25 {:e} #   mh\n".format(process.mass_H))
        f_madspin_LHE.write("      82 0.000000e+00 #   gh : 0.0\n")
        f_madspin_LHE.write("###################################\n")
        f_madspin_LHE.write("## INFORMATION FOR SMINPUTS\n")
        f_madspin_LHE.write("###################################\n")
        f_madspin_LHE.write("BLOCK SMINPUTS #\n")
        f_madspin_LHE.write("      1 {:e} #   aewm1\n".format(process.alphaem_inv))
        f_madspin_LHE.write("      2 {:e} #   gf\n".format(process.G_F))
        f_madspin_LHE.write("      3 {:e} #   as\n".format(process.alphaqcd))
        f_madspin_LHE.write("###################################\n")
        f_madspin_LHE.write("## INFORMATION FOR YUKAWA\n")
        f_madspin_LHE.write("###################################\n")
        f_madspin_LHE.write("BLOCK YUKAWA #\n")
        f_madspin_LHE.write("      5 {:e} #   ymb\n".format(process.mass_b))
        f_madspin_LHE.write("      6 {:e} #   ymt\n".format(process.mass_t))
        f_madspin_LHE.write("      15 {:e} #   ymtau\n".format(process.mass_tau))
        f_madspin_LHE.write("###################################\n")
        f_madspin_LHE.write("## INFORMATION FOR QNUMBERS 82\n")
        f_madspin_LHE.write("###################################\n")
        f_madspin_LHE.write("BLOCK QNUMBERS 82 #   gh\n")
        f_madspin_LHE.write("      1 0 #   3 times electric charge\n")
        f_madspin_LHE.write("      2 1 #   number of spin states (2s+1)\n")
        f_madspin_LHE.write("      3 8 #   colour rep (1: singlet, 3: triplet, 8: octet)\n")
        f_madspin_LHE.write("      4 1 #   particle/antiparticle distinction (0=own anti)\n")
        f_madspin_LHE.write("#\n")
        f_madspin_LHE.write("#*************************\n")
        f_madspin_LHE.write("#      Decay widths      *\n")
        f_madspin_LHE.write("#*************************\n")
        f_madspin_LHE.write("#\n")
        f_madspin_LHE.write("#      PDG        Width\n")
        f_madspin_LHE.write("DECAY  1   0.000000e+00\n")
        f_madspin_LHE.write("#\n")
        f_madspin_LHE.write("#      PDG        Width\n")
        f_madspin_LHE.write("DECAY  2   0.000000e+00\n")
        f_madspin_LHE.write("#\n")
        f_madspin_LHE.write("#      PDG        Width\n")
        f_madspin_LHE.write("DECAY  3   0.000000e+00\n")
        f_madspin_LHE.write("#\n")
        f_madspin_LHE.write("#      PDG        Width\n")
        f_madspin_LHE.write("DECAY  4   0.000000e+00\n")
        f_madspin_LHE.write("#\n")
        f_madspin_LHE.write("#      PDG        Width\n")
        f_madspin_LHE.write("DECAY  5   0.000000e+00\n")
        f_madspin_LHE.write("#\n")
        # Scale t BRs down so that sum is 1.0
        BR_t_to_any = process.BR_t_to_Wb + process.BR_t_to_Ws + process.BR_t_to_Wd
        f_madspin_LHE.write("#      PDG        Width\n")
        f_madspin_LHE.write("DECAY  6   {:e}\n".format(process.width_t))
        f_madspin_LHE.write("#  BR             NDA  ID1    ID2   ...\n")
        f_madspin_LHE.write("   {:e}   2    5  24 # 1.32\n".format(process.BR_t_to_Wb / BR_t_to_any))
        f_madspin_LHE.write("   {:e}   2    3  24 # 1.32\n".format(process.BR_t_to_Ws / BR_t_to_any))
        f_madspin_LHE.write("   {:e}   2    1  24 # 1.32\n".format(process.BR_t_to_Wd / BR_t_to_any))
        f_madspin_LHE.write("#      PDG        Width\n")
        f_madspin_LHE.write("DECAY  -6   {:e}\n".format(process.width_t))
        f_madspin_LHE.write("#  BR             NDA  ID1    ID2   ...\n")
        f_madspin_LHE.write("   {:e}   2   -5 -24 # 1.32\n".format(process.BR_t_to_Wb / BR_t_to_any))
        f_madspin_LHE.write("   {:e}   2   -3 -24 # 1.32\n".format(process.BR_t_to_Ws / BR_t_to_any))
        f_madspin_LHE.write("   {:e}   2   -1 -24 # 1.32\n".format(process.BR_t_to_Wd / BR_t_to_any))
        f_madspin_LHE.write("#\n")
        f_madspin_LHE.write("#      PDG        Width\n")
        f_madspin_LHE.write("DECAY  11   0.000000e+00\n")
        f_madspin_LHE.write("#\n")
        f_madspin_LHE.write("#      PDG        Width\n")
        f_madspin_LHE.write("DECAY  12   0.000000e+00\n")
        f_madspin_LHE.write("#\n")
        f_madspin_LHE.write("#      PDG        Width\n")
        f_madspin_LHE.write("DECAY  13   0.000000e+00\n")
        f_madspin_LHE.write("#\n")
        f_madspin_LHE.write("#      PDG        Width\n")
        f_madspin_LHE.write("DECAY  14   0.000000e+00\n")
        f_madspin_LHE.write("#\n")
        f_madspin_LHE.write("#      PDG        Width\n")
        f_madspin_LHE.write("DECAY  15   0.000000e+00\n")
        f_madspin_LHE.write("#\n")
        f_madspin_LHE.write("#      PDG        Width\n")
        f_madspin_LHE.write("DECAY  16   0.000000e+00\n")
        f_madspin_LHE.write("#\n")
        f_madspin_LHE.write("#      PDG        Width\n")
        f_madspin_LHE.write("DECAY  21   0.000000e+00\n")
        f_madspin_LHE.write("#\n")
        f_madspin_LHE.write("#      PDG        Width\n")
        f_madspin_LHE.write("DECAY  22   0.000000e+00\n")
        f_madspin_LHE.write("#\n")
        f_madspin_LHE.write("#      PDG        Width\n")
        f_madspin_LHE.write("DECAY  23   {:e}\n".format(process.width_Z))
        f_madspin_LHE.write("#\n")
        # Scale W BRs down so that sum is 1.0
        BR_W_to_any = process.BR_W_to_hadrons + process.BR_W_to_leptons
        f_madspin_LHE.write("#      PDG        Width\n")
        f_madspin_LHE.write("DECAY  -24   {:e}\n".format(process.width_W))
        f_madspin_LHE.write("#  BR             NDA  ID1    ID2   ...\n")
        f_madspin_LHE.write("      {:e}   2    1  -2\n".format((process.BR_W_to_hadrons / 2.) / BR_W_to_any))
        f_madspin_LHE.write("      {:e}   2    3  -4\n".format((process.BR_W_to_hadrons / 2.) / BR_W_to_any))
        f_madspin_LHE.write("      {:e}   2   11 -12\n".format((process.BR_W_to_leptons / 3.) / BR_W_to_any))
        f_madspin_LHE.write("      {:e}   2   13 -14\n".format((process.BR_W_to_leptons / 3.) / BR_W_to_any))
        f_madspin_LHE.write("      {:e}   2   15 -16\n".format((process.BR_W_to_leptons / 3.) / BR_W_to_any))
        f_madspin_LHE.write("#\n")
        f_madspin_LHE.write("#      PDG        Width\n")
        f_madspin_LHE.write("DECAY  24   {}\n".format(process.width_W))
        f_madspin_LHE.write("#  BR             NDA  ID1    ID2   ...\n")
        f_madspin_LHE.write("      {:e}   2   -1   2\n".format((process.BR_W_to_hadrons / 2.) / BR_W_to_any))
        f_madspin_LHE.write("      {:e}   2   -3   4\n".format((process.BR_W_to_hadrons / 2.) / BR_W_to_any))
        f_madspin_LHE.write("      {:e}   2  -11  12\n".format((process.BR_W_to_leptons / 3.) / BR_W_to_any))
        f_madspin_LHE.write("      {:e}   2  -13  14\n".format((process.BR_W_to_leptons / 3.) / BR_W_to_any))
        f_madspin_LHE.write("      {:e}   2  -15  16\n".format((process.BR_W_to_leptons / 3.) / BR_W_to_any))
        f_madspin_LHE.write("#\n")
        f_madspin_LHE.write("#      PDG        Width\n")
        f_madspin_LHE.write("DECAY  25   {:e}\n".format(process.width_H))
        f_madspin_LHE.write("#\n")
        f_madspin_LHE.write("#      PDG        Width\n")
        f_madspin_LHE.write("DECAY  82   0.000000e+00\n")
        f_madspin_LHE.write("</slha>\n")
        f_madspin_LHE.write("{}\n".format(header_contents))
        f_madspin_LHE.write("</header>\n")
        f_madspin_LHE.write("{}\n".format(post_header))
        for event in LHE.event_iterator(input_LHE_events):
            f_madspin_LHE.write(event)
        f_madspin_LHE.write(postamble)

    # Rename LHE files
    shutil.move(input_LHE_events, "{}.undecayed".format(input_LHE_events))
    shutil.move("madspin_LHE_input.lhe", input_LHE_events)

    # Write MadSpin runcard
    with open("madspin_runcard.txt", "wb") as f_madspin_runcard:
        f_madspin_runcard.write("import {}\n".format(input_LHE_events))
        f_madspin_runcard.write("set spinmode {}\n".format(process.MadSpin_mode))
        for decay in process.MadSpin_decays:
            f_madspin_runcard.write("{0}\n".format(decay))
        f_madspin_runcard.write("launch\n")
        f_madspin_runcard.write("quit\n")


@timed("MadSpin executable")
def __run_executable(executable):
    """! Run MadSpin executable.

    @author James Robinson <james.robinson@cern.ch>
    """
    if not os.path.isfile(executable):
        raise OSError("MadSpin executable {} not found!".format(executable))
    logger.info("MadSpin executable: {}".format(executable))
    with open("madspin_runcard.txt", "rb") as runcard_input:
        processes = [SingleProcessThread([executable], stdin=runcard_input, ignore_output=["INFO:", "MadSpin>"],
                                         error_output=["Command \"launch\" interrupted with error:", "MadSpinError"],
                                         info_output=["generating the production square matrix element"])]
        manager = ProcessManager(processes)
        while manager.monitor():
            pass


def __prepare_outputs(input_LHE_events):
    """! Prepare MadSpin output.

    @author James Robinson <james.robinson@cern.ch>
    """
    logger.info("Preparing MadSpin output")

    # Unzip MadSpin events
    subprocess.call("gunzip pwgevents_decayed.lhe.gz 2> /dev/null", shell=True)
    shutil.move("pwgevents_decayed.lhe", "{}.decayed".format(input_LHE_events))

    # Get appropriate header sections from Powheg and MadSpin LHE
    powheg_LHE, madspin_LHE = "{}.undecayed".format(input_LHE_events), "{}.decayed".format(input_LHE_events)
    pre_header = "\n".join([elem for elem in [LHE.opening_tag(powheg_LHE), LHE.comment_block(powheg_LHE)] if elem])
    header = LHE.header_block(madspin_LHE)
    post_header = LHE.init_block(powheg_LHE)
    postamble = LHE.postamble(powheg_LHE)

    # Write events to LHE file with MadSpin information
    with open("pwgevents_with_MadSpin.lhe", "wb") as f_combined_LHE:
        f_combined_LHE.write("{}\n".format(pre_header))
        f_combined_LHE.write("{}\n".format(header))
        f_combined_LHE.write("{}\n".format(post_header))
        for event in LHE.event_iterator("{}.decayed".format(input_LHE_events)):
            f_combined_LHE.write(event)
        f_combined_LHE.write(postamble)

    # Rename output file
    if os.path.isfile(input_LHE_events):
        shutil.move(input_LHE_events, "{}.undecayed_ready_for_madspin".format(input_LHE_events))
    shutil.move("pwgevents_with_MadSpin.lhe", input_LHE_events)

    for LHE_tarball in list(set(glob.glob("*lhe*.gz") + glob.glob("*events*.gz"))):
        logger.info("Cleaning up MadSpin tarball: {}".format(LHE_tarball))
        os.remove(LHE_tarball)
