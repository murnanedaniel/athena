# Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration

from AthenaCommon import Logging
from ..powheg_V2 import PowhegV2
from ..external import ExternalMadSpin
import os

## Get handle to Athena logging
logger = Logging.logging.getLogger("PowhegControl")


# Dictionary to convert the PowhegControl decay mode names to the appropriate
# decay mode numbers understood by Powheg
#
# The PowhegControl decay modes with MadSpin in their name use MadSpin to
# generate the top decays, the others use Powheg
_decay_mode_lookup = {
    "t t~ > all": "22222",
    "t t~ > b j j b~ j j": "00022",
    "t t~ > b l+ vl b~ l- vl~": "22200",
    "t t~ > b emu+ vemu b~ emu- vemu~": "22000",
    "t t~ > semileptonic": "11111",
    "t t~ > undecayed" : "00000",
    "t t~ > all [MadSpin]" : "00000", # switch off decays in Powheg and let MadSpin handle them!
}


class ttj(PowhegV2):
    """! Default Powheg configuration for top pair production plus one jet.

    Create a configurable object with all applicable Powheg options.

    @author James Robinson  <james.robinson@cern.ch>
    """

    def __init__(self, base_directory, **kwargs):
        """! Constructor: all process options are set here.

        @param base_directory: path to PowhegBox code.
        @param kwargs          dictionary of arguments from Generate_tf.
        """
        super(ttj, self).__init__(base_directory, "ttbarj", **kwargs)

        # hack in place to help powheg executable find all dynamic libraries
        logger.warning("Applying manual, hard-coded fixes for library paths")
        OLPath = os.path.dirname(self.executable)+"/obj-gnu"
        os.environ['OpenLoopsPath'] = OLPath
        logger.info("OpenLoopsPath defined as = {0}".format(os.getenv('OpenLoopsPath')))

        # Add algorithms to the sequence
        self.add_algorithm(ExternalMadSpin(process="generate p p > t t~ j [QCD]"))

        # Add parameter validation functions
        self.validation_functions.append("validate_decays")

        ## List of allowed decay modes
        # (The sorting of the list is just to increase readability when it's printed)
        self.allowed_decay_modes = sorted(_decay_mode_lookup.keys())

        # Add all keywords for this process, overriding defaults if required
        self.add_keyword("alphaem")
        self.add_keyword("bcut")
        self.add_keyword("bornktmin", 5.0)
        self.add_keyword("bornonly")
        self.add_keyword("bornsuppfact", 100.0)
        self.add_keyword("bornzerodamp")
        self.add_keyword("bottomthr")
        self.add_keyword("bottomthrpdf")
        self.add_keyword("charmthr")
        self.add_keyword("charmthrpdf")
        self.add_keyword("CKM_Vcb")
        self.add_keyword("CKM_Vcd")
        self.add_keyword("CKM_Vcs")
        self.add_keyword("CKM_Vtb")
        self.add_keyword("CKM_Vtd")
        self.add_keyword("CKM_Vts")
        self.add_keyword("CKM_Vub")
        self.add_keyword("CKM_Vud")
        self.add_keyword("CKM_Vus")
        self.add_keyword("clobberlhe")
        self.add_keyword('collremnsamp"')
        self.add_keyword("colltest")
        self.add_keyword("compute_rwgt")
        self.add_keyword("elbranching")
        self.add_keyword("facscfact", self.default_scales[0])
        self.add_keyword("fastbtlbound")
        self.add_keyword("ffltest")
        self.add_keyword("flg_debug")
        self.add_keyword("foldcsi", 2)
        self.add_keyword("foldphi", 2)
        self.add_keyword("foldy", 2)
        self.add_keyword("fullrwgt")
        self.add_keyword("fullrwgtmode")
        self.add_keyword("ggproc")
        self.add_keyword("gqbproc")
        self.add_keyword("hdamp")
        self.add_keyword("hfact")
        self.add_keyword("icsimax", 3)
        self.add_keyword("ih1")
        self.add_keyword("ih2")
        self.add_keyword("itmx1")
        self.add_keyword("itmx2", 4)
        self.add_keyword("iupperfsr")
        self.add_keyword("iupperisr")
        self.add_keyword("iymax", 3)
        self.add_keyword("jacsing")
        self.add_keyword("largecorrfact")
        self.add_keyword("lhans1", self.default_PDFs)
        self.add_keyword("lhans2", self.default_PDFs)
        self.add_keyword("lhfm/bmass")
        self.add_keyword("lhfm/cmass")
        self.add_keyword("lhfm/emass")
        self.add_keyword("lhfm/mumass")
        self.add_keyword("lhfm/taumass")
        self.add_keyword("manyseeds")
        self.add_keyword("masswindow_high")
        self.add_keyword("masswindow_low")
        self.add_keyword("ncall1", 51000)
        self.add_keyword("ncall1rm")
        self.add_keyword("ncall2", 91000)
        self.add_keyword("ncall2rm")
        self.add_keyword("nubound", 11000)
        self.add_keyword("par_2gsupp")
        self.add_keyword("par_diexp")
        self.add_keyword("par_dijexp")
        self.add_keyword("par_fsrtinycsi")
        self.add_keyword("par_fsrtinyy")
        self.add_keyword("par_isrtinycsi")
        self.add_keyword("par_isrtinyy")
        self.add_keyword("parallelstage")
        self.add_keyword("pdfreweight")
        self.add_keyword("psgen")
        self.add_keyword("ptmin_jet", 20.0)
        self.add_keyword("ptsqmin")
        self.add_keyword("ptsupp")
        self.add_keyword("qgproc")
        self.add_keyword("qqbproc")
        self.add_keyword("qqpproc")
        self.add_keyword("qqproc")
        self.add_keyword("R_jet")
        self.add_keyword("radregion")
        self.add_keyword("rand1")
        self.add_keyword("rand2")
        self.add_keyword("renscfact", self.default_scales[1])
        self.add_keyword("runningscale")
        self.add_keyword("rwl_add")
        self.add_keyword("rwl_file")
        self.add_keyword("rwl_format_rwgt")
        self.add_keyword("rwl_group_events")
        self.add_keyword("semileptonic", hidden=True)
        self.add_keyword("smartsig")
        self.add_keyword("softtest")
        self.add_keyword("sthw2")
        self.add_keyword("storeinfo_rwgt")
        self.add_keyword("testplots")
        self.add_keyword("testsuda")
        self.add_keyword("topdecaymode", "t t~ > all", name="decay_mode")
        self.add_keyword("topmass")
        self.add_keyword("topwidth")
        self.add_keyword("ubsigmadetails")
        self.add_keyword("use-old-grid")
        self.add_keyword("use-old-ubound")
        self.add_keyword("use-OLP-interface")
        self.add_keyword("withdamp")
        self.add_keyword("withnegweights")
        self.add_keyword("withsubtr")
        self.add_keyword("Wmass")
        self.add_keyword("Wwidth")
        self.add_keyword("xupbound", 6)
        self.add_keyword("zerowidth")
        self.add_keyword("Zmass")
        self.add_keyword("Zwidth")

    def validate_decays(self):
        """! Validate semileptonic and topdecaymode keywords."""
        self.expose()  # convenience call to simplify syntax
        if self.decay_mode not in self.allowed_decay_modes:
            error_message = "Decay mode '{given}' not recognised, valid choices are: '{choices}'!".format(given=self.decay_mode, choices="', '".join(self.allowed_decay_modes))
            logger.warning(error_message)
            raise ValueError(error_message)

        # Check if MadSpin decays are requested.
        # Accordingly, MadSpin will run or not run.
        if "MadSpin" in self.decay_mode:
            self.externals["MadSpin"].parameters_by_keyword("powheg_top_decays_enabled")[0].value = False
            self.externals["MadSpin"].parameters_by_keyword("MadSpin_model")[0].value = "loop_sm-no_b_mass"
            self.externals["MadSpin"].parameters_by_keyword("MadSpin_nFlavours")[0].value = 5

        self.parameters_by_keyword("topdecaymode")[0].value = _decay_mode_lookup[self.decay_mode]
        if self.decay_mode == "semileptonic":
            # Parameter semileptonic must be set to 1 to actually get semileptonic decays, because the topdecaymode=11111 also allows fully hadronic decays (with one up and one charm quark)
            self.parameters_by_keyword("semileptonic")[0].value = 1
