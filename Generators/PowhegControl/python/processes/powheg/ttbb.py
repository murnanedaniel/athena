# Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration

from AthenaCommon import Logging
from ..powheg_RES import PowhegRES
from ..external import ExternalMadSpin

## Get handle to Athena logging
logger = Logging.logging.getLogger("PowhegControl")


# Dictionary to convert the PowhegControl decay mode names to the appropriate
# decay mode numbers understood by Powheg
#
# The PowhegControl decay modes with MadSpin in their name use MadSpin to
# generate the top decays, the others use Powheg
_decay_mode_lookup = {
    "t t~ > all [MadSpin]" : "00000", # switch off decays in Powheg and let MadSpin handle them!
    "t t~ > all": "22222",
    "t t~ > b j j b~ j j": "00022",
    "t t~ > b l+ vl b~ l- vl~": "22200",
    "t t~ > b emu+ vemu b~ emu- vemu~": "22000",
    "t t~ > semileptonic": "11111",
    "t t~ > undecayed" : "00000",
}



class ttbb(PowhegRES):
    """
    Powheg interface for top-antitop-bottom-antibottom production.
    The top quarks may be left undecayed, or their decays (including spin correlations)
    can be generated by Powheg or using MadSpin.

    Reference for this process: https://arxiv.org/abs/1802.00426

    @author Stefan Richter  <stefan.richter@cern.ch>
    """

    def __init__(self, base_directory, **kwargs):
        """! Constructor: all process options are set here.

        @param base_directory: path to PowhegBox code.
        @param kwargs          dictionary of arguments from Generate_tf.
        """
        super(ttbb, self).__init__(base_directory, "ttbb", **kwargs)

        # This is a hacky fix that's needed at the moment...
        self.manually_set_openloops_paths()

        # This process' integration needs Athena to be set to run at least two parallel processes
        # Advise the user about this here:
        if self.cores < 2:
            info_message = """
                Due to an apparent bug in PowhegBox, the *INTEGRATION* for this Powheg process (ttbb)
                requires running at least two (computer) processes in parallel. Please configure Athena
                to do so, e.g. by setting the environment variable ATHENA_PROC_NUMBER to 2 or a higher
                number. In Bash, do e.g.: 'export ATHENA_PROC_NUMBER=4'.
                """
            logger.info(info_message)

        # Add algorithms to the sequence
        self.add_algorithm(ExternalMadSpin(process="generate p p > t t~ b b~ [QCD]"))

        # Add parameter validation functions
        self.validation_functions.append("validate_decays")

        # List of allowed decay modes
        # (The sorting of the list is just to increase readability when it's printed)
        self.allowed_decay_modes = sorted(_decay_mode_lookup.keys())

        # Add all keywords for this process, overriding defaults if required
        self.add_keyword("alphas_from_lhapdf", 1)
        self.add_keyword("bornsuppfact")
        self.add_keyword("bornzerodamp")
        self.add_keyword("bornzerodampcut", 5.)
        self.add_keyword("btlscalect", 1)
        self.add_keyword("btlscalereal", 1)
        self.add_keyword("clobberlhe")
        self.add_keyword("compress_lhe", 1)
        self.add_keyword("compress_upb", 1)
        self.add_keyword("facscfact", self.default_scales[0])
        self.add_keyword("fastbtlbound", 1)
        self.add_keyword("foldcsi", 5)
        self.add_keyword("foldphi", 1)
        self.add_keyword("foldy", 5)
        self.add_keyword("for_reweighting")
        self.add_keyword("fullrwgt")
        self.add_keyword("fullrwgtmode")
        self.add_keyword("icsimax", 1)
        self.add_keyword("ih1")
        self.add_keyword("ih2")
        self.add_keyword("itmx1", 2)
        self.add_keyword("itmx2", 3)
        self.add_keyword("iymax", 1)
        self.add_keyword("lhans1", self.default_PDFs)
        self.add_keyword("lhans2", self.default_PDFs)
        self.add_keyword("manyseeds", 1)
        self.add_keyword("maxseeds", 1000)
        self.add_keyword("ncall1", 40000)
        self.add_keyword("ncall2", 40000)
        self.add_keyword("ncall2rm", 80000)
        self.add_keyword("nubound", 10000)
        self.add_keyword("olverbose", 0)
        self.add_keyword("parallelstage")
        self.add_keyword("renscfact", self.default_scales[1])
        self.add_keyword("runningscales", 2)
        self.add_keyword("rwl_add")
        self.add_keyword("rwl_file")
        self.add_keyword("rwl_format_rwgt")
        self.add_keyword("rwl_group_events")
        self.add_keyword("semileptonic", hidden=True)
        self.add_keyword("storemintupb", 1)
        self.add_keyword("tdec/bmass")
        self.add_keyword("tdec/cmass")
        self.add_keyword("tdec/dmass")
        self.add_keyword("tdec/elbranching")
        self.add_keyword("elbranching")
        self.add_keyword("tdec/emass")
        self.add_keyword("tdec/mumass")
        self.add_keyword("tdec/sin2cabibbo")
        self.add_keyword("tdec/smass")
        self.add_keyword("tdec/taumass")
        self.add_keyword("tdec/twidth")
        self.add_keyword("tdec/umass")
        self.add_keyword("tdec/wmass")
        self.add_keyword("tdec/wwidth")
        self.add_keyword("topdecaymode", "t t~ > all [MadSpin]", name="decay_mode")
        self.add_keyword("use-old-grid", 1)
        self.add_keyword("use-old-ubound", 1)
        self.add_keyword("withdamp", 1)
        self.add_keyword("xgriditeration")
        self.add_keyword("xupbound", 2)
        self.add_keyword("hdamp", -1)
        self.add_keyword("dynhdamp", 1)
        self.add_keyword("dynhdampPF", 0.5)


    def validate_decays(self):
        """
        Validate semileptonic and topdecaymode keywords and translate them from ATLAS input to Powheg input
        """
        self.expose() # convenience call to simplify syntax
        if self.decay_mode not in self.allowed_decay_modes:
            error_message = "Decay mode '{given}' not recognised, valid choices are: '{choices}'!".format(given=self.decay_mode, choices="', '".join(self.allowed_decay_modes))
            logger.warning(error_message)
            raise ValueError(error_message)

        # Check if MadSpin decays are requested.
        # Accordingly, MadSpin will run or not run.
        if "MadSpin" in self.decay_mode:
            self.externals["MadSpin"].parameters_by_keyword("powheg_top_decays_enabled")[0].value = False

        self.parameters_by_keyword("topdecaymode")[0].value = _decay_mode_lookup[self.decay_mode]
        if "semileptonic" in self.decay_mode:
            # Parameter semileptonic must be set to 1 to actually get semileptonic decays, because the topdecaymode=11111 also allows fully hadronic decays (with one up and one charm quark)
            self.parameters_by_keyword("semileptonic")[0].value = 1

