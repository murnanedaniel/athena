# Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration

#--------------------------------------------------------------
# Powheg Zj setup starting from ATLAS defaults
#--------------------------------------------------------------
include("PowhegControl/PowhegControl_Wj_MiNNLO_Common.py")
PowhegConfig.generate()

#--------------------------------------------------------------
# Pythia8 showering with the A14 NNPDF2.3 tune
#--------------------------------------------------------------
include("Pythia8_i/Pythia8_A14_NNPDF23LO_EvtGen_Common.py")
include("Pythia8_i/Pythia8_Powheg.py")

#--------------------------------------------------------------
# EVGEN configuration
#--------------------------------------------------------------
evgenConfig.description = "POWHEG+Pythia8 W+j MiNNLO production with A14 NNPDF2.3 tune"
evgenConfig.keywords    = [ "SM", "W", "1jet" ]
evgenConfig.contact     = [ "simone.amoroso@cern.ch" ]
