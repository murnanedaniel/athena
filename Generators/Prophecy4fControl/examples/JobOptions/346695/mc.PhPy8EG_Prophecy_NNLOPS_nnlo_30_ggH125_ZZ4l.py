# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
#-------------------------------------------------------------------------------------------------
# Prophecy4f setup with input from PowHeg LHE, merge decays from Prophecy4f, output Pythia8 evtGen
#-------------------------------------------------------------------------------------------------

# need to save the value of inputGeneratorFile into powhegInputFile, since inputGeneratorFile is
# overwritten with an OTF file name in PowhegControl
powhegInputFile                          = runArgs.inputGeneratorFile    # for prod jobs

include('PowhegControl/PowhegControl_ggF_H_Common.py')


# Setup Powheg to use some variables to be the same in the Prophecy config. These
# should eventually be taken from the incoming LHE file

PowhegConfig.generate(create_run_card_only=True)

#--------------------------------------------------------------
# configure Prophecy
#--------------------------------------------------------------
include('Prophecy4fControl/ProphecyPowhegCommon.py')

# set name of LHE file - used below in merger

# print "runArgs", runArgs

#------------------------------------------------------------------------------------------------------
# configure and run ProphecyPowhegMerger to merge the prophecy4f LHE file with the incoming PowHeg file
#------------------------------------------------------------------------------------------------------

include('Prophecy4fControl/ProphecyPowhegMergeCommon.py')

#--------------------------------------------------------------
# Pythia8 showering
# note: Main31 is set in Pythia8_AZNLO_CTEQ6L1_EvtGen_Common.py
#--------------------------------------------------------------
# include('MC15JobOptions/Pythia8_AZNLO_CTEQ6L1_EvtGen_Common.py')
include('Pythia8_i/Pythia8_AZNLO_CTEQ6L1_EvtGen_Common.py')

#--------------------------------------------------------------
# Pythia8 main31 update
#--------------------------------------------------------------
# genSeq.Pythia8.UserModes += [ 'Main31:NFinal = 1' ]

if "UserHooks" in genSeq.Pythia8.__slots__.keys():

  genSeq.Pythia8.Commands  += ['Powheg:NFinal = 2']

else:

  genSeq.Pythia8.UserModes += [ 'Main31:NFinal = 2' ]

#--------------------------------------------------------------
# Higgs at Pythia8
#--------------------------------------------------------------
# genSeq.Pythia8.Commands += ["TimeShower:QEDshowerByOther = off",
#                             "TimeShower:QEDshowerByL = off",
#                             "TimeShower:QEDshowerByGamma = off",
#                             "TimeShower:QEDshowerByL = off",
#                             "SpaceShower:QEDshowerByQ = off",
#                             "SpaceShower:QEDshowerByL = off",
#                             "UncertaintyBands:doVariations=off"] 
genSeq.Pythia8.Commands += ["TimeShower:QEDshowerByL = off",
                            "TimeShower:QEDshowerByGamma = off" ]

# genSeq.Pythia8.Commands += ["TimeShower:QEDshowerByOther = off"]

#--------------------------------------------------------------
# EVGEN configuration
#--------------------------------------------------------------
evgenConfig.description = "POWHEG+PROPHECY+PYTHIA8+EVTGEN, H+jet production with NNLOPS and A14 tune, HZZ4l mh=125 GeV"
evgenConfig.keywords    = [ "Higgs", "SMHiggs", "ZZ" ]
evgenConfig.contact     = [ 'R.D.Schaffer@cern.ch' ]

evgenConfig.nEventsPerJob = 10000
