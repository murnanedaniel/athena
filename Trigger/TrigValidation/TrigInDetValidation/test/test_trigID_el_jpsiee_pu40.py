#!/usr/bin/env python
# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

# art-description: art job for el_Jpsiee_pu40
# art-type: grid
# art-include: master/Athena
# art-include: 22.0/Athena
# art-input: mc15_13TeV.129190.Pythia8_AU2CTEQ6L1_ppToJpsie3e3.recon.RDO.e3802_s2608_s2183_r7042
# art-input-nfiles: 16
# art-athena-mt: 4
# art-html: https://idtrigger-val.web.cern.ch/idtrigger-val/TIDAWeb/TIDAart/?jobdir=
# art-output: *.txt
# art-output: *.log
# art-output: log.*
# art-output: *.out
# art-output: *.err
# art-output: *.log.tar.gz
# art-output: *.new
# art-output: *.json
# art-output: *.root
# art-output: *.check*
# art-output: HLT*
# art-output: times*
# art-output: cost-perCall
# art-output: cost-perEvent
# art-output: cost-perCall-chain
# art-output: cost-perEvent-chain
# art-output: *.dat 


Slices  = ['electron']
Events  = 8000
Threads = 8 
Slots   = 8
Input   = 'Jpsiee_pu40'    # defined in TrigValTools/share/TrigValInputs.json
GridFiles=True

preinclude_file = 'all:TrigInDetValidation/TIDV_cond_fix.py' #conditions fix for ATR-23982. In future find a more recent RDO  

Jobs = [ ( "Truth",       " TIDAdata-run3.dat                    -o data-hists.root -p 11" ),
         ( "Offline",     " TIDAdata-run3-offline.dat -r Offline -o data-hists-offline.root" ) ]

Comp = [ ( "L2ele",              "L2electron",      "data-hists.root",         " -c TIDAhisto-panel.dat  -d HLTL2-plots " ),
         ( "L2ele-lowpt",        "L2electronLowpt", "data-hists.root",         " -c TIDAhisto-panel.dat  -d HLTL2-plots-lowpt " ),
         ( "L2eleoffline",       "L2electron",      "data-hists-offline.root", " -c TIDAhisto-panel.dat  -d HLTL2-plots-offline " ),
         ( "L2eleoffline-lowpt", "L2electronLowpt", "data-hists-offline.root", " -c TIDAhisto-panel.dat  -d HLTL2-plots-lowpt-offline " ),
         ( "EFele",              "EFelectron",      "data-hists.root",         " -c TIDAhisto-panel.dat  -d HLTEF-plots " ),
         ( "EFele-lowpt",        "EFelectronLowpt", "data-hists.root",         " -c TIDAhisto-panel.dat  -d HLTEF-plots-lowpt " ),
         ( "EFeleoffline",       "EFelectron",      "data-hists-offline.root", " -c TIDAhisto-panel.dat  -d HLTEF-plots-offline " ),
         ( "EFeleoffline-lowpt", "EFelectronLowpt", "data-hists-offline.root", " -c TIDAhisto-panel.dat  -d HLTEF-plots-lowpt-offline " ) ]


from AthenaCommon.Include import include 
include("TrigInDetValidation/TrigInDetValidation_Base.py")
