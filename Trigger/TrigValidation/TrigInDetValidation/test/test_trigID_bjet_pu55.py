#!/usr/bin/env python
# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

# art-description: art job for bjet_pu55
# art-type: grid
# art-include: master/Athena
# art-include: 22.0/Athena
# art-input-nfiles: 3
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

Slices = ['bjet']
Events  = 4000
Threads = 8 
Slots   = 8
Input   = 'ttbar'    # defined in TrigValTools/share/TrigValInputs.json  
Release = "current"

Jobs = [ ( "Truth",       " TIDAdata-run3.dat                    -o data-hists.root" ),
         ( "Offline",     " TIDAdata-run3-offline.dat -r Offline -o data-hists-offline.root" ) ]

Comp = [ ( "L2bjet",              "L2bjet",      "data-hists.root",         " -c TIDAhisto-panel.dat  -d HLTL2-plots " ),
         ( "L2bjetoffline",       "L2bjet",      "data-hists-offline.root", " -c TIDAhisto-panel.dat  -d HLTL2-plots-offline " ),
         ( "EFbjet",              "EFbjet",      "data-hists.root",         " -c TIDAhisto-panel.dat  -d HLTEF-plots " ),
         ( "EFbjetoffline",       "EFbjet",      "data-hists-offline.root", " -c TIDAhisto-panel.dat  -d HLTEF-plots-offline " ) ]


from AthenaCommon.Include import include 
include("TrigInDetValidation/TrigInDetValidation_Base.py")

