#!/usr/bin/env python
# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

# art-description: art job for tau_ztautau_pu46
# art-type: grid
# art-include: master/Athena
# art-include: 22.0/Athena
# art-athena-mt: 8
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


Slices  = ['tau']
Events  = 6000
Threads = 8 
Slots   = 8
Release = "current"
Input   = 'Ztautau_pu46'    # defined in TrigValTools/share/TrigValInputs.json  

ExtraAna = " -c 'parentpdgid=15' "


Jobs = [ ( "Offline",  " TIDAdata-run3-offline.dat -r Offline -o data-hists-offline.root" ),
         ( "Truth",    " TIDAdata-run3.dat                    -o data-hists.root" ) ]


Comp = [ ( "L2tau",        "L2tau",       "data-hists.root",          " -c TIDAhisto-panel.dat  -d HLTL2-plots " ),
         ( "EFtau",        "EFtau",       "data-hists.root",          " -c TIDAhisto-panel.dat  -d HLTEF-plots " ),
         ( "L2tauOff",     "L2tau",       "data-hists-offline.root",  " -c TIDAhisto-panel.dat  -d HLTL2-plots-offline " ),
         ( "EFtauOff",     "EFtau",       "data-hists-offline.root",  " -c TIDAhisto-panel.dat  -d HLTEF-plots-offline " ),
         ( "EFvtx",        "EFtauvtx",    "data-hists-offline.root",  " -c TIDAhisto-panel-vtx.dat  -d HLTL2-plots-vtx  --ncols 3" ) ]

   
from AthenaCommon.Include import include 
include("TrigInDetValidation/TrigInDetValidation_Base.py")

