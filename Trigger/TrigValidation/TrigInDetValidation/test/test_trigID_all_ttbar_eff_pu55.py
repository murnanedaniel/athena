#!/usr/bin/env python

# art-description: art job for all_ttbar_eff_pu55
# art-type: grid
# art-include: master/Athena
# art-include: 22.0/Athena
# art-athena-mt: 8
# art-memory: 4096
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


Slices  = ['muon','electron','tau','bjet','fsjet']
Events  = 4000
Threads = 8 
Slots   = 8
Release = "current"

preinclude_file = 'RDOtoRDOTrigger:TrigInDetValidation/effpreinclude.py'

Input   = 'ttbar'    # defined in TrigValTools/share/TrigValInputs.json  

Jobs = [ ( "Offline",     " TIDAdata-run3-offline.dat            -r Offline -o data-hists-offline.root" ),
         ( "Offlinerz",   " TIDAdata-run3-offline-rzMatcher.dat  -r Offline -o data-hists-offline-rz.root" ),
         ( "OfflineVtx",  " TIDAdata-run3-offline-vtx.dat        -r Offline -o data-hists-offline-vtx.root" ) ]

Comp = [ ( "L2muon",       "L2muon",      "data-hists-offline.root",      " -c TIDAhisto-panel.dat  -d HLTL2-plots-muon " ),
         ( "L2electron",   "L2electron",  "data-hists-offline.root",      " -c TIDAhisto-panel.dat  -d HLTL2-plots-electron " ),
         ( "L2tau",        "L2tau",       "data-hists-offline.root",      " -c TIDAhisto-panel.dat  -d HLTL2-plots-tau " ),
         ( "L2bjet",       "L2bjet",      "data-hists-offline.root",      " -c TIDAhisto-panel.dat  -d HLTL2-plots-bjet " ),   
         ( "FSvtx",        "L2fsjetvtx",  "data-hists-offline-vtx.root",  " -c TIDAhisto-panel-vtx.dat  -d HLTL2-plots-vtx     --ncols 3" ),
         ( "FSvtxall",     "L2fsjetvtx",  "data-hists-offline.root",      " -c TIDAhisto-panel-vtx.dat  -d HLTL2-plots-vtxall  --ncols 3" ), 
         ( "FSjetoffline",   "L2fsjet",   "data-hists-offline.root",      " -c TIDAhisto-panel.dat  -d HLTL2-plots-FS " ),
         ( "FSjetofflinerz", "L2fsjet",   "data-hists-offline-rz.root",   " -c TIDAhisto-panel.dat  -d HLTL2-plots-rz-FS " ),

         ( "EFmuon",       "EFmuon",      "data-hists-offline.root",   " -c TIDAhisto-panel.dat  -d HLTEF-plots-muon " ),
         ( "EFelectron",   "EFelectron",  "data-hists-offline.root",   " -c TIDAhisto-panel.dat  -d HLTEF-plots-electron " ),
         ( "EFtau",        "EFtau",       "data-hists-offline.root",   " -c TIDAhisto-panel.dat  -d HLTEF-plots-tau " ),
         ( "EFbjet",       "EFbjet",      "data-hists-offline.root",   " -c TIDAhisto-panel.dat  -d HLTEF-plots-bjet " ),

         ( "EFmuonrz",     "EFmuon",      "data-hists-offline-rz.root",   " -c TIDAhisto-panel.dat  -d HLTEF-plots-rz-muon " ),
         ( "EFelectronrz", "EFelectron",  "data-hists-offline-rz.root",   " -c TIDAhisto-panel.dat  -d HLTEF-plots-rz-electron " ),
         ( "EFtaurz",      "EFtau",       "data-hists-offline-rz.root",   " -c TIDAhisto-panel.dat  -d HLTEF-plots-rz-tau " ),
         ( "EFbjetrz",     "EFbjet",      "data-hists-offline-rz.root",   " -c TIDAhisto-panel.dat  -d HLTEF-plots-rz-bjet " )  ]
   
from AthenaCommon.Include import include 
include("TrigInDetValidation/TrigInDetValidation_Base.py")


