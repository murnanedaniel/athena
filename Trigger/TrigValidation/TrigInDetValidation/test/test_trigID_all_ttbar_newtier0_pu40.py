#!/usr/bin/env python

# art-description: art job for all_ttbar_newtier0_pu40
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

ExtraAna = " -c doNewTIDATier0=True "

Input   = 'ttbar'    # defined in TrigValTools/share/TrigValInputs.json  

Jobs = [] 

Comp = [ ( "L2muon",       "L2muon",      "data-hists-tier0.root",   " -b HLT/TRIDT/Muon/Expert     -s '_HLT_IDTrack' '/HLT_IDTrack' -c TIDAhisto-tier0.dat --ncols 3 -d HLTL2-plots-muon " ),
         ( "L2electron",   "L2electron",  "data-hists-tier0.root",   " -b HLT/TRIDT/Egamma/Expert   -s '_HLT_IDTrack' '/HLT_IDTrack' -c TIDAhisto-tier0.dat --ncols 3 -d HLTL2-plots-electron " ),
         ( "L2tau",        "L2tau",       "data-hists-tier0.root",   " -b HLT/TRIDT/Tau/Expert      -s '_HLT_IDTrack' '/HLT_IDTrack' -c TIDAhisto-tier0.dat --ncols 3 -d HLTL2-plots-tau " ),
         ( "L2bjet",       "L2bjet",      "data-hists-tier0.root",   " -b HLT/TRIDT/Bjet/Expert     -s '_HLT_IDTrack' '/HLT_IDTrack' -c TIDAhisto-tier0.dat --ncols 3 -d HLTL2-plots-bjet " ),   
         ( "FSjetoffline", "L2fsjet",     "data-hists-tier0.root",   " -b HLT/TRIDT/Bjet/Expert     -s '_HLT_IDTrack' '/HLT_IDTrack' -c TIDAhisto-tier0.dat --ncols 3 -d HLTL2-plots-FS " ),
         ( "FSt0vtx",      "L2fst0vtx",   "data-hists-tier0.root",   " -b HLT/TRIDT/Bjet/Expert     -s '_HLT_IDTrack' '/HLT_IDTrack' -c TIDAhisto-tier0-vtx.dat --ncols 3 -d HLTL2-plots-vtx " ),

         ( "EFmuon",       "EFmuon",      "data-hists-tier0.root",   " -b HLT/TRIDT/Muon/Expert     -s '_HLT_IDTrack' '/HLT_IDTrack' -c TIDAhisto-tier0.dat --ncols 3 -d HLTEF-plots-muon " ),
         ( "EFelectron",   "EFelectron",  "data-hists-tier0.root",   " -b HLT/TRIDT/Egamma/Expert   -s '_HLT_IDTrack' '/HLT_IDTrack' -c TIDAhisto-tier0.dat --ncols 3 -d HLTEF-plots-electron " ),
         ( "EFtau",        "EFtau",       "data-hists-tier0.root",   " -b HLT/TRIDT/Tau/Expert      -s '_HLT_IDTrack' '/HLT_IDTrack' -c TIDAhisto-tier0.dat --ncols 3 -d HLTEF-plots-tau " ),
         ( "EFbjet",       "EFbjet",      "data-hists-tier0.root",   " -b HLT/TRIDT/Bjet/Expert     -s '_HLT_IDTrack' '/HLT_IDTrack' -c TIDAhisto-tier0.dat --ncols 3 -d HLTEF-plots-bjet " ) ]
   
from AthenaCommon.Include import include 
include("TrigInDetValidation/TrigInDetValidation_Base.py")

