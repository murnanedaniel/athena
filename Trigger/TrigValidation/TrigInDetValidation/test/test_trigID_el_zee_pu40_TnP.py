#!/usr/bin/env python
# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

# art-description: art job for el_zee_pu40
# art-type: grid
# art-include: master/Athena
# art-include: 22.0/Athena
# art-input: mc15_13TeV.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee.recon.RDO.e3601_s2665_s2183_r7191
# art-input-nfiles: 8
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


Slices  = ['electron-tnp']
Events  = 16000
Threads = 8 
Slots   = 8
Input   = 'Zee_pu40'    # defined in TrigValTools/share/TrigValInputs.json
Release = "current"
GridFiles = True
preinclude_file = 'all:TrigInDetValidation/TIDV_cond_fix.py' #conditions fix for ATR-23982. In future find a more recent RDO 

Jobs = [( "Truth",       " TIDAdata-run3.dat                    -o data-hists.root -p 11" ),
        ( "Offline",     " TIDAdata-run3-offline.dat -r Offline -o data-hists-offline.root" )] 
 
Comp = [( "L2ele",              "L2electronTnP",      "data-hists.root",         " -c TIDAhisto-panel-TnP.dat -l e26_e14_idperf_tight_50invmAB130_FTF_FE_1 e26_e14_idperf_tight_50invmAB130_FTF_FE_1_probe  -d HLTL2-plots " ),
        ( "L2eleoffline",       "L2electronTnP",      "data-hists-offline.root", " -c TIDAhisto-panel-TnP.dat -l e26_e14_idperf_tight_50invmAB130_FTF_FE_1 e26_e14_idperf_tight_50invmAB130_FTF_FE_1_probe -d HLTL2-plots-offline " ),
        ( "EFele",              "EFelectronTnP",      "data-hists.root",         " -c TIDAhisto-panel-TnP.dat -l e26_e14_nogsf_idperf_tight_50invmAB130_FTF_FE_1_probe e26_e14_nogsf_idperf_tight_50invmAB130_IDTrig_1_probe e26_e14_idperf_tight_50invmAB130_GSF_1_probe e26_e14_idperf_tight_50invmAB130_GSF_1 -d HLTEF-plots "         ),
        ( "EFeleoffline",       "EFelectronTnP",      "data-hists-offline.root", " -c TIDAhisto-panel-TnP.dat -l e26_e14_nogsf_idperf_tight_50invmAB130_FTF_FE_1_probe e26_e14_nogsf_idperf_tight_50invmAB130_IDTrig_1_probe e26_e14_idperf_tight_50invmAB130_GSF_1_probe e26_e14_idperf_tight_50invmAB130_GSF_1 -d HLTEF-plots-offline " ) ]


    

from AthenaCommon.Include import include 
include("TrigInDetValidation/TrigInDetValidation_Base.py")
