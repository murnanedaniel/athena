# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

# ********************************************************************
# TriggerContent.py
# 
# Configures content on EGAM derivations needed for trigger usage and 
# trigger studies
# author: fernando.monticelli@cern.ch
# ********************************************************************


# List of single photon triggers. Including those used for Bootstrap studies
singlePhotonTriggers = {}
singlePhotonTriggers['Run3'] = [
        'HLT_g10_loose_L1EM7',
        'HLT_g10_loose_L1eEM9',
        'HLT_g15_loose_L1EM10VH',
        'HLT_g15_loose_L1eEM12L',
        'HLT_g15_tight_L1EM10VH',
        'HLT_g15_tight_L1eEM12L',
        'HLT_g20_loose_L1EM15VH',
        'HLT_g20_loose_L1eEM18L',
        'HLT_g20_tight_L1EM15VHI',
        'HLT_g20_tight_L1eEM18M',
        'HLT_g20_tight_icaloloose_L1EM15VHI',
        'HLT_g20_tight_icaloloose_L1eEM18M',
        'HLT_g22_tight_L1EM15VHI',
        'HLT_g22_tight_L1eEM18M',
        'HLT_g25_loose_L1EM20VH',
        'HLT_g25_loose_L1eEM24L',
        'HLT_g25_medium_L1EM20VH',
        'HLT_g25_medium_L1eEM24L',
        'HLT_g25_tight_icaloloose_L1EM20VH',
        'HLT_g25_tight_icalomedium_L1EM20VH',
        'HLT_g25_tight_icalotight_L1EM20VH',
        'HLT_g300_etcut_L1eEM26M',
        'HLT_g30_loose_L1EM20VH',
        'HLT_g30_loose_L1eEM24L',
        'HLT_g35_medium_L1EM20VH',
        'HLT_g35_medium_L1eEM24L',
        'HLT_g40_loose_L1EM20VH',
        'HLT_g40_loose_L1eEM24L',
        'HLT_g50_loose_L1EM20VH',
        'HLT_g50_loose_L1eEM24L',
        'HLT_g60_loose_L1EM22VHI',
        'HLT_g60_loose_L1eEM26M',
        'HLT_g80_loose_L1EM22VHI',
        'HLT_g80_loose_L1eEM26M',
        'HLT_g100_loose_L1EM22VHI',
        'HLT_g100_loose_L1eEM26M',
        'HLT_g120_loose_L1EM22VHI',
        'HLT_g120_loose_L1eEM26M',
        'HLT_g140_loose_L1EM22VHI',
        'HLT_g140_loose_L1eEM26M',
        'HLT_g250_etcut_L1EM22VHI',
        'HLT_g250_etcut_L1eEM26M',
        'HLT_g300_etcut_L1EM22VHI',
        ]

singlePhotonTriggers['Run2'] = [
                        'HLT_g6_loose',
                        'HLT_g6_tight_icalotight',
                        'HLT_g10_etcut',
                        'HLT_g10_loose',
                        'HLT_g10_loose_L1EM3',
                        'HLT_g10_medium',
                        'HLT_g10_medium_L1EM3',
                        'HLT_g12_loose',
                        'HLT_g12_medium',
                        'HLT_g15_etcut_L1EM7',
                        'HLT_g15_loose',
                        'HLT_g15_loose_L1EM3',
                        'HLT_g15_loose_L1EM7',
                        'HLT_g15_loose_L1EM8VH',
                        'HLT_g20_etcut',
                        'HLT_g20_etcut_L1EM12',
                        'HLT_g20_etcut_L1EM15',
                        'HLT_g20_loose',
                        'HLT_g20_loose_L1EM12',
                        'HLT_g20_loose_L1EM15',
                        'HLT_g20_loose_L1EM18VH',
                        'HLT_g20_medium',
                        'HLT_g20_medium_L1EM15',
                        'HLT_g20_tight',
                        'HLT_g20_tight_L1EM15',
                        'HLT_g20_tight_L1EM15VHI',
                        'HLT_g20_tight_icalotight_L1EM15VHI',
                        'HLT_g20_tight_icalovloose_L1EM15VHI',
                        'HLT_g22_tight',
                        'HLT_g22_tight_L1EM15VHI',
                        'HLT_g22_tight_icalotight',
                        'HLT_g22_tight_icalotight_L1EM15VHI',
                        'HLT_g22_tight_icalovloose',
                        'HLT_g22_tight_icalovloose_L1EM15VHI',
                        'HLT_g24_loose',
                        'HLT_g25_etcut_L1EM15',
                        'HLT_g25_loose',
                        'HLT_g25_loose_L1EM15 ',
                        'HLT_g25_loose_L1EM20VH',
                        'HLT_g25_medium',
                        'HLT_g25_medium_L1EM15',
                        'HLT_g25_medium_L1EM20VH',
                        'HLT_g25_medium_L1EM22VHI',
                        'HLT_g25_tight_L1EM15',
                        'HLT_g25_tight_L1EM20VH',
                        'HLT_g30_etcut_L1EM15'
                        'HLT_g30_loose_L1EM15',
                        'HLT_g35_loose',
                        'HLT_g35_loose_L1EM15',
                        'HLT_g35_loose_L1EM20VH',
                        'HLT_g35_loose_L1EM22VHI',
                        'HLT_g35_loose_L1EM24VHI',
                        'HLT_g35_loose_L1EM24VHIM',
                        'HLT_g35_loose_L1EM26VHI',
                        'HLT_g35_medium',
                        'HLT_g35_medium_L1EM20VH',
                        'HLT_g35_medium_L1EM22VHI',
                        'HLT_g35_medium_L1EM24VHI',
                        'HLT_g35_medium_icalotight',
                        'HLT_g35_medium_icalovloose',
                        'HLT_g35_tight_icalotight_L1EM24VHI',
                        'HLT_g35_tight_icalotight_L1EM24VHIM',
                        'HLT_g40_loose_L1EM15',
                        'HLT_g40_tight',
                        'HLT_g40_tight_icalotight_L1EM24VHI',
                        'HLT_g40_tight_icalotight_L1EM24VHIM',
                        'HLT_g45_loose_L1EM15',
                        'HLT_g45_tight',
                        'HLT_g45_tight_L1EM22VHI',
                        'HLT_g45_tight_L1EM24VHI',
                        'HLT_g50_loose',
                        'HLT_g50_loose_L1EM15',
                        'HLT_g50_loose_L1EM20VH',
                        'HLT_g60_loose',
                        'HLT_g60_loose_L1EM15VH ',
                        'HLT_g60_loose_L1EM20VH',
                        'HLT_g60_loose_L1EM24VHI',
                        'HLT_g60_loose_L1EM24VHIM',
                        'HLT_g60_loose_L1EM26VHI',
                        'HLT_g70_loose',
                        'HLT_g70_loose_L1EM24VHI',
                        'HLT_g70_loose_L1EM24VHIM',
                        'HLT_g70_loose_L1EM26VHI',
                        'HLT_g70_loose_L1EN24VHIM',
                        'HLT_g80_loose',
                        'HLT_g80_loose_L1EM24VHI',
                        'HLT_g80_loose_L1EM24VHIM',
                        'HLT_g80_loose_L1EM26VHI',
                        'HLT_g80_loose_icalovloose_L1EM24VHIM',
                        'HLT_g85_tight',
                        'HLT_g85_tight_L1EM24VHI',
                        'HLT_g85_tight_L1EM24VHIM',
                        'HLT_g85_tight_icalovloose_L1EM24VHI',
                        'HLT_g85_tight_icalovloose_L1EM24VHIM',
                        'HLT_g100_loose',
                        'HLT_g100_loose_L1EM24VHI',
                        'HLT_g100_loose_L1EM24VHIM',
                        'HLT_g100_loose_L1EM26VHI',
                        'HLT_g100_tight',
                        'HLT_g100_tight_L1EM24VHI',
                        'HLT_g100_tight_L1EM24VHIM',
                        'HLT_g100_tight_icalovloose_L1EM24',
                        'HLT_g100_tight_icalovloose_L1EM24VHI',
                        'HLT_g100_tight_icalovloose_L1EM24VHIM',
                        'HLT_g120_loose',
                        'HLT_g120_loose_L1EM24VHI',
                        'HLT_g120_loose_L1EM24VHIM',
                        'HLT_g120_loose_L1EM26VHI',
                        'HLT_g140_loose',
                        'HLT_g140_loose_HLTCalo',
                        'HLT_g140_loose_L1EM24VHIM',
                        'HLT_g140_loose_L1EM26VHI',
                        'HLT_g140_tight',
                        'HLT_g140_tight_L1EM24VHIM',
                        'HLT_g160_loose',
                        'HLT_g160_loose_L1EM24VHIM',
                        'HLT_g160_loose_L1EM26VHI',
                        'HLT_g180_loose',
                        'HLT_g180_loose_L1EM24VHIM',
                        'HLT_g180_loose_L1EM26VHI',
                        'HLT_g200_etcut',
                        'HLT_g200_loose',
                        'HLT_g200_loose_L1EM24VHIM',
                        'HLT_g200_loose_L1EM26VHI',
                        'HLT_g250_etcut',
                        'HLT_g300_etcut'
                        'HLT_g300_etcut_L1EM24VHI',
                        'HLT_g300_etcut_L1EM24VHIM',
                        ]


diPhotonTriggers = {}

diPhotonTriggers['Run3'] = [
        'HLT_2g15_loose_25dphiAA_invmAA80_L12EM7', 
        'HLT_2g15_loose_25dphiAA_invmAA80_L1DPHI-M70-2eEM15M', 
        'HLT_2g15_tight_25dphiAA_invmAA80_L12EM7', 
        'HLT_2g15_tight_25dphiAA_invmAA80_L1DPHI-M70-2eEM15M', 
        'HLT_2g15_tight_25dphiAA_L12EM7', 
        'HLT_2g15_tight_25dphiAA_L1DPHI-M70-2eEM15M', 
        'HLT_2g20_loose_L12eEM18L', 
        'HLT_2g20_loose_L12EM15VH', 
        'HLT_2g20_tight_icaloloose_L12eEM18M', 
        'HLT_2g20_tight_icaloloose_L12EM15VHI', 
        'HLT_2g20_tight_L12EM15VHI', 
        'HLT_2g22_tight_L12eEM18M', 
        'HLT_2g22_tight_L12EM15VHI', 
        'HLT_2g22_tight_L1eEM9_EMPTY', 
        'HLT_2g22_tight_L1eEM9_UNPAIRED_ISO', 
        'HLT_2g22_tight_L1EM7_EMPTY', 
        'HLT_2g22_tight_L1EM7_UNPAIRED_ISO', 
        'HLT_2g25_loose_g15_loose_L12eEM24L', 
        'HLT_2g25_loose_g15_loose_L12EM20VH', 
        'HLT_2g50_loose_L12eEM24L', 
        'HLT_2g50_loose_L12EM20VH', 
        'HLT_2g50_tight_L1eEM9_EMPTY', 
        'HLT_2g50_tight_L1eEM9_UNPAIRED_ISO', 
        'HLT_2g50_tight_L1EM7_EMPTY', 
        'HLT_2g50_tight_L1EM7_UNPAIRED_ISO', 
        'HLT_2g9_loose_25dphiAA_invmAA80_L12EM7', 
        'HLT_2g9_loose_25dphiAA_invmAA80_L1DPHI-M70-2eEM9', 
        'HLT_2g9_loose_25dphiAA_invmAA80_L1DPHI-M70-2eEM9L', 
        'HLT_g35_medium_g25_medium_L12eEM24L', 
        'HLT_g35_medium_g25_medium_L12EM20VH', 
        'HLT_g35_medium_g25_medium_L1eEM9_EMPTY', 
        'HLT_g35_medium_g25_medium_L1eEM9_UNPAIRED_ISO', 
        'HLT_g35_medium_g25_medium_L1EM7_EMPTY', 
        'HLT_g35_medium_g25_medium_L1EM7_UNPAIRED_ISO', 
        ]


diPhotonTriggers['Run2'] = [
                    'HLT_2g20_loose_L12EM15',
                    'HLT_2g20_loose',
                    'HLT_2g20_tight',
                    'HLT_2g22_tight',
                    'HLT_2g25_tight',
                    'HLT_g35_loose_g25_loose',
                    'HLT_g35_medium_HLTCalo_g25_medium_HLTCalo',
                    'HLT_g35_loose_L1EM15_g25_loose_L1EM15',
                    'HLT_g35_loose_L1EM15VH_g25_loose_L1EM15VH',
                    'HLT_g35_medium_g25_medium',
                    'HLT_2g50_loose',
                    'HLT_2g60_loose_L12EM15VH ',
                    'HLT_2g10_loose',
                    'HLT_2g50_loose_L12EM18VH',
                    'HLT_2g60_loose_L12EM18VH',
                    'HLT_2g50_loose_L12EM20VH',
                    'HLT_g50_loose_L12EM18VH',
                    'HLT_g60_loose_L12EM18VH',
                    'HLT_g50_loose_L12EM20VH',
                    'HLT_g60_loose_L12EM20VH',
                    'HLT_2g25_tight_L12EM20VH',
                    'HLT_g35_loose_g25_loose_L12EM18VH',
                    'HLT_g35_loose_g25_loose_L12EM20VH ',
                    'HLT_g35_medium_g25_medium_L12EM18VH',
                    'HLT_g35_medium_g25_medium_L12EM20VH',
                    'HLT_2g20_tight_L12EM15VHI',
                    'HLT_2g20_tight_icalovloose_L12EM15VHI',
                    'HLT_2g20_tight_icalotight_L12EM15VHI',
                    'HLT_2g22_tight_L12EM15VHI',
                    'HLT_2g22_tight_icalovloose_L12EM15VHI',
                    'HLT_2g22_tight_icalotight_L12EM15VHI',
                    'HLT_2g60_loose_L12EM20VH',
                    'HLT_2g3_loose_dPhi15_L12EM3_VTE50',
                    'HLT_2g3_loose_L12EM3_VTE50',
                    'HLT_2g3_medium_dPhi15_L12EM3_VTE50',
                    'HLT_2g22_tight_icalovloose',
                    'HLT_2g22_tight_icalotight',
                    'HLT_2g10_loose_L12EM7',
                    'HLT_2g15_loose_L12EM7']


triPhotonTriggers = {}
triPhotonTriggers['Run3'] = [
        'HLT_2g25_loose_g15_loose_L12EM20VH',
        'HLT_2g25_loose_g15_loose_L12eEM24L'
        ]
triPhotonTriggers['Run2'] = [
                     'HLT_3g15_loose',
                     'HLT_g20_loose_2g15_loose_L12EM13VH',
                     'HLT_2g20_loose_g15_loose',
                     'HLT_3g20_loose',
                     'HLT_3g20_loose_L12EM18VH',
                     'HLT_2g24_loose_g15_loose',
                     'HLT_2g24_g20_loose',
                     'HLT_3g24_loose_L12EM20VH',
                     'HLT_2g25_loose_g15_loose',
                     'HLT_2g25_loose_g20_loose',
                     'HLT_3g25_loose']


JPsiTriggers = {}

JPsiTriggers['Run3'] = [
        'HLT_e9_lhtight_e4_idperf_tight_probe_1invmAB5_L1JPSI-1M5-EM7',
        'HLT_e9_lhtight_e4_idperf_tight_nogsf_probe_1invmAB5_L1JPSI-1M5-EM7',
        'HLT_e14_lhtight_e4_idperf_tight_probe_1invmAB5_L1JPSI-1M5-EM12',
        'HLT_e14_lhtight_e4_idperf_tight_nogsf_probe_1invmAB5_L1JPSI-1M5-EM12',
        'HLT_e9_lhtight_e4_idperf_tight_probe_1invmAB5_L1JPSI-1M5-eEM9',
        'HLT_e9_lhtight_e4_idperf_tight_nogsf_probe_1invmAB5_L1JPSI-1M5-eEM9',
        'HLT_e14_lhtight_e4_idperf_tight_probe_1invmAB5_L1JPSI-1M5-eEM15',
        'HLT_e14_lhtight_e4_idperf_tight_nogsf_probe_1invmAB5_L1JPSI-1M5-eEM15',
        'HLT_e9_lhtight_e4_etcut_1invmAB5_L1JPSI-1M5-EM7',
        'HLT_e5_lhtight_e9_etcut_1invmAB5_L1JPSI-1M5-EM7',
        'HLT_e5_lhtight_e14_etcut_1invmAB5_L1JPSI-1M5-EM12',
        'HLT_e14_lhtight_e4_etcut_1invmAB5_L1JPSI-1M5-EM12',
        'HLT_e9_lhtight_e4_etcut_1invmAB5_L1JPSI-1M5-eEM9',
        'HLT_e5_lhtight_e9_etcut_1invmAB5_L1JPSI-1M5-eEM9',
        'HLT_e5_lhtight_e14_etcut_1invmAB5_L1JPSI-1M5-eEM15',
        'HLT_e14_lhtight_e4_etcut_1invmAB5_L1JPSI-1M5-eEM15',
        ]
JPsiTriggers['Run2'] = [
        'HLT_e5_lhtight_e4_etcut_Jpsiee',
        'HLT_e5_lhtight_nod0_e4_etcut_Jpsiee',
        'HLT_e5_lhtight_e4_etcut',
        'HLT_e5_lhtight_nod0_e4_etcut',

        'HLT_e9_lhtight_e4_etcut_Jpsiee',
        'HLT_e9_lhtight_nod0_e4_etcut_Jpsiee',
        'HLT_e9_etcut_e5_lhtight_nod0_Jpsiee',
        'HLT_e9_etcut_e5_lhtight_Jpsiee',

        'HLT_e14_etcut_e5_lhtight_Jpsiee',
        'HLT_e14_etcut_e5_lhtight_nod0_Jpsiee',
        'HLT_e14_lhtight_e4_etcut_Jpsiee',
        'HLT_e14_lhtight_nod0_e4_etcut_Jpsiee',

        'HLT_e5_lhtight_nod0_e4_etcut_Jpsiee_L1RD0_FILLED',
        'HLT_e5_lhtight_nod0_e9_etcut_Jpsiee',
        'HLT_e5_lhtight_nod0_e14_etcut_Jpsiee',
        'HLT_e5_lhtight_nod0_e9_etcut_Jpsiee_L1JPSI-1M5-EM7',
        'HLT_e9_lhtight_nod0_e4_etcut_Jpsiee_L1JPSI-1M5-EM7',
        'HLT_e5_lhtight_nod0_e14_etcut_Jpsiee_L1JPSI-1M5-EM12',
        'HLT_e14_lhtight_nod0_e4_etcut_Jpsiee_L1JPSI-1M5-EM12',
        ]


WTnPTriggers = {}

WTnPTriggers['Run3'] = [
        # No W TnP triggers in the menu?
        ]

WTnPTriggers['Run2'] = [
        # L1Topo W T&P 
        'HLT_e13_etcut_trkcut' ,
        'HLT_e18_etcut_trkcut' ,
        ## # Non-L1Topo W TP commissioning triggers ==> in MC, in 50 ns data
        'HLT_e13_etcut_trkcut_xs15' ,
        'HLT_e18_etcut_trkcut_xs20' ,
        ## W T&P triggers ==> not in MC, in 50 ns data
        'HLT_e13_etcut_trkcut_xs15_mt25' ,
        'HLT_e18_etcut_trkcut_xs20_mt35' ,
        ###W T&P triggers ==> not in MC, not in 50 ns data, will be in 25 ns data
        'HLT_e13_etcut_trkcut_xs15_j20_perf_xe15_2dphi05' ,
        'HLT_e13_etcut_trkcut_xs15_j20_perf_xe15_2dphi05_mt25' ,
        'HLT_e13_etcut_trkcut_j20_perf_xe15_2dphi05_mt25' ,
        'HLT_e13_etcut_trkcut_j20_perf_xe15_2dphi05' ,
        'HLT_e13_etcut_trkcut_xs15_j20_perf_xe15_6dphi05' ,
        'HLT_e13_etcut_trkcut_xs15_j20_perf_xe15_6dphi05_mt25' ,
        'HLT_e13_etcut_trkcut_j20_perf_xe15_6dphi05_mt25' ,
        'HLT_e13_etcut_trkcut_j20_perf_xe15_6dphi05' ,
        'HLT_e18_etcut_trkcut_xs20_j20_perf_xe20_6dphi15' ,
        'HLT_e18_etcut_trkcut_xs20_j20_perf_xe20_6dphi15_mt35' ,
        'HLT_e18_etcut_trkcut_j20_perf_xe20_6dphi15_mt35' ,
        'HLT_e18_etcut_trkcut_j20_perf_xe20_6dphi15' ,

        # others
        'HLT_e5_etcut_L1W-05DPHI-JXE-0',
        'HLT_e5_etcut_L1W-10DPHI-JXE-0',
        'HLT_e5_etcut_L1W-15DPHI-JXE-0',
        'HLT_e5_etcut_L1W-10DPHI-EMXE-0',
        'HLT_e5_etcut_L1W-15DPHI-EMXE-0',
        'HLT_e5_etcut_L1W-05DPHI-EMXE-1',
        'HLT_e5_etcut_L1W-05RO-XEHT-0',
        'HLT_e5_etcut_L1W-90RO2-XEHT-0',
        'HLT_e5_etcut_L1W-250RO2-XEHT-0',
        'HLT_e5_etcut_L1W-HT20-JJ15.ETA49',

        'HLT_e13_etcut_L1W-NOMATCH',
        'HLT_e13_etcut_L1W-NOMATCH_W-05RO-XEEMHT',
        'HLT_e13_etcut_L1EM10_W-MT25',
        'HLT_e13_etcut_L1EM10_W-MT30',
        'HLT_e13_etcut_trkcut_L1EM12',
        'HLT_e13_etcut_trkcut_L1EM10_W-MT25_W-15DPHI-JXE-0_W-15DPHI-EMXE',
        'HLT_e13_etcut_trkcut_j20_perf_xe15_6dphi15_mt25',
        'HLT_e13_etcut_trkcut_j20_perf_xe15_6dphi15_mt25_L1EM12_W-MT25_W-15DPHI-JXE-0_W-15DPHI-EMXE_XS20',
        'HLT_e13_etcut_trkcut_j20_perf_xe15_6dphi15_mt25_L1EM12_W-MT25_W-15DPHI-JXE-0_W-15DPHI-EMXE_W-90RO2-XEHT-0',
        'HLT_e13_etcut_trkcut_xs30_xe30_mt35',
        'HLT_e13_etcut_trkcut_xs30_j15_perf_xe30_6dphi05_mt35',
        'HLT_e13_etcut_trkcut_xs30_j15_perf_xe30_6dphi15_mt35',
        'HLT_e13_etcut_trkcut_xs30_j15_perf_xe30_2dphi05_mt35',
        'HLT_e13_etcut_trkcut_xs30_j15_perf_xe30_2dphi15_mt35',
        'HLT_e13_etcut_trkcut_xs30_j15_perf_xe30_2dphi15_mt35_L1EM12_W-MT25_W-15DPHI-JXE-0_W-15DPHI-EMXE_XS20',
        'HLT_e13_etcut_trkcut_xs30_j15_perf_xe30_6dphi15_mt35_L1EM12_W-MT25_W-15DPHI-JXE-0_W-15DPHI-EMXE_W-90RO2-XEHT-0',

        'HLT_e18_etcut_L1EM15_W-MT35',
        'HLT_e18_etcut_trkcut_L1EM15',
        'HLT_e18_etcut_trkcut_L1EM15_W-MT35_W-05DPHI-JXE-0_W-05DPHI-EMXE',
        'HLT_e18_etcut_trkcut_xs30_xe30_mt35',
        'HLT_e18_etcut_trkcut_xs30_j15_perf_xe30_6dphi05_mt35',
        'HLT_e18_etcut_trkcut_xs30_j15_perf_xe30_6dphi15_mt35',
        'HLT_e18_etcut_trkcut_xs30_j15_perf_xe30_2dphi05_mt35',
        'HLT_e18_etcut_trkcut_xs30_j15_perf_xe30_2dphi15_mt35',
        'HLT_e18_etcut_trkcut_xs30_j15_perf_xe30_6dphi05_mt35_L1EM15_W-MT35_W-05DPHI-JXE-0_W-05DPHI-EM15XE_XS30', 
        'HLT_e18_etcut_trkcut_xs30_j15_perf_xe30_6dphi15_mt35_L1EM15_W-MT35_W-05DPHI-JXE-0_W-05DPHI-EM15XE_XS30', 
        'HLT_e18_etcut_trkcut_xs30_j15_perf_xe30_2dphi05_mt35_L1EM15_W-MT35_W-250RO2-XEHT-0_W-05DPHI-JXE-0_W-05DPHI-EM15XE', 
        'HLT_e18_etcut_trkcut_xs30_j15_perf_xe30_2dphi15_mt35_L1EM15_W-MT35_W-250RO2-XEHT-0_W-05DPHI-JXE-0_W-05DPHI-EM15XE', 
        'HLT_e18_etcut_trkcut_xs30_j15_perf_xe30_2dphi15_mt35_L1EM15_W-MT35_W-05DPHI-JXE-0_W-05DPHI-EM15XE_XS30',
        'HLT_e18_etcut_trkcut_xs30_j15_perf_xe30_6dphi05_mt35_L1EM15_W-MT35_W-250RO2-XEHT-0_W-05DPHI-JXE-0_W-05DPHI-EM15XE',
        'HLT_e18_etcut_trkcut_xs30_j15_perf_xe30_6dphi15_mt35_L1EM15_W-MT35_W-250RO2-XEHT-0_W-05DPHI-JXE-0_W-05DPHI-EM15XE',
        'HLT_e18_etcut_trkcut_xs30_j15_perf_xe30_6dphi15_mt35_L1EM15_W-MT35_W-250RO2-XEHT-0_W-15DPHI-JXE-0_W-15DPHI-EM15XE',
        'HLT_e18_etcut_trkcut_xs30_xe30_mt35_L1EM15_W-MT35_W-05DPHI-JXE-0_W-05DPHI-EM15XE_XS30', 
        'HLT_e18_etcut_trkcut_xs30_xe30_mt35_L1EM15_W-MT35_W-250RO2-XEHT-0_W-05DPHI-JXE-0_W-05DPHI-EM15XE', 
        'HLT_e18_etcut_trkcut_xs30_xe30_mt35_L1EM15_W-MT35_W-250RO2-XEHT-0_W-15DPHI-JXE-0_W-15DPHI-EM15XE',
        'HLT_e18_etcut_trkcut_j20_perf_xe20_6dphi15_mt35_L1EM15_W-MT35_W-05DPHI-JXE-0_W-05DPHI-EM15XE_XS30', 
        'HLT_e18_etcut_trkcut_j20_perf_xe20_6dphi15_mt35_L1EM15_W-MT35_W-250RO2-XEHT-0_W-05DPHI-JXE-0_W-05DPHI-EM15XE',

        # added for 2017
        'HLT_e60_etcut',
        'HLT_e60_etcut_L1EM24VHIM',
        'HLT_e60_etcut_trkcut_L1EM24VHIM_j15_perf_xe60_6dphi15_mt35',
        'HLT_e60_etcut_trkcut_L1EM24VHIM_xe60_mt35',
        'HLT_e60_etcut_trkcut_L1EM24VHIM_xs30_j15_perf_xe30_6dphi15_mt35',
        'HLT_e60_etcut_trkcut_L1EM24VHIM_xs30_xe30_mt35',
        'HLT_e60_lhmedium_nod0',
        'HLT_e60_lhmedium_nod0_L1EM24VHI',
        'HLT_e60_lhmedium_nod0_L1EM24VHIM',
        'HLT_e60_lhvloose_nod0',
        'HLT_e60_etcut_trkcut_j15_perf_xe60_6dphi05_mt35',
        'HLT_e60_etcut_trkcut_xs30_j15_perf_xe30_6dphi05_mt35',
        'HLT_e70_etcut',
        'HLT_e70_etcut_L1EM24VHIM',
        'HLT_e70_lhloose_nod0_L1EM24VHIM_xe70noL1',
        'HLT_e70_lhloose_nod0_xe70noL1',
        'HLT_noalg_l1topo_L1EM15',
        'HLT_noalg_l1topo_L1EM7',
        'HLT_j80_xe80',
        'HLT_xe80_tc_lcw_L1XE50',
        'HLT_xe90_mht_L1XE50',
        'HLT_xe90_tc_lcw_wEFMu_L1XE50',
        'HLT_xe90_mht_wEFMu_L1XE50',
        'HLT_xe110_mht_L1XE50',
        'HLT_xe110_pufit_L1XE50',

        #added for low-mu data analysis, 2017 and 2018 data
        'HLT_e15_lhloose_nod0_L1EM12',
        #added for low-mu data analysis, 2018 data
        'HLT_xe35',
        'HLT_e15_etcut_trkcut_xe30noL1',
]


noalgTriggers={}
noalgTriggers['Run2'] = [
        'HLT_noalg_L1EM15VH',
        'HLT_noalg_L1EM12',
        'HLT_noalg_L1EM15',
        'HLT_noalg_L1EM18VH',
        'HLT_noalg_L1EM20VH',
        'HLT_noalg_L1EM10',
        'HLT_noalg_L1EM10VH',
        'HLT_noalg_L1EM13VH',
        'HLT_noalg_L1EM20VHI',
        'HLT_noalg_L1EM22VHI',
        'HLT_noalg_L1EM8VH',
    ]
noalgTriggers['Run3'] = [
        'HLT_noalg_L1EM22VHI',
        'HLT_noalg_L1EM3',
        'HLT_noalg_L1EM7',
        'HLT_noalg_L1EM12',
        'HLT_noalg_L1EM15',
        'HLT_noalg_L1EM8VH',
        'HLT_noalg_L1EM10VH',
        'HLT_noalg_L1EM15VH',
        'HLT_noalg_L1EM20VH',
        'HLT_noalg_L1eEM5',
        'HLT_noalg_L1eEM7',
        'HLT_noalg_L1eEM9',
        'HLT_noalg_L1eEM10L',
        'HLT_noalg_L1eEM12L',
        'HLT_noalg_L1eEM15',
        'HLT_noalg_L1eEM18',
        'HLT_noalg_L1eEM18L',
        'HLT_noalg_L1eEM18M',
        'HLT_noalg_L1eEM22M',
        'HLT_noalg_L1eEM24L',
        'HLT_noalg_L1eEM24VM',
        'HLT_noalg_L1eEM26',
        'HLT_noalg_L1eEM26L',
        'HLT_noalg_L1eEM26M',
        'HLT_noalg_L1eEM26T',
        ]



#====================================================================
# trigger-based selection
# prescaled _etcut triggers
# prescaled _loose triggers
# prescaled _lhloose triggers
#====================================================================
BkgElectronTriggers = {}
BkgElectronTriggers['Run3'] = [
        'HLT_e5_nopid_L1EM3',
        'HLT_e5_etcut_L1EM3',
        'HLT_e50_etcut_L1EM22VHI',
        'HLT_e50_etcut_L1eEM26M',
        'HLT_e120_etcut_L1EM22VHI',
        'HLT_e120_etcut_L1eEM26M',
        'HLT_e250_etcut_L1EM22VHI',
        'HLT_e250_etcut_L1eEM26M',
        'HLT_e300_etcut_L1eEM26M',
        'HLT_e300_etcut_L1EM22VHI',
        'HLT_g250_etcut_L1EM22VHI',
        'HLT_g250_etcut_L1eEM26M',
        'HLT_g300_etcut_L1EM22VHI',
        'HLT_g300_etcut_L1eEM26M',
        ]

BkgElectronTriggers['Run2'] = [
        'HLT_e4_etcut'        ,
        'HLT_e5_etcut'        ,
        'HLT_e9_etcut'        ,            
        'HLT_e10_etcut_L1EM7' ,            
        'HLT_e14_etcut'       ,            
        'HLT_e15_etcut_L1EM7' ,
        'HLT_e17_etcut_L1EM15',            
        'HLT_e20_etcut_L1EM12',            
        'HLT_e25_etcut_L1EM15',            
        'HLT_e30_etcut_L1EM15',            
        'HLT_e40_etcut_L1EM15',            
        'HLT_e50_etcut_L1EM15',            
        'HLT_e60_etcut'       ,            
        'HLT_e80_etcut'       ,            
        'HLT_e100_etcut'      ,            
        'HLT_e120_etcut'      ,            
        'HLT_g10_etcut'       ,            
        'HLT_g20_etcut_L1EM12',            
        'HLT_g200_etcut'      ,            
        'HLT_e5_lhloose'                      ,
        'HLT_e5_lhvloose'                     ,
        'HLT_e5_loose'                        ,
        'HLT_e5_vloose'                       ,
        'HLT_e10_lhvloose_L1EM7'              ,
        'HLT_e10_vloose_L1EM7'                ,
        'HLT_e12_lhloose'                     ,
        'HLT_e12_lhloose_L1EM10VH'            ,
        'HLT_e12_lhvloose_L1EM10VH'           ,
        'HLT_e12_loose'                       ,
        'HLT_e12_loose_L1EM10VH'              ,
        'HLT_e12_vloose_L1EM10VH'             ,
        'HLT_e15_lhloose_L1EM13VH'            ,
        'HLT_e15_lhvloose_L1EM13VH'           ,
        'HLT_e15_lhvloose_L1EM7'              ,
        'HLT_e15_loose_L1EM13VH'              ,
        'HLT_e15_vloose_L1EM13VH'             ,
        'HLT_e15_vloose_L1EM7'                ,
        'HLT_e17_lhloose'                     ,
        'HLT_e17_lhloose_L1EM15'              ,
        'HLT_e17_lhloose_cutd0dphideta_L1EM15',
        'HLT_e17_lhloose_nod0_L1EM15'         ,
        'HLT_e17_lhloose_nodeta_L1EM15'       ,
        'HLT_e17_lhloose_nodphires_L1EM15'    ,
        'HLT_e17_lhloose_L1EM15VHJJ1523ETA49' ,
        'HLT_e17_lhvloose'                    ,
        'HLT_e17_loose'                       ,
        'HLT_e17_loose_L1EM15'                ,
        'HLT_e17_loose_L1EM15VHJJ1523ETA49'   ,
        'HLT_e17_vloose'                      ,
        'HLT_e20_lhvloose'                    ,
        'HLT_e20_lhvloose_L1EM12'             ,
        'HLT_e20_vloose'                      ,
        'HLT_e20_vloose_L1EM12'               ,
        'HLT_e25_lhvloose_L1EM15'             ,
        'HLT_e25_vloose_L1EM15'               ,
        'HLT_e30_lhvloose_L1EM15'             ,
        'HLT_e30_vloose_L1EM15'               ,
        'HLT_e40_lhvloose'                    ,
        'HLT_e40_lhvloose_L1EM15'             ,
        'HLT_e40_vloose_L1EM15'               ,
        'HLT_e50_lhvloose_L1EM15'             ,
        'HLT_e50_vloose_L1EM15'               ,
        'HLT_e60_loose'                       ,
        'HLT_e60_vloose'                      ,
        'HLT_e60_lhvloose'                    ,
        'HLT_e70_etcut'                       ,
        'HLT_e70_lhloose'                     ,
        'HLT_e70_lhvloose'                    ,
        'HLT_e70_loose'                       ,
        'HLT_e70_vloose'                      ,
        'HLT_e80_lhvloose'                    ,
        'HLT_e80_vloose'                      ,
        'HLT_e100_lhvloose'                   ,
        'HLT_e100_vloose'                     ,
        'HLT_e120_lhvloose'                   ,
        'HLT_e120_lhloose'                    ,
        'HLT_e120_loose'                      ,
        'HLT_e120_vloose'                     ,
        'HLT_e140_etcut'                      ,
        'HLT_e160_etcut'                      ,
        'HLT_e180_etcut'                      ,
        'HLT_e200_etcut'                      ,
        'HLT_e250_etcut'                      ,
        'HLT_e300_etcut'                      ,
        'HLT_g250_etcut'                      ,
        'HLT_g300_etcut'                      ,
        ]


BootstrapPhotonTriggers = {}
BootstrapPhotonTriggers['Run3'] = [
        'HLT_g25_medium_L1EM20VH',
        'HLT_g35_medium_L1EM20VH',
        'HLT_g20_tight_icaloloose_L1EM15VHI',
        'HLT_g15_tight_L1EM10VH',
        'HLT_g20_tight_L1EM15VHI',
        'HLT_g22_tight_L1EM15VHI',
        'HLT_g25_medium_L1eEM24L',
        'HLT_g35_medium_L1eEM24L',
        'HLT_g20_tight_icaloloose_L1eEM18M',
        'HLT_g15_tight_L1eEM12L',
        'HLT_g20_tight_L1eEM18M',
        'HLT_g22_tight_L1eEM18M',
        'HLT_g250_etcut_L1EM22VHI',
        'HLT_g10_loose_L1EM7',
        'HLT_g15_loose_L1EM10VH',
        'HLT_g20_loose_L1EM15VH',
        'HLT_g25_loose_L1EM20VH',
        'HLT_g30_loose_L1EM20VH',
        'HLT_g40_loose_L1EM20VH',
        'HLT_g50_loose_L1EM20VH',
        'HLT_g60_loose_L1EM22VHI',
        'HLT_g80_loose_L1EM22VHI',
        'HLT_g100_loose_L1EM22VHI',
        'HLT_g120_loose_L1EM22VHI',
        'HLT_g250_etcut_L1eEM26M',
        'HLT_g10_loose_L1eEM9',
        'HLT_g15_loose_L1eEM12L',
        'HLT_g20_loose_L1eEM18L',
        'HLT_g25_loose_L1eEM24L',
        'HLT_g30_loose_L1eEM24L',
        'HLT_g40_loose_L1eEM24L',
        'HLT_g50_loose_L1eEM24L',
        'HLT_g60_loose_L1eEM26M',
        'HLT_g80_loose_L1eEM26M',
        'HLT_g100_loose_L1eEM26M',
        'HLT_g120_loose_L1eEM26M',
        'HLT_g25_tight_icaloloose_L1EM20VH',
        'HLT_g25_tight_icalomedium_L1EM20VH',
        'HLT_g25_tight_icalotight_L1EM20VH',
        ]

BootstrapPhotonTriggers['Run2'] = [
        # pt_cut triggers
        'HLT_g20_etcut_L1EM12'
        # Passed through triggers for bootstrapping
        'HLT_g10_loose',
        'HLT_g15_loose_L1EM7',
        'HLT_g20_loose_L1EM12',
        'HLT_g20_loose',
        'HLT_g25_loose_L1EM15',
        'HLT_g60_loose',
        'HLT_g100_loose',
        'HLT_g120_loose',
        'HLT_g160_loose',
        'HLT_g160_loose_L1EM24VHIM',
        'HLT_g180_loose',
        'HLT_g180_loose_L1EM24VHIM',
        'HLT_g35_loose_L1EM15',
        'HLT_g40_loose_L1EM15',
        'HLT_g45_loose_L1EM15',
        'HLT_g50_loose_L1EM15',
        'HLT_g70_loose',
        'HLT_g80_loose',
        'HLT_g140_loose',
        'HLT_g200_loose',
        ]



noalgTriggers['Run2'] = [
                 'HLT_noalg_L1EM12',
                 'HLT_noalg_L1EM15',
                 'HLT_noalg_L1EM18VH',
                 'HLT_noalg_L1EM20VH',
                 'HLT_noalg_L1EM10',
                 'HLT_noalg_L1EM10VH',
                 'HLT_noalg_L1EM13VH',
                 'HLT_noalg_L1EM20VHI',
                 'HLT_noalg_L1EM22VHI',
                 'HLT_noalg_L1EM8VH',
                 'HLT_noalg_L1EM15VH',
                 'HLT_noalg_L12EM7',
                 'HLT_noalg_L12EM15']

# Additional contaienrs for photon trigger studies
ExtraContainersPhotonTrigger = {}
ExtraContainersPhotonTrigger['Run3']=[
        "HLT_egamma_Photons",
        "HLT_egamma_PhotonsAux.",
        "HLT_egamma_IsoPhotons",
        "HLT_egamma_IsoPhotonsAux.",
        "HLT_FastCaloRinger",
        "HLT_FastCaloRingerAux.",
        "HLT_FastCaloEMClusters",
        "HLT_FastCaloEMClustersAux.",
        "HLT_CaloEMClusters_Photon",
        "HLT_CaloEMClusters_PhotonAux.",
        ]

ExtraContainersPhotonTrigger['Run2'] = [
        "HLT_xAOD__PhotonContainer_egamma_Photons",
        "HLT_xAOD__PhotonContainer_egamma_PhotonsAux.",
        "HLT_xAOD__PhotonContainer_egamma_Iso_Photons",
        "HLT_xAOD__PhotonContainer_egamma_Iso_PhotonsAux.",
        "HLT_xAOD__TrigPhotonContainer_L2PhotonFex",
        "HLT_xAOD__TrigPhotonContainer_L2PhotonFexAux.",
        ]


# Additional contaienrs for electron trigger studies
ExtraContainersElectronTrigger = {}
ExtraContainersElectronTrigger['Run3'] = [
        "HLT_egamma_Electrons",
        "HLT_egamma_ElectronsAux.",
        "HLT_FastCaloRinger",
        "HLT_FastCaloRingerAux.",
        "HLT_FastCaloEMClusters",
        "HLT_FastCaloEMClustersAux.",
        "HLT_IDTrack_Electron_FTF",
        "HLT_IDTrack_Electron_FTFAux.",
        "HLT_IDTrack_ElecLRT_FTF",
        "HLT_IDTrack_ElecLRT_FTFAux.",
        "HLT_FastElectrons",
        "HLT_FastElectronsAux."
        "HLT_FastElectrons_LRT",
        "HLT_FastElectrons_LRTAux.",
        "HLT_CaloEMClusters_Electron",
        "HLT_CaloEMClusters_ElectronAux.",
        "HLT_TrigEMClusters_Electrons",
        "HLT_TrigEMClusters_ElectronsAux.",
        "HLT_TrigEMClusters_Electrons_GSF",
        "HLT_TrigEMClusters_Electrons_GSFAux."
        "HLT_IDTrack_Electron_IDTrig",
        "HLT_IDTrack_Electron_IDTrigAux.",
        "HLT_IDTrack_Electron_GSF",
        "HLT_IDTrack_Electron_GSFAux.",
        ]

ExtraContainersElectronTrigger['Run2'] = [
        "HLT_xAOD__ElectronContainer_egamma_Electrons",
        "HLT_xAOD__ElectronContainer_egamma_ElectronsAux.",
        "HLT_xAOD__TrigElectronContainer_L2ElectronFex",
        "HLT_xAOD__TrigElectronContainer_L2ElectronFexAux.",
        "HLT_xAOD__CaloClusterContainer_TrigEFCaloCalibFex",
        "HLT_xAOD__CaloClusterContainer_TrigEFCaloCalibFexAux.",
        "HLT_xAOD__TrackParticleContainer_InDetTrigTrackingxAODCnv_Electron_IDTrig",
        "HLT_xAOD__TrackParticleContainer_InDetTrigTrackingxAODCnv_Electron_IDTrigAux."
        "HLT_xAOD__TrackParticleContainer_InDetTrigTrackingxAODCnv_Electron_EFID",
        "HLT_xAOD__TrackParticleContainer_InDetTrigTrackingxAODCnv_Electron_EFIDAux.",
        ]

# Containers aimed for Navigation, Trigger Matching  and L1 RoIs
ExtraContainersTrigger = {}
ExtraContainersTrigger['Run3'] = [
        #And these are Run3 containers
        "HLT_eEMRoI",
        "HLT_eEMRoIAux.",
        "HLTNav_Summary_AODSlimmed",
        "HLTNav_Summary_AODSlimmedAux.",
        ]


ExtraContainersTrigger['Run2'] = [
        #These are Run2 Containers
        "LVL1EmTauRoIs",
        "LVL1EmTauRoIsAux.",
        "HLT_EMRoIs",
        "HLT_EMRoIsAux.",
        "HLT_TrigPassBitsCollection_passbits",
        "HLT_TrigPassBitsCollection_passbitsAux.",
        "HLT_TrigPassFlagsCollection_passflags",
        "HLT_TrigPassFlagsCollection_passflagsAux.",
        "HLT_TrigRoiDescriptorCollection_initialRoI",
        "HLT_TrigRoiDescriptorCollection_initialRoIAux.",
        ]


ExtraContainersMuonTrigger = {}
ExtraContainersMuonTrigger['Run3'] = [
        'HLT_Muons_RoI',
        'HLT_Muons_RoIAux.',
        'HLT_Muons_FS',
        'HLT_Muons_FSAux.',
        'HLT_MuonsCB_RoI',
        'HLT_MuonsCB_RoIAux.'
        'HLT_MuonsCB_LRT',
        'HLT_MuonsCB_LRTAux.'
        'HLT_MuonsCB_FS',
        'HLT_MuonsCB_FSAux.',
        ]
ExtraContainersMuonTrigger['Run2'] = [
        "HLT_xAOD__MuonContainer_MuonEFInfo",
        "HLT_xAOD__MuonContainer_MuonEFInfoAux.",
        "HLT_xAOD__MuonContainer_MuonEFInfo_FullScan",
        "HLT_xAOD__MuonContainer_MuonEFInfo_FullScanAux.",
        ]



ExtraVariablesHLTPhotons = {}
ExtraVariablesHLTPhotons['Run3'] = [
        "HLT_egamma_Photons.e.pt.m.author.Rhad.Rhad1.e277.Reta.Rphi.weta2.f1.fracs1.wtots1.weta1.DeltaE.Eratio.caloClusterLinks",
        "HLT_CaloEMClusters_Photon.calE.calEta.calPhi.calM.e_sampl.eta_sampl.etaCalo.phiCalo.ETACALOFRAME.PHICALOFRAME"
]
ExtraVariablesHLTPhotons['Run2'] = [
        "HLT_xAOD__PhotonContainer_egamma_Photons.e.pt.m.author.Rhad.Rhad1.e277.Reta.Rphi.weta2.f1.fracs1.wtots1.weta1.DeltaE.Eratio.caloClusterLinks",
        "HLT_xAOD__CaloClusterContainer_TrigEFCaloCalibFex.calE.calEta.calPhi.calM.e_sampl.eta_sampl.etaCalo.phiCalo.ETACALOFRAME.PHICALOFRAME"
]


ExtraContainersTriggerDataOnly = {}
ExtraContainersTriggerDataOnly['Run3'] = [
        ]

ExtraContainersTriggerDataOnly['Run2'] = [
        "HLT_xAOD__TrigEMClusterContainer_TrigT2CaloEgamma",
        "HLT_xAOD__TrigEMClusterContainer_TrigT2CaloEgammaAux.",
        "HLT_xAOD__CaloClusterContainer_TrigCaloClusterMaker",
        "HLT_xAOD__CaloClusterContainer_TrigCaloClusterMakerAux.",
        "HLT_xAOD__TrackParticleContainer_InDetTrigTrackingxAODCnv_Electron_FTF",
        "HLT_xAOD__TrackParticleContainer_InDetTrigTrackingxAODCnv_Electron_FTFAux.",
        "HLT_xAOD__TrackParticleContainer_InDetTrigTrackingxAODCnv_Electron_L2ID",
        "HLT_xAOD__TrackParticleContainer_InDetTrigTrackingxAODCnv_Electron_L2IDAux.",
        ]
