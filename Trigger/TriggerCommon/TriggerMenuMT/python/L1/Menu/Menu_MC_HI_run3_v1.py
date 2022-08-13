# Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
#
# Run this file in order to print out the empty slots

from TriggerMenuMT.L1.Base.L1MenuFlags import L1MenuFlags
from TriggerMenuMT.L1.Base.Limits import Limits

def print_available():
    import logging
    defineMenu()
    available = list(set(range(Limits.MaxTrigItems-3)) - set(L1MenuFlags.CtpIdMap.value.values()) - set([508]))
    freeItems = Limits.MaxTrigItems - len(L1MenuFlags.items.value) # correct for ZB and CALREQ items
    floatingItems = sorted(list(set(L1MenuFlags.items.value) - set(L1MenuFlags.CtpIdMap.value.keys()))) # these items get their CTPID assigned automatically
    unusedItemsWithCTPID = set(L1MenuFlags.CtpIdMap.value.keys()) - set(L1MenuFlags.items.value) # this should be empty, otherwise remove the items from the CtpIdMap
    available.sort()
    logging.info("There are %d available CTP IDs: %s", len(available), ",".join(map(str,available)))
    logging.info("IDs >= 472 go in partition 2, IDs >= 492 go in partition 3")
    logging.info("There are %d free items", freeItems)
    logging.info("There are %d floating items: %s", len(floatingItems), ",".join(map(str,floatingItems)))
    logging.info("There are %d unused items with CTP ID: %s", len(unusedItemsWithCTPID), ",".join(map(str,unusedItemsWithCTPID)))

def defineMenu():

    L1MenuFlags.CTPVersion = 4 # new CTP

    L1MenuFlags.BunchGroupPartitioning = [1, 15, 15] # partition 1: 1-10, partition 2: empty (was 14), partition 3: 15 (note that BGRP0 is used by all items)
    L1MenuFlags.BunchGroupNames = ['BCRVeto', 'Paired', 'CalReq', 'Empty', 
                                   'IsolatedUnpaired', 'NonIsolatedUnpaired', 'EmptyAfterPaired', 'InTrain', 
                                   'AbortGapNotCalReq', 'VdM', 'ALFA', 'EmptyBeforePaired',
                                   'EmptyAndPaired']

    L1MenuFlags.MenuPartitioning = [0, 472, 492] # partition 1: ctpid 0-471, partition 2: ctpid 472-491, partition 3: ctpid 492-511



    L1MenuFlags.items = [

        ##
        # single EM
        ##
        'L1_EM3', 'L1_EM7', 'L1_EM8VH', 'L1_EM10', 'L1_EM10VH', 'L1_EM12', 'L1_EM14', 'L1_EM15', 'L1_EM16', 'L1_EM18VH', 'L1_EM20VH', 'L1_EM20VHI', 'L1_EM22',
        'L1_EM22VHI',
        'L1_EM3_EMPTY', 'L1_EM7_EMPTY', 'L1_EM7_UNPAIRED_ISO', 'L1_EM7_FIRSTEMPTY',
        'L1_EM20VH_FIRSTEMPTY',
        # new calo
        'L1_eEM5', 'L1_eEM9', 'L1_eEM18', 'L1_eEM15',
        'L1_eEM18L', 'L1_eEM26', 'L1_eEM26M',

        ## 
        # MU
        ##
        'L1_MU3V', 'L1_MU5VF', 'L1_MU8F', 'L1_MU8VF', 'L1_MU14FCH', 'L1_MU14FCHR',
        'L1_MU3VF', 'L1_MU8FC', 'L1_MU15VFCH', 'L1_MU10BOM',
        'L1_2MU3V', 'L1_2MU5VF', 'L1_2MU8F', 'L1_MU8VF_2MU5VF', 'L1_MU5VF_2MU3V',
        'L1_3MU3V', 'L1_3MU5VF', 'L1_MU5VF_3MU3V', 'L1_4MU3V',
        'L1_2MU5VF_3MU3V', 'L1_2MU8VF',

        'L1_2MU14FCH_OVERLAY',
        'L1_MU3V_EMPTY', 'L1_MU5VF_EMPTY', 'L1_MU3V_FIRSTEMPTY', 'L1_MU8VF_EMPTY',
        'L1_MU3V_UNPAIRED_ISO',

        ##
        # combined lepton (e and mu)
        ##
        'L1_2EM7', 'L1_2EM10', 'L1_2EM15', 'L1_2EM16',
        'L1_2EM20VH',
        # new calo
        #'L1_2eEM7', 'L1_2eEM9', 'L1_2eEM15',
        'L1_2eEM18L',

        
        # combined mu - jet
        'L1_MU3V_J12','L1_MU3V_J15', 
        'L1_MU3V_jJ30', 'L1_MU3V_jJ40',

        'L1_TAU8', 'L1_TAU60', 'L1_TAU12IM', 'L1_TAU20IM',
        'L1_TAU8_EMPTY',
        # new calo
        'L1_eTAU12',


        # single jet
        'L1_J12','L1_J15','L1_J20','L1_J25', 'L1_J30', 'L1_J40', 'L1_J50' ,'L1_J75','L1_J85', 'L1_J100',
        'L1_J20p31ETA49', 'L1_J30p31ETA49', 'L1_J50p31ETA49', 'L1_J75p31ETA49', 'L1_J15p31ETA49',
        'L1_J12_EMPTY','L1_J12_FIRSTEMPTY', 'L1_J12_UNPAIRED_ISO', 'L1_J12_UNPAIRED_NONISO',
        'L1_J15p31ETA49_UNPAIRED_ISO',
        'L1_J30_EMPTY', 'L1_J30_FIRSTEMPTY', 'L1_J30p31ETA49_EMPTY', 'L1_J30p31ETA49_UNPAIRED_ISO', 'L1_J30p31ETA49_UNPAIRED_NONISO',
        'L1_J50_UNPAIRED_ISO', 'L1_J50_UNPAIRED_NONISO',
        'L1_J100_FIRSTEMPTY',
        'L1_J12_BGRP12',
        'L1_J400', 'L1_J400_LAR',
        # di-jet
        'L1_2J15',
        # new calo
        'L1_jJ500', 'L1_jJ500_LAR',
        'L1_jJ20', 'L1_jJ30',
        'L1_jJ40', 'L1_jJ50', 'L1_jJ55', 'L1_jJ60', 'L1_jJ80', 'L1_jJ90',
        'L1_jJ40p31ETA49', 'L1_jJ50p31ETA49', 'L1_jJ60p31ETA49', 'L1_jJ90p31ETA49', 'L1_jJ125p31ETA49',


        # XE
        'L1_XE35', 'L1_XE40', 'L1_XE45', 'L1_XE50', 
        'L1_XE55', 'L1_XE60', 'L1_XE30', 'L1_XE300',
       
        'L1_J40_XE50', 'L1_J40_XE60', 
 
         # calo
        'L1_TE4', 'L1_TE20', 'L1_TE50',
        'L1_TE100', 'L1_TE200',
        'L1_TE10000', 'L1_TE12000',
        'L1_TE3p0ETA49', 'L1_TE7p0ETA49',
        'L1_TE600p0ETA49', 'L1_TE1500p0ETA49', 'L1_TE3000p0ETA49', 'L1_TE3500p0ETA49', 'L1_TE6500p0ETA49', 'L1_TE8000p0ETA49',
        'L1_TE50_VTE600p0ETA49',
        # calo overlay
        'L1_TE50_OVERLAY', 'L1_TE600p0ETA49_OVERLAY', 'L1_TE1500p0ETA49_OVERLAY', 'L1_TE3000p0ETA49_OVERLAY',
        'L1_TE3500p0ETA49_OVERLAY', 'L1_TE6500p0ETA49_OVERLAY', 'L1_TE8000p0ETA49_OVERLAY',
        
         # new calo
        'L1_gTE200',
        #UPC - MU
        'L1_MU3V_VTE50', 'L1_MU5VF_VTE50', 'L1_2MU3V_VTE50',
        
        #UPC - EM
        'L1_TAU1_TE4_VTE200', 'L1_2TAU1_VTE50',
        'L1_EM7_VTE200',
        
        #UPC - new EM
        #'L1_eEM1_TE4_VgTE200', 'L1_2eEM1_VgTE50'
        
        #UPC - calo, MBTS, calo  
        'L1_ZDC_XOR_VTE200', 'L1_VZDC_A_VZDC_C_TE5_VTE200',
        'L1_ZDC_A_VZDC_C_VTE200', 'L1_ZDC_C_VZDC_A_VTE200',
        'L1_MBTS_1_ZDC_A_VZDC_C_VTE200', 'L1_MBTS_1_ZDC_C_VZDC_A_VTE200',
        'L1_TE3p0ETA49_ZDC_A_VZDC_C_VTE200', 'L1_TE3p0ETA49_ZDC_C_VZDC_A_VTE200', 'L1_TE5_ZDC_A_VZDC_C_VTE200','L1_TE5_ZDC_C_VZDC_A_VTE200','L1_TE20_ZDC_A_VZDC_C_VTE200', 'L1_TE20_ZDC_C_VZDC_A_VTE200', 
        #UPC - calo, MBTS - legacy
        'L1_MBTS_1_VTE50',
        'L1_MBTS_1_1_VTE50',
        'L1_MBTS_1_VTE200',
        #UPC - calo, TRT - legacy
        'L1_TRT_VTE50',
        'L1_TRT_VTE200',
        #UPC - calo only - legacy
        'L1_VTE20',
        'L1_VTE50', 'L1_TE3_VTE50', 'L1_TE4_VTE50', 'L1_TE5_VTE50',
        'L1_VTE200', 'L1_TE3_VTE200', 'L1_TE5_VTE200', 'L1_TE20_VTE200', 'L1_TE50_VTE200',
        'L1_J12_VTE200',
        
        

                
        # RNDM
        'L1_RD0_FILLED', 'L1_RD0_UNPAIRED_ISO',  'L1_RD0_EMPTY',
        'L1_RD0_FIRSTEMPTY', 'L1_RD0_BGRP11',
        'L1_RD0_BGRP7',
        'L1_RD0_FIRSTINTRAIN',
        'L1_RD1_EMPTY',
        'L1_RD2_EMPTY',
        'L1_RD2_FILLED',
        'L1_RD3_EMPTY',
        'L1_RD3_FILLED',

        #LUCID

        # ZDC
        'L1_ZDC_A','L1_ZDC_C','L1_ZDC_A_C',
        'L1_ZDC_AND',
        
        # ZDC and calo
        'L1_ZDC_A_C_VTE50',
        
        # VDM

        # TRT
        'L1_TRT_FILLED', 'L1_TRT_EMPTY',

        # TGC
        'L1_TGC_BURST',

        # LHCF
    
        #CALREQ
        'L1_CALREQ1',
        'L1_CALREQ2',

        # ZB
        'L1_ZB',

        # BPTX
        
        # BCM
        'L1_BCM_Wide_BGRP12', 'L1_BCM_AC_CA_BGRP12', 'L1_BCM_Wide_EMPTY', 'L1_BCM_Wide_UNPAIRED_ISO', 'L1_BCM_Wide_UNPAIRED_NONISO',
        'L1_BCM_AC_UNPAIRED_ISO','L1_BCM_CA_UNPAIRED_ISO',
        'L1_BCM_AC_UNPAIRED_NONISO','L1_BCM_CA_UNPAIRED_NONISO',
        'L1_BCM_Wide_CALIB',
        'L1_J12_UNPAIREDB1', 'L1_J12_UNPAIREDB2',
        'L1_BCM_2A_EMPTY', 'L1_BCM_2C_EMPTY',
        'L1_BCM_2A_UNPAIRED_ISO', 'L1_BCM_2C_UNPAIRED_ISO', 'L1_BCM_2A_UNPAIRED_NONISO', 'L1_BCM_2C_UNPAIRED_NONISO',
        'L1_BCM_2A_FIRSTINTRAIN', 'L1_BCM_2C_FIRSTINTRAIN',
        # Expected to be needed later after commissioning of the BCM_2A,2C items in other BCIDs
        # 'L1_BCM_2A_UNPAIREDB1', 'L1_BCM_2A_UNPAIREDB2',
        # 'L1_BCM_2C_UNPAIREDB1', 'L1_BCM_2C_UNPAIREDB2',
        # 'L1_BCM_2A_CALIB', 'L1_BCM_2C_CALIB',


        # AFP
        'L1_EM7_AFP_A_OR_C', 'L1_EM7_AFP_A_AND_C',
        'L1_MU5VF_AFP_A_OR_C', 'L1_MU5VF_AFP_A_AND_C',
        'L1_eEM9_AFP_A_OR_C','L1_eEM9_AFP_A_AND_C',

        'L1_AFP_A_OR_C_J12', 'L1_AFP_A_AND_C_J12',
        'L1_AFP_A_OR_C_jJ20', 'L1_AFP_A_AND_C_jJ20',
        'L1_AFP_A_OR_C_jJ30', 'L1_AFP_A_AND_C_jJ30',
     
        'L1_AFP_A_AND_C_TOF_J20', 'L1_AFP_A_AND_C_TOF_T0T1_J20', 
        'L1_AFP_A_AND_C_TOF_J30', 'L1_AFP_A_AND_C_TOF_T0T1_J30',
        'L1_AFP_A_AND_C_TOF_J50', 'L1_AFP_A_AND_C_TOF_T0T1_J50',
        'L1_AFP_A_AND_C_TOF_J75', 'L1_AFP_A_AND_C_TOF_T0T1_J75',

        'L1_AFP_A_AND_C_TOF_jJ50', 'L1_AFP_A_AND_C_TOF_T0T1_jJ50', 
        'L1_AFP_A_AND_C_TOF_jJ60', 'L1_AFP_A_AND_C_TOF_T0T1_jJ60',
        'L1_AFP_A_AND_C_TOF_jJ90', 'L1_AFP_A_AND_C_TOF_T0T1_jJ90', 
        'L1_AFP_A_AND_C_TOF_jJ125', 'L1_AFP_A_AND_C_TOF_T0T1_jJ125',

        'L1_AFP_A_OR_C', 'L1_AFP_A_AND_C', 'L1_AFP_A', 'L1_AFP_C', 'L1_AFP_A_AND_C_TOF_T0T1',
        'L1_AFP_FSA_BGRP12', 'L1_AFP_FSC_BGRP12', 'L1_AFP_NSA_BGRP12', 'L1_AFP_NSC_BGRP12',
        'L1_AFP_FSA_TOF_T0_BGRP12', 'L1_AFP_FSA_TOF_T1_BGRP12', 'L1_AFP_FSA_TOF_T2_BGRP12', 'L1_AFP_FSA_TOF_T3_BGRP12',
        'L1_AFP_FSC_TOF_T0_BGRP12', 'L1_AFP_FSC_TOF_T1_BGRP12', 'L1_AFP_FSC_TOF_T2_BGRP12', 'L1_AFP_FSC_TOF_T3_BGRP12',
        'L1_AFP_A_OR_C_UNPAIRED_ISO', 'L1_AFP_A_OR_C_UNPAIRED_NONISO',
        'L1_AFP_A_OR_C_EMPTY', 'L1_AFP_A_OR_C_FIRSTEMPTY',
        'L1_AFP_A_OR_C_MBTS_2', 'L1_AFP_A_AND_C_MBTS_2',
       

        # MBTS (ATR-24701)
        'L1_MBTS_1', 'L1_MBTS_1_1',  'L1_MBTS_2',
        'L1_MBTS_2_2', 'L1_MBTS_3_3',  'L1_MBTS_4_4',
        'L1_MBTS_1_EMPTY', 'L1_MBTS_1_1_EMPTY', 'L1_MBTS_2_EMPTY',
        #'L1_MBTS_1_UNPAIRED', 'L1_MBTS_2_UNPAIRED',
        'L1_MBTS_1_UNPAIRED_ISO', 'L1_MBTS_1_1_UNPAIRED_ISO', 'L1_MBTS_2_UNPAIRED_ISO',
        'L1_MBTS_A', 'L1_MBTS_C',
        # extra MBTS
        #'L1_MBTSA0', 'L1_MBTSA1', 'L1_MBTSA2', 'L1_MBTSA3', 'L1_MBTSA4', 'L1_MBTSA5', 'L1_MBTSA6', 'L1_MBTSA7', 'L1_MBTSA8', 'L1_MBTSA9', 'L1_MBTSA10', 'L1_MBTSA11', 'L1_MBTSA12', 'L1_MBTSA13', 'L1_MBTSA14', 'L1_MBTSA15',
        #'L1_MBTSC0', 'L1_MBTSC1', 'L1_MBTSC2', 'L1_MBTSC3', 'L1_MBTSC4', 'L1_MBTSC5', 'L1_MBTSC6', 'L1_MBTSC7', 'L1_MBTSC8', 'L1_MBTSC9', 'L1_MBTSC10', 'L1_MBTSC11', 'L1_MBTSC12', 'L1_MBTSC13', 'L1_MBTSC14', 'L1_MBTSC15',

 
        #--------------------------------
        # TOPO items
        #--------------------------------

        'L1_LAR-ZEE', 'L1_LAR-ZEE-eEM',

        #ATR-17320
        'L1_CEP-CjJ100',
        'L1_CEP-CjJ90' ,

        #ATR-21371
        'L1_ALFA_ANY',
        'L1_ALFA_ELAST15', 'L1_ALFA_ELAST18',
        'L1_ALFA_B7L1U','L1_ALFA_B7L1L','L1_ALFA_A7L1U','L1_ALFA_A7L1L','L1_ALFA_A7R1U','L1_ALFA_A7R1L','L1_ALFA_B7R1U','L1_ALFA_B7R1L',
        'L1_ALFA_SYST9', 'L1_ALFA_SYST10', 'L1_ALFA_SYST11', 'L1_ALFA_SYST12', 'L1_ALFA_SYST17', 'L1_ALFA_SYST18',

        # ATR-23602
        'L1_MBTS_1_A_ALFA_C', 'L1_MBTS_1_C_ALFA_A', 'L1_EM3_ALFA_ANY', 'L1_J12_ALFA_ANY', 'L1_MU3V_ALFA_ANY', 'L1_TE5_ALFA_ANY',  
        'L1_MU3V_ALFA_EINE', 'L1_EM3_ALFA_EINE','L1_2EM3_ALFA_EINE', 'L1_J12_ALFA_EINE', 'L1_TE5_ALFA_EINE',

        'L1_TE3', 'L1_TE5', 'L1_TE10', 'L1_TE40'# also for HMT triggers
        ]



# Run this file as python python/L1/Menu_MC_HI_run3_v1.py to print out available IDs
# CTP IDs 509-511 are reserved for CALREQ
    
    L1MenuFlags.CtpIdMap = {
 

        # NB: 508 is reserved for the zero bias trigger, and 509-511 for the CALREQ triggers (at the moment, ATR-22654)

    }

if __name__ == "__main__": print_available()
