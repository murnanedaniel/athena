# Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
#
# Run this file in order to print out the empty slots

from TriggerMenuMT.L1.Base.L1MenuFlags import L1MenuFlags

def print_available():
    import logging
    defineMenu()
    available = list(set(range(509)) - set(L1MenuFlags.CtpIdMap.value.values()) - set([508]))
    freeItems = 512 - len(L1MenuFlags.items.value) # correct for ZB and CALREQ items
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

    L1MenuFlags.ApplyCTPLimits = False # allow to have more than 512 CTP inputs or outputs, NOT FOR P1! Only for simulation

    L1MenuFlags.items = [

        ##
        # single EM
        ##
        'L1_EM3','L1_EM7','L1_EM12', 'L1_EM8VH', 'L1_EM10VH', 'L1_EM15', 'L1_EM15VH', 'L1_EM15VHI', 'L1_EM18VHI',
        'L1_EM20VH', 'L1_EM20VHI', 'L1_EM22VH', 'L1_EM22VHI', 'L1_EM24VHI',
        'L1_EM3_EMPTY', 'L1_EM7_EMPTY', 'L1_EM7_UNPAIRED_ISO', 'L1_EM7_FIRSTEMPTY',
        'L1_EM20VH_FIRSTEMPTY',
        # new calo
        'L1_eEM3', 'L1_eEM5', 'L1_eEM7', 'L1_eEM8L',
        'L1_eEM10L', 'L1_eEM12', 'L1_eEM15',  'L1_eEM15L',  'L1_eEM15M', 
        'L1_eEM18M', 'L1_eEM20L', 'L1_eEM20VM',
        'L1_eEM7_EMPTY', 'L1_eEM7_UNPAIRED_ISO',
        'L1_eEM22', 'L1_eEM22L', 'L1_eEM22M', 'L1_eEM22T', 'L1_eEM24M',

        ## 
        # MU
        ##
        'L1_MU3V', 'L1_MU5VF', 'L1_MU8F', 'L1_MU8VF', 'L1_MU14FCH', 'L1_MU14FCHR',
        'L1_MU3VF', 'L1_MU8FC', 'L1_MU8VFC', 'L1_MU15VFCH', 'L1_MU10BOM', 'L1_MU4BOM', 'L1_MU12BOM', 'L1_MU10BO', 
        'L1_2MU3V', 'L1_2MU5VF', 'L1_2MU8F', 'L1_MU8VF_2MU5VF', 'L1_MU5VF_2MU3V',
        'L1_3MU3V', 'L1_3MU5VF', 'L1_MU5VF_3MU3V', 'L1_4MU3V',
        'L1_2MU5VF_3MU3V', 'L1_2MU8VF', 
        'L1_MU8F_2MU5VF',

        'L1_2MU14FCH_OVERLAY',
        'L1_MU3V_EMPTY', 'L1_MU5VF_EMPTY', 'L1_MU3V_FIRSTEMPTY', 'L1_MU8VF_EMPTY',
        'L1_MU3V_UNPAIRED_ISO',        

        ##
        # combined lepton (e and mu)
        ##
        'L1_2EM3', 'L1_2EM7', 'L1_2EM15', 'L1_2EM15VH', 'L1_2EM15VHI',
        'L1_2EM20VH', 'L1_EM20VH_3EM10VH',
        'L1_2EM8VH_MU8F', 'L1_EM15VH_MU8F',  'L1_EM7_MU8F',

        # Combined lepton, new calo (for ATR-24182)
        'L1_2eEM15L', 'L1_2eEM15M', 'L1_2eEM20L', 'L1_3eEM10L', 'L1_eEM20L_3eEM10L',
        'L1_eEM15L_MU8F', 'L1_2eEM8L_MU8F',

        # single tau
        'L1_TAU8', 'L1_TAU12IM',   'L1_TAU20IM',  'L1_TAU40', 'L1_TAU60', 'L1_TAU100',
        'L1_TAU8_EMPTY', 'L1_TAU8_FIRSTEMPTY', 'L1_TAU8_UNPAIRED_ISO', 'L1_TAU40_EMPTY', 'L1_TAU40_UNPAIRED_ISO',
        # new calo
        'L1_eTAU8', 
        'L1_eTAU12', 'L1_jTAU12', 'L1_jTAU20', 'L1_jTAU20M', 'L1_cTAU12M',
        'L1_eTAU12L', 'L1_eTAU12M', 
        'L1_eTAU20', 'L1_cTAU20M', 
        'L1_eTAU25', 'L1_cTAU25M',
        'L1_eTAU30HM', 
        'L1_eTAU40', 'L1_eTAU60', 'L1_eTAU100',

        # multi tau
        'L1_2TAU8', 'L1_TAU20IM_2TAU12IM', 'L1_TAU60_2TAU40',

        # multi tau, new calo
        'L1_eTAU60_2eTAU40', 
        'L1_cTAU20M_2cTAU12M_4jJ12p0ETA25',
        'L1_cTAU25M_2cTAU20M_2jJ25_3jJ20',  

        # combined tau - lepton
        'L1_EM15VHI_2TAU12IM',
        'L1_MU8F_TAU12IM',  
        'L1_MU8F_TAU12IM_J25_2J12',
        'L1_MU8F_TAU12IM_3J12',
        'L1_EM15VHI_2TAU12IM_J25_3J12',
        'L1_MU8F_TAU20IM',
        'L1_MU8F_TAU20IM_J25_2J20',

        # combined tau - jet
        'L1_TAU20IM_2TAU12IM_J25_2J20_3J12',
        'L1_TAU20IM_2TAU12IM_4J12p0ETA25',
        'L1_TAU25IM_2TAU20IM_2J25_3J20',

        # combined tau - xe
        'L1_TAU20IM_2J20_XE45',
        'L1_TAU20IM_2J20_XE50',
        'L1_EM15VHI_2TAU12IM_XE35',
        'L1_EM15VHI_2TAU12IM_4J12',
        'L1_MU8F_TAU12IM_XE35',  
        'L1_TAU20IM_2TAU12IM_XE35', 
        'L1_TAU40_2TAU12IM_XE40',  

        # combined em - jet
        'L1_EM18VHI_3J20',
        'L1_EM20VH_3J20',

        # combined mu - jet 
        'L1_MU3V_J15', 'L1_BTAG-MU3VjJ15', 'L1_BTAG-MU5VFjJ20', 

        #ATR-13743 
        'L1_MU8F_2J20','L1_MU8F_3J20', 'L1_MU8F_2J15_J20',
        'L1_MU14FCH_J40',
        'L1_MU14FCH_XE30',

        # single jet
        # L1_J12 must be kept in the menu for online monitoring of BIB (can be disabled via PS is necessary)
        'L1_J12', 'L1_J15','L1_J20','L1_J25', 'L1_J30', 'L1_J40', 'L1_J50' ,'L1_J75','L1_J85', 'L1_J100', 'L1_J400',
        'L1_J400_LAR',
        'L1_J20p31ETA49', 'L1_J30p31ETA49', 'L1_J50p31ETA49', 'L1_J75p31ETA49', 'L1_J15p31ETA49',
        'L1_J12_EMPTY','L1_J12_FIRSTEMPTY', 'L1_J12_UNPAIRED_ISO', 'L1_J12_UNPAIRED_NONISO', 'L1_J12_ABORTGAPNOTCALIB',
        'L1_J15p31ETA49_UNPAIRED_ISO',
        'L1_J30_EMPTY', 'L1_J30_FIRSTEMPTY', 'L1_J30p31ETA49_EMPTY', 'L1_J30p31ETA49_UNPAIRED_ISO', 'L1_J30p31ETA49_UNPAIRED_NONISO',
        'L1_J50_UNPAIRED_ISO', 'L1_J50_UNPAIRED_NONISO', 'L1_J50_ABORTGAPNOTCALIB',         
        'L1_J100_FIRSTEMPTY',
        'L1_J12_BGRP12',

        # jJ 
        'L1_jJ12', 'L1_jJ12p0ETA25', 'L1_jJ15', 'L1_jJ15p0ETA25', 'L1_jJ20', 'L1_jJ25',  'L1_jJ25p0ETA23',
        'L1_jJ30', 'L1_jJ30_EMPTY', 'L1_jJ30_FIRSTEMPTY', 'L1_jJ35p0ETA23', 'L1_jJ40', 'L1_jJ40p0ETA25', 'L1_jJ45p0ETA21', 'L1_jJ50',
        'L1_jJ75', 'L1_jJ85', 'L1_jJ100', 'L1_jJ120', 'L1_jJ400',
        'L1_jJ400_LAR',

        'L1_jJ15p31ETA49', 'L1_jJ20p31ETA49', 'L1_jJ30p31ETA49', 'L1_jJ50p31ETA49', 'L1_jJ75p31ETA49',

        'L1_4jJ15', 'L1_3jJ50', 'L1_4jJ20', 'L1_4jJ15p0ETA25', 'L1_5jJ15p0ETA25', 
        'L1_3jJ35p0ETA23', 'L1_jJ85_3jJ30', 
        'L1_jJ40p0ETA25_2jJ25_jJ20p31ETA49', 'L1_jJ25p0ETA23_2jJ15p31ETA49', 'L1_jJ45p0ETA21_3jJ15p0ETA25', 
        'L1_MU3V_jJ15', 

  
        # jLJ
        'L1_jLJ80', 'L1_jLJ100', 'L1_jLJ140', 'L1_jLJ160', 

        # jEM
        'L1_jEM15', 'L1_jEM15M',   

        # multi jet
        'L1_J45p0ETA21_3J15p0ETA25',
        'L1_J50_2J40p0ETA25_3J15p0ETA25',
        'L1_3J50', 'L1_4J15', 'L1_4J20',
        'L1_3J15p0ETA25_XE40',
        'L1_6J15',
        'L1_J85_3J30',

        # multi jet forward
        'L1_J25p0ETA23_2J15p31ETA49',
        'L1_J40p0ETA25_2J15p31ETA49',
        'L1_J40p0ETA25_2J25_J20p31ETA49',
        
        # multi jet central
        'L1_3J25p0ETA23',
        'L1_3J35p0ETA23',
        'L1_4J15p0ETA25',
        'L1_5J15p0ETA25', 

        # combined jet
        'L1_2J15_XE55', 'L1_J40_XE50',
        'L1_2J50_XE40', 'L1_J40_XE60',
        
        # XE
        'L1_XE35', 'L1_XE40', 'L1_XE45', 'L1_XE50', 
        'L1_XE55', 'L1_XE60', 'L1_XE30', 'L1_XE300',
        # new calo
        'L1_gXERHO30', 'L1_gXERHO50',
        'L1_gXENC30', 'L1_gXENC50',
        'L1_gXEJWOJ30', 'L1_gXEJWOJ40', 'L1_gXEJWOJ50',
        'L1_gTE200',
        'L1_gMHT500',

        'L1_jXE30', 'L1_jXE40', 'L1_jXE50', 'L1_jXE55', 'L1_jXE300', 
        'L1_jXEC50', 'L1_jTE200', 'L1_jTEC200', 'L1_jTEFWD100', 'L1_jTEFWDA100', 'L1_jTEFWDC100',
           
        # ATR-24037 
        'L1_jXEPerf50',
    
        # RNDM
        'L1_RD0_FILLED', 'L1_RD0_UNPAIRED_ISO',  'L1_RD0_EMPTY', 'L1_RD0_ABORTGAPNOTCALIB',
        'L1_RD0_FIRSTEMPTY', 'L1_RD0_BGRP11',
        'L1_RD0_BGRP7',
        'L1_RD1_EMPTY',
        'L1_RD2_EMPTY',
        'L1_RD2_FILLED',
        'L1_RD3_EMPTY',
        'L1_RD3_FILLED',

        #LUCID

        # ZDC 
        # TODO: request for pilot run (ATR-22904), to be removed for pp
        'L1_ZDC_A', 'L1_ZDC_C', 'L1_ZDC_AND',

        # VDM

        # TRT
        'L1_TRT_FILLED', 'L1_TRT_EMPTY',

        # TGC
        'L1_TGC_BURST',

        # LHCF
    
        #CALREQ
        'L1_CALREQ1',
        'L1_CALREQ2',

        # BPTX
        #'L1_BPTX0_BGRP0','L1_BPTX1_BGRP0',

        # BCM
        'L1_BCM_Wide_BGRP0', 'L1_BCM_AC_CA_BGRP0', 'L1_BCM_Wide_EMPTY', 'L1_BCM_Wide_UNPAIRED_ISO', 'L1_BCM_Wide_UNPAIRED_NONISO',
        'L1_BCM_AC_UNPAIRED_ISO','L1_BCM_CA_UNPAIRED_ISO',
        'L1_BCM_AC_UNPAIRED_NONISO','L1_BCM_CA_UNPAIRED_NONISO',
        'L1_BCM_AC_ABORTGAPNOTCALIB', 'L1_BCM_CA_ABORTGAPNOTCALIB','L1_BCM_Wide_ABORTGAPNOTCALIB',
        'L1_BCM_AC_CALIB', 'L1_BCM_CA_CALIB','L1_BCM_Wide_CALIB',
        'L1_BCM_AC_UNPAIREDB1', 'L1_BCM_CA_UNPAIREDB2',
        'L1_J12_UNPAIREDB1', 'L1_J12_UNPAIREDB2',

        # AFP
        # high-priority (all mu)
        'L1_AFP_A_AND_C_TOF_T0T1',
        'L1_AFP_FSA_BGRP0', 'L1_AFP_FSC_BGRP0',
        # med-priority (all mu)
        'L1_AFP_FSA_TOF_T0_BGRP0', 'L1_AFP_FSA_TOF_T1_BGRP0','L1_AFP_FSC_TOF_T0_BGRP0', 'L1_AFP_FSC_TOF_T1_BGRP0',
        # low-priority (all mu)
        'L1_AFP_FSA_TOF_T2_BGRP0', 'L1_AFP_FSA_TOF_T3_BGRP0','L1_AFP_FSC_TOF_T2_BGRP0', 'L1_AFP_FSC_TOF_T3_BGRP0',
        'L1_AFP_A_OR_C_UNPAIRED_ISO', 'L1_AFP_A_OR_C_UNPAIRED_NONISO', 'L1_AFP_A_OR_C_EMPTY', 'L1_AFP_A_OR_C_FIRSTEMPTY',
        'L1_AFP_A_AND_C_TOF_J20', 'L1_AFP_A_AND_C_TOF_T0T1_J20', 'L1_AFP_A_AND_C_TOF_J30', 'L1_AFP_A_AND_C_TOF_T0T1_J30', 'L1_AFP_A_AND_C_TOF_J50', 'L1_AFP_A_AND_C_TOF_T0T1_J50', 'L1_AFP_A_AND_C_TOF_J75', 'L1_AFP_A_AND_C_TOF_T0T1_J75', 

        'L1_AFP_A_AND_C_TOF_jJ20', 'L1_AFP_A_AND_C_TOF_T0T1_jJ20', 'L1_AFP_A_AND_C_TOF_jJ30', 'L1_AFP_A_AND_C_TOF_T0T1_jJ30', 'L1_AFP_A_AND_C_TOF_jJ50', 'L1_AFP_A_AND_C_TOF_T0T1_jJ50', 'L1_AFP_A_AND_C_TOF_jJ75', 'L1_AFP_A_AND_C_TOF_T0T1_jJ75',

        # high-priority (low mu)
        'L1_AFP_NSA_BGRP0', 'L1_AFP_NSC_BGRP0', 
        'L1_AFP_A','L1_AFP_C', 'L1_AFP_A_AND_C',
        #'L1_AFP_A_OR_C_J5','L1_AFP_A_AND_C_J5', # J5 not available in legacy menu. Need to update to jJ threshold for low-mu
        'L1_AFP_A_OR_C_J12','L1_AFP_A_AND_C_J12',
        'L1_MU5VF_AFP_A_OR_C','L1_MU5VF_AFP_A_AND_C','L1_EM7_AFP_A_OR_C','L1_EM7_AFP_A_AND_C',
        'L1_eEM7_AFP_A_OR_C','L1_eEM7_AFP_A_AND_C',
        # med-priority (low mu)
        'L1_AFP_A_OR_C','L1_AFP_A_OR_C_MBTS_2', 'L1_AFP_A_AND_C_MBTS_2',
          

        # MBTS
        'L1_MBTS_A', 'L1_MBTS_C',
        'L1_MBTS_1_EMPTY', 'L1_MBTS_1_1_EMPTY', 'L1_MBTS_2_EMPTY', 
        'L1_MBTS_1_UNPAIRED_ISO', 'L1_MBTS_1_1_UNPAIRED_ISO', 'L1_MBTS_2_UNPAIRED_ISO',
        'L1_MBTS_1', 'L1_MBTS_1_1', 'L1_MBTS_2',
        'L1_MBTS_4_A', 'L1_MBTS_4_C',
        'L1_MBTS_1_A', 'L1_MBTS_1_C',
        'L1_MBTS_1_A_EMPTY', 'L1_MBTS_1_C_EMPTY',

        # extra MBTS 
        # TODO: to be removed for high-mu pp        
        #'L1_MBTSA0', 'L1_MBTSA1', 'L1_MBTSA2', 'L1_MBTSA3', 'L1_MBTSA4', 'L1_MBTSA5', 'L1_MBTSA6', 'L1_MBTSA7', 'L1_MBTSA8', 'L1_MBTSA9', 'L1_MBTSA10', 'L1_MBTSA11', 'L1_MBTSA12', 'L1_MBTSA13', 'L1_MBTSA14', 'L1_MBTSA15', 'L1_MBTSC0', 'L1_MBTSC1', 'L1_MBTSC2', 'L1_MBTSC3', 'L1_MBTSC4', 'L1_MBTSC5', 'L1_MBTSC6', 'L1_MBTSC7', 'L1_MBTSC8', 'L1_MBTSC9', 'L1_MBTSC10', 'L1_MBTSC11', 'L1_MBTSC12', 'L1_MBTSC13', 'L1_MBTSC14', 'L1_MBTSC15', 

         # ZB 
        'L1_ZB',

        #ATR-21371
        # TODO: to be removed for high-mu pp
        #'L1_ALFA_ANY',
        #'L1_ALFA_ELAST15', 'L1_ALFA_ELAST18',
        #'L1_ALFA_B7L1U','L1_ALFA_B7L1L','L1_ALFA_A7L1U','L1_ALFA_A7L1L','L1_ALFA_A7R1U','L1_ALFA_A7R1L','L1_ALFA_B7R1U','L1_ALFA_B7R1L',
        #'L1_ALFA_SYST9', 'L1_ALFA_SYST10', 'L1_ALFA_SYST11', 'L1_ALFA_SYST12', 'L1_ALFA_SYST17', 'L1_ALFA_SYST18',



        #--------------------------------
        # TOPO items
        #--------------------------------

        # legacy L1Topo
        'L1_HT190-J15s5pETA21', 
        'L1_BPH-0M9-EM7-EM5_2MU3V',
        'L1_BPH-0M9-EM7-EM5_MU5VF',
        'L1_BPH-0DR3-EM7J15_2MU3V',
        'L1_BPH-0DR3-EM7J15_MU5VF',
        'L1_JPSI-1M5-EM7',
        'L1_JPSI-1M5-EM12',
        'L1_MJJ-500-NFF',
        'L1_MJJ-700',
        'L1_EM18VHI_MJJ-300',
        'L1_HT150-J20s5pETA31_MJJ-400-CF',
        'L1_LLP-RO',
        'L1_LLP-NOMATCH',
        'L1_DPHI-2EM3',
        'L1_SC111-CJ15',
        'L1_J50_DETA20-J50J',
        'L1_DR-TAU20ITAU12I',
        'L1_DR-TAU20ITAU12I-J25',
        'L1_TAU60_DR-TAU20ITAU12I', 
        'L1_BPH-0M9-EM7-EM5', 
        'L1_BPH-0DR3-EM7J15',
        'L1_LAR-ZEE',

        # phase1 L1Topo
        'L1_LAR-ZEE-eEM',
        'L1_JPSI-1M5-eEM7',
        'L1_JPSI-1M5-eEM12',
        'L1_BPH-0M9-eEM7-eEM5', 'L1_BPH-0M9-eEM7-eEM5_MU5VF', 'L1_BPH-0M9-eEM7-eEM5_2MU3V',
        'L1_BPH-0DR3-eEM7jJ15', 'L1_BPH-0DR3-eEM7jJ15_MU5VF', 'L1_BPH-0DR3-eEM7jJ15_2MU3V',
        'L1_LLP-RO-eEM', 'L1_LLP-NOMATCH-eEM',
        'L1_DPHI-2eEM3', 
        'L1_jHT150-jJ20s5pETA31_jMJJ-400-CF',
        'L1_HT190-jJ15s5pETA21',
        'L1_SC111-CjJ15',  
        'L1_jJ50_DETA20-jJ50J', 

        # tau 
        'L1_DR-eTAU20eTAU12', 'L1_DR-eTAU20eTAU12-jJ25', 
        'L1_DR-eTAU20MeTAU12M', 'L1_DR-eTAU20MeTAU12M-jJ25', 
        'L1_cTAU20M_2cTAU12M_DR-eTAU20MeTAU12M', 'L1_cTAU20M_2cTAU12M_DR-eTAU20MeTAU12M-jJ25',
        'L1_eTAU60_2cTAU12M_DR-eTAU20eTAU12', 

        #ATR-19510
        'L1_DY-BOX-2MU5VF', 'L1_DY-BOX-MU5VFMU3V', 'L1_DY-BOX-2MU3V',

        #ATR-23394
        'L1_LFV-eEM8L-MU8VF', 'L1_LFV-eEM12L-MU5VF',

        'L1_LFV-MU5VF',
        'L1_LFV-MU8VF', 

        'L1_jMJJ-700-NFF',
        'L1_jMJJ-300-NFF',
        'L1_jMJJ-500-NFF',
        'L1_jMJJ-400-CF',
        'L1_eEM18M_jMJJ-300-NFF',

        'L1_10DR-MU14FCH-MU5VF', #ATR-19376
    
        #ATR-19720, ATR-19639
        'L1_BPH-2M9-0DR15-2MU3V',
        'L1_BPH-2M9-0DR15-MU5VFMU3V',
        'L1_BPH-2M9-0DR15-C-MU5VFMU3V',
        'L1_BPH-2M9-2DR15-2MU5VF',
        'L1_BPH-8M15-0DR22-MU5VFMU3V-BO',
        'L1_BPH-8M15-0DR22-2MU5VF',
        #ATR-19355
        'L1_BPH-0M10-3MU3V',
        #ATR-19638
        'L1_BPH-0M10C-3MU3V',
        #ATR-21566 
        'L1_BPH-0DR12-2MU3V',
        'L1_BPH-7M22-0DR20-2MU3V',
        'L1_BPH-0DR12C-2MU3V',
        'L1_BPH-7M22-0DR20C-2MU3V',
        #ATR-22782
        'L1_BPH-0M16-20DR99-2MU3V', 'L1_BPH-0M16-15DR99-2MU3V',
        'L1_BPH-8M15-20DR99-C-2MU3V', 'L1_BPH-8M15-20DR99-2MU3V', 'L1_BPH-8M15-15DR99-2MU3V',          

        # INVM + DPHI 
        'L1_jMJJ-400-NFF-0DPHI22',
        'L1_jMJJ-400-NFF-0DPHI24',
        'L1_jMJJ-400-NFF-0DPHI26',

        'L1_LATE-MU8F_XE30', 'L1_LATE-MU8F_J50',

        #ATR-18824
        'L1_ZAFB-04DPHI-eEM15M',
        'L1_ZAFB-25DPHI-eEM15M',
        #ATR-22109
        #'L1_ZAFB-25DPHI-eEM18M',

        'L1_DPHI-M70-2eEM10M', 'L1_DPHI-M70-2eEM12M', #ATR-19302
        'L1_DPHI-M70-2eEM7', 'L1_DPHI-M70-2eEM7L', # ATR-21637 (no or loose shower shape cuts)
                
        #ATR-17320
        'L1_CEP-CjJ60',
        'L1_CEP-CjJ50',
        'L1_AFP_A_AND_C_TOF_CEP-CjJ60','L1_AFP_A_AND_C_TOF_T0T1_CEP-CjJ60',

        ]


# Run this file as python python/l1menu/Menu_MC_pp_v7.py to print out available IDs
# CTP ID 509-511 are reserved for CALREQ
    
    L1MenuFlags.CtpIdMap = {

        # NB: 508 is reserved for the zero bias trigger, and 509-511 for the CALREQ triggers (at the moment, ATR-22654)

    }

if __name__ == "__main__": print_available()
