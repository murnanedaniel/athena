# Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration

from ..Base.Thresholds import MuonThreshold, EMThreshold, TauThreshold, JetThreshold, XEThreshold, MBTSThreshold, MBTSSIThreshold, NimThreshold

class ThresholdDef:

    alreadyExecuted = False

    eEMVar = {
        1 : {
            "eta_bin_boundaries": [0, 0.7, 0.8, 1.1, 1.3, 1.4, 1.5, 1.7, 2.5], # 8 bins => 9 boundaries
            "shift": [ 2, 1, 0, 0,  -1, -3, -1, 1]
        }
    }

    @staticmethod
    def addVaryingThrValues(thr, pt, shift_set):
        eta_bin_boundaries = ThresholdDef.eEMVar[shift_set]["eta_bin_boundaries"]
        shift = ThresholdDef.eEMVar[shift_set]["shift"]
        thr.addThrValue(pt+shift[0], priority=1)
        for idx,sh in reversed(list(enumerate(shift))):
            eta_min = int(10 * eta_bin_boundaries[idx])
            eta_max = int(10 * eta_bin_boundaries[idx+1])
            thr.addThrValue( pt + sh, -eta_max, -eta_min, priority=2)
        for idx,sh in enumerate(shift):
            eta_min = int(10 * eta_bin_boundaries[idx])
            eta_max = int(10 * eta_bin_boundaries[idx+1])
            thr.addThrValue( pt + sh, eta_min, eta_max, priority=2)
        return thr


    @staticmethod
    def registerThresholds(tc, menuName):

        if ThresholdDef.alreadyExecuted:
            raise RuntimeError("Calling ThresholdDef.registerThresholds twice")
        ThresholdDef.alreadyExecuted = True
 
        # MU
        # from Junpei https://indico.cern.ch/event/1011425

        MuonThreshold( "MU4"      ).setThrValue( thr=4 )                           # similar to Run-2 MU4 efficiency
        MuonThreshold( "MU4F"     ).setThrValue( thr=4 ).setTGCFlags("F")          # similar to Run-2 MU4 rate
        MuonThreshold( "MU6VF"    ).setThrValue( thr=5, ba=6 ).setTGCFlags("F")    # similar to Run-2 MU6
        MuonThreshold( "MU8F"     ).setThrValue( thr=8 ).setTGCFlags("F")          # similar to Run-2 MU10
        MuonThreshold( "MU10VF"   ).setThrValue( thr=8, ba=10 ).setTGCFlags("F")   # similar to Run-2 MU11
        MuonThreshold( "MU14FCH"  ).setThrValue( thr=14 ).setTGCFlags("F & C & H") # similar to Run-2 MU20
        MuonThreshold( "MU14FCHR" ).setThrValue( thr=14 ).setTGCFlags("F & C & H").setExclusionList("rpcFeet") # similar to Run-2 MU21

        # examples for regional muons
        MuonThreshold( 'MU20BA').setThrValue( thr=14 ).setTGCFlags("F & C & H").setRegion("BA")  # barrel only 
        MuonThreshold( 'MU20EC').setThrValue( thr=14 ).setTGCFlags("F & C & H").setRegion("EC,FW")  # forward muon

        # backward compatible threshold names (definition like above)
        MuonThreshold( "MU6" ).setThrValue( thr=5, ba=6 ).setTGCFlags("F")    # similar to Run-2 MU6
        MuonThreshold( "MU10").setThrValue( thr=8 ).setTGCFlags("F")          # similar to Run-2 MU10
        MuonThreshold( "MU11").setThrValue( thr=8, ba=10 ).setTGCFlags("F")   # similar to Run-2 MU11
        MuonThreshold( "MU15").setThrValue( thr=12).setTGCFlags("F & c & H")  # guess
        MuonThreshold( "MU20").setThrValue( thr=14 ).setTGCFlags("F & C & H") # similar to Run-2 MU20
        MuonThreshold( "MU21" ).setThrValue( thr=14 ).setTGCFlags("F & C & H").setExclusionList("rpcFeet") # similar to Run-2 MU21

        # special threshold for magnet-off menu
        MuonThreshold( "MU0").setThrValue( thr=0 )



        # EM 
        for thrV in [3, 7, 8, 10, 15, 20, 22]:
            EMThreshold('eEM%i' % thrV, 'eEM').addThrValue(thrV)

        # VH section

        EMThreshold( 'eEM8VH', 'eEM').setIsolation( rhad = "Tight" )\
            .addThrValue(9, priority=1)\
            .addThrValue(9, -8, 8, priority=2)\
            .addThrValue(7, -11, -9, priority=2).addThrValue(7, 9, 11, priority=2)\
            .addThrValue(6, -14, -12, priority=2).addThrValue(6, 12, 14, priority=2)\
            .addThrValue(5, -15, -15, priority=2).addThrValue(5, 15, 15, priority=2)\
            .addThrValue(7, -18, -16, priority=2).addThrValue(7, 16, 18, priority=2)\
            .addThrValue(8, -25, -19, priority=2).addThrValue(8, 19, 25, priority=2)
        
        EMThreshold( 'eEM10VH', 'eEM' ).setIsolation( rhad = "Medium" )\
            .addThrValue(11, priority=1)\
            .addThrValue(11, -8, 8, priority=2)\
            .addThrValue(9, -11, -9, priority=2).addThrValue(9, 9, 11, priority=2)\
            .addThrValue(8, -14, -12, priority=2).addThrValue(8, 12, 14, priority=2)\
            .addThrValue(7, -15, -15, priority=2).addThrValue(7, 15, 15, priority=2)\
            .addThrValue(9, -18, -16, priority=2).addThrValue(9, 16, 18, priority=2)\
            .addThrValue(10, -25, -19, priority=2).addThrValue(10, 19, 25, priority=2)

        EMThreshold( 'eEM15VH', 'eEM').setIsolation( rhad = "Loose" )\
            .addThrValue(17, priority=1)\
            .addThrValue(17, -7, 7, priority=2)\
            .addThrValue(16, -9, -8, priority=2).addThrValue(16, 8, 9, priority=2)\
            .addThrValue(15, -12, -10, priority=2).addThrValue(15, 10, 12, priority=2)\
            .addThrValue(14, -14, -13, priority=2).addThrValue(14, 13, 14, priority=2)\
            .addThrValue(13, -15, -15, priority=2).addThrValue(13, 15, 15, priority=2)\
            .addThrValue(15, -17, -16, priority=2).addThrValue(15, 16, 17, priority=2)\
            .addThrValue(16, -25, -18, priority=2).addThrValue(16, 18, 25, priority=2)  
      
        EMThreshold( 'eEM20VH', 'eEM').setIsolation( rhad = "Loose" )\
            .addThrValue(22, priority=1)\
            .addThrValue(22, -7, 7, priority=2)\
            .addThrValue(21, -8, -8, priority=2).addThrValue(21, 8, 8, priority=2)\
            .addThrValue(20, -11, -9, priority=2).addThrValue(20, 9, 11, priority=2)\
            .addThrValue(19, -13, -12, priority=2).addThrValue(19, 12, 13, priority=2)\
            .addThrValue(18, -14, -14, priority=2).addThrValue(18, 14, 14, priority=2)\
            .addThrValue(17, -15, -15, priority=2).addThrValue(17, 15, 15, priority=2)\
            .addThrValue(19, -17, -16, priority=2).addThrValue(19, 16, 17, priority=2)\
            .addThrValue(21, -25, -18, priority=2).addThrValue(21, 18, 25, priority=2)       

        # (V)HI section

        EMThreshold( 'eEM15VHI', 'eEM').setIsolation( reta = "Loose", wstot = "Medium" )\
            .addThrValue(17, priority=1)\
            .addThrValue(17, -7, 7, priority=2)\
            .addThrValue(16, -9, -8, priority=2).addThrValue(16, 8, 9, priority=2)\
            .addThrValue(15, -12, -10, priority=2).addThrValue(15, 10, 12, priority=2)\
            .addThrValue(14, -14, -13, priority=2).addThrValue(14, 13, 14, priority=2)\
            .addThrValue(13, -15, -15, priority=2).addThrValue(13, 15, 15, priority=2)\
            .addThrValue(15, -17, -16, priority=2).addThrValue(15, 16, 17, priority=2)\
            .addThrValue(16, -25, -18, priority=2).addThrValue(16, 18, 25, priority=2)

        EMThreshold( 'eEM18VHI', 'eEM').setIsolation( reta = "Loose", wstot = "Medium" )\
            .addThrValue(20, priority=1)\
            .addThrValue(20, -7, 7, priority=2)\
            .addThrValue(19, -8, -8, priority=2).addThrValue(19, 8, 8, priority=2)\
            .addThrValue(18, -11, -9, priority=2).addThrValue(18, 9, 11, priority=2)\
            .addThrValue(17, -13, -12, priority=2).addThrValue(17, 12, 13, priority=2)\
            .addThrValue(16, -14, -14, priority=2).addThrValue(16, 14, 14, priority=2)\
            .addThrValue(15, -15, -15, priority=2).addThrValue(15, 15, 15, priority=2)\
            .addThrValue(17, -17, -16, priority=2).addThrValue(17, 16, 17, priority=2)\
            .addThrValue(19, -25, -18, priority=2).addThrValue(19, 18, 25, priority=2)

        EMThreshold( 'eEM20VHI', 'eEM').setIsolation( reta = "Loose", wstot = "Medium" )\
            .addThrValue(22, priority=1)\
            .addThrValue(22, -7, 7, priority=2)\
            .addThrValue(21, -8, -8, priority=2).addThrValue(21, 8, 8, priority=2)\
            .addThrValue(20, -11, -9, priority=2).addThrValue(20, 9, 11, priority=2)\
            .addThrValue(19, -13, -12, priority=2).addThrValue(19, 12, 13, priority=2)\
            .addThrValue(18, -14, -14, priority=2).addThrValue(18, 14, 14, priority=2)\
            .addThrValue(17, -15, -15, priority=2).addThrValue(17, 15, 15, priority=2)\
            .addThrValue(19, -17, -16, priority=2).addThrValue(19, 16, 17, priority=2)\
            .addThrValue(21, -25, -18, priority=2).addThrValue(21, 18, 25, priority=2)

        ThresholdDef.addVaryingThrValues(EMThreshold('eEM22VHI', 'eEM').setIsolation(reta="Loose", wstot="Medium"), pt=22, shift_set=1)
        # EMThreshold( 'eEM22VHI', 'eEM').setIsolation( reta = "Loose", wstot = "Medium" )\
        #     .addThrValue(24, priority=1)\
        #     .addThrValue(24, -7, 7, priority=2)\
        #     .addThrValue(23, -8, -8, priority=2).addThrValue(23, 8, 8, priority=2)\
        #     .addThrValue(22, -11, -9, priority=2).addThrValue(22, 9, 11, priority=2)\
        #     .addThrValue(21, -13, -12, priority=2).addThrValue(21, 12, 13, priority=2)\
        #     .addThrValue(20, -14, -14, priority=2).addThrValue(20, 14, 14, priority=2)\
        #     .addThrValue(19, -15, -15, priority=2).addThrValue(19, 15, 15, priority=2)\
        #     .addThrValue(21, -17, -16, priority=2).addThrValue(21, 16, 17, priority=2)\
        #     .addThrValue(23, -25, -18, priority=2).addThrValue(23, 18, 25, priority=2)
        


        # TAU
        for et in [12, 20, 40, 60, 100]:
            TauThreshold('eTAU%i' % et, 'eTAU').setEt(et)

        for et in [12,20, 25]:
            TauThreshold('eTAU%iIM' % et, 'eTAU').setEt(et)

        # JET
        for thrV in [12, 15, 20, 25, 30, 40, 50, 85, 100]:
            JetThreshold('jJ%i' % thrV, 'jJ').setPt(thrV).addRange(etamin=-31, etamax=31) # jets are between -31 and 31 -ATR-11526

        # Central jet
        for (thrV, etamax) in [(12,25), (15,25), (25,23), (35,23), (40,25)]:
            JetThreshold('jJ%ip0ETA%i'  % (thrV, etamax), 'jJ').setPt(thrV).addRange(etamin = -etamax,  etamax = etamax)  

        # Standard forward jet
        for thrV in [15, 20, 75]:
            JetThreshold('jJ%ip31ETA49' % thrV, 'jJ').setPt(thrV).addRange(etamin=31, etamax=49).addRange(etamin=-49, etamax=-31)

        # XE
        for thrV in [20, 50]:
            XEThreshold('gXEPUFIT%i' % thrV, 'gXE').setXE(thrV)

        for thrV in [20, 30, 35, 40, 45, 50]:
            XEThreshold('gXERHO%i' % thrV, 'gXE').setXE(thrV)

        for thrV in [50]:
            XEThreshold('gXE%i' % thrV, 'gXE').setXE(thrV)

        for thrV in [50]:
            XEThreshold('jXE%i' % thrV, 'jXE').setXE(thrV)




        # CALREQ
            
        for i in range(3):
            NimThreshold('CAL%i' % i, 'CALREQ', mapping=i)


        ## MBTS

        # MBTS naming scheme defined in
        # https://docs.google.com/spreadsheets/d/1R0s8Lw-0nPSjqe9YTuZBCeAdedn_Ve4Ax6bbMe_4bSk/edit#gid=1818891632

        # run 1
        thresholdA=[ 32.04, 26.98, 35.00, 33.54, 32.08, 36.46, 30.63, 32.08, 33.54, 30.63, 29.17, 33.54, 32.08, 32.08, 30.63, 26.25]
        thresholdC=[ 55.42, 31.98, 32.81, 49.48, 98.44, 32.11, 32.62, 29.90, 24.06, 25.81, 25.52, 35.00, 27.71, 36.46, 26.25, 30.63]

        # run 2 above MBTS_A08 only the even numbers are used
        thresholdA=[ 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
        thresholdC=[ 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]

        for i, (vA, vC) in enumerate(zip(thresholdA, thresholdC)):
            MBTSSIThreshold('MBTS_A%i' % i).setVoltage(vA)
            MBTSSIThreshold('MBTS_C%i' % i).setVoltage(vC)

        thr_mbtsA = MBTSThreshold('MBTS_A', mapping=0)
        thr_mbtsC = MBTSThreshold('MBTS_C', mapping=1)
        for x in range(16):
            if tc.thresholdExists('MBTS_A%i' % x):
                thr_mbtsA.addSector( tc.getDefinedThreshold('MBTS_A%i' % x) )
            if tc.thresholdExists('MBTS_C%i' % x):
                thr_mbtsC.addSector( tc.getDefinedThreshold('MBTS_C%i' % x) )



        ## ZDC
        
        NimThreshold('ZDC_A',   'ZDC')
        NimThreshold('ZDC_C',   'ZDC')
        NimThreshold('ZDC_AND', 'ZDC')


        ## BCM

        NimThreshold('BCM_AtoC', 'BCM')
        NimThreshold('BCM_CtoA', 'BCM')
        NimThreshold('BCM_Wide', 'BCM')
        NimThreshold('BCM_Comb', 'BCMCMB')
        NimThreshold('BCM_MCA',  'BCM')
        NimThreshold('BCM_MCC',  'BCM')
        NimThreshold('BCM_X',    'BCM')

        ## LUCID

        NimThreshold('LUCID_A', 'LUCID')
        NimThreshold('LUCID_C', 'LUCID')
        NimThreshold('LUCID_Coinc_AC', 'LUCID')
        NimThreshold('LUCID_COMM', 'LUCID')
        NimThreshold('LUCID_05', 'LUCID')
        NimThreshold('LUCID_06', 'LUCID')

        ## AFP
        NimThreshold('AFP_NSC', 'NIM', mapping=2)        
        NimThreshold('AFP_NSA', 'NIM', mapping=3)        
        NimThreshold('AFP_FSA_SIT', 'NIM', mapping=4)        
        NimThreshold('AFP_FSA_TOF', 'NIM', mapping=5)  
        NimThreshold('AFP_FSA_LOG', 'NIM', mapping=6)        
        NimThreshold('AFP_FSC_SIT', 'NIM', mapping=7)        
        NimThreshold('AFP_FSC_LOG', 'NIM', mapping=8)        
        NimThreshold('AFP_FSC_TOF', 'NIM', mapping=9)  


        ## BPTX
        NimThreshold('BPTX0', 'BPTX')
        NimThreshold('BPTX1', 'BPTX')


        ## Other NIMs
        NimThreshold('NIML1A',  'NIM', mapping=0)
        NimThreshold('NIMLHCF', 'NIM', mapping=1)
        NimThreshold('NIMTGC',  'NIM', mapping=12)
        NimThreshold('NIMRPC',  'NIM', mapping=13)
        NimThreshold('NIMTRT',  'NIM', mapping=14)

        # ALFA
        LUT1offset =  2 # this is needed because the first 4 direct inputs are in a LUT with 8 PITs so the OR with the next inputs would not work
        LUT2offset =  8
        LUT3offset = 14
        LUT4offset = 20
        LUT5offset = 26
        for i, alfa in enumerate( ['B7R1L', 'B7R1U', 'A7R1L', 'A7R1U', 'A7L1L', 'A7L1U', 'B7L1L', 'B7L1U'] ):
            phaseOffset = 32 * (i%2)
            NimThreshold('ALFA_%s'    % alfa, 'ALFA', mapping = LUT1offset + i/2 + phaseOffset )
            NimThreshold('ALFA2_%s'   % alfa, 'ALFA', mapping = LUT2offset + i/2 + phaseOffset )
            NimThreshold('ALFA3_%s'   % alfa, 'ALFA', mapping = LUT3offset + i/2 + phaseOffset )
            NimThreshold('ALFA4_%s'   % alfa, 'ALFA', mapping = LUT4offset + i/2 + phaseOffset )
            NimThreshold('ALFA_%s_OD' % alfa, 'ALFA', mapping = LUT5offset + i/2 + phaseOffset )
