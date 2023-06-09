# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration


class TauHypoProvider:
    # Provide a centralized way to store hypotheses for tau trigger chains
    # Avoid duplication of hypotheses
    # Avoid declaration of hypotheses which are not needed

    # The dictionary containing the hypos
    hypoStore = {}
    
    # Define it locally for simplicity
    GeV = 1000.0
    
    def GetHypo(self, level, threshold, criteria, part, strategy):

        hypoParameters = [level,str(threshold),criteria,part,strategy];
        hypoName = '_'.join(hypoParameters)
        
        # Check if the access key is in the dictionary
        if hypoName in self.hypoStore:
            return self.hypoStore[hypoName]

        # Otherwise, allocate it and return it
        if strategy == 'r1':
            
            if level == 'L2':

                currentHypoKey = 'l2'+part+'_tau'+threshold+'_'+criteria+'_r1'
                
                if part == 'calo':
                    if criteria == 'perf':
                        from TrigTauHypo.TrigTauHypoBase import T2CaloTauHypo_tauNoCut
                        currentHypo = T2CaloTauHypo_tauNoCut(currentHypoKey.replace(threshold, ''))
                    else:
                        from TrigTauHypo.TrigTauHypoBase import T2CaloTauHypo
                        theVars = ['EtRawMin', 'CoreFractionMin']
                        theThresh = self.thresholdsL2Calo[(criteria, int(threshold))]
                        currentHypo = T2CaloTauHypo(currentHypoKey, theVars, theThresh)
                    

                if part == 'id':
                    if criteria== 'perf':
                        from TrigTauHypo.TrigTauHypoBase import T2IDTauHypo_tauNoCut
                        currentHypo = T2IDTauHypo_tauNoCut(currentHypoKey.replace(threshold, ''))
                    else:
                        from TrigTauHypo.TrigTauHypoBase import T2IDTauHypo
                        theVars = ['NTrkMax', 'SumPtRatioMax1P', 'SumPtRatioMaxMP']
                        theThresh = self.thresholdsL2ID[(criteria, int(threshold))]
                        currentHypo = T2IDTauHypo(currentHypoKey, theVars, theThresh)
    
                if part == '':
                    if criteria== 'perf':
                        from TrigTauHypo.TrigTauHypoBase import T2TauHypo_tauNoCut
                        currentHypo = T2TauHypo_tauNoCut(currentHypoKey.replace(threshold, ''))
                    else:
                        from TrigTauHypo.TrigTauHypoBase import T2TauHypo
                        theVars = ['EtOverPtLeadTrkMax1P', 'EtOverPtLeadTrkMaxMP', 'TrkAvgDistMax1P', 'TrkAvgDistMaxMP']
                        theThresh = self.thresholdsL2[(criteria, int(threshold))]
                        currentHypo = T2TauHypo(currentHypoKey, theVars, theThresh)

            if level == 'EF':

                currentHypoKey = 'ef'+part+'_tau'+threshold+'_'+criteria+'_r1'

                if part == '':
                    if criteria== 'perf' or criteria== 'perf0' or criteria== 'cosmic': 
                        from TrigTauHypo.TrigTauHypoConfig2012 import EFTauMVHypo_tauNoCut
                        currentHypo = EFTauMVHypo_tauNoCut(currentHypoKey.replace(threshold, ''))

                    elif criteria=='dikaon' or criteria=='dikaontight' or criteria=='dikaonmass' or criteria=='dikaonmasstight' or criteria=='kaonpi1' or criteria=='kaonpi2' or criteria=='dipion1' or criteria=='dipion1loose' or criteria=='dipion2':
                        from TrigTauHypo.TrigTauHypoConfig2012 import EFTauDiKaonHypo
                        theVars = ['massTrkSysMin', 'massTrkSysMax', 'massTrkSysKaonMin', 'massTrkSysKaonMax', 'massTrkSysKaonPiMin', 'massTrkSysKaonPiMax', 'targetMassTrkSysKaonPi', 'leadTrkPtMin','EtCalibMin','EMPOverTrkSysPMax']
                        theThresh = self.thresholdsEF_dikaon[(criteria, int(threshold))]
                        currentHypo = EFTauDiKaonHypo(currentHypoKey, theVars, theThresh)

#                    elif criteria=='dikaonmass' or criteria=='dikaonmasstight':
#                        from TrigTauHypo.TrigTauHypoConfig2012 import EFTauDiKaonHypo
#                        theVars = ['massTrkSysKaonMin', 'massTrkSysKaonMax', 'leadTrkPtMin','EtCalibMin','EMPOverTrkSysPMax','nWideTrackMax']
#                        theThresh = self.thresholdsEF_dikaonmass[(criteria, int(threshold))]
#                        currentHypo = EFTauDiKaonHypo(currentHypoKey, theVars, theThresh)

                    elif criteria=='singlepion' or criteria=='singlepiontight':
                        from TrigTauHypo.TrigTauHypoConfig2012 import EFTauDiKaonHypo
                        theVars = ['leadTrkPtMin','EtCalibMin','nTrackMax','nWideTrackMax','dRmaxMax','etOverPtLeadTrkMin','etOverPtLeadTrkMax']
                        theThresh = self.thresholdsEF_singlepion[(criteria, int(threshold))]
                        currentHypo = EFTauDiKaonHypo(currentHypoKey, theVars, theThresh)

                    elif criteria=='medium1HighptL' or criteria=='medium1HighptM' or criteria=='medium1HighptH':
                        from TrigTauHypo.TrigTauHypoConfig2012 import EFTauMVHypo_highpt
                        theVars = ['NTrackMax', 'EtCalibMin', 'Level','HighptTrkThr','HighptIDThr','HighptJetThr']
                        theThresh = self.thresholdsEF[(criteria, int(threshold))]
                        theThresh.extend(self.thresholdsHighpt[(criteria)])
                        currentHypo = EFTauMVHypo_highpt(currentHypoKey, theVars, theThresh)

                    elif criteria=='medium0':
                        from TrigTauHypo.TrigTauHypoConfig2012 import EFTauMVHypo
                        theVars = ['NTrackMin','NTrackMax', 'EtCalibMin', 'Level','ApplyIDon0p']
                        theThresh = self.thresholdsEF_FTK[(criteria, int(threshold))]
                        currentHypo = EFTauMVHypo(currentHypoKey, theVars, theThresh)

                    else:
                        from TrigTauHypo.TrigTauHypoConfig2012 import EFTauMVHypo
                        theVars = ['NTrackMax', 'EtCalibMin', 'Level']
                        theThresh = self.thresholdsEF[(criteria, int(threshold))]
                        currentHypo = EFTauMVHypo(currentHypoKey, theVars, theThresh)

        if strategy == 'calo' or strategy =='ptonly' or strategy == 'mvonly' or strategy == 'caloonly' or strategy == 'track' or strategy == 'trackonly' or strategy == 'tracktwo' or strategy == 'trackcalo' or strategy == 'tracktwocalo' or strategy == 'tracktwo2015' or strategy == 'FTK' or strategy == 'FTKRefit' or strategy == 'FTKNoPrec':

            # Simple implementation of 2015 pre-selection
            currentHypoKey = 'l2'+part+'_tau'+threshold+'_'+criteria+'_'+strategy

            # Re-define the calo part using the generic hypo
            if part == 'calo':
                from TrigTauHypo.TrigTauHypoConf import HLTTauGenericHypo
                currentHypo = HLTTauGenericHypo(currentHypoKey)

                # pT cut always defined: ugly string-to-int-to-string conversion, sorry :)
                myThreshold = str(int(threshold)*self.GeV)
                theDetails  = [int(-1)]
                theFormulas = ['y > '+myThreshold]

                if strategy =='calo' or strategy == 'caloonly' or strategy == 'tracktwocalo' or strategy == 'trackcalo':
                    # centFrac cut (detail #24: 2nd order fit, turn-off at ~ 55 GeV, 95% efficiency)
                    theDetails += [24]
                    theFormulas += ['x > (0.945 - (1.26e-05*TMath::Min(y, 50000.)) + (1.05e-10*TMath::Min(y, 50000.)*TMath::Min(y, 50000.)))']
                    # centFrac cut (detail #24: 2nd order fit, turn-off at 50 GeV, 90% efficiency)
                    # theFormulas += ['x > TMath::Max(0.936 - (8.5e-06*y) + (6.54e-11*y*y), 0.660)']
                currentHypo.Details = theDetails
                currentHypo.Formulas = theFormulas
                
            #if part == 'calo':
            #    from TrigTauHypo.TrigTauHypoBase import HLTCaloTauHypo
            #    theVars = ['LowerPtCut', 'UseCellCut', 'CoreFractionCut', 'HadRadiusCut']
            #    if strategy == 'calo' or strategy == 'caloonly':
            #        theThresh = [int(threshold)*self.GeV, 1, 0.63, 0.8]
            #    if strategy == 'ptonly' or strategy == 'trackonly' or strategy == 'track' or strategy == 'tracktwo':
            #        theThresh = [int(threshold)*self.GeV, 0, 0.0, 0.8]
            #    if strategy == 'mvonly':
            #        theThresh = [0, 0, 0.0, 0.8]
            #    currentHypo = HLTCaloTauHypo(currentHypoKey, theVars, theThresh)

            if part == 'id':
                from TrigTauHypo.TrigTauHypoBase import HLTTrackTauHypo
                from TrigTauHypo.TrigTauHypoConfig2012 import EFTauMVHypo
                from TrigTauHypo.TrigTauHypoConfig2012 import EFTauMVHypo_highpt
                # Important Note: the pT cut here is an unused dummy
                if criteria == 'cosmic':
                    theVars = ['LowerPtCut', 'TracksInCoreCut', 'TracksInIsoCut', 'DeltaZ0Cut']
                    theThresh = [int(threshold)*self.GeV, 9999, 9999, 9999.]
                    currentHypo = HLTTrackTauHypo(currentHypoKey, theVars, theThresh)
                # accept 0p taus in FTK chains
                elif criteria == 'perf0' and (strategy == 'FTK' or strategy == 'FTKNoPrec'):
                    theVars = ['NTrackMin','NTrackMax','NWideTrackMax','EtCalibMin', 'Level','Method']
                    theThresh = [0,3,1,0.*self.GeV,-1111,0]
                    currentHypo = EFTauMVHypo(currentHypoKey, theVars, theThresh)
                else:
                    if strategy != 'tracktwo' and strategy != 'FTK' and strategy != 'FTKRefit' and strategy != 'FTKNoPrec':
                        theVars = ['LowerPtCut','LowerTrackPtCut']
                        theThresh = [int(threshold)*self.GeV,1.*self.GeV]
                        currentHypo = HLTTrackTauHypo(currentHypoKey, theVars, theThresh)
                    else:
                        theVars = ['NTrackMin','NTrackMax','NWideTrackMax','EtCalibMin', 'Level','Method']
                        theThresh = [1,3,1,0.*self.GeV,-1111,0]
                        if criteria=='medium1HighptL' or criteria=='medium1HighptM' or criteria=='medium1HighptH':
                            theVars.extend(['HighptTrkThr','HighptIDThr','HighptJetThr'])
                            theThresh.extend(self.thresholdsHighpt[(criteria)])
                            currentHypo = EFTauMVHypo_highpt(currentHypoKey, theVars, theThresh)
                        else:
                            currentHypo = EFTauMVHypo(currentHypoKey, theVars, theThresh)

        assert currentHypo, 'unable to find hypothesis algorithm: '+currentHypoKey
            
        self.hypoStore[hypoName] = currentHypo
        return currentHypo



    # Here we need a large repository of configuration values, defined as class members
    thresholdsL2Calo = {
        ('loose1', 20) : [15199.0, 0.8000],   
        ('loose1', 29) : [20696.0, 0.8000],   
        ('loose1', 70) : [50000.0, 0.8000],   
        ('loose1', 100): [80000.0, 0.8000],
        ('medium', 20): [ 15199.0, 0.7485],
        ('medium', 29): [ 20696.0, 0.7755],
        ('medium', 38): [ 29430.0, 0.8175],
        ('medium', 50): [ 40341.0, 0.8310],
        ('medium', 115): [100000.0, 0.8790],
        ('medium', 125): [109233.0, 0.8790],
        ('medium1', 20): [ 15199.0, 0.7485],
        ('medium1', 29): [ 20696.0, 0.7755],
        ('medium1', 38): [ 29430.0, 0.8175],
        ('medium1', 50): [ 40341.0, 0.8310],
        ('medium1', 115): [100000.0, 0.8790],
        ('medium1', 125): [109233.0, 0.8790] }
    
    thresholdsL2ID = {
        ('loose1', 20): [4, 0.001, 0.001],
        ('loose1', 29): [4, 0.001, 0.001],
        ('loose1', 70): [4, 0.002, 0.002],
        ('loose1', 100):[4, 0.002, 0.002],    
        ('medium', 20): [7, 0.001, 0.099],
        ('medium', 29): [7, 0.001, 0.001],
        ('medium', 38): [7, 0.001, 0.001],
        ('medium', 50): [7, 0.001, 0.001],
        ('medium', 115):[7, 0.001, 0.001],
        ('medium', 125):[7, 0.001, 0.001],
        ('medium1', 20): [4, 0.001, 0.099],
        ('medium1', 29): [4, 0.001, 0.001],
        ('medium1', 38): [4, 0.001, 0.001],
        ('medium1', 50): [4, 0.001, 0.001],
        ('medium1', 115):[4, 0.001, 0.001],
        ('medium1', 125):[4, 0.001, 0.001] }

    thresholdsL2 = {
        ('loose1', 20): [999.99, 999.99, 999.99, 999.99],
        ('loose1', 29): [999.99, 999.99, 999.99, 999.99],
        ('loose1', 70): [999.99, 999.99, 999.99, 999.99],
        ('loose1', 100): [999.99, 999.99, 999.99, 999.99],
        ('medium', 20): [  8.84, 2.84, 0.105, 0.073],
        ('medium', 29): [  9.72, 2.92, 0.086, 0.062],
        ('medium', 38): [  9.68, 3.04, 0.076, 0.059],
        ('medium', 50): [  8.22, 3.02, 0.062, 0.054],
        ('medium', 115): [999.99, 4.38, 0.057, 0.053],
        ('medium', 125): [999.99, 4.38, 0.057, 0.053],
        ('medium1', 20): [  8.84, 2.84, 0.105, 0.073],
        ('medium1', 29): [  9.72, 2.92, 0.086, 0.062],
        ('medium1', 38): [  9.68, 3.04, 0.076, 0.059],
        ('medium1', 50): [  8.22, 3.02, 0.062, 0.054],
        ('medium1', 115): [999.99, 4.38, 0.057, 0.053],
        ('medium1', 125): [999.99, 4.38, 0.057, 0.053] }

    thresholdsEF = {
        ('medium', 20): [6,  20.0*GeV, 1],
        ('medium', 25): [6,  25.0*GeV, 1],
        ('medium', 29): [6,  29.0*GeV, 1],
        ('medium', 38): [6,  38.0*GeV, 1],
        ('medium', 50): [6,  50.0*GeV, 1],
        ('medium', 115): [6, 115.0*GeV, 1],
        ('medium', 125): [6, 125.0*GeV, 1],
        ('loose1', 20): [3,  20.0*GeV, 1],
        ('loose1', 25): [3,  25.0*GeV, 1],
        ('loose1', 29): [3,  29.0*GeV, 1],
        ('loose1', 35): [3,  35.0*GeV, 1],
        ('loose1', 38): [3,  38.0*GeV, 1],
        ('loose1', 50): [3,  50.0*GeV, 1],
        ('loose1', 60): [3,  60.0*GeV, 1],
        ('loose1', 80): [3,  80.0*GeV, 1],
        ('loose1', 115): [3, 115.0*GeV, 1],
        ('loose1', 125): [3, 125.0*GeV, 1], 
        ('loose1', 160): [3, 160.0*GeV, 1],
        ('loose1', 200): [3, 200.0*GeV, 1],
        ('medium1', 0): [3,  0.0*GeV, 2], 
        ('medium1', 12): [3,  12.0*GeV, 2],
        ('medium1', 20): [3,  20.0*GeV, 2],
        ('medium1', 25): [3,  25.0*GeV, 2],
        ('medium1', 29): [3,  29.0*GeV, 2],
        ('medium1', 35): [3,  35.0*GeV, 2],
        ('medium1', 38): [3,  38.0*GeV, 2],
        ('medium1', 50): [3,  50.0*GeV, 2],
        ('medium1', 60): [3,  60.0*GeV, 2],
        ('medium1', 80): [3,  80.0*GeV, 2],
        ('medium1', 115): [3, 115.0*GeV, 2],
        ('medium1', 125): [3, 125.0*GeV, 2], 
        ('medium1', 160): [3, 160.0*GeV, 2], 
        ('medium1', 200): [3, 200.0*GeV, 2],
        ('tight1', 20): [3,  20.0*GeV, 3],
        ('tight1', 25): [3,  25.0*GeV, 3],
        ('tight1', 29): [3,  29.0*GeV, 3],
        ('tight1', 35): [3,  35.0*GeV, 3],
        ('tight1', 38): [3,  38.0*GeV, 3],
        ('tight1', 50): [3,  50.0*GeV, 3],
        ('tight1', 60): [3,  60.0*GeV, 3],
        ('tight1', 80): [3,  80.0*GeV, 3],
        ('tight1', 115): [3, 115.0*GeV, 3],
        ('tight1', 125): [3, 125.0*GeV, 3], 
        ('tight1', 160): [3, 160.0*GeV, 3],
        ('tight1', 200): [3, 200.0*GeV, 3],
        ('medium1HighptL', 20): [3,  20.0*GeV, 2],
        ('medium1HighptL', 25): [3,  25.0*GeV, 2],
        ('medium1HighptL', 29): [3,  29.0*GeV, 2],
        ('medium1HighptL', 35): [3,  35.0*GeV, 2],
        ('medium1HighptL', 38): [3,  38.0*GeV, 2],
        ('medium1HighptL', 50): [3,  50.0*GeV, 2],
        ('medium1HighptL', 60): [3,  60.0*GeV, 2],
        ('medium1HighptL', 80): [3,  80.0*GeV, 2],
        ('medium1HighptL', 115): [3, 115.0*GeV, 2],
        ('medium1HighptL', 125): [3, 125.0*GeV, 2],
        ('medium1HighptL', 160): [3, 160.0*GeV, 2],
        ('medium1HighptM', 20): [3,  20.0*GeV, 2],
        ('medium1HighptM', 25): [3,  25.0*GeV, 2],
        ('medium1HighptM', 29): [3,  29.0*GeV, 2],
        ('medium1HighptM', 35): [3,  35.0*GeV, 2],
        ('medium1HighptM', 38): [3,  38.0*GeV, 2],
        ('medium1HighptM', 50): [3,  50.0*GeV, 2],
        ('medium1HighptM', 60): [3,  60.0*GeV, 2],
        ('medium1HighptM', 80): [3,  80.0*GeV, 2],
        ('medium1HighptM', 115): [3, 115.0*GeV, 2],
        ('medium1HighptM', 125): [3, 125.0*GeV, 2],
        ('medium1HighptM', 160): [3, 160.0*GeV, 2],
        ('medium1HighptH', 20): [3,  20.0*GeV, 2],
        ('medium1HighptH', 25): [3,  25.0*GeV, 2],
        ('medium1HighptH', 29): [3,  29.0*GeV, 2],
        ('medium1HighptH', 35): [3,  35.0*GeV, 2],
        ('medium1HighptH', 38): [3,  38.0*GeV, 2],
        ('medium1HighptH', 50): [3,  50.0*GeV, 2],
        ('medium1HighptH', 60): [3,  60.0*GeV, 2],
        ('medium1HighptH', 80): [3,  80.0*GeV, 2],
        ('medium1HighptH', 115): [3, 115.0*GeV, 2],
        ('medium1HighptH', 125): [3, 125.0*GeV, 2],
        ('medium1HighptH', 160): [3, 160.0*GeV, 2], 
        ('medium1HighptH', 200): [3, 200.0*GeV, 2]
        }

    thresholdsEF_FTK = {
        ('medium0', 0):  [0,3,  0.0*GeV, 2, False],
        ('medium0', 12): [0,3,  12.0*GeV, 2, False],
        ('medium0', 20): [0,3,  20.0*GeV, 2, False],
        ('medium0', 25): [0,3,  25.0*GeV, 2, False],
        ('medium0', 29): [0,3,  29.0*GeV, 2, False],
        ('medium0', 35): [0,3,  35.0*GeV, 2, False],
        ('medium0', 38): [0,3,  38.0*GeV, 2, False],
        ('medium0', 50): [0,3,  50.0*GeV, 2, False],
        ('medium0', 60): [0,3,  60.0*GeV, 2, False],
        ('medium0', 80): [0,3,  80.0*GeV, 2, False],
        ('medium0', 115): [0,3, 115.0*GeV, 2, False],
        ('medium0', 125): [0,3, 125.0*GeV, 2, False],
        ('medium0', 160): [0,3, 160.0*GeV, 2, False],
        ('medium0', 200): [0,3, 200.0*GeV, 2, False],
        }
    thresholdsHighpt = {
        ('medium1HighptL'):[250.0*GeV,330.0*GeV,410.0*GeV], 
        ('medium1HighptM'):[200.0*GeV,330.0*GeV,410.0*GeV],
        ('medium1HighptH'):[160.0*GeV,330.0*GeV,410.0*GeV]   
        }

    thresholdsEF_dikaon = {
        ('dikaon', 25):          [0.2*GeV, 0.45*GeV,    0.0*GeV, 1000.0*GeV,  0.0*GeV, 1000.0*GeV, 0.0*GeV,   15.0*GeV, 25.0*GeV, 1.5],
        ('dikaon', 35):          [0.2*GeV, 0.45*GeV,    0.0*GeV, 1000.0*GeV,  0.0*GeV, 1000.0*GeV, 0.0*GeV,   25.0*GeV, 35.0*GeV, 1.5], 
        ('dikaontight', 25):     [0.2*GeV, 0.45*GeV,    0.0*GeV, 1000.0*GeV,  0.0*GeV, 1000.0*GeV, 0.0*GeV,   15.0*GeV, 25.0*GeV, 1.0],
        ('dikaontight', 35):     [0.2*GeV, 0.45*GeV,    0.0*GeV, 1000.0*GeV,  0.0*GeV, 1000.0*GeV, 0.0*GeV,   25.0*GeV, 35.0*GeV, 1.0],
        ('dikaonmass', 25):      [0.0*GeV, 1000.0*GeV,  0.987*GeV, 1.060*GeV, 0.0*GeV, 1000.0*GeV, 0.0*GeV,   15.0*GeV, 25.0*GeV, 1.5],
        ('dikaonmass', 35):      [0.0*GeV, 1000.0*GeV,  0.987*GeV, 1.060*GeV, 0.0*GeV, 1000.0*GeV, 0.0*GeV,   25.0*GeV, 35.0*GeV, 1.5],
        ('dikaonmasstight', 25): [0.0*GeV, 1000.0*GeV,  0.987*GeV, 1.060*GeV, 0.0*GeV, 1000.0*GeV, 0.0*GeV,   15.0*GeV, 25.0*GeV, 1.0],
        ('dikaonmasstight', 35): [0.0*GeV, 1000.0*GeV,  0.987*GeV, 1.060*GeV, 0.0*GeV, 1000.0*GeV, 0.0*GeV,   25.0*GeV, 35.0*GeV, 1.0],
        ('kaonpi1', 25):         [0.0*GeV, 1000.0*GeV,  0.0*GeV, 1000.0*GeV,  0.79*GeV, 0.99*GeV,  0.89*GeV,  15.0*GeV, 25.0*GeV, 1.0],
        ('kaonpi1', 35):         [0.0*GeV, 1000.0*GeV,  0.0*GeV, 1000.0*GeV,  0.79*GeV, 0.99*GeV,  0.89*GeV,  25.0*GeV, 35.0*GeV, 1.0],
        ('kaonpi2', 25):         [0.0*GeV, 1000.0*GeV,  0.0*GeV, 1000.0*GeV,  1.8*GeV, 1.93*GeV,   1.865*GeV, 15.0*GeV, 25.0*GeV, 1.0],
        ('kaonpi2', 35):         [0.0*GeV, 1000.0*GeV,  0.0*GeV, 1000.0*GeV,  1.8*GeV, 1.93*GeV,   1.865*GeV, 25.0*GeV, 35.0*GeV, 1.0],
        ('dipion1', 25):         [0.475*GeV, 1.075*GeV, 0.0*GeV, 1000.0*GeV,  0.0*GeV, 1000.0*GeV, 0.0*GeV,   15.0*GeV, 25.0*GeV, 1.0],
        ('dipion1', 35):         [0.475*GeV, 1.075*GeV, 0.0*GeV, 1000.0*GeV,  0.0*GeV, 1000.0*GeV, 0.0*GeV,   25.0*GeV, 35.0*GeV, 1.0],
        ('dipion1loose', 25):    [0.475*GeV, 1.075*GeV, 0.0*GeV, 1000.0*GeV,  0.0*GeV, 1000.0*GeV, 0.0*GeV,   25.0*GeV, 25.0*GeV, 1.5],
        ('dipion1loose', 35):    [0.475*GeV, 1.075*GeV, 0.0*GeV, 1000.0*GeV,  0.0*GeV, 1000.0*GeV, 0.0*GeV,   25.0*GeV, 35.0*GeV, 1.5],
        ('dipion2', 25):         [0.460*GeV, 0.538*GeV, 0.0*GeV, 1000.0*GeV,  0.0*GeV, 1000.0*GeV, 0.0*GeV,   15.0*GeV, 25.0*GeV, 1.0],
        ('dipion2', 35):         [0.460*GeV, 0.538*GeV, 0.0*GeV, 1000.0*GeV,  0.0*GeV, 1000.0*GeV, 0.0*GeV,   25.0*GeV, 35.0*GeV, 1.0]
        }


    thresholdsEF_singlepion = {
        ('singlepion', 25): [30.0*GeV, 25.0*GeV, 1, 0, 0.06, 0.4, 0.85],
        ('singlepiontight', 25): [30.0*GeV, 25.0*GeV, 1, 0, 0.06, 0.4, 0.80]
        }

