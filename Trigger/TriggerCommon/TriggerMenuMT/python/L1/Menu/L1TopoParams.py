# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
# L1TopoParams.py
#
# *** IMPORTANT ***
# This should be generated by the L1Topo group directly based on
# the firmware configuration and alg specification.
# No edits should be made to this file by anybody else, other
# than propagating any changes from the provided json

L1TopoParams = {
 'DeltaEtaIncl1': {'comment': '',
                   'parameters': ['MinET1',
                                  'MinET2',
                                  'DeltaEtaMin',
                                  'DeltaEtaMax']},
 'DeltaEtaIncl2': {'comment': '',
                   'parameters': ['MinET1',
                                  'MinET2',
                                  'DeltaEtaMin',
                                  'DeltaEtaMax']},
 'DeltaEtaPhiIncl1': {'comment': '',
                      'parameters': ['MinET1',
                                     'MinET2',
                                     'DeltaEtaMin',
                                     'DeltaEtaMax',
                                     'DeltaPhiMin',
                                     'DeltaPhiMax']},
 'DeltaEtaPhiIncl2': {'comment': '',
                      'parameters': ['MinET1',
                                     'MinET2',
                                     'DeltaEtaMin',
                                     'DeltaEtaMax',
                                     'DeltaPhiMin',
                                     'DeltaPhiMax']},
 'DeltaPhiIncl1': {'comment': '',
                   'parameters': ['MinET1',
                                  'MinET2',
                                  'DeltaPhiMin',
                                  'DeltaPhiMax']},
 'DeltaPhiIncl2': {'comment': '',
                   'parameters': ['MinET1',
                                  'MinET2',
                                  'DeltaPhiMin',
                                  'DeltaPhiMax']},
 'DeltaRSqrIncl1': {'comment': '',
                    'parameters': ['MinET1',
                                   'MinET2',
                                   'DeltaRMin',
                                   'DeltaRMax']},
 'DeltaRSqrIncl1Charge': {'comment': '',
                          'parameters': ['MinET1',
                                         'MinET2',
                                         'DeltaRMin',
                                         'DeltaRMax']},
 'DeltaRSqrIncl2': {'comment': '',
                    'parameters': ['MinET1',
                                   'MinET2',
                                   'DeltaRMin',
                                   'DeltaRMax']},
 'DeltaRSqrIncl2Charge': {'comment': '',
                          'parameters': ['MinET1',
                                         'MinET2',
                                         'DeltaRMin',
                                         'DeltaRMax']},
 'DisambiguationDRIncl2': {'comment': '',
                           'parameters': ['MinET1',
                                          'MinET2',
                                          'MinET3',
                                          'DisambDRSqrMax']},
 'DisambiguationDRIncl3': {'comment': '',
                           'parameters': ['MinET1',
                                          'MinET2',
                                          'MinET3',
                                          'DisambDRSqrMin',
                                          'DisambDRSqrMax',
                                          'DisambDRSqr']},
 'DisambiguationIncl2': {'ApplyDR = 1': {'comment': '',
                                         'parameters': ['MinET1',
                                                        'MinET2',
                                                        'DisambDRSqrMin',
                                                        'DisambDRSqrMax']},
                         'ClusterOnly = 0 and ApplyDR = 0': {'comment': '',
                                                             'parameters': ['MinET1',
                                                                            'MinET2',
                                                                            'DisambDRSqrMin']},
                         'ClusterOnly = 1 and ApplyDR = 0': {'comment': '',
                                                             'parameters': ['MinET1',
                                                                            'MinET2']}},
 'DisambiguationIncl3': {'ApplyDR = 0': {'comment': '',
                                         'parameters': ['MinET1',
                                                        'MinET2',
                                                        'MinET3',
                                                        'DisambDRSqrMax']},
                         'ApplyDR = 1': {'comment': '',
                                         'parameters': ['MinET1',
                                                        'MinET2',
                                                        'MinET3',
                                                        'DisambDRSqrMax',
                                                        'DisambDRSqrMin',
                                                        'DisambDRSqrMax']}},
 'DisambiguationInvmIncl2': {'comment': '',
                             'parameters': ['MinET1',
                                            'MinET2',
                                            'MinMSqr',
                                            'MaxMSqr']},
 'EtCut': {'comment': '', 'parameters': ['MinET']},
 'EtaPhiWindow': {'comment': '',
                  'parameters': ['MinET',
                                 'MinEta',
                                 'MaxEta',
                                 'MinPhi',
                                 'MaxPhi']},
 'ExclusiveJets': {'comment': 'simulation still hardcodes A& B',
#                   'parameters': ['A', 'B', 'MinET', 'MinXi', 'MaxXi']},
                   'parameters': ['MinET1','MinXi','MaxXi']},
 'InvariantMassDeltaPhiInclusive1': {'ApplyEtaCut = 0': {'comment': '',
                                                         'parameters': ['MinET1',
                                                                        'MinET2',
                                                                        'MinMSqr',
                                                                        'MaxMSqr',
                                                                        'DeltaPhiMin',
                                                                        'DeltaPhiMax']},
                                     'ApplyEtaCut = 1': {'comment': '',
                                                         'parameters': ['MinET1',
                                                                        'MinET2',
                                                                        'MinMSqr',
                                                                        'MaxMSqr',
                                                                        'MinEta1',
                                                                        'MaxEta1',
                                                                        'MinEta2',
                                                                        'MaxEta2',
                                                                        'DeltaPhiMin',
                                                                        'DeltaPhiMax']}},
 'InvariantMassDeltaPhiInclusive2': {'ApplyEtaCut = 0': {'comment': '',
                                                         'parameters': ['MinET1',
                                                                        'MinET2',
                                                                        'MinMSqr',
                                                                        'MaxMSqr',
                                                                        'DeltaPhiMin',
                                                                        'DeltaPhiMax']},
                                     'ApplyEtaCut = 1': {'comment': '',
                                                         'parameters': ['MinET1',
                                                                        'MinET2',
                                                                        'MinMSqr',
                                                                        'MaxMSqr',
                                                                        'MinEta1',
                                                                        'MaxEta1',
                                                                        'MinEta2',
                                                                        'MaxEta2',
                                                                        'DeltaPhiMin',
                                                                        'DeltaPhiMax']}},
 #'InvariantMassDeltaRSqrIncl1':
 'InvariantMassInclusiveDeltaRSqrIncl1': {'comment': '',
                                 'parameters': ['MinET1',
                                                'MinET2',
                                                'MinMSqr',
                                                'MaxMSqr',
                                                'DeltaRMin',
                                                'DeltaRMax']},
#'InvariantMassDeltaRSqrIncl1Charge':
 'InvariantMassInclusiveDeltaRSqrIncl1Charge': {'comment': '',
                                       'parameters': ['MinET1',
                                                      'MinET2',
                                                      'MinMSqr',
                                                      'MaxMSqr',
                                                      'DeltaRMin',
                                                      'DeltaRMax']},
 'InvariantMassInclusiveDeltaRSqrIncl2Charge': {'comment': '',
                                       'parameters': ['MinET1',
                                                      'MinET2',
                                                      'MinMSqr',
                                                      'MaxMSqr',
                                                      'DeltaRMin',
                                                      'DeltaRMax']},
 #'InvariantMassDeltaRSqrIncl2': 
 'InvariantMassInclusiveDeltaRSqrIncl2': {'ApplyEtaCut = 0': {'comment': '',
                                                     'parameters': ['MinET1',
                                                                    'MinET2',
                                                                    'MinMSqr',
                                                                    'MaxMSqr',
                                                                    'DeltaRMin',
                                                                    'DeltaRMax']},
                                 'ApplyEtaCut = 1': {'comment': '',
                                                     'parameters': ['MinET1',
                                                                    'MinET2',
                                                                    'MinMSqr',
                                                                    'MaxMSqr',
                                                                    'MinEta1',
                                                                    'MaxEta1',
                                                                    'MinEta2',
                                                                    'MaxEta2',
                                                                    'DeltaRMin',
                                                                    'DeltaRMax']}},
 'InvariantMassInclusive1': {'comment': '',
                             'parameters': ['MinET1',
                                            'MinET2',
                                            'MinMSqr',
                                            'MaxMSqr']},
 'InvariantMassInclusive2': {'ApplyEtaCut = 0': {'comment': '',
                                                 'parameters': ['MinET1',
                                                                'MinET2',
                                                                'MinMSqr',
                                                                'MaxMSqr']},
                             'ApplyEtaCut = 1': {'comment': '',
                                                 'parameters': ['MinET1',
                                                                'MinET2',
                                                                'MinMSqr',
                                                                'MaxMSqr',
                                                                'MinEta1',
                                                                'MaxEta1',
                                                                'MinEta2',
                                                                'MaxEta2']}},
 'InvariantMassThreeTOBsIncl1': {'comment': '',
                                 'parameters': ['MinET1',
                                                'MinMSqr',
                                                'MaxMSqr']},
 'InvariantMassThreeTOBsIncl1Charge': {'comment': '',
                                       'parameters': ['MinET1',
                                                      'MinMSqr',
                                                      'MaxMSqr']},
 'JetHT': {'comment': 'All following pars are MinHt',
           'parameters': ['MinET', 'MinEta', 'MaxEta', 'MinHt']},
 'KalmanMETCorrection': {'comment': '', 'parameters': ['MinET'] + 6*['KFXE']},
 'MetNoSort': {'comment': '', 'parameters': []},
 'MetSort': {'comment': '', 'parameters': []},
 'DeltaPhiMinIncl2': {'comment': '',
                      'parameters': ['MinET1', 'MinET2', 'DeltaPhiMin']},
 'MuonNoSort': {'comment': '', 'parameters': []},
 'MuonSelect': {'comment': '',
                'parameters': ['MinEtRPC',
                               'MinEtTGC',
                               'MinEta',
                               'MaxEta',
                               'InnerCoinCut',
                               'FullStationCut',
                               'GoodMFieldCut']},
 'MuonSort': {'comment': '',
              'parameters': ['MinEta',
                             'MaxEta',
                             'InnerCoinCut',
                             'FullStationCut',
                             'GoodMFieldCut']},
 'MuonSort_1BC': {'comment': '',
                  'parameters': ['MinEta',
                                 'MaxEta',
                             'InnerCoinCut',
                             'FullStationCut',
                             'GoodMFieldCut']},
 'NotMatch': {'comment': '',
              'parameters': ['MinET1',
                             'MinET2',
                             'MinEta1',
                             'MaxEta1',
                             'MinEta2',
                             'MaxEta2',
                             'DRCut']},
 'Ratio': {'comment': '',
           'parameters': ['MinET2', 'MinEta', 'MaxEta', 'MinET1', 'MinHt']},
 'RatioMatch': {'comment': '', 'parameters': ['MinET1', 'MinET2','Ratio']},
 'RatioSum': {'comment': '',
              'parameters': ['minEtHt',
                             'minEtaHt',
                             'maxEtaHt',
                             'minEt',
                             'minEta',
                             'maxEta']},
 'SimpleCone': {'comment': '',
                'parameters': ['MinET', 'MinSumET', 'MaxRSqr']},
 'TransverseMassInclusive1': {'comment': '',
                              'parameters': ['MinET1',
                                             'MinET2',
                                             'MinTransverseMassSqr']},
 'eEmNoSort': {'comment': '', 'parameters': ['REtaMin', 'RHadMin', 'WsTotMin']},
 'eEmSelect': {'comment': '',
               'parameters': ['MinET', 'REtaMin', 'RHadMin', 'WsTotMin']},
 'eEmSort': {'comment': '', 'parameters': ['REtaMin', 'RHadMin', 'WsTotMin']},
 'eTauNoSort': {'comment': '', 'parameters': ['RCoreMin', 'RHadMin']},
 'eTauSelect': {'comment': '', 'parameters': ['MinET', 'RCoreMin', 'RHadMin']},
 'eTauSort': {'comment': '', 'parameters': ['RCoreMin', 'RHadMin']},
 'gJetSort': {'comment': '', 'parameters': ['MinEta', 'MaxEta']},
 'jEmNoSort': {'comment': '',
               'parameters': ['IsolationMin',
                              'EmFraction1Min',
                              'EmFraction2Min']},
 'jEmSort': {'comment': 'misses Iso cuts', 'parameters': ['MinEta', 'MinEta']},
 'jJetSelect': {'comment': '', 'parameters': ['MinET', 'MinEta', 'MaxEta']},
 'jJetSort': {'comment': '', 'parameters': ['MinEta', 'MaxEta']},

 'JetNoSort': {'comment': '', 'parameters': []},
    # For the time being, no dedicated algs for gJetNoSort, jJetNoSort, jLJetNoSort
 'jXENoSort': {'comment': '', 'parameters': []},
 'gXENoSort': {'comment': '', 'parameters': []},
}
