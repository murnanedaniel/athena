# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration

#
# File: CaloClusterCorrection/python/CaloSwLongWeights_atlfast_v1.py
# Created: July 2008, NK, SP
# Purpose: Longitudinal weights corrections for atlfast2, v1
#
# These corrections were derived by Nicolas Kerschen and Stathes Paganis
# The geometry is ATLAS-CSC-05-00-00
# 
#

from CaloClusterCorrection.common import *

#######################################################################
# 5x5 cluster size electrons.
CaloSwLongWeights_atlfast_v1_ele55 = [
               # w0        w3       escale    eoffset
[ 0.012500,  1.15582,   2.05996,   1.05900,    97.478  ],
[ 0.037500,  1.13916,   2.15886,   1.05651,   148.866  ],  
[ 0.062500,  1.19836,   2.19800,   1.05592,   166.515  ],
[ 0.087500,  1.10783,   1.80372,   1.05899,   132.295  ],
[ 0.112500,  1.13891,   1.97621,   1.05677,   127.392  ],
[ 0.137500,  1.09163,   1.92105,   1.05715,   142.410  ], 
[ 0.162500,  1.12075,   1.91242,   1.05696,   133.289  ],
[ 0.187500,  1.11606,   1.79531,   1.05610,   185.959  ],
[ 0.212500,  1.12132,   1.77841,   1.05517,   169.402  ],
[ 0.237500,  1.13230,   1.96470,   1.05234,   161.960  ], 
[ 0.262500,  1.09136,   1.81338,   1.05404,   164.730  ], 
[ 0.287500,  1.13126,   1.74608,   1.05306,   140.633  ],
[ 0.312500,  1.13924,   1.80526,   1.05212,   179.100  ], 
[ 0.337500,  1.11728,   1.74262,   1.05107,   141.292  ],
[ 0.362500,  1.09862,   1.62397,   1.05319,   115.513  ],
[ 0.387500,  1.10853,   1.55663,   1.05341,   101.907  ],
[ 0.412500,  1.10842,   1.57455,   1.05434,   115.979  ],
[ 0.437500,  1.12054,   1.35661,   1.05318,    98.632  ],
[ 0.462500,  1.08514,   1.25478,   1.05505,   131.338  ],
[ 0.487500,  1.11564,   1.20856,   1.05296,   114.356  ],
[ 0.512500,  1.10298,   1.26542,   1.05263,   143.596  ],
[ 0.537500,  1.09626,   1.11066,   1.05354,   160.638  ],
[ 0.562500,  1.10833,   1.09301,   1.05377,   186.482  ],
[ 0.587500,  1.08126,   0.98197,   1.05344,   184.023  ], 
[ 0.612500,  1.11874,   0.95467,   1.05404,   145.009  ],
[ 0.637500,  1.08873,   0.90502,   1.05419,   180.130  ], 
[ 0.662500,  1.11698,   0.86721,   1.05569,   239.301  ],
[ 0.687500,  1.17124,   1.27586,   1.05024,   372.521  ],
[ 0.712500,  1.19670,   1.22965,   1.05220,   325.370  ], 
[ 0.737500,  1.32462,   1.69236,   1.05483,   379.715  ],
[ 0.762500,  1.05935,   2.60483,   1.06103,   796.297  ],
[ 0.787500,  0.99460,   1.13088,   1.04724,   877.263  ],
[ 0.812500,  1.00205,   1.69343,   1.04777,   768.794  ],
[ 0.837500,  1.07966,   1.35299,   1.05711,   506.078  ],
[ 0.862500,  1.19279,   1.38961,   1.05492,   440.185  ],
[ 0.887500,  1.12274,   0.76549,   1.06415,   435.548  ],
[ 0.912500,  1.13079,   0.93369,   1.06180,   531.104  ],
[ 0.937500,  1.14241,   0.74080,   1.06451,   528.717  ],
[ 0.962500,  1.13424,   0.71617,   1.06641,   479.644  ], 
[ 0.987500,  1.11276,   0.49831,   1.07085,   572.522  ],
[ 1.012500,  1.12715,   0.54889,   1.06988,   499.717  ],
[ 1.037500,  1.12116,   0.35821,   1.07431,   430.218  ],
[ 1.062500,  1.11848,   0.44614,   1.07251,   478.190  ], 
[ 1.087500,  1.08113,   0.50462,   1.07500,   469.831  ],
[ 1.112500,  1.08277,   0.53183,   1.07512,   341.536  ],
[ 1.137500,  1.07590,   0.38499,   1.07673,   374.248  ],
[ 1.162500,  1.05906,   0.30245,   1.08114,   427.511  ],
[ 1.187500,  1.08548,   0.41874,   1.07902,   569.091  ],
[ 1.212500,  1.09595,   0.36830,   1.07860,   549.295  ],
[ 1.237500,  1.05873,   0.40780,   1.08157,   414.223  ],
[ 1.262500,  1.05463,   0.21300,   1.08592,   475.051  ],
[ 1.287500,  1.04984,   0.00000,   1.09287,   515.972  ],  
[ 1.312500,  1.05735,   0.00000,   1.09522,   486.909  ],
[ 1.337500,  1.00294,   0.00000,   1.12102,   658.315  ],
[ 1.362500,  1.01132,   0.00000,   1.12722,   593.622  ],
[ 1.387500,  0.95776,   0.00000,   1.14166,   592.800  ], 
[ 1.412500,  0.94683,   2.85610,   1.14432,   572.922  ],
[ 1.437500,  0.04897,   0.00000,   1.36610,  6964.540  ],
[ 1.462500,  1.32877,   0.00000,   1.41732,  3766.470  ],
[ 1.487500,  1.61020,   0.00000,   1.15711,   877.198  ],
[ 1.512500,  1.87721,   0.00000,   1.06779,  1581.020  ],
[ 1.537500,  1.31987,   0.00000,   1.06769,  1140.000  ],  
[ 1.562500,  1.26017,   0.12484,   1.06353,  1323.400  ],   	  
[ 1.587500,  1.22728,   0.38637,   1.06567,   975.832  ],
[ 1.612500,  1.23909,   0.43694,   1.06339,  1032.210  ], 
[ 1.637500,  1.27395,   0.04805,   1.06823,   883.986  ],
[ 1.662500,  1.27565,   0.54514,   1.06329,  1158.130  ],
[ 1.687500,  1.21483,   0.32551,   1.07197,   919.775  ],  
[ 1.712500,  1.20215,   0.62569,   1.06544,   592.673  ],
[ 1.737500,  1.26415,   0.75686,   1.06422,   588.125  ],
[ 1.762500,  1.20572,   0.80570,   1.06393,   544.170  ], 
[ 1.787500,  1.34147,   0.38310,   1.07915,   668.003  ],
[ 1.812500,  0.00000,   0.19805,   1.08871,   878.546  ],
[ 1.837500,  0.00000,   0.01502,   1.09278,  1164.570  ],
[ 1.862500,  0.00000,   0.00000,   1.10147,  1120.800  ],
[ 1.887500,  0.00000,   0.00000,   1.10251,  1119.080  ],
[ 1.912500,  0.00000,   0.05039,   1.09672,   978.081  ],
[ 1.937500,  0.00000,   0.42872,   1.08848,  1005.550  ],
[ 1.962500,  0.00000,   0.47971,   1.08864,   844.663  ],
[ 1.987500,  0.00000,   0.25487,   1.09347,   950.044  ],
[ 2.012500,  0.00000,   0.22506,   1.09962,   759.410  ], 
[ 2.037500,  0.00000,   0.01070,   1.10661,   903.627  ],
[ 2.062500,  0.00000,   0.05780,   1.10798,   941.159  ],
[ 2.087500,  0.00000,   0.25594,   1.10752,   918.558  ],
[ 2.112500,  0.00000,   0.17908,   1.10787,   814.686  ],
[ 2.137500,  0.00000,   0.20956,   1.10335,   866.134  ],
[ 2.162500,  0.00000,   0.35464,   1.09899,   765.715  ],
[ 2.187500,  0.00000,   0.34045,   1.10217,   562.155  ],
[ 2.212500,  0.00000,   0.43862,   1.10084,   784.489  ],
[ 2.237500,  0.00000,   0.53924,   1.10119,   732.490  ], 
[ 2.262500,  0.00000,   0.59195,   1.10183,   577.309  ],
[ 2.287500,  0.00000,   0.55789,   1.10613,   572.924  ],
[ 2.312500,  0.00000,   0.60550,   1.10433,   624.400  ],  
[ 2.337500,  0.00000,   0.57225,   1.10188,   666.793  ],
[ 2.362500,  0.00000,   0.74692,   1.09852,   774.555  ],
[ 2.387500,  0.00000,   0.46868,   1.10233,   651.411  ],
[ 2.412500,  0.00000,   0.50674,   1.10387,   571.244  ], 
[ 2.437500,  0.00000,   0.92592,   1.09862,   653.818  ],
[ 2.462500,  0.00000,   0.81821,   1.11154,   766.515  ],
[ 2.487500,  0.00000,   2.58419,   1.14068,  1559.690  ]
] 



#######################################################################
#  3x5 cluster size electrons. (same as 3x7 weights)
CaloSwLongWeights_atlfast_v1_ele35 = [
                # w0        w3       escale    eoffset
[ 0.012500,  1.17039,   1.99402,   1.07444,   113.241  ], 
[ 0.037500,  1.15606,   2.07609,   1.07166,   184.178  ], 
[ 0.062500,  1.20302,   2.15518,   1.07079,   187.128  ],  
[ 0.087500,  1.10861,   1.82859,   1.07329,   164.160  ],
[ 0.112500,  1.14951,   1.88668,   1.07232,   122.715  ], 
[ 0.137500,  1.10047,   1.89203,   1.07171,   159.400  ], 	
[ 0.162500,  1.12619,   1.87306,   1.07172,   140.913  ], 
[ 0.187500,  1.12955,   1.77372,   1.07060,   183.628  ],  
[ 0.212500,  1.13334,   1.79868,   1.06881,   204.930  ], 	 
[ 0.237500,  1.14723,   1.91770,   1.06705,   172.542  ],  
[ 0.262500,  1.11585,   1.76266,   1.06824,   171.207  ],  
[ 0.287500,  1.12420,   1.65895,   1.06777,   150.064  ],  
[ 0.312500,  1.11603,   1.67064,   1.06736,   197.703  ], 
[ 0.337500,  1.11167,   1.69190,   1.06572,   157.148  ],  
[ 0.362500,  1.10962,   1.63716,   1.06732,   147.766  ],  
[ 0.387500,  1.10361,   1.46945,   1.06831,   132.505  ],  	
[ 0.412500,  1.11111,   1.54038,   1.06893,   129.236  ],     
[ 0.437500,  1.10642,   1.26619,   1.06841,   144.412  ], 
[ 0.462500,  1.08661,   1.21179,   1.07009,   143.624  ],  
[ 0.487500,  1.11085,   1.11156,   1.06834,   131.031  ],  
[ 0.512500,  1.10385,   1.19456,   1.06755,   162.563  ], 
[ 0.537500,  1.10683,   1.09002,   1.06812,   165.173  ],  
[ 0.562500,  1.10958,   1.05997,   1.06851,   185.132  ], 
[ 0.587500,  1.08264,   0.96275,   1.06764,   232.110  ],
[ 0.612500,  1.10523,   0.94908,   1.06849,   177.050  ],
[ 0.637500,  1.08747,   0.87083,   1.06953,   179.769  ], 
[ 0.662500,  1.10958,   0.84346,   1.07086,   203.593  ], 
[ 0.687500,  1.14963,   1.22490,   1.06534,   365.457  ], 
[ 0.712500,  1.18510,   1.21161,   1.06683,   308.674  ], 
[ 0.737500,  1.30321,   1.70146,   1.06919,   395.361  ], 
[ 0.762500,  1.04238,   2.70256,   1.07504,   775.099  ], 
[ 0.787500,  0.98215,   1.52503,   1.06121,   844.266  ], 
[ 0.812500,  0.98956,   1.84893,   1.06163,   780.106  ], 
[ 0.837500,  1.09142,   1.37793,   1.07203,   442.118  ], 
[ 0.862500,  1.15588,   1.36896,   1.07038,   391.231  ], 
[ 0.887500,  1.12147,   0.74533,   1.07890,   361.064  ], 
[ 0.912500,  1.11854,   0.86479,   1.07837,   458.709  ],  
[ 0.937500,  1.12036,   0.60487,   1.08281,   457.385  ], 
[ 0.962500,  1.13129,   0.60305,   1.08366,   467.228  ], 
[ 0.987500,  1.09764,   0.54156,   1.08647,   568.482  ], 
[ 1.012500,  1.10807,   0.51886,   1.08660,   466.142  ],	
[ 1.037500,  1.10603,   0.32934,   1.09077,   420.791  ], 
[ 1.062500,  1.11080,   0.41998,   1.08874,   448.248  ], 	
[ 1.087500,  1.07048,   0.51764,   1.09111,   448.792  ],
[ 1.112500,  1.09542,   0.48777,   1.09130,   291.311  ],	
[ 1.137500,  1.06700, 	0.35870,   1.09404,   371.073  ],
[ 1.162500,  1.05647,   0.33691,   1.09709,   371.747  ], 
[ 1.187500,  1.08067,   0.33635,   1.09747,   465.777  ], 
[ 1.212500,  1.08377,   0.30476,   1.09742,   461.657  ], 
[ 1.237500,  1.05294,   0.36360,   1.09930,   363.047  ], 
[ 1.262500,  1.05798,   0.18570,   1.10208,   466.147  ],  
[ 1.287500,  1.06031,   0.00000,   1.10756,   405.161  ], 
[ 1.312500,  1.06866,   0.00000,   1.10917,   414.499  ], 
[ 1.337500,  1.01246,   0.00000,   1.13681,   503.358  ],
[ 1.362500,  1.01476,   0.00000,   1.13567,   516.879  ], 
[ 1.387500,  0.99368,   0.00000,   1.20087,   594.180  ],
[ 1.412500,  0.96416,   0.10336,   1.15026,   541.750  ],
[ 1.437500,  0.00000,   16.8394,   1.41195,  7371.370  ],
[ 1.462500,  1.88034,   20.0000,   1.45519,  3780.570  ],
[ 1.487500,  1.15640,   0.00000,   1.23951,  1216.910  ],
[ 1.512500,  1.88653,   0.00000,   1.07859,  1153.040  ],
[ 1.537500,  1.31293,   0.00000,   1.07609,   992.383  ],
[ 1.562500,  1.23548,   0.08084,   1.07599,  1220.110  ],
[ 1.587500,  1.21889,   0.45404,   1.07692,   883.685  ],
[ 1.612500,  1.22889,   0.37477,   1.07699,   827.017  ],
[ 1.637500,  1.26288,   0.16545,   1.08028,   733.912  ],
[ 1.662500,  1.26440,   0.54997,   1.07704,   965.327  ],
[ 1.687500,  1.28363,   0.33796,   1.08175,   809.316  ],
[ 1.712500,  1.23950,   0.80043,   1.07695,   466.174  ],
[ 1.737500,  1.29797,   0.72355,   1.07909,   537.325  ],
[ 1.762500,  1.17293,   0.91190,   1.07794,   444.046  ],
[ 1.787500,  1.38069,   0.35092,   1.09419,   678.945  ],
[ 1.812500,  0.00000,   0.35651,   1.10221,   910.758  ],
[ 1.837500,  0.00000,   0.00000,   1.11176,  1178.500  ],
[ 1.862500,  0.00000,   0.00000,   1.11996,  1126.800  ],
[ 1.887500,  0.00000,   0.00000,   1.12222,  1105.680  ],
[ 1.912500,  0.00000,   0.00000,   1.11756,   911.105  ],
[ 1.937500,  0.00000,   0.37413,   1.10750,   910.053  ],	 
[ 1.962500,  0.00000,   0.45851,   1.10705,   806.992  ],
[ 1.987500,  0.00000,   0.31399,   1.11144,   806.586  ],
[ 2.012500,  0.00000,   0.18987,   1.11917,   743.685  ],
[ 2.037500,  0.00000,   0.02381,   1.12758,   828.436  ],
[ 2.062500,  0.00000,   0.01256,   1.13162,   843.278  ],
[ 2.087500,  0.00000,   0.29957,   1.12993,   819.772  ],
[ 2.112500,  0.00000,   0.31098,   1.12770,   780.175  ],
[ 2.137500,  0.00000,   0.25477,   1.12333,   848.701  ],
[ 2.162500,  0.00000,   0.38074,   1.11943,   738.869  ], 
[ 2.187500,  0.00000,   0.47069,   1.12100,   636.899  ], 
[ 2.212500,  0.00000,   0.46815,   1.12265,   740.430  ],
[ 2.237500,  0.00000,   0.53536,   1.12476,   677.919  ], 	
[ 2.262500,  0.00000,   0.62058,   1.12561,   528.987  ],
[ 2.287500,  0.00000,   0.50183,   1.13185,   554.038  ],  
[ 2.312500,  0.00000,   0.57380,   1.12881,   675.490  ],
[ 2.337500,  0.00000,   0.57348,   1.12611,   659.311  ],
[ 2.362500,  0.00000,   0.66744,   1.12436,   704.394  ],
[ 2.387500,  0.00000,   0.41550,   1.12833,   642.067  ],  
[ 2.412500,  0.00000,   0.47125,   1.13026,   586.030  ],
[ 2.437500,  0.00000,   0.81466,   1.12760,   656.386  ],
[ 2.462500,  0.00000,   0.82526,   1.13050,   703.149  ],	
[ 2.487500,  0.00000,   1.98347,   1.15770,  1301.860  ]
] 


#######################################################################
# 3x7 cluster size electrons.
CaloSwLongWeights_atlfast_v1_ele37 = [
                # w0        w3       escale    eoffset
[ 0.012500,  1.17039,   1.99402,   1.07444,   113.241  ], 
[ 0.037500,  1.15606,   2.07609,   1.07166,   184.178  ], 
[ 0.062500,  1.20302,   2.15518,   1.07079,   187.128  ],  
[ 0.087500,  1.10861,   1.82859,   1.07329,   164.160  ],
[ 0.112500,  1.14951,   1.88668,   1.07232,   122.715  ], 
[ 0.137500,  1.10047,   1.89203,   1.07171,   159.400  ], 	
[ 0.162500,  1.12619,   1.87306,   1.07172,   140.913  ], 
[ 0.187500,  1.12955,   1.77372,   1.07060,   183.628  ],  
[ 0.212500,  1.13334,   1.79868,   1.06881,   204.930  ], 	 
[ 0.237500,  1.14723,   1.91770,   1.06705,   172.542  ],  
[ 0.262500,  1.11585,   1.76266,   1.06824,   171.207  ],  
[ 0.287500,  1.12420,   1.65895,   1.06777,   150.064  ],  
[ 0.312500,  1.11603,   1.67064,   1.06736,   197.703  ], 
[ 0.337500,  1.11167,   1.69190,   1.06572,   157.148  ],  
[ 0.362500,  1.10962,   1.63716,   1.06732,   147.766  ],  
[ 0.387500,  1.10361,   1.46945,   1.06831,   132.505  ],  	
[ 0.412500,  1.11111,   1.54038,   1.06893,   129.236  ],     
[ 0.437500,  1.10642,   1.26619,   1.06841,   144.412  ], 
[ 0.462500,  1.08661,   1.21179,   1.07009,   143.624  ],  
[ 0.487500,  1.11085,   1.11156,   1.06834,   131.031  ],  
[ 0.512500,  1.10385,   1.19456,   1.06755,   162.563  ], 
[ 0.537500,  1.10683,   1.09002,   1.06812,   165.173  ],  
[ 0.562500,  1.10958,   1.05997,   1.06851,   185.132  ], 
[ 0.587500,  1.08264,   0.96275,   1.06764,   232.110  ],
[ 0.612500,  1.10523,   0.94908,   1.06849,   177.050  ],
[ 0.637500,  1.08747,   0.87083,   1.06953,   179.769  ], 
[ 0.662500,  1.10958,   0.84346,   1.07086,   203.593  ], 
[ 0.687500,  1.14963,   1.22490,   1.06534,   365.457  ], 
[ 0.712500,  1.18510,   1.21161,   1.06683,   308.674  ], 
[ 0.737500,  1.30321,   1.70146,   1.06919,   395.361  ], 
[ 0.762500,  1.04238,   2.70256,   1.07504,   775.099  ], 
[ 0.787500,  0.98215,   1.52503,   1.06121,   844.266  ], 
[ 0.812500,  0.98956,   1.84893,   1.06163,   780.106  ], 
[ 0.837500,  1.09142,   1.37793,   1.07203,   442.118  ], 
[ 0.862500,  1.15588,   1.36896,   1.07038,   391.231  ], 
[ 0.887500,  1.12147,   0.74533,   1.07890,   361.064  ], 
[ 0.912500,  1.11854,   0.86479,   1.07837,   458.709  ],  
[ 0.937500,  1.12036,   0.60487,   1.08281,   457.385  ], 
[ 0.962500,  1.13129,   0.60305,   1.08366,   467.228  ], 
[ 0.987500,  1.09764,   0.54156,   1.08647,   568.482  ], 
[ 1.012500,  1.10807,   0.51886,   1.08660,   466.142  ],	
[ 1.037500,  1.10603,   0.32934,   1.09077,   420.791  ], 
[ 1.062500,  1.11080,   0.41998,   1.08874,   448.248  ], 	
[ 1.087500,  1.07048,   0.51764,   1.09111,   448.792  ],
[ 1.112500,  1.09542,   0.48777,   1.09130,   291.311  ],	
[ 1.137500,  1.06700, 	0.35870,   1.09404,   371.073  ],
[ 1.162500,  1.05647,   0.33691,   1.09709,   371.747  ], 
[ 1.187500,  1.08067,   0.33635,   1.09747,   465.777  ], 
[ 1.212500,  1.08377,   0.30476,   1.09742,   461.657  ], 
[ 1.237500,  1.05294,   0.36360,   1.09930,   363.047  ], 
[ 1.262500,  1.05798,   0.18570,   1.10208,   466.147  ],  
[ 1.287500,  1.06031,   0.00000,   1.10756,   405.161  ], 
[ 1.312500,  1.06866,   0.00000,   1.10917,   414.499  ], 
[ 1.337500,  1.01246,   0.00000,   1.13681,   503.358  ],
[ 1.362500,  1.01476,   0.00000,   1.13567,   516.879  ], 
[ 1.387500,  0.99368,   0.00000,   1.20087,   594.180  ],
[ 1.412500,  0.96416,   0.10336,   1.15026,   541.750  ],
[ 1.437500,  0.00000,   16.8394,   1.41195,  7371.370  ],
[ 1.462500,  1.88034,   20.0000,   1.45519,  3780.570  ],
[ 1.487500,  1.15640,   0.00000,   1.23951,  1216.910  ],
[ 1.512500,  1.88653,   0.00000,   1.07859,  1153.040  ],
[ 1.537500,  1.31293,   0.00000,   1.07609,   992.383  ],
[ 1.562500,  1.23548,   0.08084,   1.07599,  1220.110  ],
[ 1.587500,  1.21889,   0.45404,   1.07692,   883.685  ],
[ 1.612500,  1.22889,   0.37477,   1.07699,   827.017  ],
[ 1.637500,  1.26288,   0.16545,   1.08028,   733.912  ],
[ 1.662500,  1.26440,   0.54997,   1.07704,   965.327  ],
[ 1.687500,  1.28363,   0.33796,   1.08175,   809.316  ],
[ 1.712500,  1.23950,   0.80043,   1.07695,   466.174  ],
[ 1.737500,  1.29797,   0.72355,   1.07909,   537.325  ],
[ 1.762500,  1.17293,   0.91190,   1.07794,   444.046  ],
[ 1.787500,  1.38069,   0.35092,   1.09419,   678.945  ],
[ 1.812500,  0.00000,   0.35651,   1.10221,   910.758  ],
[ 1.837500,  0.00000,   0.00000,   1.11176,  1178.500  ],
[ 1.862500,  0.00000,   0.00000,   1.11996,  1126.800  ],
[ 1.887500,  0.00000,   0.00000,   1.12222,  1105.680  ],
[ 1.912500,  0.00000,   0.00000,   1.11756,   911.105  ],
[ 1.937500,  0.00000,   0.37413,   1.10750,   910.053  ],	 
[ 1.962500,  0.00000,   0.45851,   1.10705,   806.992  ],
[ 1.987500,  0.00000,   0.31399,   1.11144,   806.586  ],
[ 2.012500,  0.00000,   0.18987,   1.11917,   743.685  ],
[ 2.037500,  0.00000,   0.02381,   1.12758,   828.436  ],
[ 2.062500,  0.00000,   0.01256,   1.13162,   843.278  ],
[ 2.087500,  0.00000,   0.29957,   1.12993,   819.772  ],
[ 2.112500,  0.00000,   0.31098,   1.12770,   780.175  ],
[ 2.137500,  0.00000,   0.25477,   1.12333,   848.701  ],
[ 2.162500,  0.00000,   0.38074,   1.11943,   738.869  ], 
[ 2.187500,  0.00000,   0.47069,   1.12100,   636.899  ], 
[ 2.212500,  0.00000,   0.46815,   1.12265,   740.430  ],
[ 2.237500,  0.00000,   0.53536,   1.12476,   677.919  ], 	
[ 2.262500,  0.00000,   0.62058,   1.12561,   528.987  ],
[ 2.287500,  0.00000,   0.50183,   1.13185,   554.038  ],  
[ 2.312500,  0.00000,   0.57380,   1.12881,   675.490  ],
[ 2.337500,  0.00000,   0.57348,   1.12611,   659.311  ],
[ 2.362500,  0.00000,   0.66744,   1.12436,   704.394  ],
[ 2.387500,  0.00000,   0.41550,   1.12833,   642.067  ],  
[ 2.412500,  0.00000,   0.47125,   1.13026,   586.030  ],
[ 2.437500,  0.00000,   0.81466,   1.12760,   656.386  ],
[ 2.462500,  0.00000,   0.82526,   1.13050,   703.149  ],	
[ 2.487500,  0.00000,   1.98347,   1.15770,  1301.860  ]
] 


#######################################################################
# 5x5 cluster size photons. 
CaloSwLongWeights_atlfast_v1_gam55 = [
[ 0.012500,   1.1834,   2.0521,   1.0564,   56.35  ],
[ 0.037500,   1.2141,   2.2984,   1.0520,  136.76  ],
[ 0.062500,   1.2288,   2.1245,   1.0548,   58.47  ],
[ 0.087500,   1.1333,   1.9911,   1.0547,   98.60  ],
[ 0.112500,   1.1576,   2.1096,   1.0543,   70.33  ],
[ 0.137500,   1.1317,   2.0462,   1.0532,   73.53  ],
[ 0.162500,   1.1352,   1.9724,   1.0533,  113.38  ],
[ 0.187500,   1.1289,   1.8964,   1.0534,  146.48  ],
[ 0.212500,   1.1614,   1.9695,   1.0517,   75.59  ],
[ 0.237500,   1.1605,   2.0308,   1.0504,   88.21  ],
[ 0.262500,   1.1935,   1.9299,   1.0506,   71.62  ],
[ 0.287500,   1.1503,   2.0380,   1.0481,  107.98  ],
[ 0.312500,   1.0974,   1.8786,   1.0501,  112.18  ],
[ 0.337500,   1.0966,   1.9082,   1.0490,   51.71  ],
[ 0.362500,   1.1175,   1.9185,   1.0489,   63.03  ],
[ 0.387500,   1.1122,   1.6876,   1.0497,   54.58  ],
[ 0.412500,   1.1379,   1.6298,   1.0512,   51.25  ],
[ 0.437500,   1.1322,   1.4257,   1.0510,    5.98  ],
[ 0.462500,   1.1497,   1.3447,   1.0521,    8.83  ],
[ 0.487500,   1.1360,   1.3159,   1.0497,   49.14  ],
[ 0.512500,   1.1423,   1.3673,   1.0495,   53.99  ],
[ 0.537500,   1.1531,   1.2597,   1.0488,   55.77  ],
[ 0.562500,   1.1608,   1.1927,   1.0494,   45.27  ],
[ 0.587500,   1.1353,   1.1415,   1.0489,   70.70  ],
[ 0.612500,   1.1191,   1.0781,   1.0491,   87.84  ],
[ 0.637500,   1.1209,   0.9833,   1.0500,   58.58  ],
[ 0.662500,   1.1871,   1.2135,   1.0465,   49.17  ],
[ 0.687500,   1.2586,   1.3473,   1.0439,   29.70  ],
[ 0.712500,   1.2761,   1.3179,   1.0449,   31.24  ],
[ 0.737500,   1.4385,   1.6727,   1.0440,    0.00  ],
[ 0.762500,   1.1413,   2.3042,   1.0490,  133.77  ],
[ 0.787500,   1.0275,   0.7179,   1.0301,  233.98  ],
[ 0.812500,   1.0350,   2.3544,   1.0309,  154.87  ],
[ 0.837500,   1.1156,   1.5142,   1.0484,    6.92  ],
[ 0.862500,   1.1671,   1.5382,   1.0468,    0.00  ],
[ 0.887500,   1.1794,   1.1915,   1.0498,    0.00  ],
[ 0.912500,   1.1631,   1.0907,   1.0507,    0.00  ],
[ 0.937500,   1.1612,   1.0872,   1.0502,    0.00  ],
[ 0.962500,   1.1426,   1.0101,   1.0526,    0.00  ],
[ 0.987500,   1.1357,   1.0437,   1.0526,    0.00  ],
[ 1.012500,   1.1254,   0.8780,   1.0550,    0.00  ],
[ 1.037500,   1.1378,   0.9226,   1.0548,    0.00  ],
[ 1.062500,   1.1260,   0.8556,   1.0562,    0.00  ],
[ 1.087500,   1.1341,   0.8516,   1.0569,    0.00  ],
[ 1.112500,   1.1220,   0.8567,   1.0555,   29.72  ],
[ 1.137500,   1.1098,   0.6804,   1.0606,    0.00  ],
[ 1.162500,   1.1085,   0.6981,   1.0608,   22.57  ],
[ 1.187500,   1.0956,   0.5804,   1.0653,    0.00  ],
[ 1.212500,   1.0971,   0.5012,   1.0660,    0.00  ],
[ 1.237500,   1.0735,   0.4945,   1.0688,    0.00  ],
[ 1.262500,   1.0764,   0.4962,   1.0680,    0.00  ],
[ 1.287500,   1.0549,   0.0342,   1.0803,    0.00  ],
[ 1.312500,   1.0427,   0.0017,   1.0811,    0.00  ],
[ 1.337500,   1.0445,   0.0000,   1.1063,    0.00  ],
[ 1.362500,   1.0343,   0.0000,   1.1109,    0.00  ],
[ 1.387500,   0.9382,   0.0000,   1.1344,    0.00  ],
[ 1.412500,   0.9475,   0.2234,   1.1376,    0.00  ],
[ 1.437500,   0.8723,   0.0000,   1.2561, 5848.21  ],
[ 1.462500,   1.5569,   0.0233,   1.3842, 3072.03  ],
[ 1.487500,   2.0306,   0.0000,   1.0892,    0.00  ],
[ 1.512500,   1.8586,   0.0000,   1.0532,    0.00  ],
[ 1.537500,   1.1502,   0.4063,   1.0502,    0.00  ],
[ 1.562500,   1.1686,   0.3692,   1.0517,    0.00  ],
[ 1.587500,   1.1701,   0.7160,   1.0477,    0.00  ],
[ 1.612500,   1.1726,   0.7018,   1.0493,    0.00  ],
[ 1.637500,   1.1725,   0.9298,   1.0464,    0.00  ],
[ 1.662500,   1.1691,   1.0691,   1.0474,    0.00  ],
[ 1.687500,   1.1469,   0.9871,   1.0483,    0.00  ],
[ 1.712500,   1.1076,   1.0546,   1.0481,    0.00  ],
[ 1.737500,   1.2413,   1.1551,   1.0454,    0.00  ],
[ 1.762500,   1.2015,   1.0967,   1.0471,    5.88  ],
[ 1.787500,   1.2200,   0.9552,   1.0538,  132.13  ],
[ 1.812500,   0.0000,   0.8264,   1.0590,  114.34  ],
[ 1.837500,   0.0000,   0.8106,   1.0624,  155.81  ],
[ 1.862500,   0.0000,   0.6688,   1.0662,  119.34  ],
[ 1.887500,   0.0000,   0.7479,   1.0682,  125.54  ],
[ 1.912500,   0.0000,   0.8319,   1.0653,  131.00  ],
[ 1.937500,   0.0000,   0.8327,   1.0659,  183.96  ],
[ 1.962500,   0.0000,   0.8386,   1.0655,  135.18  ],
[ 1.987500,   0.0000,   0.6959,   1.0715,  135.06  ],
[ 2.012500,   0.0000,   0.6280,   1.0760,   67.34  ],
[ 2.037500,   0.0000,   0.6943,   1.0816,   64.77  ],
[ 2.062500,   0.0000,   0.6632,   1.0840,  111.72  ],
[ 2.087500,   0.0000,   0.6749,   1.0855,  172.60  ],
[ 2.112500,   0.0000,   0.7101,   1.0860,  119.39  ],
[ 2.137500,   0.0000,   0.8788,   1.0814,  127.23  ],
[ 2.162500,   0.0000,   0.8995,   1.0793,  175.50  ],
[ 2.187500,   0.0000,   0.8078,   1.0834,  141.36  ],
[ 2.212500,   0.0000,   0.8496,   1.0844,   85.72  ],
[ 2.237500,   0.0000,   0.8349,   1.0879,  104.40  ],
[ 2.262500,   0.0000,   0.9106,   1.0884,   96.28  ],
[ 2.287500,   0.0000,   0.8413,   1.0907,  129.90  ],
[ 2.312500,   0.0000,   0.9083,   1.0895,  137.24  ],
[ 2.337500,   0.0000,   0.8813,   1.0862,  175.49  ],
[ 2.362500,   0.0000,   1.0185,   1.0826,  201.62  ],
[ 2.387500,   0.0000,   0.8034,   1.0872,  133.74  ],
[ 2.412500,   0.0000,   0.7219,   1.0898,  110.63  ],
[ 2.437500,   0.0000,   1.0294,   1.0866,  186.49  ],
[ 2.462500,   0.0000,   1.0379,   1.0974,  211.31  ],
[ 2.487500,   0.0000,   2.2830,   1.1356,  738.87  ]
]



#######################################################################
# 3x5 cluster size photons. 
CaloSwLongWeights_atlfast_v1_gam35 = [
[ 0.012500,   1.2366,   2.0684,   1.0754,   104.65  ],
[ 0.037500,   1.2496,   2.3059,   1.0707,   180.67  ],
[ 0.062500,   1.2778,   2.1103,   1.0737,    77.09  ],
[ 0.087500,   1.1605,   1.9479,   1.0744,   115.28  ],
[ 0.112500,   1.1879,   2.0967,   1.0732,   101.19  ],
[ 0.137500,   1.1832,   2.0284,   1.0724,    74.83  ],
[ 0.162500,   1.1623,   1.9560,   1.0724,   126.48  ],
[ 0.187500,   1.1639,   1.8758,   1.0721,   168.02  ],
[ 0.212500,   1.2040,   2.0096,   1.0696,   113.91  ],
[ 0.237500,   1.1886,   2.0048,   1.0691,   115.86  ],
[ 0.262500,   1.2536,   1.9001,   1.0688,    96.66  ],
[ 0.287500,   1.2149,   2.0653,   1.0655,   147.32  ],
[ 0.312500,   1.1529,   1.8813,   1.0681,   152.10  ],
[ 0.337500,   1.1354,   1.9359,   1.0672,    74.76  ],
[ 0.362500,   1.1569,   1.9701,   1.0671,   104.97  ],
[ 0.387500,   1.1513,   1.7207,   1.0683,    83.58  ],
[ 0.412500,   1.1606,   1.5623,   1.0706,    69.54  ],
[ 0.437500,   1.1568,   1.3449,   1.0706,    35.25  ],
[ 0.462500,   1.1808,   1.2517,   1.072,     18.08  ],
[ 0.487500,   1.1528,   1.2519,   1.0693,    57.37  ],
[ 0.512500,   1.1745,   1.3477,   1.0684,    62.06  ],
[ 0.537500,   1.1638,   1.1854,   1.0685,    90.53  ],
[ 0.562500,   1.1998,   1.1687,   1.0685,    61.96  ],
[ 0.587500,   1.1569,   1.1163,   1.0685,    69.02  ],
[ 0.612500,   1.1405,   1.0606,   1.0683,   113.71  ],
[ 0.637500,   1.1402,   0.9642,   1.0695,    79.39  ],
[ 0.662500,   1.2031,   1.1779,   1.0662,    80.03  ],
[ 0.687500,   1.2893,   1.4292,   1.0621,    64.08  ],
[ 0.712500,   1.2935,   1.3535,   1.0644,    58.06  ],
[ 0.737500,   1.4294,   1.9157,   1.0651,     0.00  ],
[ 0.762500,   1.1575,   2.6176,   1.0692,   136.61  ],
[ 0.787500,   1.0230,   1.2167,   1.0500,   293.81  ],
[ 0.812500,   1.0347,   2.6494,   1.0517,   183.29  ],
[ 0.837500,   1.1184,   1.4828,   1.0708,    14.44  ],
[ 0.862500,   1.1791,   1.4896,   1.0682,     0.00  ],
[ 0.887500,   1.1900,   1.1609,   1.0715,     0.00  ],
[ 0.912500,   1.1838,   1.0168,   1.0731,     0.00  ],
[ 0.937500,   1.1842,   1.0360,   1.0730,     0.00  ],
[ 0.962500,   1.1653,   0.8954,   1.0767,     0.00  ],
[ 0.987500,   1.1495,   0.9831,   1.0764,     0.00  ],
[ 1.012500,   1.1436,   0.8156,   1.0785,     0.00  ],
[ 1.037500,   1.1563,   0.8833,   1.0784,     0.00  ],
[ 1.062500,   1.1448,   0.8273,   1.0796,     0.00  ],
[ 1.087500,   1.1570,   0.8368,   1.0800,     0.00  ],
[ 1.112500,   1.1476,   0.8148,   1.0787,    56.65  ],
[ 1.137500,   1.1274,   0.6277,   1.0851,     0.00  ],
[ 1.162500,   1.1298,   0.6566,   1.0851,    24.79  ],
[ 1.187500,   1.1192,   0.5677,   1.0894,     0.00  ],
[ 1.212500,   1.1175,   0.4432,   1.0909,     0.00  ],
[ 1.237500,   1.0944,   0.4741,   1.0926,     0.00  ],
[ 1.262500,   1.1011,   0.4594,   1.0920,     0.00  ],
[ 1.287500,   1.0788,   0.0000,   1.1041,     0.00  ],
[ 1.312500,   1.0743,   0.0000,   1.1045,     0.00  ],
[ 1.337500,   1.0711,   0.0000,   1.1306,     0.00  ],
[ 1.362500,   1.0647,   0.0000,   1.1275,     0.00  ],
[ 1.387500,   0.9922,   0.0000,   1.2030,    50.83  ],
[ 1.412500,   0.9741,   0.0000,   1.1536,     0.00  ],
[ 1.437500,   0.9923,   0.0000,   1.2851,  6663.88  ],
[ 1.462500,   1.8990,   1.0000,   1.4332,  3742.08  ],
[ 1.487500,   1.8639,   0.0000,   1.1538,   887.95  ],
[ 1.512500,   1.9022,   0.0000,   1.0728,     0.00  ],
[ 1.537500,   1.1643,   0.3892,   1.0690,     0.00  ],
[ 1.562500,   1.1853,   0.3382,   1.0710,     0.00  ],
[ 1.587500,   1.1939,   0.7337,   1.0676,     0.00  ],
[ 1.612500,   1.1855,   0.7274,   1.0692,     0.00  ],
[ 1.637500,   1.2032,   0.9399,   1.0678,     0.00  ],
[ 1.662500,   1.2118,   1.0308,   1.0695,     0.00  ],
[ 1.687500,   1.1916,   1.0267,   1.0704,     0.00  ],
[ 1.712500,   1.1726,   1.0959,   1.0694,     6.73  ],
[ 1.737500,   1.3024,   1.1697,   1.0673,     1.35  ],
[ 1.762500,   1.2483,   1.1115,   1.0686,    39.34  ],
[ 1.787500,   1.3181,   0.9399,   1.0776,   153.75  ],
[ 1.812500,   0.0000,   0.7468,   1.0843,   151.04  ],
[ 1.837500,   0.0000,   0.7300,   1.0891,   192.28  ],
[ 1.862500,   0.0000,   0.6231,   1.0928,   151.78  ],
[ 1.887500,   0.0000,   0.6660,   1.0962,   164.38  ],
[ 1.912500,   0.0000,   0.7777,   1.0925,   202.87  ],
[ 1.937500,   0.0000,   0.8338,   1.0923,   229.46  ],
[ 1.962500,   0.0000,   0.8240,   1.0919,   178.22  ],
[ 1.987500,   0.0000,   0.6841,   1.0991,   175.48  ],
[ 2.012500,   0.0000,   0.5993,   1.1045,   123.92  ],
[ 2.037500,   0.0000,   0.6250,   1.1134,    96.80  ],
[ 2.062500,   0.0000,   0.6207,   1.1160,   139.54  ],
[ 2.087500,   0.0000,   0.6770,   1.1173,   242.58  ],
[ 2.112500,   0.0000,   0.7299,   1.1179,   170.76  ],
[ 2.137500,   0.0000,   0.9429,   1.1115,   161.83  ],
[ 2.162500,   0.0000,   0.9943,   1.1079,   205.78  ],
[ 2.187500,   0.0000,   0.8126,   1.1149,   138.95  ],
[ 2.212500,   0.0000,   0.8273,   1.1173,   106.89  ],
[ 2.237500,   0.0000,   0.8111,   1.1224,   101.28  ],
[ 2.262500,   0.0000,   0.9059,   1.1230,   121.63  ],
[ 2.287500,   0.0000,   0.7795,   1.1286,   140.45  ],
[ 2.312500,   0.0000,   0.8174,   1.1260,   171.33  ],
[ 2.337500,   0.0000,   0.8410,   1.1212,   176.56  ],
[ 2.362500,   0.0000,   1.0673,   1.1149,   223.40  ],
[ 2.387500,   0.0000,   0.7754,   1.1229,   156.68  ],
[ 2.412500,   0.0000,   0.6620,   1.1271,   131.10  ],
[ 2.437500,   0.0000,   0.9470,   1.1259,   176.52  ],
[ 2.462500,   0.0000,   0.9901,   1.1277,   218.71  ],
[ 2.487500,   0.0000,   1.6550,   1.1664,   634.44  ]
]


#######################################################################
# 3x7 cluster size photons. (same as 3x5 weights)
CaloSwLongWeights_atlfast_v1_gam37 = [
               # w0        w3       escale    eoffset 
[ 0.012500,   1.2366,   2.0684,   1.0754,   104.65  ],
[ 0.037500,   1.2496,   2.3059,   1.0707,   180.67  ],
[ 0.062500,   1.2778,   2.1103,   1.0737,    77.09  ],
[ 0.087500,   1.1605,   1.9479,   1.0744,   115.28  ],
[ 0.112500,   1.1879,   2.0967,   1.0732,   101.19  ],
[ 0.137500,   1.1832,   2.0284,   1.0724,    74.83  ],
[ 0.162500,   1.1623,   1.9560,   1.0724,   126.48  ],
[ 0.187500,   1.1639,   1.8758,   1.0721,   168.02  ],
[ 0.212500,   1.2040,   2.0096,   1.0696,   113.91  ],
[ 0.237500,   1.1886,   2.0048,   1.0691,   115.86  ],
[ 0.262500,   1.2536,   1.9001,   1.0688,    96.66  ],
[ 0.287500,   1.2149,   2.0653,   1.0655,   147.32  ],
[ 0.312500,   1.1529,   1.8813,   1.0681,   152.10  ],
[ 0.337500,   1.1354,   1.9359,   1.0672,    74.76  ],
[ 0.362500,   1.1569,   1.9701,   1.0671,   104.97  ],
[ 0.387500,   1.1513,   1.7207,   1.0683,    83.58  ],
[ 0.412500,   1.1606,   1.5623,   1.0706,    69.54  ],
[ 0.437500,   1.1568,   1.3449,   1.0706,    35.25  ],
[ 0.462500,   1.1808,   1.2517,   1.072,     18.08  ],
[ 0.487500,   1.1528,   1.2519,   1.0693,    57.37  ],
[ 0.512500,   1.1745,   1.3477,   1.0684,    62.06  ],
[ 0.537500,   1.1638,   1.1854,   1.0685,    90.53  ],
[ 0.562500,   1.1998,   1.1687,   1.0685,    61.96  ],
[ 0.587500,   1.1569,   1.1163,   1.0685,    69.02  ],
[ 0.612500,   1.1405,   1.0606,   1.0683,   113.71  ],
[ 0.637500,   1.1402,   0.9642,   1.0695,    79.39  ],
[ 0.662500,   1.2031,   1.1779,   1.0662,    80.03  ],
[ 0.687500,   1.2893,   1.4292,   1.0621,    64.08  ],
[ 0.712500,   1.2935,   1.3535,   1.0644,    58.06  ],
[ 0.737500,   1.4294,   1.9157,   1.0651,     0.00  ],
[ 0.762500,   1.1575,   2.6176,   1.0692,   136.61  ],
[ 0.787500,   1.0230,   1.2167,   1.0500,   293.81  ],
[ 0.812500,   1.0347,   2.6494,   1.0517,   183.29  ],
[ 0.837500,   1.1184,   1.4828,   1.0708,    14.44  ],
[ 0.862500,   1.1791,   1.4896,   1.0682,     0.00  ],
[ 0.887500,   1.1900,   1.1609,   1.0715,     0.00  ],
[ 0.912500,   1.1838,   1.0168,   1.0731,     0.00  ],
[ 0.937500,   1.1842,   1.0360,   1.0730,     0.00  ],
[ 0.962500,   1.1653,   0.8954,   1.0767,     0.00  ],
[ 0.987500,   1.1495,   0.9831,   1.0764,     0.00  ],
[ 1.012500,   1.1436,   0.8156,   1.0785,     0.00  ],
[ 1.037500,   1.1563,   0.8833,   1.0784,     0.00  ],
[ 1.062500,   1.1448,   0.8273,   1.0796,     0.00  ],
[ 1.087500,   1.1570,   0.8368,   1.0800,     0.00  ],
[ 1.112500,   1.1476,   0.8148,   1.0787,    56.65  ],
[ 1.137500,   1.1274,   0.6277,   1.0851,     0.00  ],
[ 1.162500,   1.1298,   0.6566,   1.0851,    24.79  ],
[ 1.187500,   1.1192,   0.5677,   1.0894,     0.00  ],
[ 1.212500,   1.1175,   0.4432,   1.0909,     0.00  ],
[ 1.237500,   1.0944,   0.4741,   1.0926,     0.00  ],
[ 1.262500,   1.1011,   0.4594,   1.0920,     0.00  ],
[ 1.287500,   1.0788,   0.0000,   1.1041,     0.00  ],
[ 1.312500,   1.0743,   0.0000,   1.1045,     0.00  ],
[ 1.337500,   1.0711,   0.0000,   1.1306,     0.00  ],
[ 1.362500,   1.0647,   0.0000,   1.1275,     0.00  ],
[ 1.387500,   0.9922,   0.0000,   1.2030,    50.83  ],
[ 1.412500,   0.9741,   0.0000,   1.1536,     0.00  ],
[ 1.437500,   0.9923,   0.0000,   1.2851,  6663.88  ],
[ 1.462500,   1.8990,   1.0000,   1.4332,  3742.08  ],
[ 1.487500,   1.8639,   0.0000,   1.1538,   887.95  ],
[ 1.512500,   1.9022,   0.0000,   1.0728,     0.00  ],
[ 1.537500,   1.1643,   0.3892,   1.0690,     0.00  ],
[ 1.562500,   1.1853,   0.3382,   1.0710,     0.00  ],
[ 1.587500,   1.1939,   0.7337,   1.0676,     0.00  ],
[ 1.612500,   1.1855,   0.7274,   1.0692,     0.00  ],
[ 1.637500,   1.2032,   0.9399,   1.0678,     0.00  ],
[ 1.662500,   1.2118,   1.0308,   1.0695,     0.00  ],
[ 1.687500,   1.1916,   1.0267,   1.0704,     0.00  ],
[ 1.712500,   1.1726,   1.0959,   1.0694,     6.73  ],
[ 1.737500,   1.3024,   1.1697,   1.0673,     1.35  ],
[ 1.762500,   1.2483,   1.1115,   1.0686,    39.34  ],
[ 1.787500,   1.3181,   0.9399,   1.0776,   153.75  ],
[ 1.812500,   0.0000,   0.7468,   1.0843,   151.04  ],
[ 1.837500,   0.0000,   0.7300,   1.0891,   192.28  ],
[ 1.862500,   0.0000,   0.6231,   1.0928,   151.78  ],
[ 1.887500,   0.0000,   0.6660,   1.0962,   164.38  ],
[ 1.912500,   0.0000,   0.7777,   1.0925,   202.87  ],
[ 1.937500,   0.0000,   0.8338,   1.0923,   229.46  ],
[ 1.962500,   0.0000,   0.8240,   1.0919,   178.22  ],
[ 1.987500,   0.0000,   0.6841,   1.0991,   175.48  ],
[ 2.012500,   0.0000,   0.5993,   1.1045,   123.92  ],
[ 2.037500,   0.0000,   0.6250,   1.1134,    96.80  ],
[ 2.062500,   0.0000,   0.6207,   1.1160,   139.54  ],
[ 2.087500,   0.0000,   0.6770,   1.1173,   242.58  ],
[ 2.112500,   0.0000,   0.7299,   1.1179,   170.76  ],
[ 2.137500,   0.0000,   0.9429,   1.1115,   161.83  ],
[ 2.162500,   0.0000,   0.9943,   1.1079,   205.78  ],
[ 2.187500,   0.0000,   0.8126,   1.1149,   138.95  ],
[ 2.212500,   0.0000,   0.8273,   1.1173,   106.89  ],
[ 2.237500,   0.0000,   0.8111,   1.1224,   101.28  ],
[ 2.262500,   0.0000,   0.9059,   1.1230,   121.63  ],
[ 2.287500,   0.0000,   0.7795,   1.1286,   140.45  ],
[ 2.312500,   0.0000,   0.8174,   1.1260,   171.33  ],
[ 2.337500,   0.0000,   0.8410,   1.1212,   176.56  ],
[ 2.362500,   0.0000,   1.0673,   1.1149,   223.40  ],
[ 2.387500,   0.0000,   0.7754,   1.1229,   156.68  ],
[ 2.412500,   0.0000,   0.6620,   1.1271,   131.10  ],
[ 2.437500,   0.0000,   0.9470,   1.1259,   176.52  ],
[ 2.462500,   0.0000,   0.9901,   1.1277,   218.71  ],
[ 2.487500,   0.0000,   1.6550,   1.1664,   634.44  ]
]

#######################################################################


class CaloSwLongWeights_atlfast_v1_parms:
    eta_start_crack = 1.425
    eta_end_crack = 1.55
    etamax = 2.5
    use_raw_eta = False

    # If "preserve_offset" is set to True, then any offset that has
    # been applied prior to this correction is preserved. We do not want
    # to do this here because gapCorrections and lwcorrections are exclusive.
    preserve_offset = False

    region = CALOCORR_COMBINED2
    degree = 3
    correction = {'ele55' : CaloSwLongWeights_atlfast_v1_ele55,
                  'ele35' : CaloSwLongWeights_atlfast_v1_ele35,
                  'ele37' : CaloSwLongWeights_atlfast_v1_ele37,
                  'gam55' : CaloSwLongWeights_atlfast_v1_gam55,
                  'gam35' : CaloSwLongWeights_atlfast_v1_gam35,
                  'gam37' : CaloSwLongWeights_atlfast_v1_gam37,

                  # Use 5x5 for cluster sizes that aren't explicitly derived.
                  'ele33' : CaloSwLongWeights_atlfast_v1_ele55,
                  'ele57' : CaloSwLongWeights_atlfast_v1_ele55,
                  'ele77' : CaloSwLongWeights_atlfast_v1_ele55,
                  'gam33' : CaloSwLongWeights_atlfast_v1_gam55,
                  'gam57' : CaloSwLongWeights_atlfast_v1_gam55,
                  'gam77' : CaloSwLongWeights_atlfast_v1_gam55,
                  }
    
