# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration

# File: CaloClusterCorrection/python/CaloSwLongWeights_v5.py
# Created: July 2008, NK, SP
# Purpose: Longitudinal weights corrections, v5.
#
# These corrections were derived by Nicolas Kerschen and Stathes Paganis
# The geometry is ATLAS-CSC-05-00-00

from CaloClusterCorrection.common import *

#######################################################################
# 5x5 cluster size electrons.
CaloSwLongWeights_v5_ele55 = [
[ 0.012500,   1.2390,   2.0846,   1.0476,   213.99  ],
[ 0.037500,   1.1702,   1.7575,   1.0352,   155.85  ],
[ 0.062500,   1.1259,   1.6343,   1.0361,   124.12  ],
[ 0.087500,   1.1289,   1.6161,   1.0352,   176.74  ],
[ 0.112500,   1.1280,   1.7661,   1.0351,   132.06  ],
[ 0.137500,   1.1388,   1.6654,   1.0348,   109.69  ],
[ 0.162500,   1.1703,   1.6898,   1.0346,    99.32  ],
[ 0.187500,   1.1147,   1.6902,   1.0342,   121.35  ],
[ 0.212500,   1.0699,   1.6895,   1.0321,   199.09  ],
[ 0.237500,   1.1348,   1.6140,   1.0328,   116.05  ],
[ 0.262500,   1.1195,   1.6315,   1.0325,   111.55  ],
[ 0.287500,   1.1304,   1.7732,   1.0311,   108.97  ],
[ 0.312500,   1.1174,   1.6195,   1.0314,   160.64  ],
[ 0.337500,   1.1253,   1.5807,   1.0324,    97.38  ],
[ 0.362500,   1.0996,   1.5291,   1.0320,   153.75  ],
[ 0.387500,   1.1145,   1.4401,   1.0319,   159.75  ],
[ 0.412500,   1.0762,   1.3587,   1.0343,   138.28  ],
[ 0.437500,   1.1332,   1.4135,   1.0328,   129.46  ],
[ 0.462500,   1.0926,   1.1644,   1.0335,   158.80  ],
[ 0.487500,   1.0771,   1.1340,   1.0343,   163.70  ],
[ 0.512500,   1.1285,   1.1046,   1.0329,   160.11  ],
[ 0.537500,   1.1086,   1.0676,   1.0337,   172.98  ],
[ 0.562500,   1.1260,   0.8618,   1.0348,   139.77  ],
[ 0.587500,   1.1691,   1.1005,   1.0316,   104.27  ],
[ 0.612500,   1.1181,   0.9425,   1.0326,   148.84  ],
[ 0.637500,   1.1471,   1.0532,   1.0319,   155.48  ],
[ 0.662500,   1.1960,   0.9939,   1.0335,   233.28  ],
[ 0.687500,   1.2357,   0.9273,   1.0332,   211.84  ],
[ 0.712500,   1.2059,   0.9697,   1.0343,   166.10  ],
[ 0.737500,   1.2206,   0.8630,   1.0353,   313.50  ],
[ 0.762500,   1.2694,   0.8732,   1.0384,   312.68  ],
[ 0.787500,   1.4244,   1.4695,   1.0439,   477.07  ],
[ 0.812500,   1.4607,   2.8855,   1.0380,   692.86  ],
[ 0.837500,   1.3401,   1.5164,   1.0439,   419.56  ],
[ 0.862500,   1.3096,   1.1554,   1.0470,   564.03  ],
[ 0.887500,   1.2745,   0.7866,   1.0523,   591.79  ],
[ 0.912500,   1.3070,   0.9353,   1.0502,   674.58  ],
[ 0.937500,   1.2911,   0.4064,   1.0571,   543.60  ],
[ 0.962500,   1.3258,   0.6528,   1.0527,   597.39  ],
[ 0.987500,   1.3296,   0.8651,   1.0521,   587.16  ],
[ 1.012500,   1.3056,   0.9003,   1.0524,   646.67  ],
[ 1.037500,   1.3273,   0.6194,   1.0540,   644.06  ],
[ 1.062500,   1.3168,   0.8407,   1.0546,   563.44  ],
[ 1.087500,   1.2977,   0.3042,   1.0600,   705.29  ],
[ 1.112500,   1.3221,   0.7072,   1.0552,   499.96  ],
[ 1.137500,   1.3105,   0.5291,   1.0582,   507.21  ],
[ 1.162500,   1.2803,   0.4459,   1.0621,   690.14  ],
[ 1.187500,   1.2717,   0.3023,   1.0650,   813.09  ],
[ 1.212500,   1.2292,   0.4306,   1.0652,   630.36  ],
[ 1.237500,   1.2349,   0.1100,   1.0676,   849.24  ],
[ 1.262500,   1.2499,   1.0337,   1.0575,   777.23  ],
[ 1.287500,   1.2379,   0.2648,   1.0710,   953.35  ],
[ 1.312500,   1.2559,   0.0475,   1.0768,  1035.74  ],
[ 1.337500,   1.2927,   0.0000,   1.0751,  1227.18  ],
[ 1.362500,   1.2651,   0.0000,   1.0866,  1397.72  ],
[ 1.387500,   1.2463,   0.0000,   1.1046,  1320.56  ],
[ 1.412500,   1.1413,   2.5275,   1.1783,  1346.26  ],
[ 1.437500,   1.0550,   0.6514,   1.3318,  1659.39  ],
[ 1.462500,   1.4508,   0.0000,   1.5012,  2368.37  ],
[ 1.487500,   1.0000,   1.0000,   1.0000,     0.00  ],
[ 1.512500,   1.8724,   0.0000,   1.2256,  6678.39  ],
[ 1.537500,   3.3176,   4.8188,   0.8639,  4492.20  ],
[ 1.562500,   1.6295,   1.0928,   1.0316,  1769.67  ],
[ 1.587500,   1.4969,   0.9104,   1.0411,  1663.12  ],
[ 1.612500,   1.5762,   1.0707,   1.0314,  1914.88  ],
[ 1.637500,   1.8416,   1.4885,   1.0137,  1778.83  ],
[ 1.662500,   2.0373,   2.3206,   0.9960,  2468.89  ],
[ 1.687500,   2.4065,   1.4726,   0.9712,  2834.88  ],
[ 1.712500,   1.6986,   1.4688,   1.0187,   809.69  ],
[ 1.737500,   1.0492,   0.7768,   1.0603,   755.41  ],
[ 1.762500,   1.3879,   1.3087,   1.0328,   295.66  ],
[ 1.787500,   1.2988,   1.3610,   1.0369,   376.27  ],
[ 1.812500,   0.0000,   0.9128,   1.0715,  1510.09  ],
[ 1.837500,   0.0000,   0.0000,   1.0834,  1731.48  ],
[ 1.862500,   0.0000,   0.0000,   1.0973,  1259.12  ],
[ 1.887500,   0.0000,   0.8289,   1.0658,  1285.90  ],
[ 1.912500,   0.0000,   0.3426,   1.0696,  1157.56  ],
[ 1.937500,   0.0000,   0.3491,   1.0683,   960.62  ],
[ 1.962500,   0.0000,   0.7117,   1.0639,   917.14  ],
[ 1.987500,   0.0000,   0.8786,   1.0616,  1226.03  ],
[ 2.012500,   0.0000,   0.5981,   1.0729,  1125.34  ],
[ 2.037500,   0.0000,   0.5918,   1.0775,  1126.81  ],
[ 2.062500,   0.0000,   0.1689,   1.0783,  1243.26  ],
[ 2.087500,   0.0000,   0.1577,   1.0755,  1216.73  ],
[ 2.112500,   0.0000,   0.5082,   1.0714,  1274.68  ],
[ 2.137500,   0.0000,   0.1831,   1.0793,  1213.67  ],
[ 2.162500,   0.0000,   0.5614,   1.0739,  1124.78  ],
[ 2.187500,   0.0000,   0.4306,   1.0767,   975.30  ],
[ 2.212500,   0.0000,   0.4266,   1.0784,   950.76  ],
[ 2.237500,   0.0000,   0.4476,   1.0770,   889.05  ],
[ 2.262500,   0.0000,   0.7153,   1.0699,   845.34  ],
[ 2.287500,   0.0000,   0.6677,   1.0724,   772.32  ],
[ 2.312500,   0.0000,   0.7363,   1.0756,   797.66  ],
[ 2.337500,   0.0000,   0.8298,   1.0782,   699.72  ],
[ 2.362500,   0.0000,   0.8270,   1.0769,   665.95  ],
[ 2.387500,   0.0000,   0.7566,   1.0772,   671.49  ],
[ 2.412500,   0.0000,   0.5520,   1.0795,   664.48  ],
[ 2.437500,   0.0000,   0.5365,   1.0798,   748.75  ],
[ 2.462500,   0.0000,   0.5798,   1.0883,   979.33  ],
[ 2.487500,   0.0000,   0.2019,   1.1502,  1280.57  ]
]



#######################################################################
#  3x5 cluster size electrons. (same as 3x7 weights)
CaloSwLongWeights_v5_ele35 = [
[ 0.012500,   1.1726,   2.0961,   1.0619,   279.17  ],
[ 0.037500,   1.1319,   1.6419,   1.0502,   164.82  ],
[ 0.062500,   1.1255,   1.5955,   1.0501,   123.05  ],
[ 0.087500,   1.1134,   1.5646,   1.0499,   194.94  ],
[ 0.112500,   1.1264,   1.8305,   1.0486,   149.14  ],
[ 0.137500,   1.1170,   1.6628,   1.0494,   121.07  ],
[ 0.162500,   1.1405,   1.6612,   1.0492,   125.06  ],
[ 0.187500,   1.0977,   1.7029,   1.0486,   137.85  ],
[ 0.212500,   1.0615,   1.6352,   1.0469,   175.92  ],
[ 0.237500,   1.1172,   1.6850,   1.0469,   118.93  ],
[ 0.262500,   1.1079,   1.6659,   1.0463,   134.37  ],
[ 0.287500,   1.1172,   1.6581,   1.0460,   129.53  ],
[ 0.312500,   1.1066,   1.5687,   1.0461,   156.48  ],
[ 0.337500,   1.1299,   1.6230,   1.0465,    99.22  ],
[ 0.362500,   1.0797,   1.4627,   1.0468,   145.70  ],
[ 0.387500,   1.1066,   1.3764,   1.0470,   142.64  ],
[ 0.412500,   1.0992,   1.2348,   1.0492,   103.85  ],
[ 0.437500,   1.1260,   1.3332,   1.0478,   112.43  ],
[ 0.462500,   1.0693,   1.1486,   1.0484,   154.71  ],
[ 0.487500,   1.0793,   1.0923,   1.0490,   145.42  ],
[ 0.512500,   1.1073,   1.0651,   1.0476,   170.07  ],
[ 0.537500,   1.0943,   1.0001,   1.0488,   133.84  ],
[ 0.562500,   1.1099,   0.8505,   1.0494,   126.63  ],
[ 0.587500,   1.1434,   1.0497,   1.0469,   107.43  ],
[ 0.612500,   1.1179,   0.8732,   1.0473,   161.67  ],
[ 0.637500,   1.1197,   0.9718,   1.0475,   164.70  ],
[ 0.662500,   1.1777,   0.9191,   1.0485,   206.86  ],
[ 0.687500,   1.2048,   0.9240,   1.0480,   178.66  ],
[ 0.712500,   1.1873,   0.9031,   1.0487,   152.16  ],
[ 0.737500,   1.1840,   0.8229,   1.0496,   288.60  ],
[ 0.762500,   1.2486,   0.7950,   1.0513,   272.52  ],
[ 0.787500,   1.3759,   1.5886,   1.0578,   422.71  ],
[ 0.812500,   1.4210,   2.5368,   1.0558,   478.29  ],
[ 0.837500,   1.3115,   1.7018,   1.0614,   328.60  ],
[ 0.862500,   1.3031,   1.0591,   1.0633,   447.27  ],
[ 0.887500,   1.2510,   0.5574,   1.0696,   563.35  ],
[ 0.912500,   1.2900,   0.6957,   1.0680,   494.58  ],
[ 0.937500,   1.2892,   0.4412,   1.0718,   462.20  ],
[ 0.962500,   1.3157,   0.5955,   1.0681,   513.12  ],
[ 0.987500,   1.2962,   0.8004,   1.0693,   502.74  ],
[ 1.012500,   1.3003,   0.9145,   1.0668,   483.17  ],
[ 1.037500,   1.3086,   0.6060,   1.0699,   468.88  ],
[ 1.062500,   1.2948,   0.6745,   1.0706,   542.41  ],
[ 1.087500,   1.2996,   0.4673,   1.0722,   570.74  ],
[ 1.112500,   1.2806,   0.6771,   1.0717,   409.13  ],
[ 1.137500,   1.2908,   0.3339,   1.0755,   475.39  ],
[ 1.162500,   1.2771,   0.4599,   1.0758,   622.78  ],
[ 1.187500,   1.2612,   0.3123,   1.0791,   681.79  ],
[ 1.212500,   1.2287,   0.4969,   1.0764,   516.53  ],
[ 1.237500,   1.2036,   0.1111,   1.0829,   650.66  ],
[ 1.262500,   1.2114,   0.2555,   1.0826,   558.49  ],
[ 1.287500,   1.2115,   0.4246,   1.0818,   716.62  ],
[ 1.312500,   1.2285,   0.1116,   1.0876,   765.42  ],
[ 1.337500,   1.2419,   0.0000,   1.0897,   873.41  ],
[ 1.362500,   1.2793,   0.0000,   1.0783,   931.25  ],
[ 1.387500,   1.2766,   0.0000,   1.1555,  1144.31  ],
[ 1.412500,   1.1295,   1.1427,   1.1822,   999.40  ],
[ 1.437500,   1.0466,   3.6084,   1.3332,  1261.04  ],
[ 1.462500,   1.4513,   0.0000,   1.5093,  1845.27  ],
[ 1.487500,   1.0000,   1.0000,   1.0000,     0.00  ],
[ 1.512500,   1.2990,   0.0000,   1.4697,  6232.77  ],
[ 1.537500,   3.1296,   4.0353,   0.8901,  2854.24  ],
[ 1.562500,   1.5865,   1.2907,   1.0403,  1338.39  ],
[ 1.587500,   1.4645,   0.8798,   1.0526,  1345.53  ],
[ 1.612500,   1.5739,   1.2047,   1.0407,  1422.47  ],
[ 1.637500,   1.7665,   1.3677,   1.0305,  1274.61  ],
[ 1.662500,   2.0453,   2.1407,   1.0072,  2016.56  ],
[ 1.687500,   2.5428,   2.0164,   0.9704,  2226.92  ],
[ 1.712500,   1.5628,   1.5231,   1.0405,   397.30  ],
[ 1.737500,   1.0688,   0.8388,   1.0745,   680.20  ],
[ 1.762500,   1.3221,   1.2675,   1.0497,   212.04  ],
[ 1.787500,   1.2634,   1.3564,   1.0532,   296.10  ],
[ 1.812500,   0.0000,   0.8514,   1.0881,  1484.51  ],
[ 1.837500,   0.0000,   0.0000,   1.1005,  1601.46  ],
[ 1.862500,   0.0000,   0.0000,   1.1148,  1252.55  ],
[ 1.887500,   0.0000,   0.7791,   1.0845,  1143.76  ],
[ 1.912500,   0.0000,   0.3186,   1.0894,  1040.59  ],
[ 1.937500,   0.0000,   0.4597,   1.0864,   889.15  ],
[ 1.962500,   0.0000,   0.8429,   1.0817,   888.68  ],
[ 1.987500,   0.0000,   0.8760,   1.0815,  1178.38  ],
[ 2.012500,   0.0000,   0.4807,   1.0942,   988.72  ],
[ 2.037500,   0.0000,   0.5720,   1.0986,  1024.01  ],
[ 2.062500,   0.0000,   0.2687,   1.0988,  1182.37  ],
[ 2.087500,   0.0000,   0.2406,   1.0965,  1103.27  ],
[ 2.112500,   0.0000,   0.5218,   1.0937,  1148.97  ],
[ 2.137500,   0.0000,   0.2214,   1.1018,  1084.40  ],
[ 2.162500,   0.0000,   0.5758,   1.0974,  1029.42  ],
[ 2.187500,   0.0000,   0.4799,   1.0996,   881.26  ],
[ 2.212500,   0.0000,   0.3768,   1.1032,   857.99  ],
[ 2.237500,   0.0000,   0.4834,   1.1016,   802.30  ],
[ 2.262500,   0.0000,   0.6157,   1.0965,   751.20  ],
[ 2.287500,   0.0000,   0.7436,   1.0970,   757.82  ],
[ 2.312500,   0.0000,   0.7701,   1.1022,   720.44  ],
[ 2.337500,   0.0000,   0.8407,   1.1060,   648.96  ],
[ 2.362500,   0.0000,   0.7962,   1.1057,   635.50  ],
[ 2.387500,   0.0000,   0.7612,   1.1059,   646.92  ],
[ 2.412500,   0.0000,   0.5343,   1.1097,   629.98  ],
[ 2.437500,   0.0000,   0.5670,   1.1088,   714.57  ],
[ 2.462500,   0.0000,   0.4751,   1.1128,   764.55  ],
[ 2.487500,   0.0000,   0.0000,   1.1554,  1106.84  ]
]



#######################################################################
# 3x7 cluster size electrons.
CaloSwLongWeights_v5_ele37 = [
[ 0.012500,   1.1726,   2.0961,   1.0619,   279.17  ],
[ 0.037500,   1.1319,   1.6419,   1.0502,   164.82  ],
[ 0.062500,   1.1255,   1.5955,   1.0501,   123.05  ],
[ 0.087500,   1.1134,   1.5646,   1.0499,   194.94  ],
[ 0.112500,   1.1264,   1.8305,   1.0486,   149.14  ],
[ 0.137500,   1.1170,   1.6628,   1.0494,   121.07  ],
[ 0.162500,   1.1405,   1.6612,   1.0492,   125.06  ],
[ 0.187500,   1.0977,   1.7029,   1.0486,   137.85  ],
[ 0.212500,   1.0615,   1.6352,   1.0469,   175.92  ],
[ 0.237500,   1.1172,   1.6850,   1.0469,   118.93  ],
[ 0.262500,   1.1079,   1.6659,   1.0463,   134.37  ],
[ 0.287500,   1.1172,   1.6581,   1.0460,   129.53  ],
[ 0.312500,   1.1066,   1.5687,   1.0461,   156.48  ],
[ 0.337500,   1.1299,   1.6230,   1.0465,    99.22  ],
[ 0.362500,   1.0797,   1.4627,   1.0468,   145.70  ],
[ 0.387500,   1.1066,   1.3764,   1.0470,   142.64  ],
[ 0.412500,   1.0992,   1.2348,   1.0492,   103.85  ],
[ 0.437500,   1.1260,   1.3332,   1.0478,   112.43  ],
[ 0.462500,   1.0693,   1.1486,   1.0484,   154.71  ],
[ 0.487500,   1.0793,   1.0923,   1.0490,   145.42  ],
[ 0.512500,   1.1073,   1.0651,   1.0476,   170.07  ],
[ 0.537500,   1.0943,   1.0001,   1.0488,   133.84  ],
[ 0.562500,   1.1099,   0.8505,   1.0494,   126.63  ],
[ 0.587500,   1.1434,   1.0497,   1.0469,   107.43  ],
[ 0.612500,   1.1179,   0.8732,   1.0473,   161.67  ],
[ 0.637500,   1.1197,   0.9718,   1.0475,   164.70  ],
[ 0.662500,   1.1777,   0.9191,   1.0485,   206.86  ],
[ 0.687500,   1.2048,   0.9240,   1.0480,   178.66  ],
[ 0.712500,   1.1873,   0.9031,   1.0487,   152.16  ],
[ 0.737500,   1.1840,   0.8229,   1.0496,   288.60  ],
[ 0.762500,   1.2486,   0.7950,   1.0513,   272.52  ],
[ 0.787500,   1.3759,   1.5886,   1.0578,   422.71  ],
[ 0.812500,   1.4210,   2.5368,   1.0558,   478.29  ],
[ 0.837500,   1.3115,   1.7018,   1.0614,   328.60  ],
[ 0.862500,   1.3031,   1.0591,   1.0633,   447.27  ],
[ 0.887500,   1.2510,   0.5574,   1.0696,   563.35  ],
[ 0.912500,   1.2900,   0.6957,   1.0680,   494.58  ],
[ 0.937500,   1.2892,   0.4412,   1.0718,   462.20  ],
[ 0.962500,   1.3157,   0.5955,   1.0681,   513.12  ],
[ 0.987500,   1.2962,   0.8004,   1.0693,   502.74  ],
[ 1.012500,   1.3003,   0.9145,   1.0668,   483.17  ],
[ 1.037500,   1.3086,   0.6060,   1.0699,   468.88  ],
[ 1.062500,   1.2948,   0.6745,   1.0706,   542.41  ],
[ 1.087500,   1.2996,   0.4673,   1.0722,   570.74  ],
[ 1.112500,   1.2806,   0.6771,   1.0717,   409.13  ],
[ 1.137500,   1.2908,   0.3339,   1.0755,   475.39  ],
[ 1.162500,   1.2771,   0.4599,   1.0758,   622.78  ],
[ 1.187500,   1.2612,   0.3123,   1.0791,   681.79  ],
[ 1.212500,   1.2287,   0.4969,   1.0764,   516.53  ],
[ 1.237500,   1.2036,   0.1111,   1.0829,   650.66  ],
[ 1.262500,   1.2114,   0.2555,   1.0826,   558.49  ],
[ 1.287500,   1.2115,   0.4246,   1.0818,   716.62  ],
[ 1.312500,   1.2285,   0.1116,   1.0876,   765.42  ],
[ 1.337500,   1.2419,   0.0000,   1.0897,   873.41  ],
[ 1.362500,   1.2793,   0.0000,   1.0783,   931.25  ],
[ 1.387500,   1.2766,   0.0000,   1.1555,  1144.31  ],
[ 1.412500,   1.1295,   1.1427,   1.1822,   999.40  ],
[ 1.437500,   1.0466,   3.6084,   1.3332,  1261.04  ],
[ 1.462500,   1.4513,   0.0000,   1.5093,  1845.27  ],
[ 1.487500,   1.0000,   1.0000,   1.0000,     0.00  ],
[ 1.512500,   1.2990,   0.0000,   1.4697,  6232.77  ],
[ 1.537500,   3.1296,   4.0353,   0.8901,  2854.24  ],
[ 1.562500,   1.5865,   1.2907,   1.0403,  1338.39  ],
[ 1.587500,   1.4645,   0.8798,   1.0526,  1345.53  ],
[ 1.612500,   1.5739,   1.2047,   1.0407,  1422.47  ],
[ 1.637500,   1.7665,   1.3677,   1.0305,  1274.61  ],
[ 1.662500,   2.0453,   2.1407,   1.0072,  2016.56  ],
[ 1.687500,   2.5428,   2.0164,   0.9704,  2226.92  ],
[ 1.712500,   1.5628,   1.5231,   1.0405,   397.30  ],
[ 1.737500,   1.0688,   0.8388,   1.0745,   680.20  ],
[ 1.762500,   1.3221,   1.2675,   1.0497,   212.04  ],
[ 1.787500,   1.2634,   1.3564,   1.0532,   296.10  ],
[ 1.812500,   0.0000,   0.8514,   1.0881,  1484.51  ],
[ 1.837500,   0.0000,   0.0000,   1.1005,  1601.46  ],
[ 1.862500,   0.0000,   0.0000,   1.1148,  1252.55  ],
[ 1.887500,   0.0000,   0.7791,   1.0845,  1143.76  ],
[ 1.912500,   0.0000,   0.3186,   1.0894,  1040.59  ],
[ 1.937500,   0.0000,   0.4597,   1.0864,   889.15  ],
[ 1.962500,   0.0000,   0.8429,   1.0817,   888.68  ],
[ 1.987500,   0.0000,   0.8760,   1.0815,  1178.38  ],
[ 2.012500,   0.0000,   0.4807,   1.0942,   988.72  ],
[ 2.037500,   0.0000,   0.5720,   1.0986,  1024.01  ],
[ 2.062500,   0.0000,   0.2687,   1.0988,  1182.37  ],
[ 2.087500,   0.0000,   0.2406,   1.0965,  1103.27  ],
[ 2.112500,   0.0000,   0.5218,   1.0937,  1148.97  ],
[ 2.137500,   0.0000,   0.2214,   1.1018,  1084.40  ],
[ 2.162500,   0.0000,   0.5758,   1.0974,  1029.42  ],
[ 2.187500,   0.0000,   0.4799,   1.0996,   881.26  ],
[ 2.212500,   0.0000,   0.3768,   1.1032,   857.99  ],
[ 2.237500,   0.0000,   0.4834,   1.1016,   802.30  ],
[ 2.262500,   0.0000,   0.6157,   1.0965,   751.20  ],
[ 2.287500,   0.0000,   0.7436,   1.0970,   757.82  ],
[ 2.312500,   0.0000,   0.7701,   1.1022,   720.44  ],
[ 2.337500,   0.0000,   0.8407,   1.1060,   648.96  ],
[ 2.362500,   0.0000,   0.7962,   1.1057,   635.50  ],
[ 2.387500,   0.0000,   0.7612,   1.1059,   646.92  ],
[ 2.412500,   0.0000,   0.5343,   1.1097,   629.98  ],
[ 2.437500,   0.0000,   0.5670,   1.1088,   714.57  ],
[ 2.462500,   0.0000,   0.4751,   1.1128,   764.55  ],
[ 2.487500,   0.0000,   0.0000,   1.1554,  1106.84  ]
]


#######################################################################
# 5x5 cluster size photons.
CaloSwLongWeights_v5_gam55 = [
[ 0.012500,   1.3831,   2.2642,   1.0416,   83.21  ],
[ 0.037500,   1.2342,   1.9064,   1.0334,    0.00  ],
[ 0.062500,   1.2403,   1.9413,   1.0310,    7.42  ],
[ 0.087500,   1.2359,   1.9406,   1.0316,   23.40  ],
[ 0.112500,   1.2495,   1.9983,   1.0307,   34.03  ],
[ 0.137500,   1.2573,   1.9510,   1.0296,   54.05  ],
[ 0.162500,   1.2004,   1.8324,   1.0305,   57.08  ],
[ 0.187500,   1.2140,   1.9750,   1.0302,    0.00  ],
[ 0.212500,   1.2368,   2.0824,   1.0274,   25.59  ],
[ 0.237500,   1.2199,   1.8893,   1.0287,    0.00  ],
[ 0.262500,   1.1455,   1.8157,   1.0292,   25.70  ],
[ 0.287500,   1.2067,   1.8691,   1.0277,   42.05  ],
[ 0.312500,   1.1563,   1.7892,   1.0284,   16.83  ],
[ 0.337500,   1.1839,   1.7434,   1.0287,    2.80  ],
[ 0.362500,   1.2021,   1.8119,   1.0275,   21.55  ],
[ 0.387500,   1.2502,   1.9279,   1.0265,   24.15  ],
[ 0.412500,   1.1840,   1.6614,   1.0294,   17.32  ],
[ 0.437500,   1.1751,   1.6016,   1.0281,   89.47  ],
[ 0.462500,   1.2033,   1.6899,   1.0279,   33.25  ],
[ 0.487500,   1.2252,   1.4174,   1.0291,   16.18  ],
[ 0.512500,   1.2367,   1.3829,   1.0279,   25.22  ],
[ 0.537500,   1.1918,   1.3536,   1.0272,   79.69  ],
[ 0.562500,   1.1755,   1.2680,   1.0272,   72.49  ],
[ 0.587500,   1.2249,   1.2449,   1.0267,   39.85  ],
[ 0.612500,   1.2266,   1.1638,   1.0281,   20.16  ],
[ 0.637500,   1.2212,   1.0331,   1.0280,    0.00  ],
[ 0.662500,   1.2372,   1.0612,   1.0283,    0.00  ],
[ 0.687500,   1.2644,   1.0331,   1.0272,   23.97  ],
[ 0.712500,   1.2420,   0.9015,   1.0292,    0.00  ],
[ 0.737500,   1.2450,   0.8038,   1.0311,    0.00  ],
[ 0.762500,   1.2458,   0.8879,   1.0326,   11.69  ],
[ 0.787500,   1.3788,   1.1138,   1.0405,    0.00  ],
[ 0.812500,   1.4211,   2.1600,   1.0363,   43.74  ],
[ 0.837500,   1.2923,   1.3416,   1.0379,    0.00  ],
[ 0.862500,   1.2816,   1.4698,   1.0377,    0.00  ],
[ 0.887500,   1.2831,   1.2856,   1.0380,    0.00  ],
[ 0.912500,   1.2521,   1.1569,   1.0396,    0.00  ],
[ 0.937500,   1.2605,   1.1141,   1.0385,    0.00  ],
[ 0.962500,   1.2393,   1.0503,   1.0398,    0.00  ],
[ 0.987500,   1.2466,   0.8722,   1.0416,    0.00  ],
[ 1.012500,   1.2764,   0.9630,   1.0400,    0.00  ],
[ 1.037500,   1.2592,   0.8500,   1.0424,    0.00  ],
[ 1.062500,   1.2935,   0.7922,   1.0429,    0.00  ],
[ 1.087500,   1.3003,   0.8501,   1.0423,    0.00  ],
[ 1.112500,   1.2845,   0.8956,   1.0418,    0.00  ],
[ 1.137500,   1.2967,   0.7745,   1.0441,    0.00  ],
[ 1.162500,   1.2842,   0.7224,   1.0465,    0.00  ],
[ 1.187500,   1.2948,   1.0498,   1.0413,    0.00  ],
[ 1.212500,   1.2329,   0.6979,   1.0463,    0.00  ],
[ 1.237500,   1.2298,   0.4115,   1.0522,    0.00  ],
[ 1.262500,   1.2089,   0.0131,   1.0612,    0.00  ],
[ 1.287500,   1.2612,   0.6237,   1.0510,    0.00  ],
[ 1.312500,   1.2533,   0.2761,   1.0572,    0.00  ],
[ 1.337500,   1.2450,   0.1477,   1.0621,    0.00  ],
[ 1.362500,   1.2969,   0.0000,   1.0716,    0.00  ],
[ 1.387500,   1.2010,   0.0000,   1.0951,    0.00  ],
[ 1.412500,   1.0961,   0.0000,   1.1762,    6.07  ],
[ 1.437500,   0.9653,   2.6319,   1.3562,  323.77  ],
[ 1.462500,   1.1328,   0.0000,   1.6836, 2972.03  ],
[ 1.487500,   1.0000,   1.0000,   1.0000,    0.00  ],
[ 1.512500,   2.2146,   0.0000,   1.0738, 4903.17  ],
[ 1.537500,   2.3812,   2.0610,   0.9593,  347.27  ],
[ 1.562500,   1.3405,   0.7829,   1.0263,    0.00  ],
[ 1.587500,   1.3897,   1.0295,   1.0261,    0.00  ],
[ 1.612500,   1.3203,   1.1321,   1.0257,    0.00  ],
[ 1.637500,   1.3415,   1.2541,   1.0233,    0.00  ],
[ 1.662500,   1.5808,   1.3396,   1.0182,    0.00  ],
[ 1.687500,   1.9167,   2.0134,   0.9974,    0.00  ],
[ 1.712500,   1.1974,   1.2405,   1.0290,    0.00  ],
[ 1.737500,   1.1227,   1.2916,   1.0294,    0.00  ],
[ 1.762500,   1.1194,   1.2475,   1.0286,    0.00  ],
[ 1.787500,   1.1691,   1.2898,   1.0279,    0.00  ],
[ 1.812500,   0.0000,   0.9397,   1.0457,    0.00  ],
[ 1.837500,   0.0000,   0.8743,   1.0508,    0.00  ],
[ 1.862500,   0.0000,   1.0764,   1.0427,   70.58  ],
[ 1.887500,   0.0000,   1.0562,   1.0446,   29.37  ],
[ 1.912500,   0.0000,   1.0039,   1.0456,   35.62  ],
[ 1.937500,   0.0000,   1.0070,   1.0447,   44.73  ],
[ 1.962500,   0.0000,   1.1182,   1.0439,   93.52  ],
[ 1.987500,   0.0000,   1.0496,   1.0458,    9.69  ],
[ 2.012500,   0.0000,   1.0578,   1.0498,  130.20  ],
[ 2.037500,   0.0000,   0.9513,   1.0539,  114.23  ],
[ 2.062500,   0.0000,   0.9242,   1.0520,  107.42  ],
[ 2.087500,   0.0000,   0.9420,   1.0507,   68.93  ],
[ 2.112500,   0.0000,   0.8343,   1.0534,   58.19  ],
[ 2.137500,   0.0000,   0.7497,   1.0576,   50.07  ],
[ 2.162500,   0.0000,   0.7823,   1.0589,   59.21  ],
[ 2.187500,   0.0000,   0.8397,   1.0591,  113.25  ],
[ 2.212500,   0.0000,   0.9340,   1.0593,   45.54  ],
[ 2.237500,   0.0000,   0.9282,   1.0596,   90.12  ],
[ 2.262500,   0.0000,   0.9415,   1.0596,   39.91  ],
[ 2.287500,   0.0000,   1.0516,   1.0575,   54.98  ],
[ 2.312500,   0.0000,   1.0364,   1.0619,   68.63  ],
[ 2.337500,   0.0000,   1.0171,   1.0669,   35.61  ],
[ 2.362500,   0.0000,   0.9068,   1.0679,   67.93  ],
[ 2.387500,   0.0000,   0.9778,   1.0648,   44.19  ],
[ 2.412500,   0.0000,   1.0286,   1.0638,   84.44  ],
[ 2.437500,   0.0000,   0.8938,   1.0650,   96.10  ],
[ 2.462500,   0.0000,   0.8683,   1.0743,   87.02  ],
[ 2.487500,   0.0000,   0.8005,   1.1245,  307.63  ]
]


#######################################################################
# 3x5 cluster size photons.
CaloSwLongWeights_v5_gam35 = [
[ 0.012500,   1.4186,   2.2335,   1.0614,   81.53  ],
[ 0.037500,   1.2713,   1.8493,   1.0520,    0.00  ],
[ 0.062500,   1.2846,   2.0031,   1.0486,    8.54  ],
[ 0.087500,   1.2864,   1.9154,   1.0504,   25.05  ],
[ 0.112500,   1.2777,   1.9656,   1.0494,   49.38  ],
[ 0.137500,   1.2750,   1.8522,   1.0494,   42.34  ],
[ 0.162500,   1.2373,   1.8402,   1.0494,   37.97  ],
[ 0.187500,   1.2400,   1.9029,   1.0498,    0.00  ],
[ 0.212500,   1.2697,   2.0380,   1.0463,   27.78  ],
[ 0.237500,   1.2443,   1.8284,   1.0480,    0.29  ],
[ 0.262500,   1.1982,   1.7991,   1.0477,   15.30  ],
[ 0.287500,   1.2428,   1.8294,   1.0468,   44.58  ],
[ 0.312500,   1.2236,   1.8105,   1.0468,   13.24  ],
[ 0.337500,   1.2339,   1.7390,   1.0475,    0.03  ],
[ 0.362500,   1.2310,   1.6782,   1.0469,   20.63  ],
[ 0.387500,   1.2721,   1.8400,   1.0461,   10.00  ],
[ 0.412500,   1.2112,   1.5788,   1.0486,   35.93  ],
[ 0.437500,   1.1940,   1.5520,   1.0473,   85.80  ],
[ 0.462500,   1.2167,   1.6082,   1.0474,   39.93  ],
[ 0.487500,   1.2241,   1.3337,   1.0487,   21.93  ],
[ 0.512500,   1.2637,   1.3467,   1.0472,   24.97  ],
[ 0.537500,   1.2126,   1.2365,   1.0476,   65.79  ],
[ 0.562500,   1.2075,   1.2076,   1.0467,   70.11  ],
[ 0.587500,   1.2327,   1.1944,   1.0473,    0.00  ],
[ 0.612500,   1.2363,   1.1210,   1.0475,   18.52  ],
[ 0.637500,   1.2305,   0.9843,   1.0480,    8.53  ],
[ 0.662500,   1.2510,   0.9832,   1.0484,    7.80  ],
[ 0.687500,   1.2878,   1.0163,   1.0471,   13.25  ],
[ 0.712500,   1.2629,   0.8242,   1.0493,    0.00  ],
[ 0.737500,   1.2569,   0.7574,   1.0512,    0.00  ],
[ 0.762500,   1.2782,   0.9138,   1.0503,   13.39  ],
[ 0.787500,   1.3824,   1.1470,   1.0617,    0.00  ],
[ 0.812500,   1.4473,   2.0213,   1.0580,   42.80  ],
[ 0.837500,   1.3038,   1.3784,   1.0620,    0.00  ],
[ 0.862500,   1.2911,   1.3270,   1.0615,    0.00  ],
[ 0.887500,   1.3129,   1.2210,   1.0611,    0.00  ],
[ 0.912500,   1.2777,   1.0814,   1.0626,    0.00  ],
[ 0.937500,   1.2811,   1.0702,   1.0619,    0.00  ],
[ 0.962500,   1.2612,   0.9885,   1.0626,    0.00  ],
[ 0.987500,   1.2604,   0.8087,   1.0655,    0.00  ],
[ 1.012500,   1.2983,   0.8930,   1.0636,    0.00  ],
[ 1.037500,   1.2850,   0.7972,   1.0659,    0.00  ],
[ 1.062500,   1.3138,   0.7260,   1.0663,    0.00  ],
[ 1.087500,   1.3197,   0.8153,   1.0661,    0.00  ],
[ 1.112500,   1.3175,   0.8557,   1.0645,    0.00  ],
[ 1.137500,   1.3204,   0.7447,   1.0680,    0.00  ],
[ 1.162500,   1.3159,   0.6524,   1.0696,    0.00  ],
[ 1.187500,   1.3158,   1.0359,   1.0655,    0.00  ],
[ 1.212500,   1.2604,   0.6751,   1.0688,    0.00  ],
[ 1.237500,   1.2405,   0.3713,   1.0767,    0.00  ],
[ 1.262500,   1.2311,   0.0000,   1.0838,    0.00  ],
[ 1.287500,   1.2726,   0.6314,   1.0749,    0.00  ],
[ 1.312500,   1.2633,   0.1106,   1.0824,    0.00  ],
[ 1.337500,   1.2617,   0.0062,   1.0869,    0.00  ],
[ 1.362500,   1.3370,   0.0000,   1.0863,    0.00  ],
[ 1.387500,   1.3030,   0.0000,   1.1675,    3.40  ],
[ 1.412500,   1.1261,   0.0000,   1.1924,    0.00  ],
[ 1.437500,   0.9909,   0.3812,   1.3719,  249.50  ],
[ 1.462500,   1.1236,   0.0000,   1.7456, 3018.05  ],
[ 1.487500,   1.0000,   1.0000,   1.0000,    0.00  ],
[ 1.512500,   2.0474,   0.0000,   1.1349, 6891.45  ],
[ 1.537500,   2.4966,   2.1395,   0.9680,  499.05  ],
[ 1.562500,   1.3810,   0.7912,   1.0454,    0.00  ],
[ 1.587500,   1.3942,   0.9695,   1.0477,    0.00  ],
[ 1.612500,   1.3345,   1.0945,   1.0472,    0.00  ],
[ 1.637500,   1.3676,   1.2421,   1.0452,    0.00  ],
[ 1.662500,   1.6021,   1.3375,   1.0403,    0.00  ],
[ 1.687500,   2.0264,   2.1715,   1.0147,   14.47  ],
[ 1.712500,   1.2049,   1.1964,   1.0530,    0.00  ],
[ 1.737500,   1.1249,   1.2936,   1.0536,    0.00  ],
[ 1.762500,   1.1347,   1.2246,   1.0532,    0.00  ],
[ 1.787500,   1.1701,   1.2867,   1.0533,    0.00  ],
[ 1.812500,   0.0000,   0.8555,   1.0727,    0.00  ],
[ 1.837500,   0.0000,   0.8239,   1.0782,   17.40  ],
[ 1.862500,   0.0000,   1.0117,   1.0707,   46.28  ],
[ 1.887500,   0.0000,   1.0167,   1.0732,   12.92  ],
[ 1.912500,   0.0000,   1.0164,   1.0736,   67.46  ],
[ 1.937500,   0.0000,   0.9943,   1.0741,   53.60  ],
[ 1.962500,   0.0000,   1.0622,   1.0737,   72.03  ],
[ 1.987500,   0.0000,   1.0107,   1.0767,    6.60  ],
[ 2.012500,   0.0000,   0.9800,   1.0817,  111.64  ],
[ 2.037500,   0.0000,   0.9055,   1.0865,   87.25  ],
[ 2.062500,   0.0000,   0.9485,   1.0842,  108.32  ],
[ 2.087500,   0.0000,   0.9057,   1.0838,   55.43  ],
[ 2.112500,   0.0000,   0.7463,   1.0873,   90.90  ],
[ 2.137500,   0.0000,   0.7154,   1.0923,   39.86  ],
[ 2.162500,   0.0000,   0.6399,   1.0955,   64.02  ],
[ 2.187500,   0.0000,   0.7927,   1.0948,  122.85  ],
[ 2.212500,   0.0000,   0.8830,   1.0953,   84.57  ],
[ 2.237500,   0.0000,   0.8920,   1.0963,   86.02  ],
[ 2.262500,   0.0000,   0.8954,   1.0962,   70.49  ],
[ 2.287500,   0.0000,   1.0388,   1.0957,   23.22  ],
[ 2.312500,   0.0000,   0.9922,   1.1010,   73.24  ],
[ 2.337500,   0.0000,   0.9590,   1.1072,   44.06  ],
[ 2.362500,   0.0000,   0.8677,   1.1084,   55.37  ],
[ 2.387500,   0.0000,   0.9179,   1.1066,   54.89  ],
[ 2.412500,   0.0000,   0.9489,   1.1064,   70.17  ],
[ 2.437500,   0.0000,   0.8526,   1.1068,   68.58  ],
[ 2.462500,   0.0000,   0.7528,   1.1091,   58.29  ],
[ 2.487500,   0.0000,   0.6683,   1.1425,  323.50  ]
]



#######################################################################
# 3x7 cluster size photons. (same as 3x5)
CaloSwLongWeights_v5_gam37 = [
[ 0.012500,   1.4186,   2.2335,   1.0614,   81.53  ],
[ 0.037500,   1.2713,   1.8493,   1.0520,    0.00  ],
[ 0.062500,   1.2846,   2.0031,   1.0486,    8.54  ],
[ 0.087500,   1.2864,   1.9154,   1.0504,   25.05  ],
[ 0.112500,   1.2777,   1.9656,   1.0494,   49.38  ],
[ 0.137500,   1.2750,   1.8522,   1.0494,   42.34  ],
[ 0.162500,   1.2373,   1.8402,   1.0494,   37.97  ],
[ 0.187500,   1.2400,   1.9029,   1.0498,    0.00  ],
[ 0.212500,   1.2697,   2.0380,   1.0463,   27.78  ],
[ 0.237500,   1.2443,   1.8284,   1.0480,    0.29  ],
[ 0.262500,   1.1982,   1.7991,   1.0477,   15.30  ],
[ 0.287500,   1.2428,   1.8294,   1.0468,   44.58  ],
[ 0.312500,   1.2236,   1.8105,   1.0468,   13.24  ],
[ 0.337500,   1.2339,   1.7390,   1.0475,    0.03  ],
[ 0.362500,   1.2310,   1.6782,   1.0469,   20.63  ],
[ 0.387500,   1.2721,   1.8400,   1.0461,   10.00  ],
[ 0.412500,   1.2112,   1.5788,   1.0486,   35.93  ],
[ 0.437500,   1.1940,   1.5520,   1.0473,   85.80  ],
[ 0.462500,   1.2167,   1.6082,   1.0474,   39.93  ],
[ 0.487500,   1.2241,   1.3337,   1.0487,   21.93  ],
[ 0.512500,   1.2637,   1.3467,   1.0472,   24.97  ],
[ 0.537500,   1.2126,   1.2365,   1.0476,   65.79  ],
[ 0.562500,   1.2075,   1.2076,   1.0467,   70.11  ],
[ 0.587500,   1.2327,   1.1944,   1.0473,    0.00  ],
[ 0.612500,   1.2363,   1.1210,   1.0475,   18.52  ],
[ 0.637500,   1.2305,   0.9843,   1.0480,    8.53  ],
[ 0.662500,   1.2510,   0.9832,   1.0484,    7.80  ],
[ 0.687500,   1.2878,   1.0163,   1.0471,   13.25  ],
[ 0.712500,   1.2629,   0.8242,   1.0493,    0.00  ],
[ 0.737500,   1.2569,   0.7574,   1.0512,    0.00  ],
[ 0.762500,   1.2782,   0.9138,   1.0503,   13.39  ],
[ 0.787500,   1.3824,   1.1470,   1.0617,    0.00  ],
[ 0.812500,   1.4473,   2.0213,   1.0580,   42.80  ],
[ 0.837500,   1.3038,   1.3784,   1.0620,    0.00  ],
[ 0.862500,   1.2911,   1.3270,   1.0615,    0.00  ],
[ 0.887500,   1.3129,   1.2210,   1.0611,    0.00  ],
[ 0.912500,   1.2777,   1.0814,   1.0626,    0.00  ],
[ 0.937500,   1.2811,   1.0702,   1.0619,    0.00  ],
[ 0.962500,   1.2612,   0.9885,   1.0626,    0.00  ],
[ 0.987500,   1.2604,   0.8087,   1.0655,    0.00  ],
[ 1.012500,   1.2983,   0.8930,   1.0636,    0.00  ],
[ 1.037500,   1.2850,   0.7972,   1.0659,    0.00  ],
[ 1.062500,   1.3138,   0.7260,   1.0663,    0.00  ],
[ 1.087500,   1.3197,   0.8153,   1.0661,    0.00  ],
[ 1.112500,   1.3175,   0.8557,   1.0645,    0.00  ],
[ 1.137500,   1.3204,   0.7447,   1.0680,    0.00  ],
[ 1.162500,   1.3159,   0.6524,   1.0696,    0.00  ],
[ 1.187500,   1.3158,   1.0359,   1.0655,    0.00  ],
[ 1.212500,   1.2604,   0.6751,   1.0688,    0.00  ],
[ 1.237500,   1.2405,   0.3713,   1.0767,    0.00  ],
[ 1.262500,   1.2311,   0.0000,   1.0838,    0.00  ],
[ 1.287500,   1.2726,   0.6314,   1.0749,    0.00  ],
[ 1.312500,   1.2633,   0.1106,   1.0824,    0.00  ],
[ 1.337500,   1.2617,   0.0062,   1.0869,    0.00  ],
[ 1.362500,   1.3370,   0.0000,   1.0863,    0.00  ],
[ 1.387500,   1.3030,   0.0000,   1.1675,    3.40  ],
[ 1.412500,   1.1261,   0.0000,   1.1924,    0.00  ],
[ 1.437500,   0.9909,   0.3812,   1.3719,  249.50  ],
[ 1.462500,   1.1236,   0.0000,   1.7456, 3018.05  ],
[ 1.487500,   1.0000,   1.0000,   1.0000,    0.00  ],
[ 1.512500,   2.0474,   0.0000,   1.1349, 6891.45  ],
[ 1.537500,   2.4966,   2.1395,   0.9680,  499.05  ],
[ 1.562500,   1.3810,   0.7912,   1.0454,    0.00  ],
[ 1.587500,   1.3942,   0.9695,   1.0477,    0.00  ],
[ 1.612500,   1.3345,   1.0945,   1.0472,    0.00  ],
[ 1.637500,   1.3676,   1.2421,   1.0452,    0.00  ],
[ 1.662500,   1.6021,   1.3375,   1.0403,    0.00  ],
[ 1.687500,   2.0264,   2.1715,   1.0147,   14.47  ],
[ 1.712500,   1.2049,   1.1964,   1.0530,    0.00  ],
[ 1.737500,   1.1249,   1.2936,   1.0536,    0.00  ],
[ 1.762500,   1.1347,   1.2246,   1.0532,    0.00  ],
[ 1.787500,   1.1701,   1.2867,   1.0533,    0.00  ],
[ 1.812500,   0.0000,   0.8555,   1.0727,    0.00  ],
[ 1.837500,   0.0000,   0.8239,   1.0782,   17.40  ],
[ 1.862500,   0.0000,   1.0117,   1.0707,   46.28  ],
[ 1.887500,   0.0000,   1.0167,   1.0732,   12.92  ],
[ 1.912500,   0.0000,   1.0164,   1.0736,   67.46  ],
[ 1.937500,   0.0000,   0.9943,   1.0741,   53.60  ],
[ 1.962500,   0.0000,   1.0622,   1.0737,   72.03  ],
[ 1.987500,   0.0000,   1.0107,   1.0767,    6.60  ],
[ 2.012500,   0.0000,   0.9800,   1.0817,  111.64  ],
[ 2.037500,   0.0000,   0.9055,   1.0865,   87.25  ],
[ 2.062500,   0.0000,   0.9485,   1.0842,  108.32  ],
[ 2.087500,   0.0000,   0.9057,   1.0838,   55.43  ],
[ 2.112500,   0.0000,   0.7463,   1.0873,   90.90  ],
[ 2.137500,   0.0000,   0.7154,   1.0923,   39.86  ],
[ 2.162500,   0.0000,   0.6399,   1.0955,   64.02  ],
[ 2.187500,   0.0000,   0.7927,   1.0948,  122.85  ],
[ 2.212500,   0.0000,   0.8830,   1.0953,   84.57  ],
[ 2.237500,   0.0000,   0.8920,   1.0963,   86.02  ],
[ 2.262500,   0.0000,   0.8954,   1.0962,   70.49  ],
[ 2.287500,   0.0000,   1.0388,   1.0957,   23.22  ],
[ 2.312500,   0.0000,   0.9922,   1.1010,   73.24  ],
[ 2.337500,   0.0000,   0.9590,   1.1072,   44.06  ],
[ 2.362500,   0.0000,   0.8677,   1.1084,   55.37  ],
[ 2.387500,   0.0000,   0.9179,   1.1066,   54.89  ],
[ 2.412500,   0.0000,   0.9489,   1.1064,   70.17  ],
[ 2.437500,   0.0000,   0.8526,   1.1068,   68.58  ],
[ 2.462500,   0.0000,   0.7528,   1.1091,   58.29  ],
[ 2.487500,   0.0000,   0.6683,   1.1425,  323.50  ]
]



#######################################################################


class CaloSwLongWeights_v5_parms:
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
    correction = {'ele55' : CaloSwLongWeights_v5_ele55,
                  'ele35' : CaloSwLongWeights_v5_ele35,
                  'ele37' : CaloSwLongWeights_v5_ele37,
                  'gam55' : CaloSwLongWeights_v5_gam55,
                  'gam35' : CaloSwLongWeights_v5_gam35,
                  'gam37' : CaloSwLongWeights_v5_gam37,

                  # Use 5x5 for cluster sizes that aren't explicitly derived.
                  'ele33' : CaloSwLongWeights_v5_ele55,
                  'ele57' : CaloSwLongWeights_v5_ele55,
                  'ele77' : CaloSwLongWeights_v5_ele55,
                  'gam33' : CaloSwLongWeights_v5_gam55,
                  'gam57' : CaloSwLongWeights_v5_gam55,
                  'gam77' : CaloSwLongWeights_v5_gam55,
                  }
    
