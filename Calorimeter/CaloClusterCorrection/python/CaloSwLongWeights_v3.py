# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration

#
# $Id: CaloSwLongWeights_v3.py,v 1.2 2006-11-10 03:47:27 ssnyder Exp $
#
# File: CaloClusterCorrection/python/CaloSwLongWeights_v3.py
# Created: Nov 2006, sss
# Purpose: Longitudinal weights corrections, v3 and v3_1.
#
# These corrections were derived by Stathes Paganis.
# The electron weights were derived using 11.0.41 simulation
# and 12.0.0 reconstruction; the photons were derived using
# 11.0.42 simulation and 12.0.2 reconstruction.
#
# There are two versions exported by this file.
# `v3' has electrons only (the electron weights are reused for photon
# clusters).  This was added in CaloClusterCorrection-00-02-38,
# in 12.0.0.
# `v3_1' has both electrons and photons.  This was added in
# CaloClusterCorrection-00-02-59, in 12.0.3.
#


from CaloClusterCorrection.common import *

#######################################################################
# 5x5 cluster size electrons.
CaloSwLongWeights_v3_ele55 = [
               # w0        w3       escale    eoffset
[ 0.012500,   1.1087,   1.4394,   1.0348,    438.28  ],
[ 0.037500,   1.1628,   1.7311,   1.0293,    569.12  ],
[ 0.062500,   1.0486,   1.4997,   1.0336,    312.08  ],
[ 0.087500,   1.1137,   1.7598,   1.0301,    299.17  ],
[ 0.112500,   1.1899,   1.8084,   1.0273,    323.22  ],
[ 0.137500,   1.1116,   1.5163,   1.0310,    286.38  ],
[ 0.162500,   1.1408,   1.7358,   1.0307,    225.10  ],
[ 0.187500,   1.0602,   1.5640,   1.0290,    430.02  ],
[ 0.212500,   1.1092,   1.4456,   1.0291,    383.11  ],
[ 0.237500,   1.0937,   1.6006,   1.0295,    246.89  ],
[ 0.262500,   1.0983,   1.3089,   1.0309,    225.16  ],
[ 0.287500,   1.0859,   1.4976,   1.0280,    337.51  ],
[ 0.312500,   1.0812,   1.2598,   1.0303,    337.95  ],
[ 0.337500,   1.1025,   1.1983,   1.0303,    218.11  ],
[ 0.362500,   1.0792,   1.2802,   1.0310,    269.20  ],
[ 0.387500,   1.0521,   1.0107,   1.0312,    338.67  ],
[ 0.412500,   1.0201,   0.9478,   1.0308,    454.73  ],
[ 0.437500,   1.0933,   1.1465,   1.0308,    204.20  ],
[ 0.462500,   1.1180,   1.0511,   1.0297,    282.18  ],
[ 0.487500,   1.1125,   1.0930,   1.0268,    436.26  ],
[ 0.512500,   1.1264,   0.9691,   1.0265,    397.12  ],
[ 0.537500,   1.0884,   1.0117,   1.0291,    354.48  ],
[ 0.562500,   1.1487,   1.0134,   1.0270,    328.52  ],
[ 0.587500,   1.1124,   0.9723,   1.0269,    489.05  ],
[ 0.612500,   1.1453,   1.0802,   1.0259,    583.07  ],
[ 0.637500,   1.1356,   0.7776,   1.0285,    529.18  ],
[ 0.662500,   1.1856,   1.0408,   1.0288,    280.36  ],
[ 0.687500,   1.1632,   1.0302,   1.0276,    448.45  ],
[ 0.712500,   1.1810,   1.0796,   1.0271,    635.77  ],
[ 0.737500,   1.1577,   0.7884,   1.0310,    632.10  ],
[ 0.762500,   1.1370,   1.0622,   1.0330,    501.73  ],
[ 0.787500,   1.1387,   1.1553,   1.0378,    860.87  ],
[ 0.812500,   1.1977,   1.1208,   1.0368,   1113.08  ],
[ 0.837500,   1.1539,   1.0814,   1.0421,    769.26  ],
[ 0.862500,   1.1686,   1.0709,   1.0450,    678.60  ],
[ 0.887500,   1.1637,   0.9776,   1.0431,    895.34  ],
[ 0.912500,   1.1699,   0.9769,   1.0422,   1004.46  ],
[ 0.937500,   1.1368,   0.5496,   1.0468,    850.96  ],
[ 0.962500,   1.1527,   0.9254,   1.0430,    989.12  ],
[ 0.987500,   1.1279,   0.9543,   1.0408,   1394.30  ],
[ 1.012500,   1.1399,   0.8376,   1.0443,   1275.74  ],
[ 1.037500,   1.1130,   1.0119,   1.0409,   1455.50  ],
[ 1.062500,   1.1187,   0.8448,   1.0512,    890.40  ],
[ 1.087500,   1.1184,   0.6387,   1.0504,   1192.13  ],
[ 1.112500,   1.1068,   0.7431,   1.0472,   1628.83  ],
[ 1.137500,   1.0839,   0.5367,   1.0508,   1704.68  ],
[ 1.162500,   1.1094,   0.5407,   1.0551,   1297.52  ],
[ 1.187500,   1.0830,   0.2000,   1.0552,   2100.08  ],
[ 1.212500,   1.0781,   0.5847,   1.0637,   1138.79  ],
[ 1.237500,   1.0806,   0.6480,   1.0675,    763.43  ],
[ 1.262500,   1.0541,   0.5805,   1.0631,   1576.79  ],
[ 1.287500,   1.0779,   0.6831,   1.0539,   2059.66  ],
[ 1.312500,   1.0519,   0.0622,   1.0838,   1527.45  ],
[ 1.337500,   1.0171,   0.2389,   1.0781,   1738.25  ],
[ 1.362500,   1.0390,   0.3211,   1.0985,   1552.42  ],
[ 1.387500,   1.0928,   0.0000,   1.1124,   2254.65  ],
[ 1.412500,   1.0384,   0.0000,   1.1850,   1699.25  ],
[ 1.437500,   0.9463,   0.0511,   1.3116,   1576.02  ],
[ 1.462500,   0.9577,   2.1322,   1.7316,    645.82  ],
[ 1.487500,   1.1005,   1.0682,   1.6952,   2180.34  ],
[ 1.512500,   1.3766,   0.0000,   1.1370,   5186.30  ],
[ 1.537500,   1.2736,   1.1088,   0.9693,   2739.47  ],
[ 1.562500,   1.1518,   0.8964,   0.9388,   1949.70  ],
[ 1.587500,   1.1014,   0.7691,   0.9416,   1816.98  ],
[ 1.612500,   1.0826,   0.7639,   1.0476,   1217.71  ],
[ 1.637500,   1.1698,   0.9779,   1.0298,   1646.34  ],
[ 1.662500,   1.2216,   1.1649,   1.0264,   2138.17  ],
[ 1.687500,   1.1696,   0.9328,   1.0417,    582.29  ],
[ 1.712500,   1.1329,   1.1481,   1.0268,   1199.53  ],
[ 1.737500,   0.9474,   0.9052,   1.0373,    628.78  ],
[ 1.762500,   0.9184,   1.1401,   1.0268,   1154.51  ],
[ 1.787500,   0.8678,   1.0166,   1.0291,   1026.23  ],
[ 1.812500,   0.5815,   1.0010,   1.0437,    893.80  ],
[ 1.837500,   0.0000,   0.8413,   1.0549,    820.39  ],
[ 1.862500,   0.0000,   0.9787,   1.0395,   1686.15  ],
[ 1.887500,   0.0000,   1.0671,   1.0400,   1337.49  ],
[ 1.912500,   0.0000,   1.0465,   1.0371,   1603.04  ],
[ 1.937500,   0.0000,   0.9555,   1.0443,   1090.36  ],
[ 1.962500,   0.0000,   0.9394,   1.0482,    865.37  ],
[ 1.987500,   0.0000,   0.9893,   1.0453,   1171.31  ],
[ 2.012500,   0.0000,   0.8645,   1.0467,   1279.64  ],
[ 2.037500,   0.0000,   0.8551,   1.0465,   1388.42  ],
[ 2.062500,   0.0000,   0.9752,   1.0501,    946.58  ],
[ 2.087500,   0.0000,   0.9828,   1.0479,   1175.13  ],
[ 2.112500,   0.0000,   0.8812,   1.0538,    795.14  ],
[ 2.137500,   0.0000,   0.8945,   1.0508,   1044.38  ],
[ 2.162500,   0.0000,   0.8921,   1.0474,   1220.05  ],
[ 2.187500,   0.0000,   0.9384,   1.0517,   1103.02  ],
[ 2.212500,   0.0000,   0.9697,   1.0506,   1133.90  ],
[ 2.237500,   0.0000,   0.9897,   1.0502,   1095.14  ],
[ 2.262500,   0.0000,   1.1112,   1.0513,   1030.12  ],
[ 2.287500,   0.0000,   0.9331,   1.0534,   1023.85  ],
[ 2.312500,   0.0000,   0.9844,   1.0548,   1179.54  ],
[ 2.337500,   0.0000,   0.7608,   1.0652,    840.74  ],
[ 2.362500,   0.0000,   0.8112,   1.0641,    761.19  ],
[ 2.387500,   0.0000,   0.6904,   1.0653,    893.81  ],
[ 2.412500,   0.0000,   0.8266,   1.0608,   1153.13  ],
[ 2.437500,   0.0000,   0.7724,   1.0657,   1171.53  ],
[ 2.462500,   0.0000,   0.9821,   1.0749,   1024.67  ],
[ 2.487500,   0.0000,   0.7133,   1.0882,   1776.87  ]
] 



#######################################################################
# 3x5 cluster size electrons.
CaloSwLongWeights_v3_ele35 = [
                # w0        w3       escale    eoffset
[ 0.012500,   1.1848,   1.5063,   1.0519,    408.42  ],
[ 0.037500,   1.2170,   1.8462,   1.0458,    514.07  ],
[ 0.062500,   1.0680,   1.5460,   1.0496,    414.59  ],
[ 0.087500,   1.1275,   1.8190,   1.0486,    194.68  ],
[ 0.112500,   1.2420,   1.8119,   1.0432,    420.76  ],
[ 0.137500,   1.1397,   1.5874,   1.0480,    260.73  ],
[ 0.162500,   1.1740,   1.7560,   1.0475,    264.72  ],
[ 0.187500,   1.0852,   1.5555,   1.0450,    504.79  ],
[ 0.212500,   1.1437,   1.5178,   1.0458,    383.10  ],
[ 0.237500,   1.1297,   1.6750,   1.0460,    260.88  ],
[ 0.262500,   1.1285,   1.2323,   1.0469,    395.97  ],
[ 0.287500,   1.1335,   1.6444,   1.0444,    357.80  ],
[ 0.312500,   1.1039,   1.3134,   1.0478,    318.47  ],
[ 0.337500,   1.1258,   1.5403,   1.0471,    142.26  ],
[ 0.362500,   1.1282,   1.2033,   1.0477,    285.42  ],
[ 0.387500,   1.0725,   1.0603,   1.0470,    424.73  ],
[ 0.412500,   1.0526,   0.9990,   1.0489,    399.41  ],
[ 0.437500,   1.1298,   1.2062,   1.0483,    205.07  ],
[ 0.462500,   1.1461,   1.1301,   1.0457,    367.46  ],
[ 0.487500,   1.1218,   0.9882,   1.0460,    352.06  ],
[ 0.512500,   1.1807,   1.0533,   1.0434,    431.41  ],
[ 0.537500,   1.1063,   1.1177,   1.0438,    592.56  ],
[ 0.562500,   1.1428,   1.0280,   1.0425,    541.39  ],
[ 0.587500,   1.1295,   0.9865,   1.0452,    454.59  ],
[ 0.612500,   1.1417,   0.9568,   1.0434,    624.84  ],
[ 0.637500,   1.1517,   0.8204,   1.0457,    594.59  ],
[ 0.662500,   1.1765,   0.8894,   1.0463,    405.14  ],
[ 0.687500,   1.1589,   0.9909,   1.0461,    506.73  ],
[ 0.712500,   1.1916,   1.0466,   1.0448,    619.34  ],
[ 0.737500,   1.1848,   0.6712,   1.0498,    595.24  ],
[ 0.762500,   1.1438,   1.0814,   1.0476,   1016.26  ],
[ 0.787500,   1.1412,   1.2296,   1.0617,    540.46  ],
[ 0.812500,   1.1909,   1.1473,   1.0549,   1462.53  ],
[ 0.837500,   1.1797,   1.2503,   1.0661,    500.23  ],
[ 0.862500,   1.1794,   0.9916,   1.0671,    739.56  ],
[ 0.887500,   1.1690,   1.0053,   1.0645,    832.52  ],
[ 0.912500,   1.1616,   1.0932,   1.0626,    913.82  ],
[ 0.937500,   1.1599,   0.5702,   1.0675,    971.65  ],
[ 0.962500,   1.1730,   0.8626,   1.0637,   1097.07  ],
[ 0.987500,   1.1328,   0.8980,   1.0646,   1180.28  ],
[ 1.012500,   1.1154,   0.6181,   1.0742,    731.23  ],
[ 1.037500,   1.1277,   0.9371,   1.0712,   1185.31  ],
[ 1.062500,   1.1362,   0.8273,   1.0704,   1307.14  ],
[ 1.087500,   1.1294,   0.6805,   1.0766,    814.29  ],
[ 1.112500,   1.1112,   0.6254,   1.0732,   1139.40  ],
[ 1.137500,   1.0907,   0.4058,   1.0739,   1642.35  ],
[ 1.162500,   1.1229,   0.4496,   1.0765,   1335.58  ],
[ 1.187500,   1.0969,   0.2561,   1.0794,   1706.09  ],
[ 1.212500,   1.0805,   0.5855,   1.0884,    788.80  ],
[ 1.237500,   1.0907,   0.6977,   1.0859,   1082.58  ],
[ 1.262500,   1.0557,   0.4270,   1.0905,   1476.11  ],
[ 1.287500,   1.0675,   0.3902,   1.0855,   1442.34  ],
[ 1.312500,   1.0686,   0.3362,   1.1141,   1264.18  ],
[ 1.337500,   1.0312,   0.0000,   1.1091,   2029.16  ],
[ 1.362500,   1.0634,   0.5998,   1.1165,   2152.70  ],
[ 1.387500,   1.0865,   0.0000,   1.1602,   2561.16  ],
[ 1.412500,   1.0456,   0.0000,   1.1985,   1953.92  ],
[ 1.437500,   0.9354,   0.0000,   1.3276,   1900.39  ],
[ 1.462500,   0.9251,   7.2773,   1.5419,   2322.21  ],
[ 1.487500,   1.2450,   0.9958,   1.7161,   2251.01  ],
[ 1.512500,   1.6097,   0.0000,   1.0517,   6605.28  ],
[ 1.537500,   1.6660,   1.1509,   0.9688,   2714.72  ],
[ 1.562500,   1.1976,   0.9228,   0.9606,   1315.60  ],
[ 1.587500,   1.1578,   0.9245,   0.9545,   2032.47  ],
[ 1.612500,   1.0677,   0.7646,   1.0534,   2258.41  ],
[ 1.637500,   1.2003,   0.9117,   1.0515,   1403.33  ],
[ 1.662500,   1.2385,   1.0554,   1.0527,   1632.95  ],
[ 1.687500,   1.1866,   1.1964,   1.0479,   2464.76  ],
[ 1.712500,   1.1156,   0.9799,   1.0570,    469.03  ],
[ 1.737500,   0.9744,   0.8564,   1.0590,    645.70  ],
[ 1.762500,   0.9343,   1.1836,   1.0553,    544.71  ],
[ 1.787500,   0.8757,   0.9472,   1.0551,    836.24  ],
[ 1.812500,   0.5715,   0.9879,   1.0717,    646.91  ],
[ 1.837500,   0.0000,   0.9124,   1.0666,   1589.19  ],
[ 1.862500,   0.0000,   0.9266,   1.0662,   1607.24  ],
[ 1.887500,   0.0000,   1.1198,   1.0650,   1323.03  ],
[ 1.912500,   0.0000,   0.9838,   1.0650,   1424.29  ],
[ 1.937500,   0.0000,   0.9409,   1.0692,   1165.36  ],
[ 1.962500,   0.0000,   0.8421,   1.0744,   1096.42  ],
[ 1.987500,   0.0000,   0.9942,   1.0721,   1160.02  ],
[ 2.012500,   0.0000,   0.8700,   1.0728,   1391.64  ],
[ 2.037500,   0.0000,   0.7495,   1.0773,   1230.34  ],
[ 2.062500,   0.0000,   0.9696,   1.0794,    853.50  ],
[ 2.087500,   0.0000,   0.9761,   1.0783,   1117.55  ],
[ 2.112500,   0.0000,   0.8822,   1.0777,   1009.11  ],
[ 2.137500,   0.0000,   0.8992,   1.0805,    968.45  ],
[ 2.162500,   0.0000,   0.8339,   1.0810,   1121.40  ],
[ 2.187500,   0.0000,   1.0059,   1.0839,    932.20  ],
[ 2.212500,   0.0000,   0.9789,   1.0761,   1504.28  ],
[ 2.237500,   0.0000,   0.8513,   1.0841,   1047.02  ],
[ 2.262500,   0.0000,   1.0396,   1.0863,    883.01  ],
[ 2.287500,   0.0000,   0.9051,   1.0894,    836.69  ],
[ 2.312500,   0.0000,   0.9292,   1.0911,   1077.28  ],
[ 2.337500,   0.0000,   0.8360,   1.0941,   1194.95  ],
[ 2.362500,   0.0000,   0.7917,   1.0968,    974.36  ],
[ 2.387500,   0.0000,   0.6834,   1.1039,    593.64  ],
[ 2.412500,   0.0000,   0.7681,   1.0978,   1214.30  ],
[ 2.437500,   0.0000,   0.8970,   1.1031,   1048.70  ],
[ 2.462500,   0.0000,   0.8694,   1.1128,    826.24  ],
[ 2.487500,   0.0000,   0.6584,   1.1183,   1202.27  ]
] 


#######################################################################
# 3x7 cluster size electrons.
CaloSwLongWeights_v3_ele37 = [
                # w0        w3       escale    eoffset
[ 0.012500,   1.1090,   1.4792,   1.0506,    285.86  ],
[ 0.037500,   1.1687,   1.6084,   1.0429,    548.75  ],
[ 0.062500,   1.0553,   1.5097,   1.0472,    299.69  ],
[ 0.087500,   1.1239,   1.8141,   1.0423,    322.40  ],
[ 0.112500,   1.2004,   1.7386,   1.0397,    413.90  ],
[ 0.137500,   1.1206,   1.5018,   1.0442,    233.44  ],
[ 0.162500,   1.1137,   1.6722,   1.0440,    235.29  ],
[ 0.187500,   1.0360,   1.5743,   1.0422,    426.75  ],
[ 0.212500,   1.0927,   1.4979,   1.0431,    334.21  ],
[ 0.237500,   1.1029,   1.7038,   1.0427,    196.92  ],
[ 0.262500,   1.0820,   1.2112,   1.0438,    301.53  ],
[ 0.287500,   1.1055,   1.5975,   1.0421,    292.28  ],
[ 0.312500,   1.0526,   1.3368,   1.0446,    214.57  ],
[ 0.337500,   1.1303,   1.4256,   1.0424,    281.62  ],
[ 0.362500,   1.0537,   1.1747,   1.0467,    179.05  ],
[ 0.387500,   1.0614,   1.2950,   1.0434,    326.42  ],
[ 0.412500,   1.0106,   0.9349,   1.0444,    410.44  ],
[ 0.437500,   1.0971,   1.1597,   1.0417,    340.89  ],
[ 0.462500,   1.1071,   0.9741,   1.0429,    337.30  ],
[ 0.487500,   1.0891,   1.0149,   1.0422,    339.33  ],
[ 0.512500,   1.1437,   1.0195,   1.0412,    286.40  ],
[ 0.537500,   1.0753,   1.0728,   1.0412,    461.81  ],
[ 0.562500,   1.1404,   1.0520,   1.0390,    474.87  ],
[ 0.587500,   1.1175,   1.0131,   1.0414,    313.26  ],
[ 0.612500,   1.1492,   0.9430,   1.0411,    438.99  ],
[ 0.637500,   1.1242,   0.7135,   1.0438,    366.33  ],
[ 0.662500,   1.1578,   0.8875,   1.0443,    248.44  ],
[ 0.687500,   1.1616,   0.9870,   1.0414,    465.70  ],
[ 0.712500,   1.1502,   0.8994,   1.0401,    643.69  ],
[ 0.737500,   1.1527,   0.6771,   1.0433,    779.22  ],
[ 0.762500,   1.1410,   1.0643,   1.0444,    795.55  ],
[ 0.787500,   1.1089,   1.1035,   1.0541,    777.75  ],
[ 0.812500,   1.1961,   1.1362,   1.0538,   1014.60  ],
[ 0.837500,   1.1721,   1.3049,   1.0558,    766.60  ],
[ 0.862500,   1.1504,   0.8711,   1.0571,    957.68  ],
[ 0.887500,   1.1591,   1.0365,   1.0575,    841.71  ],
[ 0.912500,   1.1555,   1.0357,   1.0574,    897.30  ],
[ 0.937500,   1.1474,   0.5860,   1.0610,    932.18  ],
[ 0.962500,   1.1622,   0.7804,   1.0644,    640.52  ],
[ 0.987500,   1.1370,   0.8225,   1.0592,   1012.62  ],
[ 1.012500,   1.1092,   0.5653,   1.0674,    724.37  ],
[ 1.037500,   1.1237,   0.9220,   1.0638,    876.92  ],
[ 1.062500,   1.1178,   0.7389,   1.0690,    716.11  ],
[ 1.087500,   1.1261,   0.5657,   1.0653,   1198.13  ],
[ 1.112500,   1.1011,   0.6321,   1.0664,   1186.29  ],
[ 1.137500,   1.0791,   0.5562,   1.0706,   1027.87  ],
[ 1.162500,   1.1148,   0.3668,   1.0752,    916.40  ],
[ 1.187500,   1.0829,   0.1989,   1.0784,   1059.64  ],
[ 1.212500,   1.0755,   0.6655,   1.0793,    698.29  ],
[ 1.237500,   1.0732,   0.4654,   1.0790,   1228.19  ],
[ 1.262500,   1.0682,   0.3837,   1.0780,   1613.14  ],
[ 1.287500,   1.0568,   0.5968,   1.0746,   1500.14  ],
[ 1.312500,   1.0321,   0.2134,   1.1031,   1002.42  ],
[ 1.337500,   1.0285,   0.0000,   1.1169,    873.83  ],
[ 1.362500,   1.0615,   0.6869,   1.1140,   1271.76  ],
[ 1.387500,   1.0720,   0.0000,   1.1559,   1824.99  ],
[ 1.412500,   1.0198,   4.4531,   1.1911,   1526.96  ],
[ 1.437500,   0.9107,   2.8316,   1.3170,   1525.35  ],
[ 1.462500,   0.9538,   0.0000,   1.5629,   1835.07  ],
[ 1.487500,   1.2071,   0.3759,   1.6634,   1969.46  ],
[ 1.512500,   1.3950,   0.0000,   1.1142,   5641.55  ],
[ 1.537500,   1.3007,   0.9268,   0.9993,   1159.79  ],
[ 1.562500,   1.1635,   0.9845,   0.9521,   1281.32  ],
[ 1.587500,   1.1175,   0.6901,   0.9521,   1512.79  ],
[ 1.612500,   1.0434,   1.0425,   1.0482,   1950.68  ],
[ 1.637500,   1.2087,   0.9221,   1.0494,    844.34  ],
[ 1.662500,   1.2503,   1.0595,   1.0458,   1381.54  ],
[ 1.687500,   1.1831,   1.0579,   1.0548,   1200.58  ],
[ 1.712500,   1.1048,   0.9879,   1.0496,    466.85  ],
[ 1.737500,   0.9515,   0.8083,   1.0514,    437.15  ],
[ 1.762500,   0.8917,   1.0884,   1.0515,    432.15  ],
[ 1.787500,   0.8498,   0.8942,   1.0527,    524.90  ],
[ 1.812500,   0.5684,   0.7764,   1.0461,   2177.41  ],
[ 1.837500,   0.0000,   0.8913,   1.0612,   1498.82  ],
[ 1.862500,   0.0000,   0.8869,   1.0640,   1197.62  ],
[ 1.887500,   0.0000,   1.0172,   1.0581,   1401.70  ],
[ 1.912500,   0.0000,   0.9239,   1.0601,   1267.57  ],
[ 1.937500,   0.0000,   0.8474,   1.0633,   1106.15  ],
[ 1.962500,   0.0000,   0.8869,   1.0654,   1060.37  ],
[ 1.987500,   0.0000,   0.9183,   1.0672,    903.41  ],
[ 2.012500,   0.0000,   0.7920,   1.0673,   1254.36  ],
[ 2.037500,   0.0000,   0.6475,   1.0700,   1292.34  ],
[ 2.062500,   0.0000,   0.9030,   1.0735,    768.58  ],
[ 2.087500,   0.0000,   0.9405,   1.0695,   1024.00  ],
[ 2.112500,   0.0000,   0.7275,   1.0763,    877.14  ],
[ 2.137500,   0.0000,   0.8283,   1.0736,    963.74  ],
[ 2.162500,   0.0000,   0.8598,   1.0706,   1204.03  ],
[ 2.187500,   0.0000,   0.9662,   1.0719,   1106.23  ],
[ 2.212500,   0.0000,   0.9476,   1.0723,   1128.05  ],
[ 2.237500,   0.0000,   0.8364,   1.0781,    756.40  ],
[ 2.262500,   0.0000,   0.9468,   1.0749,   1008.45  ],
[ 2.287500,   0.0000,   0.9312,   1.0784,    902.82  ],
[ 2.312500,   0.0000,   0.9785,   1.0801,    982.63  ],
[ 2.337500,   0.0000,   0.8044,   1.0876,    894.29  ],
[ 2.362500,   0.0000,   0.8222,   1.0860,    857.14  ],
[ 2.387500,   0.0000,   0.6816,   1.0921,    763.41  ],
[ 2.412500,   0.0000,   0.9225,   1.0859,    999.92  ],
[ 2.437500,   0.0000,   0.8284,   1.0899,   1136.84  ],
[ 2.462500,   0.0000,   0.8836,   1.0965,    923.35  ],
[ 2.487500,   0.0000,   0.6019,   1.1049,   1461.76  ]
] 


#######################################################################
# 5x5 cluster size photons.
CaloSwLongWeights_v3_gam55 = [
               # w0        w3       escale    eoffset 
[ 0.012500,   1.2320,   1.9277,   1.0310,    228.35  ],
[ 0.037500,   1.1892,   1.8346,   1.0314,     97.29  ],
[ 0.062500,   1.2456,   1.8723,   1.0286,    103.99  ],
[ 0.087500,   1.1545,   1.8897,   1.0295,     17.22  ],
[ 0.112500,   1.1801,   1.8602,   1.0292,     45.47  ],
[ 0.137500,   1.2507,   1.8928,   1.0283,     24.00  ],
[ 0.162500,   1.2151,   1.9118,   1.0269,     91.76  ],
[ 0.187500,   1.1940,   1.7987,   1.0282,     69.23  ],
[ 0.212500,   1.2100,   1.9840,   1.0260,    101.54  ],
[ 0.237500,   1.1584,   1.8433,   1.0268,    105.20  ],
[ 0.262500,   1.1577,   1.8074,   1.0268,    102.36  ],
[ 0.287500,   1.1787,   1.9542,   1.0274,     46.56  ],
[ 0.312500,   1.1678,   1.7777,   1.0275,     48.33  ],
[ 0.337500,   1.2158,   1.8381,   1.0262,     59.03  ],
[ 0.362500,   1.2155,   1.7207,   1.0259,     77.24  ],
[ 0.387500,   1.1987,   1.7742,   1.0245,    144.14  ],
[ 0.412500,   1.1654,   1.5289,   1.0263,    106.75  ],
[ 0.437500,   1.1850,   1.4924,   1.0246,    185.31  ],
[ 0.462500,   1.1640,   1.3619,   1.0269,     86.27  ],
[ 0.487500,   1.1713,   1.3380,   1.0246,    199.93  ],
[ 0.512500,   1.1682,   1.3077,   1.0246,    139.57  ],
[ 0.537500,   1.1925,   1.2715,   1.0248,    102.41  ],
[ 0.562500,   1.2106,   1.2073,   1.0238,    174.60  ],
[ 0.587500,   1.2381,   1.2154,   1.0243,     92.79  ],
[ 0.612500,   1.1882,   1.1280,   1.0253,     64.99  ],
[ 0.637500,   1.2229,   1.1188,   1.0244,     85.90  ],
[ 0.662500,   1.2439,   1.0424,   1.0240,    177.24  ],
[ 0.687500,   1.2289,   1.0791,   1.0245,    158.20  ],
[ 0.712500,   1.2545,   1.0254,   1.0246,    153.80  ],
[ 0.737500,   1.2356,   0.9995,   1.0253,    166.04  ],
[ 0.762500,   1.2306,   1.0731,   1.0235,    364.64  ],
[ 0.787500,   1.2155,   1.0998,   1.0305,    401.44  ],
[ 0.812500,   1.2150,   1.1269,   1.0382,    174.33  ],
[ 0.837500,   1.2706,   1.3962,   1.0344,    156.12  ],
[ 0.862500,   1.2774,   1.3859,   1.0335,    200.99  ],
[ 0.887500,   1.2586,   1.2318,   1.0351,    172.04  ],
[ 0.912500,   1.2392,   1.1421,   1.0357,    156.72  ],
[ 0.937500,   1.2529,   1.0809,   1.0359,    219.08  ],
[ 0.962500,   1.2369,   1.0386,   1.0384,     87.19  ],
[ 0.987500,   1.2264,   0.9833,   1.0387,     89.79  ],
[ 1.012500,   1.2271,   0.9905,   1.0377,    218.97  ],
[ 1.037500,   1.2152,   0.8746,   1.0387,    215.28  ],
[ 1.062500,   1.1925,   0.8669,   1.0395,    204.76  ],
[ 1.087500,   1.2187,   0.9063,   1.0388,    360.30  ],
[ 1.112500,   1.1946,   0.9010,   1.0406,    253.95  ],
[ 1.137500,   1.1797,   0.7121,   1.0443,    243.85  ],
[ 1.162500,   1.1700,   0.6631,   1.0466,    250.54  ],
[ 1.187500,   1.1808,   0.8557,   1.0422,    374.21  ],
[ 1.212500,   1.1563,   0.8603,   1.0431,    364.90  ],
[ 1.237500,   1.1382,   0.7282,   1.0494,    205.60  ],
[ 1.262500,   1.1512,   0.7386,   1.0453,    554.32  ],
[ 1.287500,   1.1445,   0.8152,   1.0522,    101.58  ],
[ 1.312500,   1.1317,   0.4315,   1.0573,    441.83  ],
[ 1.337500,   1.1179,   0.2830,   1.0673,    174.06  ],
[ 1.362500,   1.0962,   0.0000,   1.0756,    609.80  ],
[ 1.387500,   1.1699,   0.0000,   1.0968,    614.77  ],
[ 1.412500,   1.1390,   0.0000,   1.1317,   1116.73  ],
[ 1.437500,   1.0969,   0.0000,   1.2353,    535.74  ],
[ 1.462500,   1.0558,   0.0000,   1.7886,      0.00  ],
[ 1.487500,   1.0997,   0.0000,   1.7768,    874.06  ],
[ 1.512500,   1.5748,   0.8580,   0.9889,   4702.83  ],
[ 1.537500,   1.2267,   0.9432,   1.0368,    221.62  ],
[ 1.562500,   1.1449,   0.9437,   1.0392,      0.00  ],
[ 1.587500,   1.0602,   1.0411,   1.0385,    116.73  ],
[ 1.612500,   1.0757,   1.0346,   1.0380,      0.00  ],
[ 1.637500,   1.1733,   1.1059,   1.0332,     66.85  ],
[ 1.662500,   1.2189,   1.0528,   1.0317,    119.40  ],
[ 1.687500,   1.1219,   1.2963,   1.0339,    123.09  ],
[ 1.712500,   1.0434,   1.1079,   1.0349,      0.00  ],
[ 1.737500,   0.9118,   1.1520,   1.0342,      0.00  ],
[ 1.762500,   0.7633,   1.0967,   1.0332,     51.59  ],
[ 1.787500,   0.7331,   1.0885,   1.0325,     94.08  ],
[ 1.812500,   0.5372,   1.0779,   1.0353,    131.38  ],
[ 1.837500,   0.0000,   1.0950,   1.0386,    154.79  ],
[ 1.862500,   0.0000,   1.1222,   1.0361,    361.39  ],
[ 1.887500,   0.0000,   1.1855,   1.0375,    203.20  ],
[ 1.912500,   0.0000,   1.1702,   1.0377,    239.68  ],
[ 1.937500,   0.0000,   1.1683,   1.0387,    248.40  ],
[ 1.962500,   0.0000,   1.1595,   1.0409,    156.40  ],
[ 1.987500,   0.0000,   1.1080,   1.0414,    248.86  ],
[ 2.012500,   0.0000,   1.1460,   1.0421,    263.80  ],
[ 2.037500,   0.0000,   1.0857,   1.0425,    260.08  ],
[ 2.062500,   0.0000,   1.1225,   1.0438,    220.40  ],
[ 2.087500,   0.0000,   1.0762,   1.0438,    218.94  ],
[ 2.112500,   0.0000,   1.0444,   1.0450,    215.91  ],
[ 2.137500,   0.0000,   1.0672,   1.0454,    232.68  ],
[ 2.162500,   0.0000,   1.0965,   1.0461,    192.64  ],
[ 2.187500,   0.0000,   1.0476,   1.0480,    172.52  ],
[ 2.212500,   0.0000,   1.0469,   1.0485,    201.31  ],
[ 2.237500,   0.0000,   1.0614,   1.0490,    165.11  ],
[ 2.262500,   0.0000,   1.0659,   1.0505,     71.97  ],
[ 2.287500,   0.0000,   1.0606,   1.0507,    225.80  ],
[ 2.312500,   0.0000,   1.0871,   1.0528,    247.59  ],
[ 2.337500,   0.0000,   1.0781,   1.0547,    245.43  ],
[ 2.362500,   0.0000,   1.0633,   1.0560,    159.95  ],
[ 2.387500,   0.0000,   1.0035,   1.0579,    246.27  ],
[ 2.412500,   0.0000,   0.8924,   1.0589,    350.22  ],
[ 2.437500,   0.0000,   0.9605,   1.0606,    269.22  ],
[ 2.462500,   0.0000,   0.9384,   1.0718,    229.01  ],
[ 2.487500,   0.0000,   0.9212,   1.1000,     11.85  ]
] 


#######################################################################
# 3x5 cluster size photons.
CaloSwLongWeights_v3_gam35 = [
               # w0        w3       escale    eoffset 
[ 0.012500,   1.2810,   1.9386,   1.0485,    182.27  ],
[ 0.037500,   1.2409,   1.8327,   1.0487,    122.47  ],
[ 0.062500,   1.3068,   1.9660,   1.0450,    131.43  ],
[ 0.087500,   1.2393,   1.8602,   1.0450,    125.97  ],
[ 0.112500,   1.2232,   1.8189,   1.0460,     53.42  ],
[ 0.137500,   1.2931,   1.8720,   1.0442,     91.49  ],
[ 0.162500,   1.2618,   1.9226,   1.0443,     67.27  ],
[ 0.187500,   1.2386,   1.7504,   1.0449,     83.15  ],
[ 0.212500,   1.2572,   1.9441,   1.0437,     80.22  ],
[ 0.237500,   1.2359,   1.8837,   1.0438,    101.54  ],
[ 0.262500,   1.1978,   1.7575,   1.0437,    139.06  ],
[ 0.287500,   1.2426,   1.9438,   1.0443,     69.64  ],
[ 0.312500,   1.2298,   1.8009,   1.0446,     53.41  ],
[ 0.337500,   1.2569,   1.7960,   1.0443,     15.21  ],
[ 0.362500,   1.2641,   1.7433,   1.0437,     41.53  ],
[ 0.387500,   1.2446,   1.7239,   1.0427,    113.20  ],
[ 0.412500,   1.2183,   1.6355,   1.0430,    126.49  ],
[ 0.437500,   1.2176,   1.4400,   1.0423,    180.79  ],
[ 0.462500,   1.2088,   1.2790,   1.0438,    149.48  ],
[ 0.487500,   1.2198,   1.3516,   1.0443,     59.23  ],
[ 0.512500,   1.2052,   1.3042,   1.0434,     86.25  ],
[ 0.537500,   1.2053,   1.2513,   1.0439,     27.14  ],
[ 0.562500,   1.2377,   1.1652,   1.0421,    132.79  ],
[ 0.587500,   1.2494,   1.1331,   1.0433,     27.44  ],
[ 0.612500,   1.2174,   1.0778,   1.0428,    107.71  ],
[ 0.637500,   1.2463,   1.0602,   1.0425,    120.14  ],
[ 0.662500,   1.2569,   0.9978,   1.0422,    199.89  ],
[ 0.687500,   1.2475,   1.0413,   1.0425,    156.62  ],
[ 0.712500,   1.2380,   0.9618,   1.0440,     75.42  ],
[ 0.737500,   1.2295,   0.9521,   1.0422,    197.71  ],
[ 0.762500,   1.2422,   1.0501,   1.0416,    370.73  ],
[ 0.787500,   1.2165,   1.0484,   1.0482,    500.31  ],
[ 0.812500,   1.2080,   1.0357,   1.0580,    266.48  ],
[ 0.837500,   1.2782,   1.3377,   1.0560,    179.51  ],
[ 0.862500,   1.2816,   1.3135,   1.0561,    188.57  ],
[ 0.887500,   1.2632,   1.1713,   1.0567,    163.96  ],
[ 0.912500,   1.2562,   1.0926,   1.0575,    138.64  ],
[ 0.937500,   1.2685,   1.0378,   1.0573,    248.38  ],
[ 0.962500,   1.2452,   0.9589,   1.0594,     68.07  ],
[ 0.987500,   1.2285,   0.9095,   1.0613,     70.65  ],
[ 1.012500,   1.2284,   0.9310,   1.0596,    220.28  ],
[ 1.037500,   1.2217,   0.8344,   1.0628,     73.08  ],
[ 1.062500,   1.2171,   0.9261,   1.0604,    198.13  ],
[ 1.087500,   1.2219,   0.8660,   1.0615,    254.65  ],
[ 1.112500,   1.2074,   0.8793,   1.0641,    188.01  ],
[ 1.137500,   1.1843,   0.6720,   1.0673,    224.32  ],
[ 1.162500,   1.1743,   0.6212,   1.0684,    358.16  ],
[ 1.187500,   1.1976,   0.8424,   1.0645,    417.62  ],
[ 1.212500,   1.1605,   0.8389,   1.0641,    443.69  ],
[ 1.237500,   1.1537,   0.7108,   1.0723,    257.05  ],
[ 1.262500,   1.1477,   0.7091,   1.0690,    590.52  ],
[ 1.287500,   1.1569,   0.8127,   1.0717,    339.16  ],
[ 1.312500,   1.1338,   0.4488,   1.0796,    437.01  ],
[ 1.337500,   1.1255,   0.5058,   1.0850,    153.02  ],
[ 1.362500,   1.1108,   0.1751,   1.0896,    532.50  ],
[ 1.387500,   1.1568,   0.0000,   1.1333,    718.11  ],
[ 1.412500,   1.0906,   0.0000,   1.1759,   1185.25  ],
[ 1.437500,   0.9638,   0.0000,   1.3378,      0.00  ],
[ 1.462500,   0.9846,   1.4185,   1.8720,      0.00  ],
[ 1.487500,   1.1697,   0.0000,   1.8597,      0.00  ],
[ 1.512500,   1.4049,   0.3605,   1.0027,   8308.79  ],
[ 1.537500,   1.2294,   0.9858,   1.0560,     47.09  ],
[ 1.562500,   1.1765,   1.0022,   1.0583,      0.28  ],
[ 1.587500,   1.0855,   1.0401,   1.0582,     94.91  ],
[ 1.612500,   1.0960,   1.0336,   1.0575,      0.00  ],
[ 1.637500,   1.2077,   1.0923,   1.0535,     79.35  ],
[ 1.662500,   1.2380,   1.0007,   1.0530,     57.24  ],
[ 1.687500,   1.1170,   1.2089,   1.0547,    109.69  ],
[ 1.712500,   1.0433,   1.1063,   1.0569,      0.00  ],
[ 1.737500,   0.9510,   1.1354,   1.0560,      6.77  ],
[ 1.762500,   0.7669,   1.1635,   1.0563,      0.00  ],
[ 1.787500,   0.7564,   1.0737,   1.0557,    105.55  ],
[ 1.812500,   0.5387,   1.0607,   1.0594,    161.26  ],
[ 1.837500,   0.0000,   1.0463,   1.0607,    312.65  ],
[ 1.862500,   0.0000,   1.1331,   1.0621,    250.99  ],
[ 1.887500,   0.0000,   1.1367,   1.0629,    211.31  ],
[ 1.912500,   0.0000,   1.1555,   1.0641,    200.52  ],
[ 1.937500,   0.0000,   1.1185,   1.0647,    261.43  ],
[ 1.962500,   0.0000,   1.1228,   1.0683,    112.91  ],
[ 1.987500,   0.0000,   1.1018,   1.0685,    255.58  ],
[ 2.012500,   0.0000,   1.1166,   1.0697,    264.13  ],
[ 2.037500,   0.0000,   1.0625,   1.0714,    222.45  ],
[ 2.062500,   0.0000,   1.0698,   1.0716,    265.84  ],
[ 2.087500,   0.0000,   1.0519,   1.0728,    222.23  ],
[ 2.112500,   0.0000,   1.0154,   1.0756,    160.33  ],
[ 2.137500,   0.0000,   1.0404,   1.0750,    232.38  ],
[ 2.162500,   0.0000,   1.0628,   1.0765,    182.32  ],
[ 2.187500,   0.0000,   1.0412,   1.0785,    161.09  ],
[ 2.212500,   0.0000,   1.0224,   1.0803,    154.38  ],
[ 2.237500,   0.0000,   1.0394,   1.0802,    205.70  ],
[ 2.262500,   0.0000,   1.0331,   1.0826,    131.50  ],
[ 2.287500,   0.0000,   1.0506,   1.0831,    212.38  ],
[ 2.312500,   0.0000,   1.0555,   1.0857,    330.70  ],
[ 2.337500,   0.0000,   1.0676,   1.0886,    223.58  ],
[ 2.362500,   0.0000,   1.0357,   1.0898,    254.36  ],
[ 2.387500,   0.0000,   0.9806,   1.0936,    171.62  ],
[ 2.412500,   0.0000,   0.8819,   1.0934,    228.62  ],
[ 2.437500,   0.0000,   0.9529,   1.0963,    199.14  ],
[ 2.462500,   0.0000,   0.8650,   1.1052,    136.78  ],
[ 2.487500,   0.0000,   0.8026,   1.1147,    626.61  ]
] 


#######################################################################
# 3x7 cluster size photons.
CaloSwLongWeights_v3_gam37 = [
               # w0        w3       escale    eoffset 
[ 0.012500,   1.2304,   1.9246,   1.0470,    102.76  ],
[ 0.037500,   1.1978,   1.8114,   1.0444,    168.34  ],
[ 0.062500,   1.2574,   1.8459,   1.0408,    188.64  ],
[ 0.087500,   1.1927,   1.8010,   1.0417,    117.64  ],
[ 0.112500,   1.1776,   1.7897,   1.0430,     13.34  ],
[ 0.137500,   1.2314,   1.8290,   1.0418,     33.88  ],
[ 0.162500,   1.2252,   1.8831,   1.0408,     66.22  ],
[ 0.187500,   1.1762,   1.7752,   1.0419,     53.97  ],
[ 0.212500,   1.2102,   1.8851,   1.0407,     41.77  ],
[ 0.237500,   1.1664,   1.8504,   1.0413,     35.33  ],
[ 0.262500,   1.1436,   1.7707,   1.0416,     35.52  ],
[ 0.287500,   1.1601,   1.8561,   1.0403,     94.48  ],
[ 0.312500,   1.1633,   1.7289,   1.0402,     97.14  ],
[ 0.337500,   1.2095,   1.7828,   1.0401,     66.67  ],
[ 0.362500,   1.1962,   1.6445,   1.0398,     68.94  ],
[ 0.387500,   1.1841,   1.6796,   1.0395,     65.22  ],
[ 0.412500,   1.1720,   1.5542,   1.0395,    136.30  ],
[ 0.437500,   1.1700,   1.4439,   1.0401,     73.97  ],
[ 0.462500,   1.1634,   1.3136,   1.0398,    116.15  ],
[ 0.487500,   1.1690,   1.3353,   1.0402,     52.28  ],
[ 0.512500,   1.1646,   1.2932,   1.0393,     74.65  ],
[ 0.537500,   1.1770,   1.1693,   1.0384,    116.11  ],
[ 0.562500,   1.2046,   1.1813,   1.0378,    161.79  ],
[ 0.587500,   1.2179,   1.1498,   1.0377,    126.43  ],
[ 0.612500,   1.1916,   1.0801,   1.0383,    154.68  ],
[ 0.637500,   1.2182,   1.0787,   1.0382,    138.61  ],
[ 0.662500,   1.2189,   0.9796,   1.0387,    115.31  ],
[ 0.687500,   1.2143,   1.0189,   1.0377,    163.35  ],
[ 0.712500,   1.2303,   0.9666,   1.0390,    122.29  ],
[ 0.737500,   1.2316,   0.9874,   1.0389,    171.34  ],
[ 0.762500,   1.2158,   1.0543,   1.0392,    278.14  ],
[ 0.787500,   1.1928,   1.0372,   1.0440,    467.65  ],
[ 0.812500,   1.1914,   1.0309,   1.0528,    267.31  ],
[ 0.837500,   1.2540,   1.2857,   1.0512,    182.68  ],
[ 0.862500,   1.2449,   1.2400,   1.0496,    250.00  ],
[ 0.887500,   1.2335,   1.1727,   1.0510,    159.40  ],
[ 0.912500,   1.2384,   1.1335,   1.0535,     69.44  ],
[ 0.937500,   1.2527,   1.0530,   1.0519,    238.22  ],
[ 0.962500,   1.2204,   0.9382,   1.0551,     33.89  ],
[ 0.987500,   1.2224,   0.9533,   1.0524,    321.12  ],
[ 1.012500,   1.2131,   0.9360,   1.0534,    235.09  ],
[ 1.037500,   1.2191,   0.8442,   1.0558,    215.43  ],
[ 1.062500,   1.2238,   0.9432,   1.0548,    311.07  ],
[ 1.087500,   1.2039,   0.8682,   1.0567,    179.49  ],
[ 1.112500,   1.1933,   0.8551,   1.0584,    176.86  ],
[ 1.137500,   1.1699,   0.6586,   1.0621,    206.57  ],
[ 1.162500,   1.1701,   0.6121,   1.0643,    249.65  ],
[ 1.187500,   1.1877,   0.8368,   1.0596,    347.87  ],
[ 1.212500,   1.1425,   0.8250,   1.0628,    227.87  ],
[ 1.237500,   1.1282,   0.6655,   1.0694,    148.49  ],
[ 1.262500,   1.1284,   0.6587,   1.0660,    295.51  ],
[ 1.287500,   1.1351,   0.7793,   1.0672,    279.79  ],
[ 1.312500,   1.1190,   0.4040,   1.0759,    372.52  ],
[ 1.337500,   1.1193,   0.4031,   1.0793,    234.30  ],
[ 1.362500,   1.0990,   0.1830,   1.0787,    840.36  ],
[ 1.387500,   1.0991,   0.0000,   1.1350,    544.18  ],
[ 1.412500,   1.0247,   0.0000,   1.1862,    415.93  ],
[ 1.437500,   0.9833,   0.0000,   1.3040,      0.00  ],
[ 1.462500,   0.9878,   0.0000,   1.8607,      0.00  ],
[ 1.487500,   1.1283,   0.0000,   1.8557,      0.00  ],
[ 1.512500,   1.5991,   0.9377,   1.0001,   4577.19  ],
[ 1.537500,   1.2030,   0.9377,   1.0501,    107.35  ],
[ 1.562500,   1.1337,   0.9592,   1.0529,     15.09  ],
[ 1.587500,   1.0144,   1.0272,   1.0535,     77.68  ],
[ 1.612500,   1.0451,   1.0422,   1.0521,      0.00  ],
[ 1.637500,   1.1628,   1.1221,   1.0488,      0.00  ],
[ 1.662500,   1.2170,   1.0474,   1.0458,    156.98  ],
[ 1.687500,   1.1131,   1.2690,   1.0503,      0.00  ],
[ 1.712500,   1.0343,   1.1254,   1.0514,      0.00  ],
[ 1.737500,   0.9055,   1.1059,   1.0505,     26.15  ],
[ 1.762500,   0.7359,   1.1607,   1.0505,      0.00  ],
[ 1.787500,   0.7086,   1.0649,   1.0496,     97.77  ],
[ 1.812500,   0.4937,   1.0645,   1.0532,    146.53  ],
[ 1.837500,   0.0000,   1.0482,   1.0555,    297.94  ],
[ 1.862500,   0.0000,   1.1225,   1.0530,    451.33  ],
[ 1.887500,   0.0000,   1.1323,   1.0561,    248.39  ],
[ 1.912500,   0.0000,   1.1735,   1.0576,    173.38  ],
[ 1.937500,   0.0000,   1.1278,   1.0580,    250.13  ],
[ 1.962500,   0.0000,   1.1122,   1.0604,    207.52  ],
[ 1.987500,   0.0000,   1.1055,   1.0613,    248.12  ],
[ 2.012500,   0.0000,   1.1042,   1.0644,    159.01  ],
[ 2.037500,   0.0000,   1.0566,   1.0637,    233.18  ],
[ 2.062500,   0.0000,   1.0811,   1.0644,    218.35  ],
[ 2.087500,   0.0000,   1.0604,   1.0642,    260.02  ],
[ 2.112500,   0.0000,   1.0198,   1.0667,    183.28  ],
[ 2.137500,   0.0000,   1.0409,   1.0666,    213.60  ],
[ 2.162500,   0.0000,   1.0615,   1.0685,    147.14  ],
[ 2.187500,   0.0000,   1.0377,   1.0707,    143.10  ],
[ 2.212500,   0.0000,   1.0306,   1.0713,    180.12  ],
[ 2.237500,   0.0000,   1.0486,   1.0711,    211.72  ],
[ 2.262500,   0.0000,   1.0354,   1.0735,    137.95  ],
[ 2.287500,   0.0000,   1.0472,   1.0746,    175.50  ],
[ 2.312500,   0.0000,   1.0614,   1.0777,    196.48  ],
[ 2.337500,   0.0000,   1.0581,   1.0798,    231.88  ],
[ 2.362500,   0.0000,   1.0309,   1.0810,    215.53  ],
[ 2.387500,   0.0000,   0.9822,   1.0843,    105.21  ],
[ 2.412500,   0.0000,   0.8744,   1.0881,     47.17  ],
[ 2.437500,   0.0000,   0.9230,   1.0860,    272.58  ],
[ 2.462500,   0.0000,   0.8689,   1.0959,    110.56  ],
[ 2.487500,   0.0000,   0.7638,   1.1095,      0.00  ]
] 


#######################################################################


# For this version, the photon corrections are the same as the electrons.
class CaloSwLongWeights_v3_parms:
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
    correction = {'ele55' : CaloSwLongWeights_v3_ele55,
                  'ele35' : CaloSwLongWeights_v3_ele35,
                  'ele37' : CaloSwLongWeights_v3_ele37,
                  'gam55' : CaloSwLongWeights_v3_ele55,
                  'gam35' : CaloSwLongWeights_v3_ele35,
                  'gam37' : CaloSwLongWeights_v3_ele37}



# This version has separate electron and photon corrections.
class CaloSwLongWeights_v3_1_parms:
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
    correction = {'ele55' : CaloSwLongWeights_v3_ele55,
                  'ele35' : CaloSwLongWeights_v3_ele35,
                  'ele37' : CaloSwLongWeights_v3_ele37,
                  'gam55' : CaloSwLongWeights_v3_gam55,
                  'gam35' : CaloSwLongWeights_v3_gam35,
                  'gam37' : CaloSwLongWeights_v3_gam37,

                  # Use 5x5 for cluster sizes that aren't explicitly derived.
                  'ele33' : CaloSwLongWeights_v3_ele55,
                  'ele57' : CaloSwLongWeights_v3_ele55,
                  'ele77' : CaloSwLongWeights_v3_ele55,
                  'gam33' : CaloSwLongWeights_v3_gam55,
                  'gam57' : CaloSwLongWeights_v3_gam55,
                  'gam77' : CaloSwLongWeights_v3_gam55,
                  }
    
