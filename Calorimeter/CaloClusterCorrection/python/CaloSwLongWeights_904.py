# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration

#
# $Id: CaloSwLongWeights_904.py,v 1.1 2006-11-10 03:47:27 ssnyder Exp $
#
# File: CaloClusterCorrection/python/CaloSwLongWeights_904.py
# Created: Nov 2006, sss
# Purpose: Longitudinal weights corrections, from 9.0.4 MC.
#
# This version of the weight corrections was derived by Stathes Paganis
# from 9.0.4 MC.  Added in LArClusterRec-02-05-32, in 10.0.1.
#


from CaloClusterCorrection.common import *

#######################################################################
# Stathes's weights from 9.0.4 MC, 5x5 cluster size.
CaloSwLongWeights_lwc904_55 = [
               # w0        w3       escale    eoffset
[ 0.012500,   0.9343,   1.1867,   1.0448,    394.03  ],
[ 0.037500,   0.9268,   1.0577,   1.0455,    231.82  ],
[ 0.062500,   0.9368,   1.0342,   1.0415,    279.60  ],
[ 0.087500,   0.9260,   1.3891,   1.0370,    347.31  ],
[ 0.112500,   0.9073,   0.9917,   1.0405,    274.33  ],
[ 0.137500,   0.9652,   1.1285,   1.0395,    241.92  ],
[ 0.162500,   0.9346,   0.9778,   1.0406,    266.90  ],
[ 0.187500,   0.9558,   1.0634,   1.0397,    270.84  ],
[ 0.212500,   0.8699,   0.9729,   1.0409,    323.17  ],
[ 0.237500,   0.9303,   1.3831,   1.0372,    326.01  ],
[ 0.262500,   0.9623,   1.2417,   1.0376,    263.50  ],
[ 0.287500,   0.8749,   0.7676,   1.0406,    298.02  ],
[ 0.312500,   0.9207,   0.9720,   1.0386,    330.80  ],
[ 0.337500,   0.9028,   0.5574,   1.0425,    260.96  ],
[ 0.362500,   0.8960,   0.6516,   1.0412,    302.13  ],
[ 0.387500,   0.9316,   0.5173,   1.0416,    254.88  ],
[ 0.412500,   0.9792,   0.6303,   1.0389,    302.57  ],
[ 0.437500,   0.9797,   0.5677,   1.0399,    257.54  ],
[ 0.462500,   0.9501,   0.7712,   1.0383,    305.92  ],
[ 0.487500,   0.9449,   0.8346,   1.0387,    354.37  ],
[ 0.512500,   0.9501,   0.7885,   1.0377,    341.62  ],
[ 0.537500,   0.9446,   0.7755,   1.0383,    298.63  ],
[ 0.562500,   0.9491,   0.7907,   1.0371,    369.21  ],
[ 0.587500,   0.9974,   0.5588,   1.0416,    224.43  ],
[ 0.612500,   0.9686,   0.6554,   1.0410,    314.70  ],
[ 0.637500,   1.0264,   0.9027,   1.0348,    496.14  ],
[ 0.662500,   1.0352,   0.8275,   1.0388,    350.94  ],
[ 0.687500,   1.0646,   0.8034,   1.0396,    377.58  ],
[ 0.712500,   1.0122,   0.6641,   1.0379,    661.55  ],
[ 0.737500,   1.0706,   0.7990,   1.0376,    653.28  ],
[ 0.762500,   1.0728,   0.9123,   1.0400,    764.24  ],
[ 0.787500,   1.0794,   1.1530,   1.0457,    732.03  ],
[ 0.812500,   1.0562,   1.1189,   1.0442,    926.27  ],
[ 0.837500,   1.0791,   1.2835,   1.0401,    801.79  ],
[ 0.862500,   1.0268,   0.8689,   1.0448,    734.06  ],
[ 0.887500,   1.0247,   0.5932,   1.0496,    550.66  ],
[ 0.912500,   0.9820,   0.5473,   1.0492,    807.07  ],
[ 0.937500,   0.9957,   0.9838,   1.0452,    974.74  ],
[ 0.962500,   1.0279,   0.8905,   1.0466,    738.45  ],
[ 0.987500,   1.0229,   0.8480,   1.0441,   1199.41  ],
[ 1.012500,   1.0552,   0.5227,   1.0546,    591.82  ],
[ 1.037500,   1.0258,   0.5747,   1.0554,    545.29  ],
[ 1.062500,   0.9943,   0.7615,   1.0566,    741.43  ],
[ 1.087500,   0.9748,   0.5806,   1.0522,   1196.85  ],
[ 1.112500,   1.0228,   0.7716,   1.0519,    888.87  ],
[ 1.137500,   1.0049,   0.7174,   1.0541,   1190.06  ],
[ 1.162500,   1.0075,   0.6331,   1.0558,   1101.85  ],
[ 1.187500,   0.9952,   0.5912,   1.0598,    985.82  ],
[ 1.212500,   0.9969,   0.5459,   1.0624,   1295.31  ],
[ 1.237500,   0.9866,   0.5321,   1.0690,   1124.31  ],
[ 1.262500,   1.0228,   0.6740,   1.0672,   1207.39  ],
[ 1.287500,   1.0223,   0.5405,   1.0747,    869.47  ],
[ 1.312500,   1.0390,   0.6780,   1.0602,   1080.48  ],
[ 1.337500,   1.0082,   0.7914,   1.0665,    902.94  ],
[ 1.362500,   1.0171,   0.8837,   1.0625,   1131.72  ],
[ 1.387500,   1.0021,   0.1942,   1.0823,   1005.70  ],
[ 1.412500,   1.0041,   0.0000,   1.0913,   1485.55  ],
[ 1.437500,   1.0099,   3.0936,   1.2246,    414.66  ],
[ 1.462500,   0.9669,   1.5888,   1.2401,   4045.49  ],
[ 1.487500,   0.9836,   0.0000,   1.8763,    532.63  ],
[ 1.512500,   1.2030,   0.0000,   1.3768,   3497.43  ],
[ 1.537500,   1.3200,   1.6057,   0.9326,   8438.36  ],
[ 1.562500,   1.1712,   0.8722,   1.0258,   1207.75  ],
[ 1.587500,   1.0714,   0.7994,   1.0255,   1376.02  ],
[ 1.612500,   1.1066,   0.6338,   1.0308,   1231.62  ],
[ 1.637500,   1.1452,   1.0362,   1.0243,   1241.25  ],
[ 1.662500,   1.2470,   0.7528,   1.0281,   1447.92  ],
[ 1.687500,   1.0919,   1.1062,   1.0502,   1046.75  ],
[ 1.712500,   1.0502,   1.0055,   1.0234,    891.85  ],
[ 1.737500,   0.9732,   0.8125,   1.0288,    581.87  ],
[ 1.762500,   0.9570,   1.0825,   1.0260,    693.45  ],
[ 1.787500,   0.8658,   0.8170,   1.0280,    826.53  ],
[ 1.812500,   0.7129,   0.8628,   1.0349,   1153.20  ],
[ 1.837500,   0.0000,   0.7408,   1.0408,   1289.23  ],
[ 1.862500,   0.0000,   0.7918,   1.0449,    791.21  ],
[ 1.887500,   0.0000,   0.8012,   1.0352,   1049.38  ],
[ 1.912500,   0.0000,   0.9658,   1.0385,    799.23  ],
[ 1.937500,   0.0000,   0.8229,   1.0404,    786.25  ],
[ 1.962500,   0.0000,   0.9296,   1.0402,    802.84  ],
[ 1.987500,   0.0000,   0.9209,   1.0377,   1096.35  ],
[ 2.012500,   0.0000,   0.8103,   1.0423,    923.86  ],
[ 2.037500,   0.0000,   0.6936,   1.0531,    569.91  ],
[ 2.062500,   0.0000,   0.8557,   1.0463,    753.87  ],
[ 2.087500,   0.0000,   0.8838,   1.0453,    791.14  ],
[ 2.112500,   0.0000,   0.9726,   1.0469,    819.39  ],
[ 2.137500,   0.0000,   0.8607,   1.0444,   1167.60  ],
[ 2.162500,   0.0000,   0.9665,   1.0490,    797.37  ],
[ 2.187500,   0.0000,   0.9452,   1.0518,    708.08  ],
[ 2.212500,   0.0000,   0.9224,   1.0546,    657.72  ],
[ 2.237500,   0.0000,   1.0537,   1.0463,   1132.76  ],
[ 2.262500,   0.0000,   0.8840,   1.0539,    720.57  ],
[ 2.287500,   0.0000,   0.8347,   1.0584,    485.20  ],
[ 2.312500,   0.0000,   0.8347,   1.0584,    485.20  ],
[ 2.337500,   0.0000,   0.8347,   1.0584,    485.20  ],
[ 2.362500,   0.0000,   0.9633,   1.0525,    776.85  ],
[ 2.387500,   0.0000,   1.0539,   1.0459,   1351.36  ],
[ 2.412500,   0.0000,   0.9633,   1.0525,    776.85  ],
[ 2.437500,   0.0000,   0.9633,   1.0625,    776.85  ],
[ 2.462500,   0.0000,   0.9633,   1.0725,    776.85  ],
[ 2.487500,   0.0000,   0.9633,   1.0925,    776.85  ]
    ]


# Stathes's weights from 9.0.4 MC, 3x7 cluster size.
CaloSwLongWeights_lwc904_37 = [
               # w0        w3       escale    eoffset
[ 0.012500,   0.8892,   0.9032,   1.0630,    291.57  ],
[ 0.037500,   0.9300,   1.1016,   1.0586,    285.16  ],
[ 0.062500,   0.9517,   1.1392,   1.0564,    222.00  ],
[ 0.087500,   0.9465,   1.2823,   1.0540,    230.89  ],
[ 0.112500,   0.9129,   0.9628,   1.0568,    216.19  ],
[ 0.137500,   0.9706,   1.1101,   1.0535,    262.42  ],
[ 0.162500,   0.9372,   0.9033,   1.0575,    134.48  ],
[ 0.187500,   0.9540,   0.9241,   1.0562,    188.00  ],
[ 0.212500,   0.8626,   0.9307,   1.0570,    241.32  ],
[ 0.237500,   0.9257,   1.1502,   1.0544,    241.92  ],
[ 0.262500,   0.9712,   1.2229,   1.0531,    220.07  ],
[ 0.287500,   0.8730,   0.6762,   1.0572,    237.34  ],
[ 0.312500,   0.9270,   0.8581,   1.0555,    251.88  ],
[ 0.337500,   0.8918,   0.6355,   1.0565,    239.41  ],
[ 0.362500,   0.8916,   0.5569,   1.0566,    253.06  ],
[ 0.387500,   0.9194,   0.5919,   1.0567,    233.31  ],
[ 0.412500,   0.9526,   0.5690,   1.0559,    209.50  ],
[ 0.437500,   0.9948,   0.5503,   1.0552,    206.71  ],
[ 0.462500,   0.9205,   0.6765,   1.0547,    247.26  ],
[ 0.487500,   0.9394,   0.7648,   1.0547,    276.76  ],
[ 0.512500,   0.9241,   0.7079,   1.0547,    281.51  ],
[ 0.537500,   0.9317,   0.7267,   1.0552,    223.95  ],
[ 0.562500,   0.9535,   0.6460,   1.0554,    271.93  ],
[ 0.587500,   0.9715,   0.4642,   1.0580,    166.56  ],
[ 0.612500,   0.9769,   0.5155,   1.0577,    240.41  ],
[ 0.637500,   1.0141,   0.8449,   1.0507,    447.44  ],
[ 0.662500,   1.0311,   0.8733,   1.0522,    404.67  ],
[ 0.687500,   1.0565,   0.7707,   1.0554,    315.85  ],
[ 0.712500,   1.0013,   0.6132,   1.0578,    430.26  ],
[ 0.737500,   1.0541,   0.7684,   1.0539,    578.92  ],
[ 0.762500,   1.0528,   0.8208,   1.0544,    668.04  ],
[ 0.787500,   1.0505,   0.9249,   1.0624,    613.96  ],
[ 0.812500,   1.0509,   1.3189,   1.0572,    970.35  ],
[ 0.837500,   1.0767,   1.2807,   1.0620,    461.66  ],
[ 0.862500,   1.0196,   0.8965,   1.0624,    677.53  ],
[ 0.887500,   1.0265,   0.4980,   1.0652,    708.89  ],
[ 0.912500,   1.0355,   0.4591,   1.0685,    604.19  ],
[ 0.937500,   1.0106,   0.9454,   1.0620,    833.61  ],
[ 0.962500,   1.0012,   0.7778,   1.0677,    497.80  ],
[ 0.987500,   1.0085,   0.6660,   1.0696,    706.91  ],
[ 1.012500,   1.0041,   0.4212,   1.0765,    518.46  ],
[ 1.037500,   1.0241,   0.4562,   1.0723,    479.77  ],
[ 1.062500,   1.0119,   0.7416,   1.0651,    956.73  ],
[ 1.087500,   0.9551,   0.4185,   1.0718,   1062.27  ],
[ 1.112500,   1.0220,   0.6132,   1.0744,    672.87  ],
[ 1.137500,   1.0105,   0.5200,   1.0720,   1057.30  ],
[ 1.162500,   1.0189,   0.5727,   1.0757,    740.99  ],
[ 1.187500,   1.0002,   0.4458,   1.0810,    669.58  ],
[ 1.212500,   1.0029,   0.5328,   1.0799,   1031.32  ],
[ 1.237500,   1.0063,   0.4499,   1.0813,    906.67  ],
[ 1.262500,   0.9931,   0.6120,   1.0812,   1012.75  ],
[ 1.287500,   1.0041,   0.3659,   1.0927,    934.92  ],
[ 1.312500,   1.0076,   0.6030,   1.0742,   1105.42  ],
[ 1.337500,   1.0232,   0.9876,   1.0601,   1123.00  ],
[ 1.362500,   1.0197,   0.7940,   1.0567,   1057.89  ],
[ 1.387500,   0.9978,   0.0000,   1.0929,   1085.28  ],
[ 1.412500,   0.9913,   0.0000,   1.0925,    994.85  ],
[ 1.437500,   1.0000,   1.0000,   1.0500,    200.00  ],
[ 1.462500,   1.0000,   1.0000,   1.0500,    200.00  ],
[ 1.487500,   0.8357,   0.0000,   1.9602,    729.94  ],
[ 1.512500,   0.8134,   0.3616,   1.6136,   2946.56  ],
[ 1.537500,   1.3166,   1.1273,   1.0197,   1950.41  ],
[ 1.562500,   1.2000,   0.8741,   1.0293,   1386.19  ],
[ 1.587500,   1.1108,   0.7276,   1.0362,   1254.45  ],
[ 1.612500,   1.1264,   0.6335,   1.0454,    783.84  ],
[ 1.637500,   1.1338,   0.8817,   1.0378,   1190.53  ],
[ 1.662500,   1.2454,   0.5257,   1.0417,   1298.18  ],
[ 1.687500,   1.1051,   0.8468,   1.0468,   1107.36  ],
[ 1.712500,   1.0817,   0.8677,   1.0454,    625.64  ],
[ 1.737500,   0.9662,   0.7460,   1.0466,    592.16  ],
[ 1.762500,   0.9385,   1.0220,   1.0453,    450.67  ],
[ 1.787500,   0.8529,   0.7996,   1.0532,    422.59  ],
[ 1.812500,   0.6995,   0.7658,   1.0522,   1232.79  ],
[ 1.837500,   0.0000,   0.6848,   1.0642,    956.80  ],
[ 1.862500,   0.0000,   0.8379,   1.0610,    887.55  ],
[ 1.887500,   0.0000,   0.7217,   1.0589,    815.04  ],
[ 1.912500,   0.0000,   0.9994,   1.0577,    736.42  ],
[ 1.937500,   0.0000,   0.8090,   1.0611,    694.32  ],
[ 1.962500,   0.0000,   0.8933,   1.0615,    720.61  ],
[ 1.987500,   0.0000,   0.8778,   1.0619,    787.72  ],
[ 2.012500,   0.0000,   0.7399,   1.0671,    719.25  ],
[ 2.037500,   0.0000,   0.4216,   1.0714,    687.65  ],
[ 2.062500,   0.0000,   0.8238,   1.0717,    603.10  ],
[ 2.087500,   0.0000,   0.9046,   1.0689,    698.10  ],
[ 2.112500,   0.0000,   0.6974,   1.0735,    688.43  ],
[ 2.137500,   0.0000,   0.9092,   1.0713,    736.23  ],
[ 2.162500,   0.0000,   0.9138,   1.0725,    697.00  ],
[ 2.187500,   0.0000,   0.8419,   1.0738,    781.59  ],
[ 2.212500,   0.0000,   0.9369,   1.0791,    485.75  ],
[ 2.237500,   0.0000,   0.9528,   1.0751,    836.69  ],
[ 2.262500,   0.0000,   0.8920,   1.0766,    730.50  ],
[ 2.287500,   0.0000,   0.8646,   1.0806,    560.12  ],
[ 2.312500,   0.0000,   0.8646,   1.0806,    560.12  ],
[ 2.337500,   0.0000,   0.8646,   1.0806,    560.12  ],
[ 2.362500,   0.0000,   1.0179,   1.0732,   1106.06  ],
[ 2.387500,   0.0000,   0.9999,   1.0756,    954.93  ],
[ 2.412500,   0.0000,   0.9999,   1.0756,    954.93  ],
[ 2.437500,   0.0000,   0.9999,   1.0756,    954.93  ],
[ 2.462500,   0.0000,   0.9999,   1.0756,    954.93  ],
[ 2.487500,   0.0000,   0.3876,   1.1039,     51.26  ]
    ]

# Stathes's weights from 9.0.4 MC, 3x5 cluster size.
CaloSwLongWeights_lwc904_35 = [
               # w0        w3       escale    eoffset
[ 0.012500,   0.9228,   1.1606,   1.0651,    332.16  ],
[ 0.037500,   0.9531,   1.2393,   1.0619,    299.04  ],
[ 0.062500,   0.9640,   1.2332,   1.0602,    230.84  ],
[ 0.087500,   0.9865,   1.4762,   1.0566,    244.50  ],
[ 0.112500,   0.9438,   1.1244,   1.0589,    254.87  ],
[ 0.137500,   0.9992,   1.2007,   1.0578,    240.17  ],
[ 0.162500,   0.9776,   1.0558,   1.0594,    230.90  ],
[ 0.187500,   0.9931,   1.0823,   1.0596,    196.89  ],
[ 0.212500,   0.8926,   0.9984,   1.0602,    295.85  ],
[ 0.237500,   0.9528,   1.4422,   1.0567,    274.14  ],
[ 0.262500,   1.0026,   1.4036,   1.0564,    259.46  ],
[ 0.287500,   0.9154,   0.8086,   1.0610,    237.02  ],
[ 0.312500,   0.9347,   0.9250,   1.0595,    249.07  ],
[ 0.337500,   0.9515,   0.6759,   1.0590,    321.99  ],
[ 0.362500,   0.9019,   0.5605,   1.0607,    286.29  ],
[ 0.387500,   0.9697,   0.7745,   1.0599,    235.99  ],
[ 0.412500,   1.0000,   0.6023,   1.0593,    244.45  ],
[ 0.437500,   1.0268,   0.7170,   1.0580,    263.82  ],
[ 0.462500,   0.9829,   0.7058,   1.0578,    288.63  ],
[ 0.487500,   0.9769,   0.8456,   1.0581,    295.45  ],
[ 0.512500,   0.9572,   0.8988,   1.0561,    411.85  ],
[ 0.537500,   0.9767,   0.8343,   1.0578,    268.96  ],
[ 0.562500,   0.9651,   0.8121,   1.0565,    395.77  ],
[ 0.587500,   1.0088,   0.5408,   1.0613,    228.50  ],
[ 0.612500,   0.9755,   0.7164,   1.0600,    333.25  ],
[ 0.637500,   1.0420,   1.0077,   1.0531,    518.80  ],
[ 0.662500,   1.0472,   0.9101,   1.0574,    383.17  ],
[ 0.687500,   1.0591,   0.8601,   1.0566,    591.02  ],
[ 0.712500,   1.0209,   0.6403,   1.0606,    589.34  ],
[ 0.737500,   1.0725,   0.8405,   1.0578,    650.82  ],
[ 0.762500,   1.0600,   1.1002,   1.0561,    876.04  ],
[ 0.787500,   1.0779,   1.1511,   1.0718,    501.47  ],
[ 0.812500,   1.0873,   1.3652,   1.0589,   1240.53  ],
[ 0.837500,   1.0907,   1.5369,   1.0602,    907.09  ],
[ 0.862500,   1.0362,   0.9148,   1.0695,    633.31  ],
[ 0.887500,   1.0461,   0.5562,   1.0683,    830.54  ],
[ 0.912500,   1.0255,   0.5571,   1.0724,    742.31  ],
[ 0.937500,   0.9997,   1.0666,   1.0673,   1049.66  ],
[ 0.962500,   1.0127,   0.9345,   1.0709,    711.72  ],
[ 0.987500,   1.0409,   0.8904,   1.0681,   1065.23  ],
[ 1.012500,   1.0339,   0.4266,   1.0816,    634.21  ],
[ 1.037500,   1.0553,   0.6045,   1.0768,    763.40  ],
[ 1.062500,   0.9827,   0.6141,   1.0830,    812.46  ],
[ 1.087500,   1.0072,   0.5523,   1.0796,    963.27  ],
[ 1.112500,   1.0200,   0.6091,   1.0813,    770.35  ],
[ 1.137500,   0.9959,   0.6207,   1.0805,   1049.84  ],
[ 1.162500,   1.0223,   0.5104,   1.0903,    760.52  ],
[ 1.187500,   0.9921,   0.5236,   1.0879,    963.03  ],
[ 1.212500,   1.0375,   0.4606,   1.0826,   1443.04  ],
[ 1.237500,   1.0120,   0.3027,   1.0898,   1048.63  ],
[ 1.262500,   1.0201,   0.6342,   1.0914,   1169.70  ],
[ 1.287500,   1.0129,   0.5491,   1.0970,   1017.87  ],
[ 1.312500,   1.0494,   0.6029,   1.0879,   1077.24  ],
[ 1.337500,   1.0264,   0.9500,   1.0803,   1115.53  ],
[ 1.362500,   1.0198,   1.1023,   1.0713,   1705.04  ],
[ 1.387500,   1.0069,   0.0000,   1.0915,   1888.33  ],
[ 1.412500,   1.0152,   0.0000,   1.0945,   1856.34  ],
[ 1.437500,   1.0152,   0.0000,   1.0945,   1856.34  ],
[ 1.462500,   1.0152,   0.0000,   1.0945,   1856.34  ],
[ 1.487500,   0.8527,   0.0000,   1.8586,   1258.57  ],
[ 1.512500,   0.8518,   0.0000,   1.8816,   1099.45  ],
[ 1.537500,   1.3276,   1.6994,   0.9402,   8856.68  ],
[ 1.562500,   1.2440,   0.8959,   1.0428,   1248.09  ],
[ 1.587500,   1.1168,   0.6786,   1.0460,   1269.71  ],
[ 1.612500,   1.1430,   0.6154,   1.0465,   1478.23  ],
[ 1.637500,   1.1687,   1.0108,   1.0397,   1564.35  ],
[ 1.662500,   1.2987,   0.9178,   1.0462,   1670.90  ],
[ 1.687500,   1.1213,   1.0466,   1.0525,   1214.16  ],
[ 1.712500,   1.0855,   1.0627,   1.0493,    749.46  ],
[ 1.737500,   0.9834,   0.8117,   1.0527,    575.04  ],
[ 1.762500,   0.9471,   1.0871,   1.0492,    633.21  ],
[ 1.787500,   0.8620,   0.8148,   1.0598,    498.94  ],
[ 1.812500,   0.7145,   0.8325,   1.0592,   1224.71  ],
[ 1.837500,   0.0000,   0.6673,   1.0714,   1027.39  ],
[ 1.862500,   0.0000,   0.8234,   1.0687,    847.57  ],
[ 1.887500,   0.0000,   0.7228,   1.0649,    943.13  ],
[ 1.912500,   0.0000,   1.0028,   1.0652,    776.91  ],
[ 1.937500,   0.0000,   0.8586,   1.0680,    787.55  ],
[ 1.962500,   0.0000,   0.9345,   1.0679,    804.07  ],
[ 1.987500,   0.0000,   0.9014,   1.0703,    789.12  ],
[ 2.012500,   0.0000,   0.8104,   1.0755,    711.20  ],
[ 2.037500,   0.0000,   0.5367,   1.0780,    761.96  ],
[ 2.062500,   0.0000,   0.8259,   1.0778,    705.46  ],
[ 2.087500,   0.0000,   0.9065,   1.0760,    807.11  ],
[ 2.112500,   0.0000,   0.8529,   1.0788,    853.30  ],
[ 2.137500,   0.0000,   0.8397,   1.0758,   1208.62  ],
[ 2.162500,   0.0000,   0.9673,   1.0819,    700.92  ],
[ 2.187500,   0.0000,   0.9317,   1.0790,    950.69  ],
[ 2.212500,   0.0000,   0.9120,   1.0845,    825.22  ],
[ 2.237500,   0.0000,   1.0158,   1.0858,    688.01  ],
[ 2.262500,   0.0000,   0.7985,   1.0863,    856.60  ],
[ 2.287500,   0.0000,   0.7544,   1.0956,    469.72  ],
[ 2.312500,   0.0000,   0.7544,   1.0956,    469.72  ],
[ 2.337500,   0.0000,   0.7544,   1.0956,    469.72  ],
[ 2.362500,   0.0000,   0.9612,   1.0849,   1063.86  ],
[ 2.387500,   0.0000,   1.0352,   1.0847,   1001.93  ],
[ 2.412500,   0.0000,   1.0352,   1.0847,   1001.93  ],
[ 2.437500,   0.0000,   1.0352,   1.0847,   1001.93  ],
[ 2.462500,   0.0000,   1.0352,   1.0847,   1001.93  ],
[ 2.487500,   0.0000,   0.3921,   1.1125,    195.27  ]
    ]


#######################################################################

# The 9.0.4 corrections, with the gap region excluded.
# This was used for the `DC2new' corrections.
class CaloSwLongWeights_904_parms:
    eta_start_crack = 1.375
    eta_end_crack = 1.525
    etamax = 2.5
    use_raw_eta = False

    # If "preserve_offset" is set to True, then any offset that has
    # been applied prior to this correction is preserved. We do not want
    # to do this here because gapCorrections and lwcorrections are exclusive.
    preserve_offset = False

    region = CALOCORR_COMBINED2
    degree = 3
    correction = {'ele55' : CaloSwLongWeights_lwc904_55,
                  'ele35' : CaloSwLongWeights_lwc904_35,
                  'ele37' : CaloSwLongWeights_lwc904_37,
                  'gam55' : CaloSwLongWeights_lwc904_55,
                  'gam35' : CaloSwLongWeights_lwc904_35,
                  'gam37' : CaloSwLongWeights_lwc904_37}


# The 9.0.4 corrections, applying them in the gap too.
# This was used for the `Rome' corrections.
class CaloSwLongWeights_904gap_parms:
    eta_start_crack = 1.5
    eta_end_crack = 1.5
    etamax = 2.5
    use_raw_eta = False

    # We don't want to throw away scintillator energy, so set this to true.
    # The scintillator contribution will still get scaled by the
    # multiplicative factor though...
    preserve_offset = True

    region = CALOCORR_COMBINED2
    degree = 3
    correction = {'ele55' : CaloSwLongWeights_lwc904_55,
                  'ele35' : CaloSwLongWeights_lwc904_35,
                  'ele37' : CaloSwLongWeights_lwc904_37,
                  'gam55' : CaloSwLongWeights_lwc904_55,
                  'gam35' : CaloSwLongWeights_lwc904_35,
                  'gam37' : CaloSwLongWeights_lwc904_37}
