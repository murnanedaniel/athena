# Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration

def setupProfile(flags, scaleTaskLength=1):

  def _evts(x):
    return int(scaleTaskLength * x)

  return [
    {'run':284500, 'lb':1, 'starttstamp':1446539185, 'dt':0.000, 'evts':_evts(1), 'mu':0.500, 'force_new':False},
    {'run':284500, 'lb':2, 'starttstamp':1446539245, 'dt':0.000, 'evts':_evts(1), 'mu':1.500, 'force_new':False},
    {'run':284500, 'lb':3, 'starttstamp':1446539305, 'dt':0.000, 'evts':_evts(1), 'mu':2.500, 'force_new':False},
    {'run':284500, 'lb':4, 'starttstamp':1446539365, 'dt':0.000, 'evts':_evts(1), 'mu':3.500, 'force_new':False},
    {'run':284500, 'lb':5, 'starttstamp':1446539425, 'dt':0.000, 'evts':_evts(1), 'mu':4.500, 'force_new':False},
    {'run':284500, 'lb':6, 'starttstamp':1446539485, 'dt':0.000, 'evts':_evts(1), 'mu':5.500, 'force_new':False},
    {'run':284500, 'lb':7, 'starttstamp':1446539545, 'dt':0.000, 'evts':_evts(1), 'mu':6.500, 'force_new':False},
    {'run':284500, 'lb':8, 'starttstamp':1446539605, 'dt':0.000, 'evts':_evts(1), 'mu':7.500, 'force_new':False},
    {'run':284500, 'lb':9, 'starttstamp':1446539665, 'dt':0.000, 'evts':_evts(6), 'mu':8.500, 'force_new':False},
    {'run':284500, 'lb':10, 'starttstamp':1446539725, 'dt':0.000, 'evts':_evts(16), 'mu':9.500, 'force_new':False},
    {'run':284500, 'lb':11, 'starttstamp':1446539785, 'dt':0.000, 'evts':_evts(29), 'mu':10.500, 'force_new':False},
    {'run':284500, 'lb':12, 'starttstamp':1446539845, 'dt':0.000, 'evts':_evts(43), 'mu':11.500, 'force_new':False},
    {'run':284500, 'lb':13, 'starttstamp':1446539905, 'dt':0.000, 'evts':_evts(55), 'mu':12.500, 'force_new':False},
    {'run':284500, 'lb':14, 'starttstamp':1446539965, 'dt':0.000, 'evts':_evts(62), 'mu':13.500, 'force_new':False},
    {'run':284500, 'lb':15, 'starttstamp':1446540025, 'dt':0.000, 'evts':_evts(68), 'mu':14.500, 'force_new':False},
    {'run':284500, 'lb':16, 'starttstamp':1446540085, 'dt':0.000, 'evts':_evts(74), 'mu':15.500, 'force_new':False},
    {'run':284500, 'lb':17, 'starttstamp':1446540145, 'dt':0.000, 'evts':_evts(79), 'mu':16.500, 'force_new':False},
    {'run':284500, 'lb':18, 'starttstamp':1446540205, 'dt':0.000, 'evts':_evts(83), 'mu':17.500, 'force_new':False},
    {'run':284500, 'lb':19, 'starttstamp':1446540265, 'dt':0.000, 'evts':_evts(86), 'mu':18.500, 'force_new':False},
    {'run':284500, 'lb':20, 'starttstamp':1446540325, 'dt':0.000, 'evts':_evts(88), 'mu':19.500, 'force_new':False},
    {'run':284500, 'lb':21, 'starttstamp':1446540385, 'dt':0.000, 'evts':_evts(89), 'mu':20.500, 'force_new':False},
    {'run':284500, 'lb':22, 'starttstamp':1446540445, 'dt':0.000, 'evts':_evts(90), 'mu':21.500, 'force_new':False},
    {'run':284500, 'lb':23, 'starttstamp':1446540505, 'dt':0.000, 'evts':_evts(92), 'mu':22.500, 'force_new':False},
    {'run':284500, 'lb':24, 'starttstamp':1446540565, 'dt':0.000, 'evts':_evts(95), 'mu':23.500, 'force_new':False},
    {'run':284500, 'lb':25, 'starttstamp':1446540625, 'dt':0.000, 'evts':_evts(95), 'mu':24.500, 'force_new':False},
    {'run':284500, 'lb':26, 'starttstamp':1446540685, 'dt':0.000, 'evts':_evts(94), 'mu':25.500, 'force_new':False},
    {'run':284500, 'lb':27, 'starttstamp':1446540745, 'dt':0.000, 'evts':_evts(91), 'mu':26.500, 'force_new':False},
    {'run':284500, 'lb':28, 'starttstamp':1446540805, 'dt':0.000, 'evts':_evts(86), 'mu':27.500, 'force_new':False},
    {'run':284500, 'lb':29, 'starttstamp':1446540865, 'dt':0.000, 'evts':_evts(81), 'mu':28.500, 'force_new':False},
    {'run':284500, 'lb':30, 'starttstamp':1446540925, 'dt':0.000, 'evts':_evts(75), 'mu':29.500, 'force_new':False},
    {'run':284500, 'lb':31, 'starttstamp':1446540985, 'dt':0.000, 'evts':_evts(68), 'mu':30.500, 'force_new':False},
    {'run':284500, 'lb':32, 'starttstamp':1446541045, 'dt':0.000, 'evts':_evts(60), 'mu':31.500, 'force_new':False},
    {'run':284500, 'lb':33, 'starttstamp':1446541105, 'dt':0.000, 'evts':_evts(53), 'mu':32.500, 'force_new':False},
    {'run':284500, 'lb':34, 'starttstamp':1446541165, 'dt':0.000, 'evts':_evts(45), 'mu':33.500, 'force_new':False},
    {'run':284500, 'lb':35, 'starttstamp':1446541225, 'dt':0.000, 'evts':_evts(38), 'mu':34.500, 'force_new':False},
    {'run':284500, 'lb':36, 'starttstamp':1446541285, 'dt':0.000, 'evts':_evts(31), 'mu':35.500, 'force_new':False},
    {'run':284500, 'lb':37, 'starttstamp':1446541345, 'dt':0.000, 'evts':_evts(26), 'mu':36.500, 'force_new':False},
    {'run':284500, 'lb':38, 'starttstamp':1446541405, 'dt':0.000, 'evts':_evts(21), 'mu':37.500, 'force_new':False},
    {'run':284500, 'lb':39, 'starttstamp':1446541465, 'dt':0.000, 'evts':_evts(17), 'mu':38.500, 'force_new':False},
    {'run':284500, 'lb':40, 'starttstamp':1446541525, 'dt':0.000, 'evts':_evts(13), 'mu':39.500, 'force_new':False},
    {'run':284500, 'lb':41, 'starttstamp':1446541585, 'dt':0.000, 'evts':_evts(10), 'mu':40.500, 'force_new':False},
    {'run':284500, 'lb':42, 'starttstamp':1446541645, 'dt':0.000, 'evts':_evts(8), 'mu':41.500, 'force_new':False},
    {'run':284500, 'lb':43, 'starttstamp':1446541705, 'dt':0.000, 'evts':_evts(6), 'mu':42.500, 'force_new':False},
    {'run':284500, 'lb':44, 'starttstamp':1446541765, 'dt':0.000, 'evts':_evts(4), 'mu':43.500, 'force_new':False},
    {'run':284500, 'lb':45, 'starttstamp':1446541825, 'dt':0.000, 'evts':_evts(3), 'mu':44.500, 'force_new':False},
    {'run':284500, 'lb':46, 'starttstamp':1446541885, 'dt':0.000, 'evts':_evts(1), 'mu':45.500, 'force_new':False},
    {'run':284500, 'lb':47, 'starttstamp':1446541945, 'dt':0.000, 'evts':_evts(1), 'mu':46.500, 'force_new':False},
    {'run':284500, 'lb':48, 'starttstamp':1446542005, 'dt':0.000, 'evts':_evts(1), 'mu':47.500, 'force_new':False},
    {'run':284500, 'lb':49, 'starttstamp':1446542065, 'dt':0.000, 'evts':_evts(1), 'mu':48.500, 'force_new':False},
    {'run':284500, 'lb':50, 'starttstamp':1446542125, 'dt':0.000, 'evts':_evts(1), 'mu':49.500, 'force_new':False},
    {'run':284500, 'lb':51, 'starttstamp':1446542185, 'dt':0.000, 'evts':_evts(1), 'mu':50.500, 'force_new':False},
    {'run':284500, 'lb':52, 'starttstamp':1446542245, 'dt':0.000, 'evts':_evts(1), 'mu':51.500, 'force_new':False},
    {'run':284500, 'lb':53, 'starttstamp':1446542305, 'dt':0.000, 'evts':_evts(1), 'mu':52.500, 'force_new':False},
    {'run':284500, 'lb':54, 'starttstamp':1446542365, 'dt':0.000, 'evts':_evts(1), 'mu':53.500, 'force_new':False},
    {'run':284500, 'lb':55, 'starttstamp':1446542425, 'dt':0.000, 'evts':_evts(1), 'mu':54.500, 'force_new':False},
    {'run':284500, 'lb':56, 'starttstamp':1446542485, 'dt':0.000, 'evts':_evts(1), 'mu':55.500, 'force_new':False},
    {'run':284500, 'lb':57, 'starttstamp':1446542545, 'dt':0.000, 'evts':_evts(1), 'mu':56.500, 'force_new':False},
]
