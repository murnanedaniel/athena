#
#  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
#

'''
@file PixelAthMVAMonAlgCfg.py
@brief Configuration of Pixel MVA Monitoring Histograms for Run 3
'''

from PixelMonitoring.PixelAthMonitoringBase import define2DProfHist


def PixelAthMVAMonAlgCfg(helper, alg, **kwargs):
    '''
    Arguments:
         helper   -- AthMonitorCfgHelper(Old) instance
         alg      -- algorithm Configurable object returned from addAlgorithm
         kwargs   -- jo agruments
    '''

    dumpTree = kwargs.get('dumpTree', False)
    path     = '/Pixel/MVA/'

    if dumpTree:
        mvaGroup = helper.addGroup(alg, 'MVA')
        mvaGroup.defineTree('pixmvamontool_lb,status_vec,\
holes_vec,outliers_vec,meas_vec,\
clsontrkf_vec,clsontrk_vec,clsofftrk_vec,clsall_vec,\
clsontrksize_vec,clsontrkrowsize_vec,clsontrkcolsize_vec,\
clsofftrksize_vec,clsofftrkrowsize_vec,clsofftrkcolsize_vec,\
trkalpha_vec,trkchi2byndf_vec,trknpixdead_vec,trknblayerhits_vec,\
clsontrktot_vec,clsofftrktot_vec,mva_vec;\
vartree', path=path, treedef='pixmvamontool_lb/i:status_vec/vector<float>\
:h_vec/vector<float>:o_vec/vector<float>:m_vec/vector<float>\
:clsontrkf_vec/vector<float>:clsontrk_vec/vector<float>:clsofftrk_vec/vector<float>:clsall_vec/vector<float>\
:clsontrksize_vec/vector<float>:clsontrkrowsize_vec/vector<float>:clsontrkcolsize_vec/vector<float>\
:clsofftrksize_vec/vector<float>:clsofftrkrowsize_vec/vector<float>:clsofftrkcolsize_vec/vector<float>\
:trkalpha_vec/vector<float>:trkchi2byndf_vec/vector<float>\
:trknpixdead_vec/vector<float>:trknblayerhits_vec/vector<float>\
:clsontrktot_vec/vector<float>:clsofftrktot_vec/vector<float>:mva_vec/vector<float>')

    histoGroupName = 'BDTWeights'
    title = 'BDT score'
    define2DProfHist(helper, alg, histoGroupName, title, path, type='TProfile2D')

    histoGroupName = 'BDTWeightsReset'
    title = 'BDT score reset every 10 LB'
    define2DProfHist(helper, alg, 'BDTWeights', title, path, type='TProfile2D', opt='kLBNHistoryDepth=10', histname=histoGroupName)
