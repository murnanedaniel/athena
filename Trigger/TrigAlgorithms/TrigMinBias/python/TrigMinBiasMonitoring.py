# Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration

def SpCountMonitoring():

    from AthenaMonitoringKernel.GenericMonitoringTool import GenericMonitoringTool
    monTool = GenericMonitoringTool('MonTool')
    monTool.defineHistogram('pixCL_1', path='EXPERT', type='TH1I', title='pixCL_1',xbins=100, xmin=0, xmax=6000)
    monTool.defineHistogram('pixCLBeforeCuts', path='EXPERT', type='TH1I', title='totPixBeforeCuts', xbins = 250, xmin=0, xmax=6000)
    monTool.defineHistogram('pixCL_2', path='EXPERT', type='TH1I', title='pixCL_2',xbins=100, xmin=0, xmax=6000)
    monTool.defineHistogram('pixCLmin3', path='EXPERT', type='TH1I', title='pixCLmin3',xbins=100, xmin=0, xmax=6000)
    monTool.defineHistogram('pixCL', path='EXPERT', type='TH1I', title='pixCL', xbins = 500, xmin=0, xmax=6000)
    monTool.defineHistogram('pixCLBarrel', path='EXPERT', type='TH1I', title='pixClBarrel', xbins = 500, xmin=0, xmax=6000)
    monTool.defineHistogram('pixCLEndcapA', path='EXPERT', type='TH1I', title='pixClEndcapA', xbins = 500, xmin=0, xmax=3000)
    monTool.defineHistogram('pixCLEndcapC', path='EXPERT', type='TH1I', title='pixClEndcapC', xbins = 500, xmin=0, xmax=3000)
    monTool.defineHistogram('pixModulesOverThreshold', path='EXPERT', type='TH1I', title='Pixels ModulesOverThreshold', xbins = 100, xmin=1, xmax=200)
    monTool.defineHistogram('sctSP', path='EXPERT', type='TH1I', title='sctSP', xbins = 500, xmin=0, xmax=6000)
    monTool.defineHistogram('sctSPEndcapC', path='EXPERT', type='TH1I', title='sctSPEndcapC', xbins = 500, xmin=0, xmax=3000)
    monTool.defineHistogram('sctSPBarrel', path='EXPERT', type='TH1I', title='sctSPBarrel', xbins = 500, xmin=0, xmax=6000)
    monTool.defineHistogram('sctSPEndcapA', path='EXPERT', type='TH1I', title='sctSPEndcapA', xbins = 500, xmin=0, xmax=3000)
    monTool.defineHistogram('pixCL, sctSP', path='EXPERT', type='TH2I', title='SP ; pix SP after ToT cut; sct SP', xbins = 50, xmin=0, xmax=6000, ybins = 50, ymin=0, ymax=6000)
    monTool.defineHistogram('sctModulesOverThreshold', path='EXPERT', type='TH1I', title='SCT ModulesOverThreshold', xbins = 100, xmin=1, xmax=200)

    return monTool


def MbtsFexMonitoring():

    from AthenaMonitoringKernel.GenericMonitoringTool import GenericMonitoringTool
    monTool = GenericMonitoringTool('MonTool')
    monTool.defineHistogram('triggerEnergies', path='EXPERT', type='TH1D', title='triggerEnergies',xbins=100, xmin=-0, xmax=60999.5)
    monTool.defineHistogram('triggerEta', path='EXPERT', type='TH1D', title='triggerEta', xbins = 250, xmin=-5, xmax=5)
    monTool.defineHistogram('triggerPhi', path='EXPERT', type='TH1I', title='triggerPhi',xbins=100, xmin=-3.2, xmax=3.2)
    monTool.defineHistogram('triggerTimes', path='EXPERT', type='TH1I', title='triggerTimes',xbins=100, xmin=-0, xmax=60999.5)
    return monTool


def TrackCountMonitoring(hypoAlg):
    from AthenaMonitoringKernel.GenericMonitoringTool import GenericMonitoringTool
    monTool = GenericMonitoringTool('MonTool')
    monTool.defineHistogram('ntrks', path='EXPERT', type='TH1I', title='ntrks', xbins=200, xmin=0, xmax=200)
    for i in range(len(hypoAlg.minPt)):
        monTool.defineHistogram('countsSelection{}'.format(i),
                                path='EXPERT', type='TH1I', title='counts for min pT and max z0 cut',
                                xbins=200, xmin=0, xmax=200)
    monTool.defineHistogram( "trkPt", path='EXPERT', type='TH1I', title="Tracks pt (low pt part);p_{T} [GeV]", xbins=100, xmin=0, xmax=10)
    monTool.defineHistogram( "trkEta", path='EXPERT', type='TH1I', title="Tracks eta;#eta", xbins=50, xmin=-2.5, xmax=2.5)

    return monTool
