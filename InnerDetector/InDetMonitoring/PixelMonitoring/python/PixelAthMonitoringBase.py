#
#  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
#

'''
@file PixelAthMonitoringBase.py
@brief Helper functions for Run 3 Pixel monitoring algorithm configuration
'''
runtext = ''

NumLayersDisk = 3
NumStavesIBL  = 14
NumStavesL0   = 22
NumStavesL1   = 38
NumStavesL2   = 52

NumModulesBarrel= 13
NumModulesDisk  = 48
NumPP0sEC       = 24

btxt     = ";eta index of module"
sbtxt    = ";shifted eta index of module"
ditxt    = ";disk number"
layers   = ["ECA","ECC","BLayer","Layer1","Layer2","IBL2D","IBL3D"]
baselayers=["ECA","ECC","BLayer","Layer1","Layer2","IBL"]
xbinsl   = [    3,    3,  13,  13,  13,   32]
xminsl   = [ -0.5, -0.5,-6.5,-6.5,-6.5,-16.5]
ybinsl   = [   48,   48,  22,  38,  52,   14]
etatxt   = [ditxt,ditxt,btxt,btxt,btxt,sbtxt]
phitext  = ";phi index of module"
xbinsfel = [   24,   24, 104, 104, 104,   32]
xminsfel = [ -0.5, -0.5, -52, -52, -52,-16.5]
ybinsfel = [   96,   96,  44,  76, 104,   14]

totcuts  = [   15,   15,  15,  15,  15,    4,    4]
xbinsem  = [    3,    3,  13,  13,  13,   20]
xminsem  = [  0.5,  0.5,-6.5,-6.5,-6.5,  -10]
xbinstotz= [    3,    3,  13,  13,  13,   20]
xminstotz= [ -0.5, -0.5,-0.5,-0.5,-0.5, -0.5]
ztotbinsy= [   20,   20,  20,  20,  20,   20,   20]
ztotminsy= [ 19.5, 19.5, 7.5,19.5,19.5, -0.5, -0.5]
errbbinsy= [   17,   17,  17,  17,  17,   27,   27]
errbminsy= [    0,    0,   0,   0,   0,    0,    0]
errbbsizy= [    1,    1,   1,   1,   1,    1,    1]
errtbinsy= [    7,    7,   7,   7,   7,    7,    7]

pp0layers= ["ECA","ECC","BLayer","Layer1","Layer2","IBLA","IBLC"]
pp0xbins = [   24,   24,  22,  38,  52,   14,     14]

fei3layers= ["ECA","ECC","BLayer","Layer1","Layer2"]
fei4layers= ["IBL2D","IBL3D"]

lumitext  = ";lumi block"
lumibinsx = 3000

bcidbinsx = 3600

#labels
LayersDisk = ["Disk 1", "Disk 2", "Disk 3"] #x EC
#xBarrel
ModulesBarrel = [
    "M6C", "M5C", "M4C", "M3C", "M2C", "M1C", "M0","M1A", "M2A", "M3A", "M4A",
    "M5A", "M6A"]
#xIBL
ModulesIBL = [
    "C8_0", "C8_1", "C7_0", "C7_1", "C6_0", "C6_1", "C5_0", "C5_1", "C4_0", "C4_1", "C3_0", "C3_1", "C2_0", "C2_1",
    "C1_0", "C1_1", "A1_0", "A1_1", "A2_0", "A2_1", "A3_0", "A3_1", "A4_0", "A4_1", "A5_0", "A5_1", "A6_0", "A6_1",
    "A7_0", "A7_1", "A8_0", "A8_1"]

LabelX = [LayersDisk, LayersDisk, ModulesBarrel, ModulesBarrel, ModulesBarrel, ModulesIBL]

#same for FE
LayersDiskFE = ["Disk1_FE_8_7","Disk1_FE_9_6","Disk1_FE_10_5","Disk1_FE_11_4","Disk1_FE_12_3","Disk1_FE_13_2","Disk1_FE_14_1","Disk1_FE_15_0",
                "Disk2_FE_7_8","Disk2_FE_6_9","Disk2_FE_5_10","Disk2_FE_4_11","Disk2_FE_3_12","Disk2_FE_2_13","Disk2_FE_1_14","Disk2_FE_0_15",
                "Disk3_FE_8_7","Disk3_FE_9_6","Disk3_FE_10_5","Disk3_FE_11_4","Disk3_FE_12_3","Disk3_FE_13_2","Disk3_FE_14_1","Disk3_FE_15_0"]

ModulesBarrelFE = [
    "M6C_FE_8_7","M6C_FE_9_6","M6C_FE_10_5","M6C_FE_11_4","M6C_FE_12_3","M6C_FE_13_2","M6C_FE_14_1","M6C_FE_15_0",
    "M5C_FE_8_7","M5C_FE_9_6","M5C_FE_10_5","M5C_FE_11_4","M5C_FE_12_3","M5C_FE_13_2","M5C_FE_14_1","M5C_FE_15_0",
    "M4C_FE_8_7","M4C_FE_9_6","M4C_FE_10_5","M4C_FE_11_4","M4C_FE_12_3","M4C_FE_13_2","M4C_FE_14_1","M4C_FE_15_0",
    "M3C_FE_8_7","M3C_FE_9_6","M3C_FE_10_5","M3C_FE_11_4","M3C_FE_12_3","M3C_FE_13_2","M3C_FE_14_1","M3C_FE_15_0",
    "M2C_FE_8_7","M2C_FE_9_6","M2C_FE_10_5","M2C_FE_11_4","M2C_FE_12_3","M2C_FE_13_2","M2C_FE_14_1","M2C_FE_15_0",
    "M1C_FE_8_7","M1C_FE_9_6","M1C_FE_10_5","M1C_FE_11_4","M1C_FE_12_3","M1C_FE_13_2","M1C_FE_14_1","M1C_FE_15_0",
    "M0_FE_8_7","M0_FE_9_6","M0_FE_10_5","M0_FE_11_4","M0_FE_12_3","M0_FE_13_2","M0_FE_14_1","M0_FE_15_0",
    "M1A_FE_8_7","M1A_FE_9_6","M1A_FE_10_5","M1A_FE_11_4","M1A_FE_12_3","M1A_FE_13_2","M1A_FE_14_1","M1A_FE_15_0",
    "M2A_FE_8_7","M2A_FE_9_6","M2A_FE_10_5","M2A_FE_11_4","M2A_FE_12_3","M2A_FE_13_2","M2A_FE_14_1","M2A_FE_15_0",
    "M3A_FE_8_7","M3A_FE_9_6","M3A_FE_10_5","M3A_FE_11_4","M3A_FE_12_3","M3A_FE_13_2","M3A_FE_14_1","M3A_FE_15_0",
    "M4A_FE_8_7","M4A_FE_9_6","M4A_FE_10_5","M4A_FE_11_4","M4A_FE_12_3","M4A_FE_13_2","M4A_FE_14_1","M4A_FE_15_0",
    "M5A_FE_8_7","M5A_FE_9_6","M5A_FE_10_5","M5A_FE_11_4","M5A_FE_12_3","M5A_FE_13_2","M5A_FE_14_1","M5A_FE_15_0",
    "M6A_FE_8_7","M6A_FE_9_6","M6A_FE_10_5","M6A_FE_11_4","M6A_FE_12_3","M6A_FE_13_2","M6A_FE_14_1","M6A_FE_15_0"]

LabelFEX = [LayersDiskFE, LayersDiskFE, ModulesBarrelFE, ModulesBarrelFE, ModulesBarrelFE, ModulesIBL]

ModulesECA = [
    "B01_S2_M1", "B01_S2_M6", "B01_S2_M2", "B01_S2_M5", "B01_S2_M3",
    "B01_S2_M4", "B02_S1_M1", "B02_S1_M6", "B02_S1_M2", "B02_S1_M5",
    "B02_S1_M3", "B02_S1_M4", "B02_S2_M1", "B02_S2_M6", "B02_S2_M2",
    "B02_S2_M5", "B02_S2_M3", "B02_S2_M4", "B03_S1_M1", "B03_S1_M6",
    "B03_S1_M2", "B03_S1_M5", "B03_S1_M3", "B03_S1_M4", "B03_S2_M1",
    "B03_S2_M6", "B03_S2_M2", "B03_S2_M5", "B03_S2_M3", "B03_S2_M4",
    "B04_S1_M1", "B04_S1_M6", "B04_S1_M2", "B04_S1_M5", "B04_S1_M3",
    "B04_S1_M4", "B04_S2_M1", "B04_S2_M6", "B04_S2_M2", "B04_S2_M5",
    "B04_S2_M3", "B04_S2_M4", "B01_S1_M1", "B01_S1_M6", "B01_S1_M2",
    "B01_S1_M5", "B01_S1_M3", "B01_S1_M4"]
ModulesECC = [
    "B01_S2_M4", "B01_S2_M3", "B01_S2_M5", "B01_S2_M2", "B01_S2_M6",
    "B01_S2_M1", "B02_S1_M4", "B02_S1_M3", "B02_S1_M5", "B02_S1_M2",
    "B02_S1_M6", "B02_S1_M1", "B02_S2_M4", "B02_S2_M3", "B02_S2_M5",
    "B02_S2_M2", "B02_S2_M6", "B02_S2_M1", "B03_S1_M4", "B03_S1_M3",
    "B03_S1_M5", "B03_S1_M2", "B03_S1_M6", "B03_S1_M1", "B03_S2_M4",
    "B03_S2_M3", "B03_S2_M5", "B03_S2_M2", "B03_S2_M6", "B03_S2_M1",
    "B04_S1_M4", "B04_S1_M3", "B04_S1_M5", "B04_S1_M2", "B04_S1_M6",
    "B04_S1_M1", "B04_S2_M4", "B04_S2_M3", "B04_S2_M5", "B04_S2_M2",
    "B04_S2_M6", "B04_S2_M1", "B01_S1_M4", "B01_S1_M3", "B01_S1_M5",
    "B01_S1_M2", "B01_S1_M6", "B01_S1_M1"]


#yBarrel
StavesL0 = [
    "B11_S2", "B01_S1", "B01_S2", "B02_S1", "B02_S2", "B03_S1", "B03_S2",
    "B04_S1", "B04_S2", "B05_S1", "B05_S2", "B06_S1", "B06_S2", "B07_S1",
    "B07_S2", "B08_S1", "B08_S2", "B09_S1", "B09_S2", "B10_S1", "B10_S2",
    "B11_S1"]


StavesL1 = [
    "B01_S1", "B01_S2", "B02_S1", "B02_S2", "B03_S1", "B03_S2", "B04_S1",
    "B04_S2", "B05_S1", "B05_S2", "B06_S1", "B06_S2", "B07_S1", "B07_S2",
    "B08_S1", "B08_S2", "B09_S1", "B09_S2", "B10_S1", "B10_S2", "B11_S1",
    "B11_S2", "B12_S1", "B12_S2", "B13_S1", "B13_S2", "B14_S1", "B14_S2",
    "B15_S1", "B15_S2", "B16_S1", "B16_S2", "B17_S1", "B17_S2", "B18_S1",
    "B18_S2", "B19_S1", "B19_S2"]

StavesL2 = [
    "B01_S2", "B02_S1", "B02_S2", "B03_S1", "B03_S2", "B04_S1", "B04_S2",
    "B05_S1", "B05_S2", "B06_S1", "B06_S2", "B07_S1", "B07_S2", "B08_S1",
    "B08_S2", "B09_S1", "B09_S2", "B10_S1", "B10_S2", "B11_S1", "B11_S2",
    "B12_S1", "B12_S2", "B13_S1", "B13_S2", "B14_S1", "B14_S2", "B15_S1",
    "B15_S2", "B16_S1", "B16_S2", "B17_S1", "B17_S2", "B18_S1", "B18_S2",
    "B19_S1", "B19_S2", "B20_S1", "B20_S2", "B21_S1", "B21_S2", "B22_S1",
    "B22_S2", "B23_S1", "B23_S2", "B24_S1", "B24_S2", "B25_S1", "B25_S2",
    "B26_S1", "B26_S2", "B01_S1"]
#yIBL
StavesIBL = [
    "S01", "S02", "S03", "S04", "S05", "S06", "S07", "S08", "S09", "S10", "S11",
    "S12", "S13","S14"]

LabelY = [ModulesECA, ModulesECC, StavesL0, StavesL1, StavesL2, StavesIBL]

#same for per-FE
ModulesFEECA = []
for i in ModulesECA: 
    ModulesFEECA.extend([i+'_FE#1', i+'_FE#2'])
ModulesFEECC = []
for i in ModulesECC: 
    ModulesFEECC.extend([i+'_FE#1', i+'_FE#2'])
StavesFEL0 = []
for i in StavesL0:
    StavesFEL0.extend([i+'_FE#1', i+'_FE#2'])
StavesFEL1 = []
for i in StavesL1:
    StavesFEL1.extend([i+'_FE#1', i+'_FE#2'])
StavesFEL2 = []
for i in StavesL2:
    StavesFEL2.extend([i+'_FE#1', i+'_FE#2'])
LabelFEY = [ModulesFEECA, ModulesFEECC, StavesFEL0, StavesFEL1, StavesFEL2, StavesIBL]

#PP0
PP0sEC = [
    "D1_B01_S2", "D1_B02_S1", "D1_B02_S2", "D1_B03_S1", "D1_B03_S2", "D1_B04_S1", "D1_B04_S2", "D1_B01_S1",
    "D2_B01_S2", "D2_B02_S1", "D2_B02_S2", "D2_B03_S1", "D2_B03_S2", "D2_B04_S1", "D2_B04_S2", "D2_B01_S1",
    "D3_B01_S2", "D3_B02_S1", "D3_B02_S2", "D3_B03_S1", "D3_B03_S2", "D3_B04_S1", "D3_B04_S2", "D3_B01_S1"]

PP0LabelX = [PP0sEC, PP0sEC, StavesL0, StavesL1, StavesL2, StavesIBL, StavesIBL]

#Errors
ErrStateLabelsFEI3 = [
      ("RODTimeout", "ROD Formatter Timeout"),
      ("RODSyncBCIDErrors", "ROD BCID Sync"),
      ("RODSyncLVL1IDErrors", "ROD LVL1ID Sync"),
      ("OpticalErrors", "Preamble/Header"),
      ("ModSyncLVL1IDErrors", "FE/MCC LVL1ID Sync"),
      ("ModSyncBCID2Errors", "FE/MCC BCID2 Sync"),
      ("ModSyncBCID1Errors", "FE/MCC BCID1 Sync"),
      ("ModTruncEoEOverflowErrors", "FE/MCC EoE Overflow Trunc"),
      ("ModTruncHitOverflowErrors", "FE/MCC Hit Overflow Trunc"),
      ("FEWarning", "FE Warning"),
      ("SEUHitParity", "SEU Hit Parity"),
      ("SEURegisterParity", "SEU Register Parity"),
      ("SEUHamming", "SEU Hamming Code"),
      ("ModTruncEOCErrors", "FE/MCC EoC Trunc"),
      ("RODTrailerBitErrors", "ROD Trailer Bit"),
      ("RODTruncHTLimitErrors", "ROD H/T Limit Trunc"),
      ("RODTruncRODOFErrors", "ROD Overflow Trunc")
]
ErrStateLabelsFEI4 = [
      ("RODTimeoutErrors", "ROD Timeout"),
      ("RODBCIDErrors", "ROD BCID synchronization"),
      ("RODLVL1IDErrors", "ROD LVL1ID synchronization"),
      ("RODPreambleErrors", "ROD Preamble"),
      ("RODTrailerErrors", "ROD Trailer"),
      ("RODRowColumnErrors", "ROD row-column"),
      ("RODMaskedLink", "ROD Masked link"),
      ("RODLimitErrors", "ROD Header Trailer limit"),
      ("SR0BCIDCounterErrors", "SR0 BCID counter"),
      ("SR1HammingCode0Errors", "SR1 Hamming code in word 0"),
      ("SR2HammingCode1Errors", "SR2 Hamming code in word 1"),
      ("SR3HammingCode2Errors", "SR3 Hamming code in word 2"),
      ("SR4L1InCounterErrors", "SR4 L1 in counter"),
      ("SR5L1RequestCounterErrors", "SR5 L1 request counter"),
      ("SR6L1RegisterErrors", "SR6 L1 register"),
      ("SR7L1TriggerIDErrors", "SR7 L1 trigger ID"),
      ("SR8FEReadoutProcessErrors", "SR8 FE readout process"),
      ("SR15SkipppedTrigCountErrors", "SR15 Skipped trigger counter"),
      ("SR16TruncatedEventFlagErrors", "SR16 Truncated event"),
      ("SR24TripleRedundantErrorsCNFGMEM", "SR24 Triple redundant CNFGMEM"),
      ("SR25WriteRegDataErrors", "SR25 Write register data"),
      ("SR26AddressErrors", "SR26 Address"),
      ("SR27OtherCMDDecoderErrors", "SR27 CMD decoder"),
      ("SR28CMDDecoderBitflipErrors", "SR28 CMD decoder bit flip"),
      ("SR29TripleRedundantErrorsCMD", "SR29 Triple redundant CMD"),
      ("SR30DataBusAddressErrors", "SR30 Data bus address"),
      ("SR31TripleRedundantErrorsEFUSE", "SR31 Triple redundant EFUSE")
]

ErrCatRODModLabels = [
      ("SyncErrorsMod", "FE/MCC Sync Errors"),
      ("SyncErrorsROD", "ROD Sync Errors"),
      ("TruncErrorsMod", "FE/MCC Trunc Errors"),
      ("TruncErrorsROD", "ROD Trunc Errors"),
      ("OpticalErrorsRODMod", "Preamble/Header Errors"),
      ("SEUErrorsRODMod", "SEU Errors"),
      ("TimeoutErrorsRODMod", "Timeout Errors")
]

ErrCatRODModLabelsNorm = [
      "SyncErrorsModFracPerEvent",
      "SyncErrorsRODFracPerEvent",
      "TruncErrorsModFracPerEvent",
      "TruncErrorsRODFracPerEvent"
]

ErrCatLabels = [
      ("SyncErrors", "Sync Errors (FE/MCC & ROD)"),
      ("TruncErrors", "Trunc Errors (FE/MCC & ROD)"),
      ("OpticalErrors", "Preamble/Header Errors"),
      ("SEUErrors", "SEU Errors"),
      ("TimeoutErrors", "Timeout Errors")
]

ErrCatLabelsNorm = [
      "SyncErrorsFracPerEvent", 
      "TruncationErrorsFracPerEvent", 
      "OpticalErrorsFracPerEvent", 
      "SEUErrorsFracPerEvent", 
      "TimeoutErrorsFracPerEvent"
]

ReadingDataErrLabels = [
      "Invalid container",
      "Invalid collection", 
      "Container empty"
]

layergroups = {}
def getLayerGroup(helper, alg, layer):
    if alg not in layergroups:
        layergroups[alg] = {}
    if layer not in layergroups[alg]:
        layergroups[alg][layer] = helper.addGroup(alg, layer)
    return layergroups[alg][layer]


def define2DProfHist(helper, alg, name, title, path, type='TProfile2D', doWeight=False, lifecycle='run', zmin=None, zmax=None, opt='', histname=None, onlylayers=baselayers):
    '''
    This function configures 2D (Profile) histograms (or maps) for Pixel base (i.e. w/o split in 2D/3D in IBL) layers.

    Arguments:
         helper     -- AthMonitorCfgHelper(Old) instance
         alg        -- algorithm - Configurable object returned from addAlgorithm
         name       -- Name of histogram (Name = name_layer)
         title      -- Title of histogram (Title = title +' '+layer)
         path       -- Path in output file for histogram
         type       -- Type of histogram (TH2D, TProfile2D)
         lifecycle  -- global life duration of histograms (run, lowStat [i.e. 20 LB], lumiblock) - APPLIES to MonGroup only
         zmin(zmax) -- fix the displayed range - simply chopping the range!!!
         opt        -- history depth of a histogram e.g. 'kLBNHistoryDepth=10'
         histname   -- another way of naming the histogram(s), useful when multiple histograms are filled from exactly the same variables, but in a different way
         onlylayers -- sublist of Pixel layers
    '''
    assert(set(onlylayers).issubset(baselayers))
    if histname is None:
        histname = name
    for i, layer in enumerate(baselayers):
        if layer not in onlylayers: 
            continue
        fulltitle   = title + ', {0}'.format(layer) + runtext + etatxt[i] + phitext
        layerGroup = getLayerGroup(helper, alg, layer)

        fullvarstring = '{0}_{1},{0}_{2}'.format(name, 'em', 'pm')
        weightvar = ''
        if 'Profile' in type: 
            fullvarstring += ',{0}_{1}'.format(name, 'val')
        elif doWeight:
            weightvar = '{0}_{1}'.format(name, 'val') # re-use the same variable for TH2 for now.
            
        fullvarstring += ';' + histname + '_{0}'.format(layer)
        layerGroup.defineHistogram(fullvarstring, 
                                    type=type, path=path, title=fulltitle, weight=weightvar,
                                    xbins=xbinsl[i], xmin=xminsl[i], xmax=xminsl[i]+xbinsl[i], 
                                    ybins=ybinsl[i], ymin=-0.5, ymax=-0.5+ybinsl[i],
                                    zmin=zmin, zmax=zmax,
                                    duration=lifecycle,
                                    opt=opt, xlabels=LabelX[i], ylabels=LabelY[i])

def define2DProfPerFEHist(helper, alg, name, title, path, type='TProfile2D', doWeight=False, lifecycle='run', zmin=None, zmax=None, opt='', histname=None, onlylayers=baselayers):
    '''
    This function configures 2D (Profile) histograms (or maps) for Pixel base (i.e. w/o split in 2D/3D in IBL) layers per FE.

    Arguments:
         helper     -- AthMonitorCfgHelper(Old) instance
         alg        -- algorithm - Configurable object returned from addAlgorithm
         name       -- Name of histogram (Name = name_layer)
         title      -- Title of histogram (Title = title +' '+layer)
         path       -- Path in output file for histogram
         type       -- Type of histogram (TH2D, TProfile2D)
         lifecycle  -- global life duration of histograms (run, lowStat [i.e. 20 LB], lumiblock) - APPLIES to MonGroup only
         zmin(zmax) -- fix the displayed range - simply chopping the range!!!
         opt        -- history depth of a histogram e.g. 'kLBNHistoryDepth=10'
         histname   -- another way of naming the histogram(s), useful when multiple histograms are filled from exactly the same variables, but in a different way
         onlylayers -- sublist of Pixel layers
    '''
    assert(set(onlylayers).issubset(baselayers))
    if histname is None:
        histname = name
    for i, layer in enumerate(baselayers):
        if layer not in onlylayers: 
            continue
        fulltitle   = title + ', {0}'.format(layer) + runtext + etatxt[i] + phitext
        layerGroup = getLayerGroup(helper, alg, layer)

        fullvarstring = '{0}_{1},{0}_{2}'.format(name, 'em', 'pm')
        weightvar = ''
        if 'Profile' in type: 
            fullvarstring += ',{0}_{1}'.format(name, 'val')
        elif doWeight:
            weightvar = '{0}_{1}'.format(name, 'val') # re-use the same variable for TH2 for now.
            
        fullvarstring += ';' + histname + '_{0}'.format(layer)
        layerGroup.defineHistogram(fullvarstring, 
                                    type=type, path=path, title=fulltitle, weight=weightvar,
                                    xbins=xbinsfel[i], xmin=xminsfel[i], xmax=xminsfel[i]+xbinsfel[i], 
                                    ybins=ybinsfel[i], ymin=-0.5, ymax=-0.5+ybinsfel[i],
                                    zmin=zmin, zmax=zmax,
                                    duration=lifecycle,
                                    opt=opt, xlabels=LabelFEX[i], ylabels=LabelFEY[i])


def definePP0Histos(helper, alg, name, title, path, opt=''):
    '''
    This function configures TProfile histograms for Pixel layers per pp0 or ROD.

    Arguments:
         helper  -- AthMonitorCfgHelper(Old) instance
         alg     -- algorithm Configurable object returned from addAlgorithm
         name    -- Name of histogram (Name = name_layer)
         title   -- Title of histogram (Title = title +' '+layer)
         path    -- Path in ouput file for histogram
         opt     -- history depth of a histogram e.g. 'kLBNHistoryDepth=10'
    '''

    for i, layer in enumerate(pp0layers):
        if i<2:
            xatxt = ';sector'
        else :
            xatxt = ';stave'
        yatxt = ';average # clusters per event per '
        if i>4:
            yatxt += 'FE'
        else :
            yatxt += 'module'
        fulltitle   = title + ', {0}'.format(layer) + runtext + xatxt + yatxt
        groupname   = name  + '_{0}'.format(layer)
        layerGroup = getLayerGroup(helper, alg, layer)
        fullvarstring = '{0}_{1},{0}_{2}'.format(name, 'pospp0x', 'val')
        fullvarstring += ';' + groupname
        layerGroup.defineHistogram(fullvarstring, 
                                    type='TProfile', path=path, title=fulltitle,
                                    xbins=pp0xbins[i], xmin=-0.5, xmax=-0.5+pp0xbins[i],
                                    opt=opt, xlabels=PP0LabelX[i])



def define1DProfLumiLayers(helper, alg, name, title, path, yaxistext, type='TProfile', opt='', histname=None, onlylayers=layers):
    '''
    This function configures 1D (Profile) vs lumi histograms for Pixel layers.

    Arguments:
         helper    -- AthMonitorCfgHelper(Old) instance
         alg       -- algorithm Configurable object returned from addAlgorithm
         name      -- Name of histogram (Name = name_layer)
         title     -- Title of histogram (Title = title +' '+layer)
         path      -- Path in ouput file for histogram
         yaxistext -- Text on the y-axis
         type      -- Type of histogram (TH1D, TProfile)
         opt       -- history depth of a histogram e.g. 'kLBNHistoryDepth=10' or 'kLive=100'
         histname  -- another way of naming the histogram(s), useful when multiple histograms are filled from exactly the same variables, but in a different way 
         onlylayers -- sublist of layers
    '''
    assert(set(onlylayers).issubset(layers))
    if histname is None:
        histname = name
    if "kLive" in opt:
        xbins = int(opt.split('=')[1])
    else:
        xbins = lumibinsx
    for layer in layers:
        if layer not in onlylayers: 
            continue
        fulltitle   = title + ', {0}'.format(layer) + runtext + lumitext + yaxistext
        fullvarstring = '{0}_{1}'.format(name,'lb')
        if 'Profile' in type: fullvarstring += ',{0}_{1}'.format(name, 'val')
        fullvarstring += ';' + histname + '_{0}'.format(layer)
        layerGroup = getLayerGroup(helper, alg, layer)
        layerGroup.defineHistogram(fullvarstring, 
                                    type=type, path=path, title=fulltitle,
                                    xbins=xbins, xmin=-0.5, xmax=-0.5+xbins, opt=opt)

def defineMapVsLumiLayers(helper, alg, name, title, path, xaxistext, yaxistext, ybins, ymins, binsizes=[1.0], ylabels=None, opt='', type='TH2F', histname=None, onlylayers=layers):
    '''
    This function configures 2D histograms vs lumi for Pixel base (i.e. w/o split in 2D/3D in IBL) layers.

    Arguments:
         helper  -- AthMonitorCfgHelper(Old) instance
         alg     -- algorithm Configurable object returned from addAlgorithm
         name    -- Name of histogram (Name = name_layer)
         title   -- Title of histogram (Title = title +' '+layer)
         path    -- Path in ouput file for histogram
         ybins, ymin, ymax, yaxistext
                 -- Configure Y-axis
         opt     -- history depth of a histogram e.g. 'kLBNHistoryDepth=10' or 'kLive=100'
         type    -- Type of histogram (TH2I, TH2F, TProfile2D)
         histname-- alternative root name of the histogram (to be filled with the same variables defined by 'name' above)
         onlylayers -- sublist of Pixel layers
    '''
    assert(set(onlylayers).issubset(layers))
    if histname is None:
        histname = name
    if "kLive" in opt:
        xbins = int(opt.split('=')[1])
    else:
        xbins = lumibinsx
    for idx,layer in enumerate(layers):
        if layer not in onlylayers: 
            continue
        fulltitle   = title + ', {0}'.format(layer) + runtext + lumitext + yaxistext
        layerGroup = getLayerGroup(helper, alg, layer)
        fullvarstring = '{0}_{1}'.format(name,'lb')
        if 'Profile' in type: fullvarstring += ',{0}_{1}'.format(name, 'cat')
        fullvarstring += ',{0}_{1}'.format(name, 'val')
        fullvarstring += ';' + histname  + '_{0}'.format(layer)
        if ( len(ybins)==1 and len(ymins)==1 and len(binsizes)==1):
            layerGroup.defineHistogram(fullvarstring, 
                                       type=type, path=path, title=fulltitle,
                                       xbins=xbins, xmin=-0.5, xmax=-0.5+xbins,
                                       ybins=ybins[0], ymin=ymins[0], ymax=ymins[0]+binsizes[0]*ybins[0], opt=opt)
        elif (len(ybins)==len(layers) and len(ymins)==len(layers) and len(binsizes)==len(layers) and len(ylabels)==len(layers)):
            layerGroup.defineHistogram(fullvarstring, 
                                       type=type, path=path, title=fulltitle,
                                       xbins=xbins, xmin=-0.5, xmax=-0.5+xbins,
                                       ybins=ybins[idx], ymin=ymins[idx], ymax=ymins[idx]+binsizes[idx]*ybins[idx], ylabels=ylabels[idx], opt=opt)
        elif (len(ybins)==len(layers) and len(ymins)==len(layers) and len(binsizes)==1 and ylabels is None):
            layerGroup.defineHistogram(fullvarstring, 
                                       type=type, path=path, title=fulltitle,
                                       xbins=xbins, xmin=-0.5, xmax=-0.5+xbins,
                                       ybins=ybins[idx], ymin=ymins[idx], ymax=ymins[idx]+ybins[idx], opt=opt)


def define1DLayers(helper, alg, name, title, path, xaxistext, yaxistext, xbins, xmins, binsizes=[1.0], type='TH1F', histname=None, onlylayers=layers):
    '''
    This function configures 1D (Profile) histograms for Pixel layers.

    Arguments:
         helper     -- AthMonitorCfgHelper(Old) instance
         alg        -- algorithm Configurable object returned from addAlgorithm
         name       -- Root name of variables and histogram (Name = name+'_'+layer)
         title      -- Title of histogram (Title = title +' '+layer)
         path       -- Path in ouput file for histogram
         xaxistext  -- X-axis title
         yaxistext  -- Y-axis title
         type       -- Type of histogram
         histname   -- alternative root name of the histogram (to be filled with the same variables defined by 'name' above)
         onlylayers -- sublist of layers
    '''
    assert(set(onlylayers).issubset(layers))
    if histname is None:
        histname = name
    for idx,layer in enumerate(layers):
        if layer not in onlylayers: 
            continue
        fulltitle   = title + ', {0}'.format(layer) + runtext + xaxistext + yaxistext
        layerGroup = getLayerGroup(helper, alg, layer)
        fullvarstring = '{0}_{1}'.format(name,'val')
        fullvarstring += ';' + histname  + '_{0}'.format(layer)
        if ( len(xbins)==1 and len(xmins)==1 and len(binsizes)==1):
            layerGroup.defineHistogram(fullvarstring, 
                                        type=type, path=path, title=fulltitle,
                                        xbins=xbins[0], xmin=xmins[0], xmax=xmins[0]+binsizes[0]*xbins[0])
        elif (len(xbins)==len(layers) and len(xmins)==len(layers) and len(binsizes)==len(layers)):
            layerGroup.defineHistogram(fullvarstring, 
                                        type=type, path=path, title=fulltitle,
                                        xbins=xbins[idx], xmin=xmins[idx], xmax=xmins[idx]+binsizes[idx]*xbins[idx])

def addOnTrackTxt(name, ontrack, wSpace=False):
    if ontrack:
        if wSpace:
            name += ' OnTrack'
        else:
            name += 'OnTrack'
    return name

def addOnTrackToPath(name, ontrack):
    if ontrack:
        name += 'OnTrack/'
    else:
        name += '/'
    return name

def fullDressTitle(name, ontrack, xaxistext, yaxistext):
    if ontrack:
        name += ' OnTrack'
    name += runtext + xaxistext + yaxistext
    return name
