#
# All needed variables are set below
#
if 'EvtNo' not in dir():
    EvtNo=10000000 #  Number of events to process
    
if 'SubDet' not in dir():
    SubDet = 'Barrel'
    
if 'InputDir' not in dir():
    InputDir = "/castor/cern.ch/grid/atlas/DAQ/lar/ElecCalib/2009/00102562"
    
if 'FullFileName' not in dir():
    FullFileName = "data09_calib.00102562.calibration_LArElec-Pedestal-7s-Low-HecFcal.daq.RAW._lb0000._EB-FCAL._0001.data"
    
if 'RunNumber' not in dir():
    RunNumber = int(FullFileName.strip().split('.')[1])
    
if 'Type' not in dir():
    Type = str(FullFileName.strip().split('.')[2].strip().split('-')[1])
    print Type

if 'Partition' not in dir():
    Partition = str(FullFileName.strip().split('.')[6].split('-')[1])
    print Partition

if not 'online' in dir():
    online = False

if not 'OutputDir' in dir():
    OutputDir="rootFiles/"

if not 'OutputNtupleDir' in dir():
    OutputNtupleDir = "rootFiles/"

if not 'runAccumulator' in dir(): 
    runAccumulator = False # :average mode, = True:transparent mode

if not 'RefRunNumber' in dir():
    RefRunNumber = 2095     # runnumber of reference file

if not 'LArDigitKey' in dir():
    LArMonFlags.LArDigitKey = 'HIGH'   # for LArRawChannelBuilder

if not 'LArRawChannelKey' in dir():
    LArRawChannelKey="LArRawChannels"

if not 'FullFileNameTab' in dir():
    FullFileNameTab = [ InputDir+"/"+FullFileName ]

if not 'DelayNtuple' in dir():
    DelayNtuple = False

if not 'PeakOF' in dir():
    PeakOF = False

if not 'simpleBuilder' in dir():
    simpleBuilder = False

#
# Output file naming
#
include('LArMonTools/LArMonOutputFileName.py')
