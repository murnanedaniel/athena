#!/usr/bin/env python
#
#  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

from AthenaConfiguration.ComponentFactory import CompFactory
from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator

if __name__=='__main__':

  import os,sys
  import argparse
  import subprocess
  from AthenaCommon import Logging
  log = Logging.logging.getLogger( 'LArSC2Ntuple' )
  
  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i','--indir', dest='indir', default="/eos/atlas/atlastier0/rucio/data_test/calibration_pulseall/00414414/data_test.00414414.calibration_pulseall.daq.RAW/", help='input files dir', type=str)
  parser.add_argument('-p','--inprefix', dest='inpref', default="data_test", help='Input filenames prefix', type=str)
  parser.add_argument('-y','--inppatt', dest='inppatt', default="lb3512", help='Input filenames pattern', type=str)
  parser.add_argument('-f','--infile', dest='infile', default="", help='Input filename (if given indir and inprefix are ignored', type=str)
  parser.add_argument('-r','--run', dest='run', default=0, help='Run number (if not given trying to judge from input file name)', type=int)
  parser.add_argument('-m','--maxev', dest='maxev', default=-1, help='Max number of events to dump', type=int)
  parser.add_argument('-x','--outlevel', dest='olevel', default=5, help='OuputLevel for dumping algo', type=int)
  parser.add_argument('-o','--outfile', dest='outfile', default="Digits.root", help='Output root filename', type=str)
  parser.add_argument('-s','--addSamples', dest='samples', default=True, help='Add Samples to output ntuple', type=bool)
  parser.add_argument('-a','--addSampBas', dest='samplesBas', default=False, help='Add ADC_BAS to output ntuple', type=bool)
  parser.add_argument('-z','--addEt', dest='Et', default=False, help='Add ET to output ntuple', type=bool)
  parser.add_argument('-g','--addEtId', dest='EtId', default=False, help='Add ET_ID to output ntuple', type=bool)
  parser.add_argument('-l','--addLatHeader', dest='lheader', default=True, help='Add LATOME Header to output ntuple', type=bool)
  parser.add_argument('-b','--addBCID', dest='bcid', default=True, help='Add BCID info to output ntuple', type=bool)
  parser.add_argument('-e','--expandId', dest='expid', default=False, help='Expand online Id to fields', type=bool)
  parser.add_argument('-n','--nsamp', dest='nsamp', default=32, help='Number of samples to dump', type=int)
  parser.add_argument('-c','--overEvNumber', dest='overEvN', default=False, help='Overwrite event number', type=bool)
  parser.add_argument('-d','--addHash', dest='ahash', default=False, help='Add hash number to output ntuple', type=bool)
  parser.add_argument('-j','--addOffline', dest='offline', default=False, help='Add offline Id to output ntuple', type=bool)
  parser.add_argument('-k','--addCalib', dest='calib', default=False, help='Add calib. info to output ntuple', type=bool)
  parser.add_argument('-t','--addGeom', dest='geom', default=False, help='Add real geom info to output ntuple', type=bool)
  parser.add_argument('-u','--addBC', dest='bc', default=False, help='Add Bad. chan info to output ntuple', type=bool)
  parser.add_argument('-w','--addROD', dest='rod', default=False, help='Add ROD energies sum to output ntuple', type=bool)

  args = parser.parse_args()
  if help in args and args.help is not None and args.help:
     parser.print_help()
     sys.exit(0)

  for _, value in args._get_kwargs():
     if value is not None:
       log.debug(value)

  #Import the flag-container that is the arguemnt to the configuration methods
  from AthenaConfiguration.AllConfigFlags import ConfigFlags
  #add SC dumping specific flags
  from LArCafJobs.LArSCDumperFlags import addSCDumpFlags
  addSCDumpFlags(ConfigFlags)

  if len(args.infile) > 0:
     ConfigFlags.Input.Files = [args.infile]
  elif len(args.inppatt) > 0:
     from LArCalibProcessing.GetInputFiles import GetInputFilesFromPattern
     ConfigFlags.Input.Files = GetInputFilesFromPattern(args.indir,args.inppatt)
  else:   
     from LArCalibProcessing.GetInputFiles import GetInputFilesFromPrefix
     ConfigFlags.Input.Files = GetInputFilesFromPrefix(args.indir,args.inpref)

  # first autoconfig
  from LArConditionsCommon.LArRunFormat import getLArDTInfoForRun
  try:
     runinfo=getLArDTInfoForRun(ConfigFlags.Input.RunNumber[0], connstring="COOLONL_LAR/CONDBR2")
  except Exception:
     log.warning("Could not get DT run info, using defaults !")
     ConfigFlags.LArSCDump.doEt=True
     if args.nsamp > 0:
        ConfigFlags.LArSCDump.nSamples=args.nsamp
     else:   
        ConfigFlags.LArSCDump.nSamples=5
     ConfigFlags.LArSCDump.nEt=1
     ConfigFlags.LArSCDump.digitsKey="SC"
     CKeys=["SC_ET"]
  else:
     CKeys=[]
     ConfigFlags.LArSCDump.digitsKey=""
     for i in range(0,len(runinfo.streamTypes())):
        if runinfo.streamTypes()[i] ==  "SelectedEnergy":
              CKeys += ["SC_ET_ID"]
              ConfigFlags.LArSCDump.doEt=True
              ConfigFlags.LArSCDump.nEt=runinfo.streamLengths()[i]
        elif runinfo.streamTypes()[i] ==  "Energy":
              CKeys += ["SC_ET"]
              ConfigFlags.LArSCDump.doEt=True
              ConfigFlags.LArSCDump.nEt=runinfo.streamLengths()[i]
        elif runinfo.streamTypes()[i] ==  "RawADC":
              ConfigFlags.LArSCDump.digitsKey="SC"
              ConfigFlags.LArSCDump.nSamples=runinfo.streamLengths()[i]
        elif runinfo.streamTypes()[i] ==  "ADC":
              CKeys += ["SC_ADC_BAS"]
              ConfigFlags.LArSCDump.nSamples=runinfo.streamLengths()[i]
     if args.nsamp < ConfigFlags.LArSCDump.nSamples:
        ConfigFlags.LArSCDump.nSamples=args.nsamp
  
  log.info("Autoconfigured: ")
  log.info("nSamples: %d nEt: %d digitsKey %s",ConfigFlags.LArSCDump.nSamples, ConfigFlags.LArSCDump.nEt, ConfigFlags.LArSCDump.digitsKey)
  log.info(CKeys)

  # now set flags according parsed options
  #if args.samples and not ("SC" in CKeys or ConfigFlags.LArSCDump.digitsKey=="SC"):
  #   log.warning("Samples asked, but they are not in RunLogger, no output !!!!")

  if args.samplesBas and "SC_ADC_BAS" not in CKeys:
     CKeys += ["SC_ADC_BAS"]
  if args.Et and "SC_ET" not in CKeys:
     CKeys += ["SC_ET"]
  if args.EtId and "SC_ET_ID" not in CKeys:
     CKeys += ["SC_ET_ID"]
  if args.lheader and "SC_LATOME_HEADER" not in CKeys:
     CKeys += ["SC_LATOME_HEADER"]

  if args.rod:
     ConfigFlags.LArSCDump.doRawChan=True  
     CKeys += ["LArRawChannels"]
     log.info("Adding ROD energies")

  # now construct the job
  ConfigFlags.LAr.doAlign=False

  ConfigFlags.lock()

  #Import the MainServices (boilerplate)
  from AthenaConfiguration.MainServicesConfig import MainServicesCfg
  from LArGeoAlgsNV.LArGMConfig import LArGMCfg

  acc = MainServicesCfg(ConfigFlags)
  acc.merge(LArGMCfg(ConfigFlags))

  if args.bc:
     # FIXME should be SC version
     from LArBadChannelTool.LArBadChannelConfig import  LArBadFebCfg, LArBadChannelCfg
     acc.merge(LArBadChannelCfg(ConfigFlags))
     acc.merge(LArBadFebCfg(ConfigFlags))

  if args.geom:
      log.warning("Adding real geometry is not working yet")
      args.geom=False
      # FIXME
      #acc.addCondAlgo(CompFactory.CaloAlignCondAlg(LArAlignmentStore="",CaloCellPositionShiftFolder=""))
      #acc.addCondAlgo(CompFactory.CaloSuperCellAlignCondAlg())
      #AthReadAlg_ExtraInputs.append(('CaloSuperCellDetDescrManager', 'ConditionStore+CaloSuperCellDetDescrManager')) 

  from LArCalibProcessing.LArSC2NtupleConfig import LArSC2NtupleCfg
  acc.merge(LArSC2NtupleCfg(ConfigFlags, AddBadChannelInfo=args.bc, AddFEBTempInfo=False, isSC=True, isFlat=False, 
                            OffId=args.offline, AddHash=args.ahash, AddCalib=args.calib, RealGeometry=args.geom, ExpandId=args.expid, # from LArCond2NtupleBase 
                            NSamples=ConfigFlags.LArSCDump.nSamples, FTlist={}, FillBCID=args.bcid, ContainerKey=ConfigFlags.LArSCDump.digitsKey,  # from LArDigits2Ntuple
                            SCContainerKeys=CKeys, OverwriteEventNumber = args.overEvN,                        # from LArSC2Ntuple
                            FillRODEnergy = ConfigFlags.LArSCDump.doRawChan,
                            OutputLevel=args.olevel
                           ))
  # ROOT file writing
  if os.path.exists(args.outfile):
      os.remove(args.outfile)
  acc.addService(CompFactory.NTupleSvc(Output = [ "FILE1 DATAFILE='"+args.outfile+"' OPT='NEW'" ]))
  acc.setAppProperty("HistogramPersistency","ROOT")

  # some logging
  log.info("Input files to be processed:")
  for f in ConfigFlags.Input.Files:
      log.info(f)
  log.info("Output file: ")
  log.info(args.outfile)


  # and run
  acc.run(args.maxev)
