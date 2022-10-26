# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

from __future__ import print_function

from AthenaCommon.Logging import logging
logOverlayReadMetadata = logging.getLogger( 'OverlayReadMetadata' )

def hitColls2SimulatedDetectors(inputlist):
    """Build a dictionary from the list of containers in the metadata"""
    simulatedDetectors = []
    simulatedDictionary = {'PixelHits': 'pixel', 'SCT_Hits': 'SCT', 'TRTUncompressedHits': 'TRT',
                           'BCMHits': 'BCM', 'LucidSimHitsVector': 'Lucid', 'LArHitEMB': 'LAr',
                           'LArHitEMEC': 'LAr', 'LArHitFCAL': 'LAr', 'LArHitHEC': 'LAr',
                           'MBTSHits': 'Tile', 'TileHitVec': 'Tile', 'MDT_Hits': 'MDT',
                           'MicromegasSensitiveDetector': 'MM',  # TODO: remove after MC21a
                           'sTGCSensitiveDetector': 'sTGC',  # TODO: remove after MC21a
                           'CSC_Hits': 'CSC', 'TGC_Hits': 'TGC', 'RPC_Hits': 'RPC',
                           'sTGC_Hits': 'sTGC', 'MM_Hits': 'MM',
                           'TruthEvent': 'Truth'} #'': 'ALFA', '': 'ZDC',
    for entry in inputlist:
        if entry[1] in simulatedDictionary.keys():
            if simulatedDictionary[entry[1]] not in simulatedDetectors:
                simulatedDetectors += [simulatedDictionary[entry[1]]]
    return simulatedDetectors

def checkLegacyEventInfo(inputlist):
    """ Check for legacy EventInfo """
    present = False
    for entry in inputlist:
        if entry[0] != 'EventInfo':
            continue

        print(entry)
        present = True

    return present

def checkLegacyNSWContainers(inputlist):
    """ Check for legacy NSW containers """
    present = False
    for entry in inputlist:
        if not (entry[0] == 'MMSimHitCollection' and entry[1] == 'MicromegasSensitiveDetector') \
            and not (entry[0] == 'sTGCSimHitCollection' and entry[1] == 'sTGCSensitiveDetector'):
            continue

        present = True
        break

    return present

def checkTileCalibrationHitFormat(inputlist):
    """Check the Tile CaloCalibrationHit format"""
    oldnames = ["TileCalibrationCellHitCnt","TileCalibrationDMHitCnt"]
    newnames = ["TileCalibHitActiveCell","TileCalibHitInactiveCell","TileCalibHitDeadMaterial"]
    nold = 0
    nnew = 0
    for entry in inputlist:
        if entry[1] in oldnames:
            logOverlayReadMetadata.debug("found %s in oldnames", entry[1])
            nold+=1
        if entry[1] in newnames:
            logOverlayReadMetadata.debug("found %s in newnames", entry[1])
            nnew+=1
    if nold > 0 and nnew > 0:
        raise SystemExit("Input file contains both old and new style TileCaloCalibrationHit containers, please check your g4sim job.")
    elif nold > 0:
        from Digitization.DigitizationFlags import digitizationFlags
        digitizationFlags.experimentalDigi += ['OldTileCalibHitContainers']
        logOverlayReadMetadata.info("Input file uses old TileCalibHitContainers names: %s", oldnames)
    elif nnew > 0:
        logOverlayReadMetadata.info("Input file uses new TileCalibHitContainers names: %s", newnames)
    return

def listOptionalContainers(inputlist):
    """Generate a list of optional containers"""
    supported = ['TrackRecordCollection', 'CaloCalibrationHitContainer', 'HijingEventParams',
                 "xAOD::TruthParticleContainer", "xAOD::JetContainer"]

    containers = {}
    for entry in inputlist:
        if entry[0] not in supported:
            continue

        if entry[0] not in containers:
            containers[entry[0]] = set()

        containers[entry[0]].add(entry[1])

    return containers

## Helper functions
def skipCheck(key):
    """This check is not required"""
    from Digitization.DigitizationFlags import digitizationFlags
    if key in digitizationFlags.overrideMetadata.get_Value():
        return True
    return False

def skipPileUpCheck(key, pileuptype):
    """This check is not required"""
    if skipCheck(key):
        return True
    pileupkey='%s_%s' %(pileuptype,key)
    return skipCheck(pileupkey)

def doMC_channel_number(f,pileUpType):
    print ("doMC_channel_number for %s", pileUpType)
    if "mc_channel_number" in f.infos.keys():
        params = dict()
        from Digitization.DigitizationFlags import digitizationFlags
        if digitizationFlags.pileupDSID.statusOn:
            params = digitizationFlags.pileupDSID.get_Value()
        print ("MC channel number from AthFile %s", f.infos["mc_channel_number"])
        params[pileUpType]= f.infos["mc_channel_number"]
        digitizationFlags.pileupDSID = params
        del params
    return

def doSpecialConfiguration(f):
    #safety checks before trying to access metadata
    if "tag_info" in f.infos.keys():
        if "specialConfiguration" in f.infos["tag_info"]:
            item = f.infos["tag_info"]["specialConfiguration"]
            logOverlayReadMetadata.info("specialConfiguration directive: %s", item)
            spcitems = item.split(";")
            preIncludes=[]
            params = {}
            from Digitization.DigitizationFlags import digitizationFlags
            if digitizationFlags.specialConfiguration.statusOn:
                logOverlayReadMetadata.info("some spcialConfiguration metadata already exists: %s", str(params))
                params = digitizationFlags.specialConfiguration.get_Value()
            for spcitem in spcitems:
                ## Ignore empty (e.g. from consecutive or trailing semicolons) or "NONE" substrings
                if not spcitem or spcitem == "NONE":
                    continue
                ## If not in key=value format, treat as v, with k="preInclude"
                if "=" not in spcitem:
                    spcitem = "preInclude=" + spcitem
                ## Handle k=v directives
                k, v = spcitem.split("=")
                logOverlayReadMetadata.info("specialConfiguration metadata item: %s => %s", k, v)
                ## Store preIncludes for including later.
                if k == "preInclude":
                    incfiles = v.split(",")
                    preIncludes+=incfiles
                else:
                    params[k] = v
            digitizationFlags.specialConfiguration = params
            ## Now that we've looked at and stored all the evgen metadata keys, we should do any requested preIncludes
            from AthenaCommon.Include import include
            for incfile in preIncludes:
                logOverlayReadMetadata.info("Including %s as instructed by specialConfiguration metadata", incfile)
                include(incfile)
            del preIncludes
            del params

def buildDict(inputtype, inputfile):
    """Build a dictionary of KEY:VALUE pairs"""
    import re
    import PyUtils.AthFile as af
    try:
        f = af.fopen(inputfile)
    except AssertionError:
        logOverlayReadMetadata.error("Failed to open input file: %s", inputfile)
        return dict(),False
    #check evt_type of input file
    if 'evt_type' in f.infos.keys():
        if not re.match(str(f.infos['evt_type'][0]), 'IS_SIMULATION') :
            logOverlayReadMetadata.error('This input file has incorrect evt_type: %s',str(f.infos['evt_type']))
            logOverlayReadMetadata.info('Please make sure you have set input file metadata correctly.')
            logOverlayReadMetadata.info('Consider using the job transforms for earlier steps if you aren\'t already.')
            #then exit gracefully
            raise SystemExit("Input file evt_type is incorrect, please check your g4sim and evgen jobs.")
    else :
        logOverlayReadMetadata.warning('Could not find \'evt_type\' key in athfile.infos. Unable to that check evt_type is correct.')

    ## Not part of building the metadata dictionary, but this is the
    ## most convenient time to access this information.
    doSpecialConfiguration(f)
    #doMC_channel_number(f,inputtype) #FIXME commented out for now until mc_channel_number is filled properly by AthFile.

    metadatadict = dict()
    digimetadatadict = dict()
    taginfometadata = dict()
    #safety checks before trying to access metadata
    if 'metadata' in f.infos.keys():
        ##if '/TagInfo' in f.infos['metadata'].keys():
        ##    taginfometadata=f.infos['metadata']['/TagInfo']
        ##    assert taginfometadata['beam_energy'] is not None
        ##    print ("beamEnergy=%s"%taginfometadata['beam_energy'])
        if '/Simulation/Parameters' in f.infos['metadata'].keys():
            metadatadict = f.infos['metadata']['/Simulation/Parameters']
            if isinstance(metadatadict, list):
                logOverlayReadMetadata.warning("%s inputfile: %s contained %s sets of Simulation Metadata. Using the final set in the list.",inputtype,inputfile,len(metadatadict))
                metadatadict=metadatadict[-1]
        ##Get IOVDbGlobalTag
        if 'IOVDbGlobalTag' not in metadatadict.keys():
            try:
                assert f.fileinfos['metadata']['/TagInfo']['IOVDbGlobalTag'] is not None
                metadatadict['IOVDbGlobalTag'] = f.fileinfos['metadata']['/TagInfo']['IOVDbGlobalTag']
            except Exception:
                try:
                    assert f.fileinfos['conditions_tag'] is not None
                    metadatadict['IOVDbGlobalTag'] = f.fileinfos['conditions_tag']
                except Exception:
                    logOverlayReadMetadata.warning("Failed to find IOVDbGlobalTag.")
                    return metadatadict,taginfometadata,digimetadatadict,False
        if '/Digitization/Parameters' in f.infos['metadata'].keys():
            ## We have the RDO file here:
            if len(f.run_numbers)>0 :
                logOverlayReadMetadata.info("Setting digitizationFlags.dataRunNumber to %i, to match the presampled RDO File", f.run_numbers[0])
                from Digitization.DigitizationFlags import digitizationFlags
                digitizationFlags.dataRunNumber = f.run_numbers[0]
            ##Do useful things with Digi MetaData
            digimetadatadict = f.infos['metadata']['/Digitization/Parameters']
            logOverlayReadMetadata.debug("Keys in /Digitization/Parameters MetaData:")
            logOverlayReadMetadata.debug(digimetadatadict.keys())

            if 'IOVDbGlobalTag' in digimetadatadict.keys():
                if 'IOVDbGlobalTag' in metadatadict.keys():
                    if not re.match(metadatadict['IOVDbGlobalTag'], digimetadatadict['IOVDbGlobalTag']):
                        logOverlayReadMetadata.info("Overriding original SimDict['IOVDbGlobalTag']", metadatadict['IOVDbGlobalTag'])
                        logOverlayReadMetadata.debug("with RDODict['IOVDbGlobalTag']", digimetadatadict['IOVDbGlobalTag'])
                        metadatadict['IOVDbGlobalTag'] = digimetadatadict['IOVDbGlobalTag']
                        logOverlayReadMetadata.debug("Updated SimDict['IOVDbGlobalTag']", metadatadict['IOVDbGlobalTag'])
            if 'DigitizedDetectors' in digimetadatadict.keys():
                if 'DigitizedDetectors' in metadatadict.keys():
                    logOverlayReadMetadata.debug("Overriding original SimDict['SimulatedDetectors']", metadatadict['SimulatedDetectors'])
                    logOverlayReadMetadata.debug("with RDODict['DigitizedDetectors']", digimetadatadict['DigitizedDetectors'])
                    metadatadict['SimulatedDetectors'] = digimetadatadict['DigitizedDetectors']
                    logOverlayReadMetadata.debug("Updated SimDict['SimulatedDetectors']", metadatadict['SimulatedDetectors'])
            pass
        if '/TagInfo' in f.infos['metadata'].keys():
            taginfometadata=f.infos['metadata']['/TagInfo']
            pass
    Nkvp = len(metadatadict)
    ## Dictionary should not be empty
    if Nkvp==0 :
        logOverlayReadMetadata.error("Found %s KEY:VALUE pairs in %s Simulation MetaData." , Nkvp,inputtype)
        return metadatadict,taginfometadata,digimetadatadict,False
    else:
        ##Patch for older hit files
        if 'hitFileMagicNumber' not in metadatadict.keys():
            metadatadict['hitFileMagicNumber'] = 0
            logOverlayReadMetadata.debug("hitFileMagicNumber key missing from %s Simulation MetaData Dictionary. Adding dummy entry.",inputtype)
        if 'SimulatedDetectors' not in metadatadict.keys():
            if 'eventdata_items' in f.infos.keys():
                metadatadict['SimulatedDetectors'] = hitColls2SimulatedDetectors(f.infos['eventdata_items'])
            else :
                metadatadict['SimulatedDetectors'] = ['pixel','SCT','TRT','BCM','Lucid','LAr','Tile','MDT','CSC','TGC','RPC','Truth']
        ## Check whether we should use the old names for the Tile CaloCalibrationHit containers
        if 'eventdata_items' in f.infos.keys():
            checkTileCalibrationHitFormat(f.infos['eventdata_items'])
        else :
            from Digitization.DigitizationFlags import digitizationFlags
            digitizationFlags.experimentalDigi += ['OldTileCalibHitContainers']

        ## Generate optional containers list
        if 'eventdata_items' in f.infos.keys():
            metadatadict['OptionalContainers'] = listOptionalContainers(f.infos['eventdata_items'])
        else:
            metadatadict['OptionalContainers'] = {}

        ## Check for legacy EventInfo
        if 'eventdata_items' in f.infos.keys():
            metadatadict['LegacyEventInfo'] = checkLegacyEventInfo(f.infos['eventdata_items'])
        else:
            metadatadict['LegacyEventInfo'] = False

        ## Check for legacy NSW containers
        if 'eventdata_items' in f.infos.keys():
            metadatadict['LegacyNSWContainers'] = checkLegacyNSWContainers(f.infos['eventdata_items'])
        else:
            metadatadict['LegacyNSWContainers'] = False

        ##End of Patch for older hit files
        logOverlayReadMetadata.debug("%s Simulation MetaData Dictionary Successfully Created.",inputtype)
        logOverlayReadMetadata.debug("Found %s KEY:VALUE pairs in %s Simulation MetaData." , Nkvp,inputtype)
        logOverlayReadMetadata.debug("KEYS FOUND: %s", metadatadict.keys())
        return metadatadict,taginfometadata,digimetadatadict,True

def signalMetaDataCheck(metadatadict):
    import re
    simkeys = metadatadict.keys()
    logOverlayReadMetadata.info("Checking Digitization properties against Signal Simulation MetaData...")
    ## Check the PhysicsList set agrees with that used in the simulation
    if not skipCheck('PhysicsList'):
        if 'PhysicsList' in simkeys:
            from Digitization.DigitizationFlags import digitizationFlags
            if re.match(metadatadict['PhysicsList'], digitizationFlags.physicsList.get_Value()):
                logOverlayReadMetadata.debug("Digitization properties matches Signal Simulation MetaData. [PhysicsList = %s]", metadatadict['PhysicsList'])
            else:
                logOverlayReadMetadata.warning("Digitization properties PhysicsList does not match the PhysicsList used in the Simulation step! Assume the PhysicsList from the Simulation step is correct!")
                digitizationFlags.physicsList = metadatadict['PhysicsList']
                logOverlayReadMetadata.info("Set digitizationFlags.physicsList = %s", digitizationFlags.physicsList.get_Value())
        else:
            logOverlayReadMetadata.error("PhysicsList key not found in Simulation MetaData!")

    ## Check the DetDescrVersion set agrees with that used in the simulation
    if not skipCheck('SimLayout'):
        if 'SimLayout' in simkeys:
            from AthenaCommon.GlobalFlags import globalflags
            if re.match(metadatadict['SimLayout'], globalflags.DetDescrVersion.get_Value()):
                logOverlayReadMetadata.debug("Digitization properties matches Signal Simulation MetaData. [DetDescrVersion = %s]",
                                            globalflags.DetDescrVersion.get_Value())
            else:
                logOverlayReadMetadata.warning("Input DetDescrVersion does not match the value used in the Simulation step!")
                from AthenaCommon.AppMgr import ServiceMgr
                ## FIXME - should not be relying on GeoModelSvc being initialized at this point.
                if hasattr( ServiceMgr, "GeoModelSvc") and ServiceMgr.GeoModelSvc.IgnoreTagDifference is True:
                    logOverlayReadMetadata.warning("Global jobproperties: [DetDescrVersion = %s], Signal Simulation MetaData: [SimLayout = %s]",
                                                  globalflags.DetDescrVersion.get_Value(), metadatadict['SimLayout'])
                    logOverlayReadMetadata.warning("Ignore Tag Difference Requested - doing nothing.")
                else:
                    logOverlayReadMetadata.warning("Assume the value from the Simulation step is correct!")
                    ## needs to be done this way as Digi_tf locks it
                    if globalflags.DetDescrVersion.is_locked() :
                        globalflags.DetDescrVersion.unlock()
                    globalflags.DetDescrVersion.set_Value_and_Lock( metadatadict['SimLayout'] )
                    logOverlayReadMetadata.warning("Set globalflags.DetDescrVersion = %s",globalflags.DetDescrVersion.get_Value())
        else:
            logOverlayReadMetadata.error("SimLayout key not found in Simulation MetaData!")

    ## Check the Conditions Tag set against that used in the simulation
    if not skipCheck('IOVDbGlobalTag'):
        if 'IOVDbGlobalTag' in simkeys:
            from Digitization.DigitizationFlags import digitizationFlags
            if (digitizationFlags.IOVDbGlobalTag.statusOn):
                logOverlayReadMetadata.info("Digitization properties: [IOVDbGlobalTag = %s], Signal Simulation MetaData: [IOVDbGlobalTag = %s]",
                                           digitizationFlags.IOVDbGlobalTag.get_Value(), metadatadict['IOVDbGlobalTag'])
            else:
                digitizationFlags.IOVDbGlobalTag = metadatadict['IOVDbGlobalTag']
                logOverlayReadMetadata.debug("Set Digitization properties to match Signal Simulation Metadata: [IOVDbGlobalTag = %s]",
                                            digitizationFlags.IOVDbGlobalTag.get_Value())
        else:
            logOverlayReadMetadata.error("IOVDbGlobalTag key not found in Simulation MetaData!")

    ## Set the TRTRangeCut digitizationFlag based on what was used during the simulation.
    if not skipCheck('TRTRangeCut'):
        if 'TRTRangeCut' in simkeys:
            from Digitization.DigitizationFlags import digitizationFlags
            if hasattr( digitizationFlags, 'TRTRangeCut'):
                digitizationFlags.TRTRangeCut = metadatadict['TRTRangeCut']
                logOverlayReadMetadata.debug("Set Digitization properties to match Signal Simulation Metadata: [TRTRangeCut = %s]",
                                            digitizationFlags.TRTRangeCut.get_Value())
        else:
            logOverlayReadMetadata.warning("TRTRangeCut key not found in Simulation MetaData!")

    ## Record the G4Version used in the simulation, so that Digitization Algorithms can use this information
    if not skipCheck('G4Version'):
        if 'G4Version' in simkeys:
            from Digitization.DigitizationFlags import digitizationFlags
            digitizationFlags.SimG4VersionUsed = metadatadict['G4Version']
            logOverlayReadMetadata.debug("digitizationFlags.SimG4VersionUsed = Value from Sim Metadata = %s ", digitizationFlags.SimG4VersionUsed.get_Value())
        else:
            logOverlayReadMetadata.error("G4Version key not found in Simulation MetaData!")

    ## Check which sub-detectors were simulated

    ## Digitization will only digitize detectors which have been simulated.
    ## If users want to digitize an un simulated detector this will be an expert
    ## action which will require hacking the python code.
    if not skipCheck('SimulatedDetectors'):
        if 'SimulatedDetectors' in simkeys:
            if isinstance(metadatadict['SimulatedDetectors'], str):
                simulatedDetectors = eval(metadatadict['SimulatedDetectors']) # convert from str to list of str
            else:
                simulatedDetectors = metadatadict['SimulatedDetectors']
            simulatedDetectors[:] = [x.lower() if x == 'Pixel' else x for x in simulatedDetectors] # to cope with CA-based inputs where Pixel rather than pixel is used
            simulatedDetectors[:] = ['MM' if x == 'Micromegas' else x for x in simulatedDetectors] # to cope with legacy NSW naming
            from AthenaCommon.DetFlags import DetFlags
            logOverlayReadMetadata.debug("Switching off subdetectors which were not simulated")
            possibleSubDetectors=['pixel','SCT','TRT','BCM','Lucid','ZDC','ALFA','AFP','FwdRegion','LAr','HGTD','Tile','MDT','CSC','TGC','RPC','MM','sTGC','Truth']
            switchedOffSubDetectors=[]
            for subdet in possibleSubDetectors:
                if subdet not in simulatedDetectors:
                    attrname = subdet+"_setOff"
                    checkfn = getattr(DetFlags, attrname, None)
                    if checkfn is not None:
                        cmd='DetFlags.%s_setOff()' % subdet
                        logOverlayReadMetadata.debug(cmd)
                        checkfn()
                    switchedOffSubDetectors+=[subdet]
            if switchedOffSubDetectors:
                logOverlayReadMetadata.info("Ensured %s sub-detectors which were not simulated were switched off: %s", len(switchedOffSubDetectors), switchedOffSubDetectors)
            else:
                logOverlayReadMetadata.info("All sub-detectors were simulated, so none needed to be switched off in digitization.")
            DetFlags.Print()

    # Check for optional containers presence
    if not skipCheck('OptionalContainers'):
        from OverlayCommonAlgs.OverlayFlags import overlayFlags
        overlayFlags.optionalContainerMap = metadatadict['OptionalContainers']

    # Check for legacy EventInfo presence
    if not skipCheck('LegacyEventInfo'):
        from OverlayCommonAlgs.OverlayFlags import overlayFlags
        overlayFlags.processLegacyEventInfo.set_Value_and_Lock(metadatadict['LegacyEventInfo'])

    # Check for legacy NSW containers
    if not skipCheck('LegacyNSWContainers'):
        if metadatadict['LegacyNSWContainers']:
            from Digitization.DigitizationFlags import digitizationFlags
            digitizationFlags.experimentalDigi += ['LegacyNSWContainers']

    ## Any other checks here
    logOverlayReadMetadata.info("Completed checks of Digitization properties against Signal Simulation MetaData.")

def pileupMetaDataCheck(sigsimdict,pileupsimdict):
    """Check the metadata for presampled pileup RDO File"""
    result = True
    import re
    pileupkeys = pileupsimdict.keys()
    sigkeys = sigsimdict.keys()
    pileuptype = "Presampled"
    longpileuptype = "Presampled Pile-up RDO File"
    ##Loop over MetaData keys which must have matching values
    SigkeysToCheck = [ 'PhysicsList', 'SimLayout', 'MagneticField','hitFileMagicNumber' ]#, 'WorldZRange' ]
    for o in SigkeysToCheck:
        if skipPileUpCheck(o, pileuptype):
            continue
        try:
            assert o in pileupkeys
        except AssertionError:
            logOverlayReadMetadata.error("%s key missing from %s Simulation MetaData!", o, longpileuptype)
            raise AssertionError("Simulation MetaData key not found")
        try:
            assert o in sigkeys
        except AssertionError:
            logOverlayReadMetadata.error("%s key missing from Signal Simulation MetaData!", o)
            raise AssertionError("Simulation MetaData key not found")
        try:
            if not isinstance(pileupsimdict[o],type(sigsimdict[o])):
                assert re.match(str(pileupsimdict[o]), str(sigsimdict[o]))
            else:
                if isinstance(pileupsimdict[o],str):
                    assert re.match(pileupsimdict[o], sigsimdict[o])
                elif isinstance(pileupsimdict[o],int):
                    assert (pileupsimdict[o]==sigsimdict[o])
                else:
                    assert re.match(str(pileupsimdict[o]), str(sigsimdict[o]))
        except AssertionError:
            logOverlayReadMetadata.error("Simulation MetaData mismatch! %s: [%s = %s] Signal: [%s = %s]", longpileuptype, o, pileupsimdict[o], o, sigsimdict[o])
            raise AssertionError("Simulation MetaData mismatch")
        logOverlayReadMetadata.debug("%s Sim MetaData matches Signal Sim MetaData. [%s = %s]", longpileuptype, o, sigsimdict[o])
    ##Ideally these keys would have matching values, but it should be OK if not.
    WarningKeys = [ 'IOVDbGlobalTag', 'G4Version' ]
    for o in WarningKeys:
        if skipPileUpCheck(o, pileuptype):
            continue
        try:
            assert o in pileupkeys
        except AssertionError:
            logOverlayReadMetadata.error("%s key missing from %s Simulation MetaData!", o, longpileuptype)
            raise AssertionError("Simulation MetaData key not found")
        try:
            assert o in sigkeys
        except AssertionError:
            logOverlayReadMetadata.error("%s key missing from Signal Simulation MetaData!", o, longpileuptype)
            raise AssertionError("Simulation MetaData key not found")
        if not re.match(pileupsimdict[o], sigsimdict[o]):
            logOverlayReadMetadata.warning("Simulation MetaData mismatch! %s: [%s = %s] Signal: [%s = %s]", longpileuptype, o, pileupsimdict[o], o, sigsimdict[o])
        else:
            logOverlayReadMetadata.debug("%s Sim MetaData matches Signal Sim MetaData. [%s = %s]", longpileuptype, o, sigsimdict[o])

    ## Check that the same sub-detectors were simulated in signal and background inputs
    if (not skipPileUpCheck('SimulatedDetectors', pileuptype)) and ('SimulatedDetectors' in sigkeys):
        switchedOffSubDetectors=[]
        for subdet in sigsimdict['SimulatedDetectors']:
            if subdet not in pileupsimdict['SimulatedDetectors']:
                if subdet == 'MM':
                    if 'Micromegas' not in pileupsimdict['SimulatedDetectors']:
                        switchedOffSubDetectors+=[subdet]    
                else:
                    switchedOffSubDetectors+=[subdet]
        if switchedOffSubDetectors:
            logOverlayReadMetadata.error("%s sub-detectors were sinmulated in the signal sample, but not in the %s background sample: %s", len(switchedOffSubDetectors), longpileuptype, switchedOffSubDetectors)
            raise AssertionError("Some simulated sub-detectors from signal sample are missing in the background samples.")
        else:
            logOverlayReadMetadata.debug("All sub-detectors simulated in the signal sample were also simulated in the %s background sample.", longpileuptype)

    # Check for optional containers presence
    from AthenaCommon.AthenaCommonFlags import athenaCommonFlags
    if not skipCheck('OptionalContainers'):
        # Combine the two dictionaries
        optionalContainers = sigsimdict['OptionalContainers']
        for key, value in pileupsimdict['OptionalContainers'].items():
            if  key not in ['TrackRecordCollection']:
                if key in optionalContainers:
                    optionalContainers[key] |= value #Here the expectation is that the values are sets
                else:
                    optionalContainers[key] = value

        from OverlayCommonAlgs.OverlayFlags import overlayFlags
        overlayFlags.optionalContainerMap.set_Value_and_Lock(optionalContainers)
    elif athenaCommonFlags.DoFullChain():
        optionalContainers = pileupsimdict['OptionalContainers']

        from OverlayCommonAlgs.OverlayFlags import overlayFlags
        overlayFlags.optionalContainerMap.set_Value_and_Lock(optionalContainers)

    return result


def tagInfoMetaDataCheck(sigtaginfodict,pileuptaginfodict):
    result = True
    """Check the metadata for presampled pileup RDO File"""
    pileupkeys = pileuptaginfodict.keys()
    logOverlayReadMetadata.debug("Signal /TagInfo ", sigtaginfodict)
    logOverlayReadMetadata.debug("Pileup /TagInfo ", pileuptaginfodict)
    sigkeys = sigtaginfodict.keys()
    sigOnlyDict = dict()
    sigOnlyKeySet = set(sigkeys).difference(set(pileupkeys))
    logOverlayReadMetadata.debug("The following keys only appear in Signal /TagInfo MetaData:")
    logOverlayReadMetadata.debug(sigOnlyKeySet)
    for key in sigOnlyKeySet:
        sigOnlyDict[key] = str(sigtaginfodict[key])
        logOverlayReadMetadata.debug("key: ", key, "value: ", sigtaginfodict[key])
        pass
    from OverlayCommonAlgs.OverlayFlags import overlayFlags
    overlayFlags.extraTagInfoPairs = sigOnlyDict
    keysToCompareSet = set(sigkeys).intersection(set(pileupkeys))
    logOverlayReadMetadata.debug("The following keys appear in Signal and Presampled /TagInfo MetaData:")
    logOverlayReadMetadata.debug(keysToCompareSet)
    return result

def readInputFileMetadata():
    logOverlayReadMetadata.info("Checking for Signal Simulation MetaData...")
    import PyUtils.AthFile as af
    af.server.load_cache('digitization-afcache.ascii')

    #--------------------------------------------------
    # Check for the Run Number in the first Input file
    #--------------------------------------------------
    from AthenaCommon.AthenaCommonFlags import athenaCommonFlags
    from Digitization.DigitizationFlags import digitizationFlags
    from OverlayCommonAlgs.OverlayFlags import overlayFlags

    digitizationFlags.simRunNumber = int(digitizationFlags.getHitFileRunNumber(athenaCommonFlags.PoolHitsInput.get_Value()[0]))

    sigsimdict,sigtaginfodict,_,result = buildDict("Signal", athenaCommonFlags.PoolHitsInput.get_Value()[0])
    if result :
        signalMetaDataCheck(sigsimdict)

        if overlayFlags.isDataOverlay():
            if 'RunNumber' in sigsimdict.keys():
                year = sigsimdict['RunNumber'] % 100
                logOverlayReadMetadata.info("Found Year = %s", year)
                from RecExConfig.RecFlags import rec
                rec.projectName = 'data' + str(year)
        else:
            ## Check Pileup Simulation Parameters match those used for signal files
            result = True
            longpileuptype= "presampled pile-up"
            logOverlayReadMetadata.info("Checking %s MetaData against Signal Simulation MetaData...", longpileuptype)
            pileupsimdict,pileuptaginfodict,pileupdigidict,result = buildDict(longpileuptype, athenaCommonFlags.PoolRDOInput.get_Value()[0])
            if not result:
                logOverlayReadMetadata.warning("Failed to Create %s Simulation MetaData Dictionary from file %s.", longpileuptype, athenaCommonFlags.PoolHitsInput.get_Value()[0])
            else:
                if pileupMetaDataCheck(sigsimdict,pileupsimdict):
                    logOverlayReadMetadata.info("Presampled RDO File Simulation MetaData matches Signal Simulation MetaData.")
                if tagInfoMetaDataCheck(sigtaginfodict,pileuptaginfodict):
                    logOverlayReadMetadata.info("Presampled RDO File TagInfo MetaData matches Signal TagInfo MetaData.")
            ## All checks completed here
            logOverlayReadMetadata.info("Completed all checks against Signal Simulation MetaData.")

            if pileupdigidict:
                from .OverlayWriteMetaData import writeOverlayDigitizationMetadata
                writeOverlayDigitizationMetadata(pileupdigidict)

            del pileupsimdict
            del pileuptaginfodict
            del pileupdigidict
    else:
        logOverlayReadMetadata.info("Failed to Create Signal MetaData Dictionaries from file %s", athenaCommonFlags.PoolHitsInput.get_Value()[0])

    ## control where metadata can be used
    del sigsimdict
    del sigtaginfodict

def readPresampledMetadata():
    logOverlayReadMetadata.info("Checking for Presampled Pile-up MetaData...")
    import os
    from AthenaCommon.AthenaCommonFlags import athenaCommonFlags
    from AthenaCommon.GlobalFlags import globalflags
    from G4AtlasApps.SimFlags import simFlags
    longpileuptype= "presampled pile-up"
    pileupsimdict,pileuptaginfodict,pileupdigidict,result = buildDict(longpileuptype, athenaCommonFlags.PoolRDOInput.get_Value()[0])

    from AthenaCommon.JobProperties import JobProperty
    simParams = [sf for sf in dir(simFlags) if isinstance(getattr(simFlags, sf), JobProperty)]
    sigsimdict = dict()
    KeysToCheck = [ 'PhysicsList', 'SimLayout', 'MagneticField', 'hitFileMagicNumber' , 'IOVDbGlobalTag', 'G4Version' ]
    for sp in simParams:
        if sp in KeysToCheck:
            sigsimdict.setdefault(sp, [])
            sigsimdict[sp] = getattr(simFlags, sp).get_Value()
    sigsimdict['hitFileMagicNumber'] = '0' # Hard-coded simulation hit (as in ISF_Metadata.py file)
    sigsimdict['G4Version'] = str(os.environ['G4VERS'])
    sigsimdict['IOVDbGlobalTag'] = globalflags.ConditionsTag()

    if pileupMetaDataCheck(sigsimdict,pileupsimdict):
        logOverlayReadMetadata.info("Presampled RDO File Simulation MetaData matches Signal Simulation MetaData.")
    if pileupdigidict:
        from EventOverlayJobTransforms.OverlayWriteMetaData import writeOverlayDigitizationMetadata
        writeOverlayDigitizationMetadata(pileupdigidict)

    del pileupsimdict
    del pileuptaginfodict
    del pileupdigidict
