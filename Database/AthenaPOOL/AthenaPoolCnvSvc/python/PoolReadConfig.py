# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from AthenaConfiguration.ComponentFactory import CompFactory
from AthenaKernel.EventIdOverrideConfig import EvtIdModifierSvcCfg, getMinMaxRunNumbers, getFirstLumiBlock
from AthenaConfiguration.Enums import ProductionStep

def EventSelectorAthenaPoolCfg(configFlags):
    result=ComponentAccumulator()
    EventSelectorAthenaPool=CompFactory.EventSelectorAthenaPool
    evSel=EventSelectorAthenaPool("EventSelector",
                                  InputCollections = configFlags.Input.Files,
                                  SkipEvents=configFlags.Exec.SkipEvents)
    if configFlags.Input.OverrideRunNumber:
        if not configFlags.Input.RunAndLumiOverrideList:
            DataRunNumber = -1
            FirstLB = 1
            InitialTimeStamp = 1
            OldRunNumber = -1
            if configFlags.Input.ConditionsRunNumber>0:
                # Behaviour for Digitization jobs using DataRunNumber
                DataRunNumber = configFlags.Input.ConditionsRunNumber
                FirstLB = 1
                InitialTimeStamp = configFlags.IOVDb.RunToTimestampDict.get(DataRunNumber, 1) # TODO fix repeated configuration
                if not configFlags.Sim.DoFullChain:
                    OldRunNumber = configFlags.Input.RunNumber[0] # CHECK this should be the Run Number from the HITS file
            elif configFlags.Input.RunNumber:
                # Behaviour for Simulation jobs
                DataRunNumber = configFlags.Input.RunNumber[0]
                FirstLB = configFlags.Input.LumiBlockNumber[0]
                InitialTimeStamp = configFlags.Input.TimeStamp[0]
            assert DataRunNumber >= 0, (
                "configFlags.Input.OverrideRunNumber was True, but provided DataRunNumber (%d) is negative. "
                "Use a real run number from data." % DataRunNumber)
            evSel.OverrideRunNumber = configFlags.Input.OverrideRunNumber
            evSel.RunNumber = DataRunNumber
            evSel.FirstLB = FirstLB
            evSel.InitialTimeStamp = InitialTimeStamp # Necessary to avoid a crash
            if hasattr(evSel, "OverrideRunNumberFromInput"):
                evSel.OverrideRunNumberFromInput = configFlags.Input.OverrideRunNumber
            if OldRunNumber > 0:
                evSel.OldRunNumber = OldRunNumber
        elif configFlags.Common.ProductionStep in [ProductionStep.Simulation, ProductionStep.FastChain]:
            # Behaviour for Simulation and FastChain jobs using RunAndLumiOverrideList
            minMax = getMinMaxRunNumbers(configFlags)
            evSel.OverrideRunNumber = configFlags.Input.OverrideRunNumber
            evSel.RunNumber = minMax[0]
            evSel.FirstLB = getFirstLumiBlock(configFlags, minMax[0])
            evSel.InitialTimeStamp = configFlags.IOVDb.RunToTimestampDict.get(minMax[0], 1) # TODO fix repeated configuration
            if hasattr(evSel, "OverrideRunNumberFromInput"):
                evSel.OverrideRunNumberFromInput = configFlags.Input.OverrideRunNumber
        else:
            # Behaviour for Digitization jobs using RunAndLumiOverrideList
            pass
        result.merge(EvtIdModifierSvcCfg(configFlags))
    result.addService(evSel)
    return result


def PoolReadCfg(configFlags):
    """
    Creates a ComponentAccumulator instance containing the 
    athena services required for POOL file reading
    """

    filenames=configFlags.Input.Files
    filenamesSecondary=configFlags.Input.SecondaryFiles

    result=ComponentAccumulator()

    PoolSvc=CompFactory.PoolSvc
    ProxyProviderSvc=CompFactory.ProxyProviderSvc
    AthenaPoolCnvSvc=CompFactory.AthenaPoolCnvSvc
    AthenaPoolAddressProviderSvc, EventSelectorAthenaPool, DoubleEventSelectorAthenaPool=CompFactory.getComps("AthenaPoolAddressProviderSvc","EventSelectorAthenaPool","DoubleEventSelectorAthenaPool",)
    EvtPersistencySvc=CompFactory.EvtPersistencySvc
    
    StoreGateSvc=CompFactory.StoreGateSvc

    result.addService(PoolSvc(MaxFilesOpen=configFlags.PoolSvc.MaxFilesOpen))
    apcs=AthenaPoolCnvSvc()
    apcs.InputPoolAttributes += ["DatabaseName = '*'; ContainerName = 'CollectionTree'; TREE_CACHE = '-1'"]
    result.addService(apcs)
    result.addService(EvtPersistencySvc("EventPersistencySvc",CnvServices=[apcs.getFullJobOptName(),])) #No service handle yet???


    result.addService(StoreGateSvc("MetaDataStore"))

    if filenamesSecondary:
        skipEventsPrimary=configFlags.Exec.SkipEvents
        skipEventsSecondary=configFlags.Exec.SkipEvents
        if configFlags.Overlay.SkipSecondaryEvents >= 0:
            skipEventsSecondary = configFlags.Overlay.SkipSecondaryEvents

        # Create DoubleEventSelector (universal for any seconday input type)
        evSel = DoubleEventSelectorAthenaPool("EventSelector",
                                              InputCollections=filenames)

        if configFlags.Overlay.DataOverlay:
            # In case of data overlay HITS are primary input
            evSel.SkipEvents = skipEventsPrimary

            # We have to check if we're running data overlay - BS is needed in this case
            from ByteStreamCnvSvc.ByteStreamConfig import ByteStreamReadCfg
            result.merge(ByteStreamReadCfg(configFlags))

            # We still have to add primary address provider
            apapsPrimary = AthenaPoolAddressProviderSvc("AthenaPoolAddressProviderSvcPrimary")
            apapsPrimary.DataHeaderKey = "EventSelector"
            result.addService(apapsPrimary)

            result.addService(ProxyProviderSvc(ProviderNames = [
                apapsPrimary.getFullJobOptName(),
            ])) #No service handle yet???
        else:
            # In case of MC overlay RDOs are primary input
            evSel.SkipEvents = skipEventsSecondary
            # Do not process secondary input metadata
            evSel.ProcessMetadata = False

            # We have primary and secondary pool inputs, create two address providers
            apapsPrimary = AthenaPoolAddressProviderSvc("AthenaPoolAddressProviderSvcPrimary")
            apapsPrimary.DataHeaderKey = "EventSelector"
            apapsPrimary.AttributeListKey = "Input"
            result.addService(apapsPrimary)
            apapsSecondary = AthenaPoolAddressProviderSvc("AthenaPoolAddressProviderSvcSecondary")
            apapsSecondary.DataHeaderKey = "SecondaryEventSelector"
            result.addService(apapsSecondary)

            result.addService(ProxyProviderSvc(ProviderNames = [
                apapsPrimary.getFullJobOptName(),
                apapsSecondary.getFullJobOptName()
            ])) #No service handle yet???

            secondarySel = EventSelectorAthenaPool("SecondaryEventSelector",
                                                   IsSecondary=True,
                                                   InputCollections=filenamesSecondary,
                                                   SkipEvents=skipEventsPrimary)
            result.addService(secondarySel)
        result.addService(evSel)
    else:
        # We have only primary inputs
        apaps=AthenaPoolAddressProviderSvc()
        result.addService(apaps)
        result.addService(ProxyProviderSvc(ProviderNames=[apaps.getFullJobOptName(),])) #No service handle yet???
        result.merge(EventSelectorAthenaPoolCfg(configFlags))
        evSel = result.getService("EventSelector")

    result.setAppProperty("EvtSel",evSel.getFullJobOptName())

    #(possibly) missing: MetaDataSvc

    return result
