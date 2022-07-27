/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef PILEUPEVENTLOOPMGR_H
#define PILEUPEVENTLOOPMGR_H
/** @file PileUpEventLoopMgr.h
    @brief The ATLAS event loop for pile-up applications.
    @author Paolo Calafiura
*/

// Base class headers
#include "AthenaKernel/IEventSeek.h"
#include "GaudiKernel/MinimalEventLoopMgr.h"

// Athena headers
#include "AthenaBaseComps/AthMessaging.h"
#include "PileUpTools/PileUpStream.h"

// Gaudi headers
#include "Gaudi/Property.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/IAlgExecStateSvc.h"
#include <string>

// Forward declarations
class IBeamIntensity;
class IBeamLuminosity;
class IBkgStreamsCache;
class IEvtSelector;
class IIncidentSvc;
class PileUpMergeSvc;
class StoreGateSvc;
class EventContext;
class EventID;
class IEvtIdModifierSvc;

/** @class PileUpEventLoopMgr
    @brief The ATLAS event loop for pile-up applications.
*/

class PileUpEventLoopMgr : virtual public IEventSeek,
                           public MinimalEventLoopMgr,
                           public AthMessaging
{
public:

  /// Standard Constructor
  PileUpEventLoopMgr(const std::string& nam, ISvcLocator* svcLoc);
  /// Standard Destructor
  virtual ~PileUpEventLoopMgr();

public:
  /// implementation of IAppMgrUI::initialize
  virtual StatusCode initialize();
  /// implementation of IAppMgrUI::finalize
  virtual StatusCode finalize();
  /// implementation of IAppMgreUI::terminate
  //  virtual StatusCode terminate();
  /// implementation of IAppMgrUI::nextEvent
  virtual StatusCode nextEvent(int maxevt);
  /// implementation of IEventProcessor::executeEvent(void* par)
  virtual StatusCode executeEvent( EventContext &&ctx );

  /// Seek to a given event
  virtual StatusCode seek(int evt);
  /// Return the current event count
  virtual int curEvent() const;

  virtual void modifyEventContext(EventContext& ctx, const EventID& eID, bool consume_modifier_stream);

  virtual StatusCode queryInterface(const InterfaceID& riid,
                                    void** ppvInterface);

  using AthMessaging::msg;
  using AthMessaging::msgLvl;


private:
  /// Reference to the Algorithm Execution State Svc
  SmartIF<IAlgExecStateSvc>  m_aess;

  /// setup input and overlay selectors and iters
  StatusCode setupStreams();

  /// Run the algorithms for the current event
  virtual StatusCode executeAlgorithms(const EventContext& ctx);

  ///return the 'fake BCID' corresponding to bunchXing
  inline unsigned int getBCID(int bunchXing, unsigned int centralBCID) const {
    //FIXME to be completely safe this should should probably depend on the bunch spacing too. Perhaps that concept should be deprecated though?
    return static_cast<unsigned int>((((bunchXing + static_cast<int>(centralBCID)) % static_cast<int>(m_maxBunchCrossingPerOrbit)) + static_cast<int>(m_maxBunchCrossingPerOrbit) )  % static_cast<int>(m_maxBunchCrossingPerOrbit));
  }

  /// Incident Service
  ServiceHandle<IIncidentSvc> m_incidentSvc;

  /// PileUp Merge Service
  ServiceHandle<PileUpMergeSvc> m_mergeSvc;

  /// Input Stream
  PileUpStream m_origStream;

  /// output store
  ServiceHandle<StoreGateSvc> m_evtStore;              // overlaid (output) event store
  
  typedef ServiceHandle<IEvtIdModifierSvc> IEvtIdModifierSvc_t;
  /// @property Reference to the EventID modifier Service
  IEvtIdModifierSvc_t m_evtIdModSvc;

  //unsigned int m_nInputs;
  //unsigned int m_nStores;

  /// @name Properties
  //@{
  /// Original (Physics) Event selector (background for overlay).
  ServiceHandle<IEvtSelector> m_origSel;
  /// Signal Event selector (for overlay).
  ServiceHandle<IEvtSelector> m_signalSel;
  /// BkgStreamsCaches managing background events
  ToolHandleArray<IBkgStreamsCache> m_caches;
  /// (max) minBias interactions per Xing, for setting MC luminosity
  Gaudi::Property<float> m_maxCollPerXing;

  /// Xing frequency(ns);
  Gaudi::Property<float> m_xingFreq;
  /// first xing to be simulated (0th xing is 1st after trigger)
  Gaudi::Property<int> m_firstXing;
  /// last xing to be simulated (0th xing is 1st after trigger)
  Gaudi::Property<int> m_lastXing;

  /// property: allow sub evts EOF condition when maxevt==-1
  Gaudi::Property<bool> m_allowSubEvtsEOF;

  /// property: process bkg events xing by xing without caching them
  Gaudi::Property<bool> m_xingByXing;

  /// property: control behaviour of event loop on algorithm failure
  Gaudi::Property<int> m_failureMode;

  /// SG key for the EventInfoContainer
  Gaudi::Property<std::string> m_evinfName;

  /// SG key for the EventInfoContainer
  Gaudi::Property<std::string> m_evinfContName;

  /// property: beam intensity service handle for beam profile in local time
  ServiceHandle<IBeamIntensity> m_beamInt;
  /// property: beam intensity service handle for luminosity profile in iovtime
  ServiceHandle<IBeamLuminosity> m_beamLumi;
  //@}

  /// current run number
  uint32_t m_currentRun;
  bool m_firstRun;

  /// max bunch crossings per orbit
  unsigned int m_maxBunchCrossingPerOrbit;

  int m_nevt;

  int m_ncurevt;
  bool m_skipExecAlgs;
  bool m_loadProxies;

  /// property: Default true. When set to false, this will allow the
  /// code to reproduce serial output in an AthenaMP job, albeit with
  /// a significant performance penalty.
  Gaudi::Property<bool> m_allowSerialAndMPToDiffer;

  Gaudi::Property<uint32_t> m_mcChannelNumber{ this, "MCChannelNumber", 0, "sample MC channel number" };
};
#endif // PILEUPTOOLS_PILEUPEVENTLOOPMGR_H
