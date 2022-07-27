/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TRIGCOSTMONITOR_TRIGCOSTSVC_H
#define TRIGCOSTMONITOR_TRIGCOSTSVC_H

#include <atomic>
#include <shared_mutex>
#include <thread>
#include <vector>

#include "GaudiKernel/ToolHandle.h"

#include "AthenaBaseComps/AthService.h"
#include "AthContainers/ConstDataVector.h"
#include "StoreGate/WriteHandle.h"

#include "xAODTrigger/TrigCompositeContainer.h"
#include "TrigSteeringEvent/TrigRoiDescriptorCollection.h"

#include "TrigCostMonitor/ITrigCostSvc.h"
#include "AlgorithmPayload.h"
#include "TrigCostDataStore.h"

/**
 * @class TrigCostSvc
 * @brief AthenaMT service to collect trigger cost data from all threads and summarise it at the end of the event
 *
 * The main hooks into this service are: 
 *  HLTSeeding - To clear the internal storage and flag the event for monitoring.
 *  TrigCostAuditor - To inform the service when algorithms start and stop executing
 *  HLTROBDataProviderSvc - To inform the service about requests for data ROBs 
 *  HLTSummaryAlg - To inform the service when the HLT has finished, and to receive the persistent payload 
 */
class TrigCostSvc : public extends <AthService, ITrigCostSvc> {
  public:

  /**
   * @brief Standard ATLAS Service constructor
   * @param[in] name The service's name
   * @param[in] svcloc A pointer to a service location service
   */
  TrigCostSvc(const std::string& name, ISvcLocator* pSvcLocator);

  /**
   * @brief Destructor. Currently nothing to delete.
   */
  virtual ~TrigCostSvc();

  /**
   * @brief Initialise, create enough storage to store m_eventSlots.
   */
  virtual StatusCode initialize() override;

  /**
   * @brief Finalize, act on m_saveHashes.
   */
  virtual StatusCode finalize() override;

  /**
   * @brief Implementation of ITrigCostSvc::startEvent.
   * @param[in] context The event context
   * @param[in] enableMonitoring Sets if the event should be monitored or not. Not monitoring will save CPU
   * @return Success unless monitoring is enabled and the service's data stores can not be cleared for some reason
   */
  virtual StatusCode startEvent(const EventContext& context, const bool enableMonitoring = true) override; 

  /**
   * @brief Implementation of ITrigCostSvc::processAlg.
   * @param[in] context The event context
   * @param[in] caller Name of the algorithm to audit CPU usage for
   * @param[in] type If we are Before or After the algorithm's execution
   */
  virtual StatusCode processAlg(const EventContext& context, const std::string& caller, const AuditType type) override; 

  /**
   * @brief Implementation of ITrigCostSvc::endEvent.
   * @param[in] context The event context
   * @param[out] costOutputHandle Write handle to fill with execution summary if the event was monitored
   * @param[out] rosOutputHandle Write handle to fill with ROS requests summary if the event was monitored
   */
  virtual StatusCode endEvent(const EventContext& context, SG::WriteHandle<xAOD::TrigCompositeContainer>& costOutputHandle, SG::WriteHandle<xAOD::TrigCompositeContainer>& rosOutputHandle) override; 

  /**
   * @return If the current context is flagged as being monitored. 
   * @param[in] context The event context
   */
  virtual bool isMonitoredEvent(const EventContext& context, const bool includeMultiSlot = true) const override;

  /**
   * @brief Implementation of ITrigCostSvc::monitorROS.
   * @param[in] context The event context
   * @param[in] payload ROB data to be associated with ROS
   */
  virtual StatusCode monitorROS(const EventContext& context, robmonitor::ROBDataMonitorStruct payload) override;

  /**
   * @return Generate timeout report with the most time consuming algorithms
   * @param[in] context The event context
   * @param[out] report Created report with algorithms and times (in ms)
   */
  virtual StatusCode generateTimeoutReport(const EventContext& context, std::string& report) override;


  /**
   * @brief Discard a cost monitored event
   * @param[in] context The event context
   */
  virtual StatusCode discardEvent(const EventContext& context) override;

  private: 

  /**
   * Internal call to save monitoring data for a given AlgorithmIdentifier
   * @param[in] context The event context
   * @param[in] ai The AlgorithmIdentifier key to store
   * @param[in] now The timestamp to store (amoung other values)
   * @param[in] type The type of the audit event to store
   * @return Success if the data are saved
   */
  StatusCode monitor(const EventContext& context, const AlgorithmIdentifier& ai, const TrigTimeStamp& now, const AuditType type);

  /**
   * Sanity check that the job is respecting the number of slots which were declared at config time
   * @param[in] context The event context
   * @return Success if the m_eventMonitored array is range, Failure if access request would overflow
   */
  StatusCode checkSlot(const EventContext& context) const;

  /**
   * @breif Internal function to return a RoI from an extended event context context
   * @param[in] context The event context
   * @return RoIId from the ATLAS extended event context. Or, AlgorithmIdentifier::s_noView = -1 for no RoIIdentifier
   */
  int32_t getROIID(const EventContext& context);

  /** 
   * @class ThreadHashCompare
   * @brief Static hash and equal members as required by tbb::concurrent_hash_map
   */
  struct ThreadHashCompare {
    static size_t hash(const std::thread::id& thread);
    static bool equal(const std::thread::id& x, const std::thread::id& y);
  };

  size_t  m_eventSlots; //!< Number of concurrent processing slots. Cached from Gaudi
  std::unique_ptr< std::atomic<bool>[] > m_eventMonitored; //!< Used to cache if the event in a given slot is being monitored.
  std::unique_ptr< std::shared_mutex[] > m_slotMutex; //!< Used to control and protect whole-table operations.
  std::mutex m_globalMutex; //!< Used to protect all-slot modifications.
  TrigCostDataStore<AlgorithmPayload> m_algStartInfo; //!< Thread-safe store of algorithm start payload.
  TrigCostDataStore<TrigTimeStamp> m_algStopTime; //!< Thread-safe store of algorithm stop times.
  TrigCostDataStore<std::vector<robmonitor::ROBDataMonitorStruct>> m_rosData; //!< Thread-safe store of ROS data

  tbb::concurrent_hash_map<std::thread::id, AlgorithmIdentifier, ThreadHashCompare> m_threadToAlgMap; //!< Keeps track of what is running right now in each thread.

  std::unordered_map<uint32_t, uint32_t> m_threadToCounterMap; //!< Map thread's hash ID to a counting numeral
  size_t m_threadCounter; //!< Count how many unique thread ID we have seen 


  Gaudi::Property<bool>        m_monitorAllEvents{this, "MonitorAllEvents", false, "Monitor every HLT event, e.g. for offline validation."};
  Gaudi::Property<bool>        m_enableMultiSlot{this, "EnableMultiSlot", false, "Monitored events in the MasterSlot collect data from events running in other slots."};
  Gaudi::Property<bool>        m_saveHashes{this, "SaveHashes", false, "Store a copy of the hash dictionary for easier debugging"};
  Gaudi::Property<size_t>      m_masterSlot{this, "MasterSlot", 0, "The slot responsible for saving MultiSlot data"};
  Gaudi::Property<std::string> m_hltSeedingName{this, "HLTSeedingName", "HLTSeeding", "The name of the Gaudi Configurable of type HLTSeeding"};
  Gaudi::Property<std::string> m_decisionSummaryMakerAlgName{this, "DecisionSummaryMakerAlgName", "DecisionSummaryMakerAlg", "The name of the Gaudi Configurable of type DecisionSummaryMakerAlg"};



};

#endif // TRIGCOSTMONITOR_TRIGCOSTSVC_H
