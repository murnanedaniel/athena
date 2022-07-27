/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/
#include "HltROBDataProviderSvc.h"
#include "TrigKernel/HltExceptions.h"
#include "TrigTimeAlgs/TrigTimeStamp.h"

// Gaudi
#include "Gaudi/Interfaces/IOptionsSvc.h"

// hltinterface / data collector
#include "hltinterface/DataCollector.h"

// eformat
#include "eformat/Status.h"

// Athena

// STL includes
#include <algorithm>    // std::find

// Maximum number of ROB fragments in ROB buffer
static const size_t MAX_ROBFRAGMENTS = 4096;

HltROBDataProviderSvc::HltROBDataProviderSvc(const std::string& name, ISvcLocator* pSvcLocator) :
  base_class(name, pSvcLocator)
{
}

StatusCode HltROBDataProviderSvc::initialize()
{
  ATH_MSG_INFO("HltROBDataProviderSvc::" << __FUNCTION__ << ": name = " << name());
//===================================================================  
//      The filtering of ROBs can be configured with job options as:
//
//      for individual ROBs as :
//      ------------------------
//      ROBDataProviderSvc.filterRobWithStatus = [ (ROB SourceId, StatusCode to remove),
//                                                 (ROB SourceId 2, StatusCode to remove 2), ... ]
//      and:
//      ROBDataProviderSvc.filterRobWithStatus += [ (ROB SourceId n, StatusCode to remove n) ]
//
//      Example:
//      ROBDataProviderSvc.filterRobWithStatus  = [ (0x42002a,0x0000000f), (0x42002e,0x00000008) ]
//      ROBDataProviderSvc.filterRobWithStatus += [ (0x42002b,0x00000000) ]
//
//      for all ROBs of a given sub detector as :
//      -----------------------------------------
//      ROBDataProviderSvc.filterSubDetWithStatus = [ (Sub Det Id, StatusCode to remove),
//                                                    (Sub Det Id 2, StatusCode to remove 2), ... ]
//      and:
//      ROBDataProviderSvc.filterSubDetWithStatus += [ (Sub Det Id n, StatusCode to remove n) ]
//
//      Example:
//      ROBDataProviderSvc.filterSubDetWithStatus  = [ (0x41,0x00000000), (0x42,0x00000000) ]
//      ROBDataProviderSvc.filterSubDetWithStatus += [ (0x41,0xcb0002) ]
//
//      For valid ROB Source Ids, Sub Det Ids and ROB Status elements see the event format
//      document ATL-D-ES-0019 (EDMS)
//===================================================================   
  // get list of ROBs to filter out by status code
  for (unsigned int i = 0; i < m_filterRobWithStatus.value().size(); i++) {
    eformat::helper::SourceIdentifier tmpsrc(m_filterRobWithStatus.value()[i].first);
    if (tmpsrc.human_detector() != "UNKNOWN") {
      m_filterRobMap[tmpsrc.code()].push_back(m_filterRobWithStatus.value()[i].second);
    }
  }
   
  // get list of subdetectors to filter out by status code
  for (unsigned int i = 0; i < m_filterSubDetWithStatus.value().size(); i++) {
    eformat::helper::SourceIdentifier tmpsrc((eformat::SubDetector)m_filterSubDetWithStatus.value()[i].first, 0);
    if (tmpsrc.human_detector() != "UNKNOWN") {
      m_filterSubDetMap[tmpsrc.subdetector_id()].push_back(m_filterSubDetWithStatus.value()[i].second);
    }
  }
  ATH_MSG_INFO(" ---> Filter out empty ROB fragments                               = " << m_filterEmptyROB);

  // print list of ROBs to filter out by status code
  ATH_MSG_INFO(" ---> Filter out specific ROBs by Status Code: # ROBs              = " << m_filterRobMap.size());
  for (const auto& it : m_filterRobMap) {
    eformat::helper::SourceIdentifier tmpsrc(it.first);
    ATH_MSG_INFO("      RobId=0x" << MSG::hex << it.first << " -> in Sub Det = " << tmpsrc.human_detector());

    for (auto it_status: it.second) { 
      eformat::helper::Status tmpstatus(it_status);
      ATH_MSG_INFO("         Status Code=0x"
		   << MSG::hex << std::setfill( '0' ) << std::setw(8) << tmpstatus.code()
		   << " Generic Part=0x" << std::setw(4) << tmpstatus.generic()
		   << " Specific Part=0x" << std::setw(4) << tmpstatus.specific() << MSG::dec);
    }
  }

  // print list of subdetectors to filter out by status code
  ATH_MSG_INFO(" ---> Filter out Sub Detector ROBs by Status Code: # Sub Detectors = " << m_filterSubDetMap.size());
  for (const auto& it : m_filterSubDetMap) {
    eformat::helper::SourceIdentifier tmpsrc(it.first, 0);
    ATH_MSG_INFO("      SubDetId=0x" << MSG::hex << it.first << " -> " << tmpsrc.human_detector());
    for (auto it_status : it.second) {
      eformat::helper::Status tmpstatus(it_status);
      ATH_MSG_INFO("         Status Code=0x"
		   << MSG::hex << std::setfill( '0' ) << std::setw(8) << tmpstatus.code()
		   << " Generic Part=0x" << std::setw(4) << tmpstatus.generic()
		   << " Specific Part=0x" << std::setw(4) << tmpstatus.specific() << MSG::dec);
    }
  }

  // get the list of enabled ROBs from OKS
  bool robOKSconfigFound = false;

  if ( m_readROBfromOKS.value() ) {
    ServiceHandle<Gaudi::Interfaces::IOptionsSvc> jobOptionsSvc("JobOptionsSvc", name());
    if ((jobOptionsSvc.retrieve()).isFailure()) {
      ATH_MSG_ERROR("Could not find JobOptionsSvc");
    } else {
      if (jobOptionsSvc->has("DataFlowConfig.DF_Enabled_ROB_IDs") &&
          m_enabledROBs.fromString(jobOptionsSvc->get("DataFlowConfig.DF_Enabled_ROB_IDs")).isSuccess()) {
        robOKSconfigFound = true;
        ATH_MSG_INFO(" ---> Read from OKS                                                = "
                     << MSG::dec << m_enabledROBs.value().size() << " enabled ROB IDs.");
      } else {
        ATH_MSG_WARNING("Could not set Property 'enabledROBs' from OKS.");
      }
    }
  }

  // print list of enabled ROBs, read from OKS
  ATH_MSG_INFO(" ---> Read list of enabled ROBs from OKS                           = " << m_readROBfromOKS);
  if (m_enabledROBs.value().empty()) {
    ATH_MSG_INFO(" ---> The list of enabled ROBs has size                            = 0. No check will be performed ");
  } else {
    if (m_readROBfromOKS.value() && robOKSconfigFound) {
      ATH_MSG_INFO(" ---> The list of enabled ROBs has size                            = " << MSG::dec << m_enabledROBs.value().size() 
                   << ". It was read from the partition database." );
    } else {
      ATH_MSG_INFO(" ---> The list of enabled ROBs has size                            = " << MSG::dec << m_enabledROBs.value().size() 
                   << ". It was read from job options." );
    }
  }

  // prefetch all ROBs in a ROS on a first retrieval of ROBs from this ROS
  ATH_MSG_INFO(" ---> Prefetch all ROBs in a ROS on first retrieval                = " << m_prefetchAllROBsfromROS);

  if ( m_prefetchAllROBsfromROS.value() && !m_enabledROBs.value().empty() ) {
    m_prefetchWholeROSList.reserve( m_enabledROBs.value().size() );
    if ( m_prefetchSubDetROS.value().empty() || 
	 (std::find(m_prefetchSubDetROS.value().begin(),m_prefetchSubDetROS.value().end(),0) != m_prefetchSubDetROS.value().end())
       )
    {
      m_prefetchWholeROSList = m_enabledROBs.value() ;
      ATH_MSG_INFO("      All enabled ROBs are used for prefetching. ROB list size     = " << m_prefetchWholeROSList.size() );
    } else {
      ATH_MSG_INFO("      The following sub-detectors or sub-detector groups are configured for prefetching : Number = "
		   << m_prefetchSubDetROS.value().size() );
      for (uint8_t it : m_prefetchSubDetROS.value() ) {
	if (it > 0xf) {
	  eformat::SubDetector sd(static_cast<eformat::SubDetector>(it));
	  eformat::helper::SourceIdentifier tmpsrc(sd, 0);
	  ATH_MSG_INFO("         SubDetId    = 0x" << MSG::hex << std::setw(2) << (int)it << MSG::dec << " -> Group : " 
		       << tmpsrc.human_group() 
		       << " -> SubDetector : " << tmpsrc.human_detector()  );
	} else {
	  eformat::SubDetectorGroup sg(static_cast<eformat::SubDetectorGroup>(it));
	  ATH_MSG_INFO("         SubDetGroup = 0x" << MSG::hex << std::setw(1) << (int)it << MSG::dec << "  -> Group : " 
		       << eformat::helper::SubDetectorGroupDictionary.string(sg));
	}
      }
      // establish prefetch list
      m_prefetchWholeROSList.clear();
      for (auto rob : m_enabledROBs.value() ) {
	eformat::helper::SourceIdentifier tmpsrc(rob);
	if( (std::find(m_prefetchSubDetROS.value().begin(),m_prefetchSubDetROS.value().end(),tmpsrc.subdetector_id()) 
	     != m_prefetchSubDetROS.value().end()) ||
	    (std::find(m_prefetchSubDetROS.value().begin(),m_prefetchSubDetROS.value().end(),tmpsrc.subdetector_group()) 
	     != m_prefetchSubDetROS.value().end())
	    ) m_prefetchWholeROSList.push_back(rob);
      }
      ATH_MSG_INFO("      Number of ROBs which are used for prefetching = " << m_prefetchWholeROSList.size() );
    }
  }

  // Setup the slot specific cache
  m_eventsCache = SG::SlotSpecificObj<EventCache>( SG::getNSlots() );

  // Retrieve the monitoring tool
  if (!m_monTool.empty()) ATH_CHECK(m_monTool.retrieve());

  if (m_doCostMonitoring) ATH_CHECK(m_trigCostSvcHandle.retrieve());

  return(StatusCode::SUCCESS);
}

StatusCode HltROBDataProviderSvc::finalize()
{
  ATH_CHECK(m_monTool.release());
  return StatusCode::SUCCESS;
}

/// --- Implementation of IROBDataProviderSvc interface ---
/// --- Legacy interface (deprecated) ---

/// Signal ROB fragments which should be considered for prefetching in online running
void HltROBDataProviderSvc::addROBData(const std::vector<uint32_t>& robIds, 
					 const std::string_view callerName)
{
  const EventContext& context{ Gaudi::Hive::currentContext() };
  return addROBData( context, robIds, callerName );
}

/// Start a new event with a set of ROB fragments, e.g. from LVL1 result, in online and add the fragments to the ROB cache
void HltROBDataProviderSvc::setNextEvent(const std::vector<ROBF>& result)
{
  const EventContext& context{ Gaudi::Hive::currentContext() };
  return setNextEvent( context, result );
}

/// Start a new event with a full event fragment and add all ROB fragments in to the ROB cache
void HltROBDataProviderSvc::setNextEvent(const RawEvent* re)
{
  const EventContext& context{ Gaudi::Hive::currentContext() };
  return setNextEvent( context, re );
}

/// Retrieve ROB fragments for given ROB ids from the ROB cache
void HltROBDataProviderSvc::getROBData(const std::vector<uint32_t>& robIds, std::vector<const ROBF*>& robFragments, 
				       const std::string_view callerName)
{
  const EventContext& context{ Gaudi::Hive::currentContext() };
  return getROBData( context, robIds, robFragments, callerName );
}

/// Retrieve the full event fragment
const RawEvent* HltROBDataProviderSvc::getEvent() 
{
  const EventContext& context{ Gaudi::Hive::currentContext() };
  return getEvent( context );
}

/// Store the status for the event.
void HltROBDataProviderSvc::setEventStatus(uint32_t status)
{
  const EventContext& context{ Gaudi::Hive::currentContext() };
  setEventStatus( context, status );
}

/// Retrieve the status for the event.
uint32_t HltROBDataProviderSvc::getEventStatus() 
{
  const EventContext& context{ Gaudi::Hive::currentContext() };
  return getEventStatus( context );
}

/// --- Implementation of IROBDataProviderSvc interface ---
/// --- Context aware interface for MT ---

/// Signal ROB fragments which should be considered for prefetching in online running
void HltROBDataProviderSvc::addROBData(const EventContext& context, const std::vector<uint32_t>& robIds, 
				       const std::string_view callerName)
{
  ATH_MSG_VERBOSE("start of " << __FUNCTION__ << ": Number of ROB Ids to add = " << robIds.size() 
		  << " caller name = " << callerName);
  EventCache* cache = m_eventsCache.get( context );

  // allocate vector of missing ROB Ids
  std::vector<uint32_t> robIds_missing ;

  // allocate vector with existing ROB fragments in cache
  std::vector<const ROBF*> robFragments_inCache ;

  // check input ROB list against cache
  eventCache_checkRobListToCache(cache,robIds, robFragments_inCache, robIds_missing ) ;

  // call data collector
  if (!robIds_missing.empty()) {
    ATH_MSG_DEBUG( __FUNCTION__ << ": Number of ROB Ids to reserve with DCM = " << robIds_missing.size()); 
    // reserve the ROBs in the DCM
    try {  
      auto mon_robres_t = Monitored::Timer("TIME_ROBReserveData");
      hltinterface::DataCollector::instance()->reserveROBData(cache->globalEventNumber, robIds_missing);
      mon_robres_t.stop();
      // Fill monitoring histograms
      auto mon_robres_nROBs = Monitored::Scalar("NUMBER_ROBReserveData",robIds_missing.size());
      auto mon = Monitored::Group(m_monTool, mon_robres_t, mon_robres_nROBs);
    } catch (const std::exception& ex) {
      ATH_MSG_ERROR( __FUNCTION__ << ":" << __LINE__ 
		     << "Failed to reserve ROB data, caught an unexpected exception: " << ex.what());
    } catch (...) {
      ATH_MSG_ERROR( __FUNCTION__ << ":" << __LINE__ 
		     << "Failed to reserve ROB data, caught an unexpected exception."); 
    }
  }
}

/// Start a new event with a set of ROB fragments, e.g. from LVL1 result, in online and add the fragments to the ROB cache
void HltROBDataProviderSvc::setNextEvent(const EventContext& context, const std::vector<ROBF>& result)
{
  ATH_MSG_FATAL("Obsolete method HltROBDataProviderSvc::setNextEvent(const EventContext& context, const std::vector<ROBF>& result) called "
		<< "\n context = " << context << " number of ROB fragments = " << result.size() );
}

/// Start a new event with a full event fragment and add all ROB fragments in to the ROB cache
void HltROBDataProviderSvc::setNextEvent(const EventContext& context, const RawEvent* re)
{
  ATH_MSG_VERBOSE("start of " << __FUNCTION__);
  EventCache* cache = m_eventsCache.get( context );
  
  // clear the event cache of the previous event
  eventCache_clear( cache );

  //----------------------------------------------+
  // Fill the event cache with the new event data |
  //----------------------------------------------+

  // store the pointer to the event
  cache->event = re;
  // set the LVL1 id
  cache->currentLvl1ID = re->lvl1_id();
  // set the global event number
  cache->globalEventNumber = re->global_id();
  // set flag for masking L2/EF module ID, this is only necessary for the separate L2 and EF systems from Run 1 
  m_maskL2EFModuleID = (re->nlvl2_trigger_info() != 0);

  //--------------------+
  // Fill the ROB cache |
  //--------------------+

  // get all the ROBFragments
  OFFLINE_FRAGMENTS_NAMESPACE::PointerType robF[MAX_ROBFRAGMENTS];
  size_t number_robs = re->children(robF,MAX_ROBFRAGMENTS);
  if (number_robs == MAX_ROBFRAGMENTS) {
    ATH_MSG_ERROR("ROB buffer overflow: ROBs found = " << number_robs 
		  << " Max. number of ROBs allowed = " << MAX_ROBFRAGMENTS);
  }
  std::vector<ROBF> rob_fragments;
  rob_fragments.reserve(number_robs);
  // loop over all ROBs
  for (size_t irob = 0; irob < number_robs; irob++) {
    rob_fragments.push_back(ROBF(robF[irob]));
  }
  // add the ROBs to the cache/rob map
  eventCache_addRobData(cache, std::move(rob_fragments)) ;

  ATH_MSG_DEBUG(" ---> setNextEvent for                " << name() );
  ATH_MSG_DEBUG("      current [global id, LVL1 id] = [" << cache->globalEventNumber << "," << cache->currentLvl1ID << "]" );
  ATH_MSG_DEBUG("      number of received ROBs      =  " << rob_fragments.size() );
  ATH_MSG_DEBUG("      size of ROB cache            =  " << cache->robmap.size() );

  //------------------------------+
  // Initiate whole ROS retrieval |
  //------------------------------+
  if ( m_prefetchAllROBsfromROS.value() && !m_prefetchWholeROSList.empty() ) {
    addROBData( context, m_prefetchWholeROSList, "prefetch_HLTROBDataProviderSvc" );
    ATH_MSG_DEBUG("      ROS prefetch init. size      =  " << m_prefetchWholeROSList.size() );
  }

  return;
}

/// Retrieve ROB fragments for given ROB ids from the ROB cache
void HltROBDataProviderSvc::getROBData(const EventContext& context, 
				       const std::vector<uint32_t>& robIds, std::vector<const ROBF*>& robFragments, 
				       const std::string_view callerName)
{
  EventCache* cache = m_eventsCache.get( context );

  // lock for event cache update with DCM
  std::lock_guard<std::mutex> lock( cache->eventCache_mtx );

  TrigTimeStamp rosStartTime;

  ATH_MSG_VERBOSE("start of " << __FUNCTION__ << ": Number of ROB Ids to get = " << robIds.size() 
		  << " caller name = " << callerName);

  // allocate vector of missing ROB Ids
  std::vector<uint32_t> robIds_missing ;

  // check input ROB list against cache
  eventCache_checkRobListToCache(cache, robIds, robFragments, robIds_missing) ;

  // no missing ROB fragments, return the found ROB fragments 
  if (robIds_missing.empty()) {
    ATH_MSG_DEBUG( __FUNCTION__ << ": All requested ROB Ids were found in the cache. "); 
    return;
  }

  // There were missing ROB fragments retrieve them from the DCM and add them to the cache
  ATH_MSG_DEBUG( __FUNCTION__ << ": Number of ROB Ids to retrieve with DCM = " << robIds_missing.size()); 

  typedef std::vector<hltinterface::DCM_ROBInfo> ROBInfoVec;
  ROBInfoVec vRobInfos;

  // Get ROB Fragments with DataCollector
  vRobInfos.reserve( robIds_missing.size() ) ;
  try {
    auto mon_rob_t = Monitored::Timer("TIME_ROBRequest");
    hltinterface::DataCollector::instance()->collect(vRobInfos, cache->globalEventNumber, robIds_missing);
    mon_rob_t.stop();
    // Fill monitoring histograms
    auto mon_rob_nROBs = Monitored::Scalar("NUMBER_ROBRequest",vRobInfos.size());
    auto mon = Monitored::Group(m_monTool, mon_rob_t, mon_rob_nROBs);
  } catch (const std::exception& ex) {
    ATH_MSG_ERROR( __FUNCTION__ << ":" << __LINE__ 
		   << "Failed to collect ROBs, caught an unexpected exception: " << ex.what()
		   << ". Throwing hltonl::Exception::EventSourceCorrupted" );
    throw hltonl::Exception::EventSourceCorrupted();
  } catch (...) {
    ATH_MSG_ERROR( __FUNCTION__ << ":" << __LINE__ 
		   << "Failed to collect ROBs, caught an unexpected exception. "
		   << "Throwing hltonl::Exception::EventSourceCorrupted" );
    throw hltonl::Exception::EventSourceCorrupted();
  }

  // Store retrieved ROB data in the cache 
  std::vector<ROBF> robFragments_missing;
  robFragments_missing.reserve( vRobInfos.size() );
  for(ROBInfoVec::const_iterator it=vRobInfos.begin(); it!=vRobInfos.end(); ++it) {
    ATH_MSG_DEBUG(__FUNCTION__ << " ROB Id = 0x" << MSG::hex << it->robFragment.source_id() << MSG::dec
		    << " retrieved from DCM for (global Id, L1 Id) = (" << cache->globalEventNumber << "," << cache->currentLvl1ID <<")" );
    robFragments_missing.push_back( it->robFragment );
  }

  // Check if event should be monitored and monitor ROB data
  auto monitorData = robmonitor::ROBDataMonitorStruct(cache->currentLvl1ID, std::string(callerName));
  if (m_doCostMonitoring && m_trigCostSvcHandle->isMonitoredEvent(context, /*includeMultiSlot =*/ false)) {
    // Monitor HLT cached ROBs
    for (const ROBF* robFrag : robFragments) {
      monitorData.requested_ROBs[robFrag->source_id()] = robmap_getRobData(*robFrag, robmonitor::HLT_CACHED);
    }

    // Add the ROBs to the cache/rob map and collect ignored robs
    std::set<uint32_t> robIds_ignored;
    eventCache_addRobData(cache, std::move(robFragments_missing), robIds_ignored);

    // Monitor DCM ROBs
    for (const hltinterface::DCM_ROBInfo& robInfo : vRobInfos) {
      robmonitor::ROBHistory status;

      // Check ROB history
      if (robIds_ignored.find(robInfo.robFragment.source_id()) != robIds_ignored.end()) {
        status = robmonitor::IGNORED;
      }
      else {
        status = robInfo.robIsCached ? robmonitor::DCM_CACHED : robmonitor::RETRIEVED;
      }

      monitorData.requested_ROBs[robInfo.robFragment.source_id()] = robmap_getRobData(robInfo.robFragment, status);
    }

    // Return all the requested ROB fragments from the cache and collect disabled ROBs
    std::set<uint32_t> robIds_disabled;
    eventCache_checkRobListToCache(cache, robIds, robFragments, robIds_missing, robIds_disabled);

    // Fill undefined (not enabled) ROBs
    for (uint32_t robId : robIds_disabled) {
        monitorData.requested_ROBs[robId] = robmonitor::ROBDataStruct(robId);
        monitorData.requested_ROBs[robId].rob_history = robmonitor::UNDEFINED;
    }
  }
  else {
    // add the ROBs to the cache/rob map
    eventCache_addRobData(cache, std::move(robFragments_missing)) ;

    // return all the requested ROB fragments from the cache
    eventCache_checkRobListToCache(cache, robIds, robFragments, robIds_missing) ;
  }

  // Save ROS processing time and pass ROS data to CostMonitor
  if (m_doCostMonitoring && m_trigCostSvcHandle->isMonitoredEvent(context, /*includeMultiSlot =*/ false)) {
    TrigTimeStamp rosEndTime;

    monitorData.start_time = rosStartTime.microsecondsSinceEpoch();
    monitorData.end_time = rosEndTime.microsecondsSinceEpoch();

    // Check if ROBMonitorDataStruct is move-constructible
    static_assert(std::is_nothrow_move_constructible<robmonitor::ROBDataMonitorStruct>::value);
    if (m_trigCostSvcHandle->monitorROS(context, std::move(monitorData)).isFailure()) {
      ATH_MSG_WARNING("TrigCost ROS monitoring failed!");
    }  
  }
}

robmonitor::ROBDataStruct HltROBDataProviderSvc::robmap_getRobData(const ROBF& robFrag, robmonitor::ROBHistory robStatus)
{
  auto robData = robmonitor::ROBDataStruct(robFrag.source_id());
  robData.rob_size = robFrag.fragment_size_word();
  robData.rob_status_word = robFrag.nstatus() ? robFrag.status()[0] : 0;
  robData.rob_history = robStatus;

  return robData;
}

/// Retrieve the full event fragment
const RawEvent* HltROBDataProviderSvc::getEvent(const EventContext& context)
{
  ATH_MSG_VERBOSE("start of " << __FUNCTION__);
  return m_eventsCache.get( context )->event;
}

/// Store the status for the event.
void HltROBDataProviderSvc::setEventStatus(const EventContext& context, uint32_t status)
{
  ATH_MSG_VERBOSE("start of " << __FUNCTION__);
  m_eventsCache.get(context)->eventStatus = status;
}

/// Retrieve the status for the event.
uint32_t HltROBDataProviderSvc::getEventStatus(const EventContext& context)
{
  ATH_MSG_VERBOSE("start of " << __FUNCTION__);
  return m_eventsCache.get( context )->eventStatus;
}

/// Apply a function to all ROBs in the cache
void HltROBDataProviderSvc::processCachedROBs(const EventContext& context, const std::function< void(const ROBF* )>& fn) const
{
  ATH_MSG_VERBOSE("start of " << __FUNCTION__);
  for ( const auto& el : m_eventsCache.get( context )->robmap )
    {
      fn( &el.second );
    }
}

/// Flag to check if all event data have been retrieved
bool HltROBDataProviderSvc::isEventComplete(const EventContext& context) const
{
  ATH_MSG_VERBOSE("start of " << __FUNCTION__);
  return m_eventsCache.get( context )->isEventComplete;
}

/// retrieve in online running all ROBs for the event from the readout system. Only those ROBs are retrieved which are not already in the cache
int HltROBDataProviderSvc::collectCompleteEventData(const EventContext& context, const std::string_view callerName)
{
  ATH_MSG_VERBOSE("start of " << __FUNCTION__  << " caller name = " << callerName);

  EventCache* cache = m_eventsCache.get( context );

  // lock for event cache update with DCM
  std::lock_guard<std::mutex> lock( cache->eventCache_mtx );

  // return if event is already complete 
  if (cache->isEventComplete) return 0;

  typedef std::vector<hltinterface::DCM_ROBInfo> ROBInfoVec;
  ROBInfoVec vRobInfos ;
  if (!m_enabledROBs.value().empty()) {
    vRobInfos.reserve( m_enabledROBs.value().size() ) ;
  } else {
    vRobInfos.reserve( MAX_ROBFRAGMENTS ) ;
  }

  // Get ROB Fragments for complete event with DataCollector
  try {
    auto mon_col_t = Monitored::Timer("TIME_CollectAllROBs");
    hltinterface::DataCollector::instance()->collect(vRobInfos, cache->globalEventNumber);
    mon_col_t.stop();
    ATH_MSG_DEBUG( __FUNCTION__ << ": Number of received ROB Ids = " << vRobInfos.size() );
    // Fill monitoring histograms
    auto mon_col_nROBs = Monitored::Scalar("NUMBER_CollectAllROBs",vRobInfos.size());
    auto mon = Monitored::Group(m_monTool, mon_col_t, mon_col_nROBs);
  } catch (const std::exception& ex) {
    ATH_MSG_ERROR( __FUNCTION__ << ":" << __LINE__ 
		   << "Failed to collect complete event, caught an unexpected exception: " << ex.what()
		   << ". Throwing hltonl::Exception::EventSourceCorrupted" );
    throw hltonl::Exception::EventSourceCorrupted();
  } catch (...) {
    ATH_MSG_ERROR( __FUNCTION__ << ":" << __LINE__ 
		   << "Failed to collect complete event, caught an unexpected exception. "
		   << "Throwing hltonl::Exception::EventSourceCorrupted" );
    throw hltonl::Exception::EventSourceCorrupted();
  }

  // Store retrieved ROB data in the cache 
  std::vector<ROBF> robFragments_missing;
  robFragments_missing.reserve( vRobInfos.size() );
  for(ROBInfoVec::const_iterator it=vRobInfos.begin(); it!=vRobInfos.end(); ++it) {
    ATH_MSG_VERBOSE(__FUNCTION__ << " ROB Id = 0x" << MSG::hex << it->robFragment.source_id() << MSG::dec
		    << " retrieved from DCM for (global Id, L1 Id) = (" << cache->globalEventNumber << "," << cache->currentLvl1ID <<")" );
    robFragments_missing.push_back( it->robFragment );
  }
  // add the ROBs to the cache/rob map
  eventCache_addRobData(cache, std::move(robFragments_missing)) ;

  // update event complete flag
  cache->isEventComplete = true;

  return vRobInfos.size();
}

/// method to filter ROBs with given Status code
bool HltROBDataProviderSvc::robmap_filterRobWithStatus(const ROBF* rob)
{
  // No filter criteria defined
  if (m_filterRobMap.empty() && m_filterSubDetMap.empty()) {
    return(false);
  }

  // There should be at least one status element if there was an error
  // in case there are 0 status elements then there was no known error
  // (see event format document ATL-D-ES-0019 (EDMS))
  const uint32_t* rob_it_status;
  const uint32_t null_status(0);
  // The ROB has no status elements
  if (rob->nstatus() == 0) {
    rob_it_status = &null_status;
  } else {
  // The ROB has at least one status element, access it via an iterator
    rob->status(rob_it_status);
  }  

  // Build the full ROB Sourceidentifier
  eformat::helper::SourceIdentifier tmpsrc(rob->rob_source_id());
  
  // Check if there is a ROB specific filter rule defined for this ROB Id and match the status code
  FilterRobMap::iterator map_it_rob = m_filterRobMap.find(tmpsrc.code());
  if (map_it_rob != m_filterRobMap.end()) {
    for (auto it_status: (*map_it_rob).second) {
      if (*rob_it_status == it_status) {
	ATH_MSG_VERBOSE(__FUNCTION__ << " ROB Id = 0x" << MSG::hex << tmpsrc.code() 
			<< " with status = 0x" << *rob_it_status  << MSG::dec
			<< " removed due to ROB filter rule.");
	return(true);
      }
    }
  }

  // Check if there is a sub detector specific filter rule defined for this ROB Id and match the status code
  FilterSubDetMap::iterator map_it_subdet = m_filterSubDetMap.find(tmpsrc.subdetector_id());
  if (map_it_subdet != m_filterSubDetMap.end()) {
    for (auto it_status: (*map_it_subdet).second) {
      if (*rob_it_status == it_status) {
	ATH_MSG_VERBOSE(__FUNCTION__ << " ROB Id = 0x" << MSG::hex << tmpsrc.code() 
			<< " with status = 0x" << *rob_it_status  << MSG::dec
			<< " removed due to SubDet filter rule.");
	return(true);
      }
    }
  }
  return(false);
}

void HltROBDataProviderSvc::eventCache_clear(EventCache* cache)
{
  ATH_MSG_VERBOSE("start of " << __FUNCTION__);
  cache->event             = nullptr;
  cache->currentLvl1ID     = 0; 
  cache->globalEventNumber = 0;
  cache->eventStatus       = 0;    
  cache->isEventComplete   = false;    
  { cache->robmap.clear(); }
}

void HltROBDataProviderSvc::eventCache_checkRobListToCache(EventCache* cache, const std::vector<uint32_t>& robIds_toCheck, 
							     std::vector<const ROBF*>& robFragments_inCache, 
							     std::vector<uint32_t>& robIds_missing,
                   std::optional<std::reference_wrapper<std::set<uint32_t>>> robIds_disabled )
{
  ATH_MSG_VERBOSE("start of " << __FUNCTION__ << " number of ROB Ids to check = " << robIds_toCheck.size());

  // clear output arrays
  robFragments_inCache.clear();
  robIds_missing.clear();

  // allocate sufficient space for output arrays
  robFragments_inCache.reserve( robIds_toCheck.size() );
  robIds_missing.reserve( robIds_toCheck.size() );

  // check input ROB ids
  for (uint32_t id : robIds_toCheck) {

    // check for duplicate IDs on the list of missing ROBs
    std::vector<uint32_t>::iterator missing_it = std::find(robIds_missing.begin(), robIds_missing.end(), id);
    if (missing_it != robIds_missing.end()) {
      ATH_MSG_VERBOSE(__FUNCTION__ << " ROB Id : 0x" << MSG::hex << id << MSG::dec <<" is already on the list of missing IDs.");
      continue;
    }

    // check if ROB is already in cache
    { ROBMAP::const_iterator map_it = cache->robmap.find(id);
      if (map_it != cache->robmap.end()) {
	ATH_MSG_VERBOSE(__FUNCTION__ << " ROB Id 0x" << MSG::hex << id << MSG::dec
			<< " found for (global Id, L1 Id) = (" << cache->globalEventNumber << "," << cache->currentLvl1ID <<")" );
	robFragments_inCache.push_back( &(map_it->second) );
	continue;
      }
    }

    // check if ROB is actually enabled for readout
    if (!m_enabledROBs.value().empty()) {
      std::vector<uint32_t>::const_iterator rob_enabled_it = 
        std::find(m_enabledROBs.value().begin(), m_enabledROBs.value().end(),id);
      if(rob_enabled_it == m_enabledROBs.value().end()) {
        ATH_MSG_VERBOSE(__FUNCTION__ << " ROB Id : 0x" << MSG::hex << id << MSG::dec
                << " will be not added, since it is not on the list of enabled ROBs.");
        if (robIds_disabled) {
          robIds_disabled->get().insert(id);
        }
        continue;
      }
    }

    // the ROB is not in the cache and should be eventually added
    ATH_MSG_VERBOSE(__FUNCTION__ << " ROB Id : 0x" << MSG::hex << id << MSG::dec <<" is missing ");
    robIds_missing.push_back( id ) ;
  } // end loop over input ROB Ids to check
}

void HltROBDataProviderSvc::eventCache_addRobData(EventCache* cache, std::vector<ROBF>&& robFragments,
              std::optional<std::reference_wrapper<std::set<uint32_t>>> robIds_ignored)
{
  ATH_MSG_VERBOSE("start of " << __FUNCTION__ << " number of ROB fragments to add = " << robFragments.size());

  for (const ROBF& rob : robFragments) {

    // Source ID
    uint32_t id = rob.source_id();
    ATH_MSG_VERBOSE(__FUNCTION__ << " Id = 0x" << std::hex << id << std::dec );

    // mask off the module ID for L2 and EF result for Run 1 data
    if ( (eformat::helper::SourceIdentifier(id).module_id() != 0) &&
	 (eformat::helper::SourceIdentifier(id).subdetector_id() == eformat::TDAQ_LVL2) ) {
      id = eformat::helper::SourceIdentifier(eformat::helper::SourceIdentifier(id).subdetector_id(),0).code();
      if (!m_maskL2EFModuleID) {
	ATH_MSG_ERROR(__FUNCTION__ << " Inconsistent flag for masking L2/EF module IDs");
	m_maskL2EFModuleID=true;
      }
    } else if ( (eformat::helper::SourceIdentifier(id).module_id() != 0) && 
		(eformat::helper::SourceIdentifier(id).subdetector_id() == eformat::TDAQ_EVENT_FILTER) &&
		(m_maskL2EFModuleID) ) {
      id = eformat::helper::SourceIdentifier(eformat::helper::SourceIdentifier(id).subdetector_id(),0).code();
    }

    // check if ROB is already in cache
    { ROBMAP::const_iterator it = cache->robmap.find(id);
      if (it != cache->robmap.end()) {
	ATH_MSG_VERBOSE(__FUNCTION__ << " Duplicate ROB Id 0x" << MSG::hex << id << MSG::dec
			<< " found for (global Id, L1 Id) = (" << cache->globalEventNumber << "," << cache->currentLvl1ID <<")" );
	continue;
      }
    }

    // filter ROBs with external criteria 
    if (robmap_filterRobWithStatus(&rob)) {
      if (rob.nstatus() > 0) {
	const uint32_t* it_status;
	rob.status(it_status);
	eformat::helper::Status tmpstatus(*it_status);
	ATH_MSG_VERBOSE(__FUNCTION__ << " ROB Id = 0x" << MSG::hex << id << std::setfill('0')
			<< " with Generic Status Code = 0x" << std::setw(4) << tmpstatus.generic()
			<< " and Specific Status Code = 0x" << std::setw(4) << tmpstatus.specific() << MSG::dec
			<< " removed for (global Id, L1 Id) = (" << cache->globalEventNumber << "," << cache->currentLvl1ID <<")" );
      }
      if (robIds_ignored) {
        robIds_ignored->get().insert(id);
      }
      continue;
    }

    // check for ROBs with no data 
    if ((rob.rod_ndata() == 0) && (m_filterEmptyROB)) {
      ATH_MSG_VERBOSE(__FUNCTION__ << " Empty ROB Id = 0x" << MSG::hex << id << MSG::dec
		      << " removed for (global Id, L1 Id) = (" << cache->globalEventNumber << "," << cache->currentLvl1ID <<")" );
      continue;
    } 

    // add ROB to map
    { cache->robmap.insert(std::make_pair(id,std::move(rob))); }
  }
}

HltROBDataProviderSvc::EventCache::~EventCache()
{
  //  delete event;
  { robmap.clear(); }
}
