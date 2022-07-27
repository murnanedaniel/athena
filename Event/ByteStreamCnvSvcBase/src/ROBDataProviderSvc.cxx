/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

//===================================================================
//  Implementation of ROBDataProviderSvc
//  Revision: November 2017
//      MT readiness
//  Revision:  July 11, 2002
//      Modified for eformat
//  Revision:  Aug 18, 2003
//      Modified to use ROBFragments directly and include methods
//      for online
//  Revision:  Apr 21, 2005
//      Remove dependency on Level-2 Data Collector, create special
//      version for online
//  Revision:  Oct 29, 2006
//      Increase MAX_ROBFRAGMENTS to 2048 to cover the special case
//      when in a "localhost" partition the complete event is given
//      to the Event Filter as one single ROS fragment (this case
//      should not happen for normal running when several ROSes are
//      used and the ROB fragments are grouped under different ROS
//      fragments)
//  Revision:  Nov 10, 2008
//      Mask off the module ID from the ROB source identifier of the
//      L2 and EF result when storing it in the ROB map. This is necessary
//      when the L2/PT node ID is stored in the source identifier as
//      module ID. With this modification the L2 and EF result can still be
//      found as 0x7b0000 and 0x7c0000
//  Revision:  Jan 12, 2009
//      Allow removal of individual ROBs and ROBs from given subdetectors
//      from the internal ROB map according to a given status code.
//      This may be necessary when corrupted and incomplete ROB fragments
//      are forwarded to the algorithms and the converters are not yet
//      prepared to handle the specific cases.
//      The filtering can be configured with job options as:
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
//  Revision:  Jan 28, 2014
//      For Run 1 the module ID from the ROB source identifier of the
//      L2 and EF result needed to be masked off before storing the ROB fragment 
//      in the ROB map. The module ID for these fragments contained an identifier 
//      of the machine on which they were produced. This produced as many different ROB IDs
//      for these fragments as HLT processors were used. The module IDs were not useful 
//      for analysis and needed to be masked off from these fragments in the ROB map in order 
//      to allow the access to the L2 or the EF result with the generic identifiers  
//      0x7b0000 and 0x7c0000. Also an event contained only one L2 and EF result.
//      From Run 2 on (eformat version 5) the HLT processor identifier is not anymore
//      stored in the module ID of the HLT result. Further there can be in one event several HLT result
//      records with the source identifier 0x7c. The different HLT results are distinguished
//      now with the module ID. A module ID 0 indicates a physiscs HLT result as before, while
//      HLT results with module IDs different from zero are produced by data scouting chains. 
//      In Run 2 the module ID should be therefore not any more masked.
//      The masking of the moduleID is switched on when a L2 result is found in the event or the
//      event header contains L2 trigger info words. This means the data were produced with run 1 HLT system.
//
//===================================================================

// Include files.
#include "ByteStreamCnvSvcBase/ROBDataProviderSvc.h"
#include "eformat/Status.h"

// Constructor.
ROBDataProviderSvc::ROBDataProviderSvc(const std::string& name, ISvcLocator* svcloc) 
  : base_class(name, svcloc) {

  declareProperty("filterRobWithStatus", m_filterRobWithStatus);
  declareProperty("filterSubDetWithStatus", m_filterSubDetWithStatus);
  declareProperty("filterEmptyROB", m_filterEmptyROB = false);
}


// Initialization
StatusCode ROBDataProviderSvc::initialize() {
   ATH_MSG_INFO("Initializing");
   m_eventsCache = SG::SlotSpecificObj<EventCache>( SG::getNSlots() );
   
   for (unsigned int i = 0; i < m_filterRobWithStatus.value().size(); i++) {
      eformat::helper::SourceIdentifier tmpsrc(m_filterRobWithStatus.value()[i].first);
      if (tmpsrc.human_detector() != "UNKNOWN") {
         m_filterRobMap[tmpsrc.code()].push_back(m_filterRobWithStatus.value()[i].second);
      }
   }
   for (unsigned int i = 0; i < m_filterSubDetWithStatus.value().size(); i++) {
      eformat::helper::SourceIdentifier tmpsrc((eformat::SubDetector)m_filterSubDetWithStatus.value()[i].first, 0);
      if (tmpsrc.human_detector() != "UNKNOWN") {
         m_filterSubDetMap[tmpsrc.subdetector_id()].push_back(m_filterSubDetWithStatus.value()[i].second);
      }
   }
   ATH_MSG_INFO(" ---> Filter out empty ROB fragments                               = " << m_filterEmptyROB);
   ATH_MSG_INFO(" ---> Filter out specific ROBs by Status Code: # ROBs = " << m_filterRobMap.size());
   for (const auto& p : m_filterRobMap) {
      eformat::helper::SourceIdentifier tmpsrc(p.first);
      ATH_MSG_INFO("      RobId=0x" << MSG::hex << p.first << " -> in Sub Det = " << tmpsrc.human_detector());
      for (uint32_t status : p.second) {
         eformat::helper::Status tmpstatus(status);
         ATH_MSG_INFO("         Status Code=0x"
	         << MSG::hex << std::setfill( '0' ) << std::setw(8) << tmpstatus.code()
	         << " Generic Part=0x" << std::setw(4) << tmpstatus.generic()
	         << " Specific Part=0x" << std::setw(4) << tmpstatus.specific());
      }
   }

   ATH_MSG_INFO(" ---> Filter out Sub Detector ROBs by Status Code: # Sub Detectors = " << m_filterSubDetMap.size());
   for (const auto& p : m_filterSubDetMap) {
      eformat::helper::SourceIdentifier tmpsrc(p.first, 0);
      ATH_MSG_INFO("      SubDetId=0x" << MSG::hex << p.first << " -> " << tmpsrc.human_detector());
      for (uint32_t status : p.second) {
         eformat::helper::Status tmpstatus(status);
         ATH_MSG_INFO("         Status Code=0x"
	         << MSG::hex << std::setfill( '0' ) << std::setw(8) << tmpstatus.code()
	         << " Generic Part=0x" << std::setw(4) << tmpstatus.generic()
	         << " Specific Part=0x" << std::setw(4) << tmpstatus.specific());
      }
   }
   return(StatusCode::SUCCESS);
}


/**
    - in offline only check that given ROB ids are in the map, issue an
      error if not
*/
void ROBDataProviderSvc::addROBData(const std::vector<uint32_t>& robIds, const std::string_view callerName) {
  const EventContext context{ Gaudi::Hive::currentContext() };
  return addROBData( context, robIds, callerName );
}

void ROBDataProviderSvc::addROBData(const EventContext& context, const std::vector<uint32_t>& robIds, const std::string_view callerName) {
    EventCache* cache = m_eventsCache.get( context );

   // Copy missing ROB ids to vector with pthread allocator
   ATH_MSG_DEBUG(" ---> Number of ROB Id s requested : " << robIds.size() << ", Caller Name = " << callerName);
   // for offline running all requested ROBs should be found in cache
   // if not issue error
   for (uint32_t id : robIds) {
      // mask off the module ID for L2 and EF result for Run 1 data
      if ( (eformat::helper::SourceIdentifier(id).module_id() != 0) &&
	   (eformat::helper::SourceIdentifier(id).subdetector_id() == eformat::TDAQ_LVL2) ) {
	 id = eformat::helper::SourceIdentifier(eformat::helper::SourceIdentifier(id).subdetector_id(),0).code();
	 // TB if it is inconsistent we should not continue like this?
	 if ( !m_maskL2EFModuleID ) {
	   ATH_MSG_ERROR("Inconsistent flag for masking L2/EF module IDs");
	   m_maskL2EFModuleID=true;
	 }
      } else if ( (eformat::helper::SourceIdentifier(id).module_id() != 0) && 
		  (eformat::helper::SourceIdentifier(id).subdetector_id() == eformat::TDAQ_EVENT_FILTER) &&
		  ( m_maskL2EFModuleID ) ) {
	 id = eformat::helper::SourceIdentifier(eformat::helper::SourceIdentifier(id).subdetector_id(),0).code();
      }
      ROBMAP& robmap( cache->robmap );
      ROBMAP::iterator map_it = robmap.find(id) ;
      if (map_it != robmap.end()) {
         ATH_MSG_DEBUG(" ---> Found   ROB Id : 0x" << MSG::hex << (*map_it).second->source_id()
	         << MSG::dec << " in cache");
      } else {
         ATH_MSG_DEBUG(" ---> ROB Id : 0x" << MSG::hex << id
	         << MSG::dec << " not found in cache for running mode OFFLINE (method addROBData),");
	 ATH_MSG_DEBUG("      Lvl1 id  = " << cache->currentLvl1ID);
    }
  }
  return;
}
/** - this is the online method to add the LVL1/LVL2 result
    - this version of ROBDataProviderSvc does not support it
    - this version is for offline use only
*/
void ROBDataProviderSvc::setNextEvent(const std::vector<ROBF>& result) {
  const EventContext context{ Gaudi::Hive::currentContext() };
  return setNextEvent( context, result );
}
void ROBDataProviderSvc::setNextEvent(const EventContext& /*context*/, const std::vector<ROBF>& result) { 
  // clear the old map
  // TB honestly, why do any action if this is FATAL mistake
  //  robmapClear( m_eventsCache.get(context)->robmap );

   // This method should never be used by offline
   ATH_MSG_FATAL(" +-----------------------------------------------------------------+ ");
   ATH_MSG_FATAL(" | The method ROBDataProviderSvc::setNextEvent(const ROBF* result) | ");
   ATH_MSG_FATAL(" |    is not implemented for this version of ROBDataProviderSvc    | ");
   ATH_MSG_FATAL(" |      Use the version from the HLT repository if you need it.    | ");
   ATH_MSG_FATAL(" +-----------------------------------------------------------------+ ");
   ATH_MSG_FATAL(" ---> The " << result.size() << " ROB fragments in the call will not be used.");
   return;
}




/** - add a new Raw event
    - rebuild the map
*/
void ROBDataProviderSvc::setNextEvent(const RawEvent* re) {
  // obtain context and redirect to the real implementation
  const EventContext context{ Gaudi::Hive::currentContext() };
  return setNextEvent( context, re );
}

void ROBDataProviderSvc::setNextEvent( const EventContext& context, const RawEvent* re ) {
  EventCache* cache = m_eventsCache.get( context );
  
   cache->event=re;
   // clear the old map
   robmapClear( cache->robmap );
   // set the LVL1 id
   cache->currentLvl1ID = re->lvl1_id();
   // set flag for masking L2/EF module ID, this is only necessary for the separate L2 and EF systems from Run 1 
   m_maskL2EFModuleID = (re->nlvl2_trigger_info() != 0);

   // get all the ROBFragments
   const size_t MAX_ROBFRAGMENTS = 4096;
   OFFLINE_FRAGMENTS_NAMESPACE::PointerType robF[MAX_ROBFRAGMENTS];
   OFFLINE_FRAGMENTS_NAMESPACE::PointerType rePointer;
   re->start(rePointer);
   size_t robcount = re->children(robF, MAX_ROBFRAGMENTS);
   if (robcount == MAX_ROBFRAGMENTS) {
      ATH_MSG_ERROR("ROB buffer overflow");
   }
   // loop over all ROBs
   for (size_t irob = 0; irob < robcount; irob++) {
      // add to the map
      //const ROBF* rob = new ROBF(robF[irob]);
      std::unique_ptr<const ROBF> rob=std::make_unique<const ROBF>(robF[irob]);
      uint32_t id =  rob->source_id();
      // mask off the module ID for L2 and EF result for Run 1 data
      if ( (eformat::helper::SourceIdentifier(id).module_id() != 0) &&
	   (eformat::helper::SourceIdentifier(id).subdetector_id() == eformat::TDAQ_LVL2) ) {
	 id = eformat::helper::SourceIdentifier(eformat::helper::SourceIdentifier(id).subdetector_id(),0).code();
	 if (!m_maskL2EFModuleID) {
	   ATH_MSG_ERROR("Inconsistent flag for masking L2/EF module IDs");
	   m_maskL2EFModuleID=true;
	 }
      } else if ( (eformat::helper::SourceIdentifier(id).module_id() != 0) && 
		  (eformat::helper::SourceIdentifier(id).subdetector_id() == eformat::TDAQ_EVENT_FILTER) &&
		  (m_maskL2EFModuleID) ) {
	 id = eformat::helper::SourceIdentifier(eformat::helper::SourceIdentifier(id).subdetector_id(),0).code();
      }
      if (filterRobWithStatus(rob.get())) {
         if (rob->nstatus() > 0) {
            const uint32_t* it_status;
            rob->status(it_status);
            eformat::helper::Status tmpstatus(*it_status);
            ATH_MSG_DEBUG(" ---> ROB Id = 0x" << MSG::hex << id << std::setfill('0')
	            << " with Generic Status Code = 0x" << std::setw(4) << tmpstatus.generic()
	            << " and Specific Status Code = 0x" << std::setw(4) << tmpstatus.specific() << MSG::dec
	            << " removed for L1 Id = " << cache->currentLvl1ID);
         }
         rob.reset();
      } else if ((rob->rod_ndata() == 0) && (m_filterEmptyROB)) {
         ATH_MSG_DEBUG( " ---> Empty ROB Id = 0x" << MSG::hex << id << MSG::dec
	         << " removed for L1 Id = " << cache->currentLvl1ID);
          rob.reset();
      } else {
         ROBMAP::const_iterator it = cache->robmap.find(id);
         if (it != cache->robmap.end()) {
            ATH_MSG_WARNING(" ROBDataProviderSvc:: Duplicate ROBID 0x" << MSG::hex << id
	            << " found. " << MSG::dec << " Overwriting the previous one ");
         } 
         cache->robmap[id]=std::move(rob);
         
      }
   }
   ATH_MSG_DEBUG(" ---> setNextEvent offline for " << name() );
   ATH_MSG_DEBUG("      current LVL1 id   = " << cache->currentLvl1ID );
   ATH_MSG_DEBUG("      size of ROB cache = " << cache->robmap.size() );
   return;
}
/** return ROBData for ROBID
 */
void ROBDataProviderSvc::getROBData(const std::vector<uint32_t>& ids, std::vector<const ROBF*>& v, const std::string_view callerName) {
  const EventContext context{ Gaudi::Hive::currentContext() };
  return getROBData( context, ids, v, callerName );
}

void ROBDataProviderSvc::getROBData(const EventContext& context, const std::vector<uint32_t>& ids, std::vector<const ROBF*>& v, 
				    const std::string_view callerName) {
  EventCache* cache = m_eventsCache.get( context );

   for (uint32_t id : ids) {
      // mask off the module ID for L2 and EF result for Run 1 data
      if ( (eformat::helper::SourceIdentifier(id).module_id() != 0) &&
	   (eformat::helper::SourceIdentifier(id).subdetector_id() == eformat::TDAQ_LVL2) ) {
	 id = eformat::helper::SourceIdentifier(eformat::helper::SourceIdentifier(id).subdetector_id(),0).code();
	 if (!m_maskL2EFModuleID) {
	   ATH_MSG_ERROR("Inconsistent flag for masking L2/EF module IDs");
	   m_maskL2EFModuleID=true;
	 }
      } else if ( (eformat::helper::SourceIdentifier(id).module_id() != 0) && 
		  (eformat::helper::SourceIdentifier(id).subdetector_id() == eformat::TDAQ_EVENT_FILTER) &&
		  (m_maskL2EFModuleID) ) {
	 id = eformat::helper::SourceIdentifier(eformat::helper::SourceIdentifier(id).subdetector_id(),0).code();
      }
      ROBMAP::iterator map_it = cache->robmap.find(id);
      if (map_it != cache->robmap.end()) {
         v.push_back((*map_it).second.get());
      } else {
	ATH_MSG_DEBUG("Failed to find ROB for id 0x" << MSG::hex << id << MSG::dec << ", Caller Name = " << callerName);
#ifndef NDEBUG
         int nrob = 0;
         ATH_MSG_VERBOSE(" --- Dump of ROB cache ids --- total size = " << cache->robmap.size());
         for (const auto& p : cache->robmap) {
	    ++nrob;
	    ATH_MSG_VERBOSE(" # = " << nrob << "  id = 0x" << MSG::hex << p.second->source_id() << MSG::dec);
         }
#endif
      }
   }
   return;
}
/** - clear ROB map
 */
void ROBDataProviderSvc::robmapClear( ROBMAP& toclear) {
   for (auto& it : toclear) {
     it.second.reset();  
  }
  toclear.clear();
}
/// Retrieve the whole event.
const RawEvent* ROBDataProviderSvc::getEvent() {
  const EventContext context{ Gaudi::Hive::currentContext() };
  return getEvent( context );
}
const RawEvent* ROBDataProviderSvc::getEvent( const EventContext& context ) {
  
  return m_eventsCache.get( context )->event;
}


/// Set the status for the event.
void ROBDataProviderSvc::setEventStatus(uint32_t status) {
  const EventContext context{ Gaudi::Hive::currentContext() };
  setEventStatus( context, status );
}

void ROBDataProviderSvc::setEventStatus(const EventContext& context, uint32_t status) {
  m_eventsCache.get(context)->eventStatus = status;
}
/// Retrieve the status for the event.
uint32_t ROBDataProviderSvc::getEventStatus() {
  const EventContext context{ Gaudi::Hive::currentContext() };
  return getEventStatus( context );
}

uint32_t ROBDataProviderSvc::getEventStatus( const EventContext& context ) {
  return m_eventsCache.get( context )->eventStatus;
}

void ROBDataProviderSvc::processCachedROBs(const EventContext& context, 
					   const std::function< void(const ROBF* )>& fn ) const {
  for ( const auto&  el : m_eventsCache.get( context )->robmap ) {
    fn( el.second.get() );
  }
}



/** - filter ROB with Sub Detector Id and Status Code
*/
bool ROBDataProviderSvc::filterRobWithStatus(const ROBF* rob) {
   // No filter criteria defined
   if ((m_filterRobMap.size() == 0) && (m_filterSubDetMap.size() == 0)) {
      return(false);
   }
   // There should be at least one status element if there was an error
   // in case there are 0 status elements then there was no known error
   // (see event format document ATL-D-ES-0019 (EDMS))
   if (rob->nstatus() == 0) {
      return(false);
   }
   // The ROB has at least one status element, access it via an iterator
   const uint32_t* rob_it_status;
   rob->status(rob_it_status);
   // Build the full ROB Sourceidentifier
   eformat::helper::SourceIdentifier tmpsrc(rob->rob_source_id());
   // Check if there is a ROB specific filter rule defined for this ROB Id and match the status code
   FilterRobMap::iterator map_it_rob = m_filterRobMap.find(tmpsrc.code());
   if (map_it_rob != m_filterRobMap.end()) {
      for (uint32_t status : map_it_rob->second) {
         if (*rob_it_status == status) {
            return(true);
         }
      }
   }
   // Check if there is a sub detector specific filter rule defined for this ROB Id and match the status code
   FilterSubDetMap::iterator map_it_subdet = m_filterSubDetMap.find(tmpsrc.subdetector_id());
   if (map_it_subdet != m_filterSubDetMap.end()) {
      for (uint32_t status : map_it_subdet->second) {
         if (*rob_it_status == status) {
            return(true);
         }
      }
   }
   return(false);
}
