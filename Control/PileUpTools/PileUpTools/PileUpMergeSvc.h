/* -*- C++ -*- */

/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

/** @file PileUpMergeSvc.h
  @brief the preferred mechanism to access information from the different event
  stores in a pileup job.
  @author Paolo Calafiura
  $Id: PileUpMergeSvc.h,v 1.17 2008-04-19 00:31:09 calaf Exp $
*/
#ifndef PILEUPTOOLS_PILEUPMERGESVC_H
#define PILEUPTOOLS_PILEUPMERGESVC_H

#include "EventInfo/PileUpTimeEventIndex.h"
#include "xAODEventInfo/EventInfo.h"
#include "PileUpTools/IPileUpTool.h"    /* subEventIterator*/
#include "xAODCnvInterfaces/IEventInfoCnvTool.h"

#include "AthenaBaseComps/AthService.h"
#include "GaudiKernel/ClassID.h"
#include "Gaudi/Property.h"   /*StringArrayProperty*/
#include "GaudiKernel/StatusCode.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"

#include "AthLinks/DataLink.h"
#include "AthenaKernel/SlotSpecificObj.h"

#include "PileUpTools/IPileUpXingFolder.h"
#include "CxxUtils/checker_macros.h"

#include <cassert>

#include <list>
#include <map>
#include <utility> /*std::pair*/
#include <mutex>

class ISvcLocator;
class StoreGateSvc;
class ITriggerTime;
class IToolSvc;


/** @class PileUpMergeSvc
  @brief the preferred mechanism to access information from the different event
  stores in a pileup job.

  @details access PileUpEventInfo to locate data objects (typically hit 
  collections) to be merged attaching sub-evt time offsets 

  @author Paolo Calafiura
  $Id: PileUpMergeSvc.h,v 1.17 2008-04-19 00:31:09 calaf Exp $
*/

class PileUpMergeSvc : public AthService {

public:
  /// Standard Gaudi Constructor
  PileUpMergeSvc(const std::string& name, ISvcLocator* svc);
        
  virtual ~PileUpMergeSvc() {}

  virtual StatusCode initialize();

  ///generate the types of the timed data objects
  template <typename DATA>
  struct TimedList {
    typedef DATA data_t;
    ///type of timed data object
    typedef std::pair<PileUpTimeEventIndex, DataLink<DATA> > value_t;
    ///type of the collection of timed data object
    typedef std::list<value_t> type;
  };

  ///retrieve keyed DATA objs for the original event only
  template <typename KEY, typename DATA>
  StatusCode 
  retrieveOriginal(const KEY& dataKey, const DATA*& data);

  template <typename KEY, typename DATA>
  StatusCode
  retrieveSingleSubEvtData(const KEY& dataKey, const DATA*& data, int bunchXing, SubEventIterator iEvt);

  ///retrieve keyed DATA objs for all sub-events and attach a time to them
  template <typename KEY, typename TIMEDDATA>
  StatusCode 
  retrieveSubEvtsData(const KEY& dataKey, //orig evt key
 		      TIMEDDATA& timedData);

  template <typename KEY, typename TIMEDDATA>
  StatusCode 
  retrieveSubSetEvtData(const KEY& dataKey, //orig evt key
			TIMEDDATA& timedData, int bunchXing,
			SubEventIterator bSubEvents, SubEventIterator eSubEvents);

  ///retrieve keyed DATA objs for all sub-events and attach a time to them
  template <typename KEY, typename TIMEDDATA>
  StatusCode 
  retrieveSubEvtsData(const KEY& dataKey, //orig evt key
		      TIMEDDATA& timedData, unsigned int& numberOfSimHits);

  
//   ///retrieve default DATA objs for all sub-events and attach a time to them
//   template <typename TIMEDDATA>
//   StatusCode 
//   retrieveSubEvtsData(TIMEDDATA& timedData);

  ///clear bkg event caches from unneeded data objects (as configured using 
  /// PileUpXingFolder CacheRefreshFrequency property)
  StatusCode clearDataCaches();

  /// Retrieve interface ID
  static const InterfaceID& interfaceID();
  virtual StatusCode queryInterface( const InterfaceID& riid, void** ppvInterface );

  /// get EventInfo from SG, by default using p_overStore
  const xAOD::EventInfo* getPileUpEvent( StoreGateSvc* sg,
                                         const std::string& einame ) const;

private:
  ServiceHandle<StoreGateSvc> p_overStore;      ///< overlaid SG (default)
  ToolHandleArray<IPileUpXingFolder> m_intervals; ///< Property: bunch xing intervals

  // Protect against multiple threads trying to make EventInfo
  // for the same slot.
  mutable SG::SlotSpecificObj<std::mutex> m_slotMutex ATLAS_THREAD_SAFE;

  template <typename DATA, typename KEY>
  bool isLive(const KEY& key, int iXing);   ///< is iXing live for DATA/key?
  bool isLive(CLID id, const std::string& dataKey, int iXing);

  void decodeIntervals();

  /// the active crossing range for a data object (CLID/key combination)
  class Range {
  public:
    Range() : m_noLimits(true), m_first(0), m_last(0), m_cacheRefreshFrequency(0.0) {}
    Range(int first, int last, double cacheRefreshFrequency) : 
      m_noLimits(false), m_first(first), m_last(last), 
      m_cacheRefreshFrequency(cacheRefreshFrequency) {
      assert(0.0<=m_cacheRefreshFrequency && m_cacheRefreshFrequency<=1.0);
    }
    Range(const Range& rhs) : 
      m_noLimits(rhs.m_noLimits),
      m_first(rhs.m_first), m_last(rhs.m_last),
      m_cacheRefreshFrequency(rhs.m_cacheRefreshFrequency) {}
    Range& operator =(const Range& rhs)
    {
      if (this == &rhs) {return *this;}   // Handle self assignment
      m_noLimits = rhs.m_noLimits;
      m_first = rhs.m_first; 
      m_last = rhs.m_last;
      m_cacheRefreshFrequency = rhs.m_cacheRefreshFrequency;
      return *this;
    }
    bool contains(int xing) const {
      return m_noLimits || (m_first<=xing && xing<=m_last);
    }
    bool doRefresh(float random) const { 
      assert(0.0 <= random && random <= 1.0);
      //      bool result(random < m_cacheRefreshFrequency);
      //      std::cerr << "Range::doRefresh: random " << random 
      //		<< " frequency " << m_cacheRefreshFrequency 
      //		<< " result " << result << std::endl;
      return ( random < m_cacheRefreshFrequency );
    }
  private:
    bool m_noLimits;
    int m_first, m_last;
    double m_cacheRefreshFrequency;
  };
  typedef std::map<std::pair<CLID, std::string>, Range> RangeContainer;
  RangeContainer m_ranges;

  ToolHandle<ITriggerTime> m_pITriggerTime; ///< allows to apply a trigger time offset
  ///< controls PileUpTimedEventIndex for TimedData returned by retrieveSubEvts
  BooleanProperty m_returnTimedData; 

  bool doRefresh(const Range& r, int iXing);

  //Default name for EventInfo
  std::string m_EventInfoKeyName;

  /// property: Handle to the EventInfo -> xAOD::EventInfo converter tool
  ToolHandle< xAODMaker::IEventInfoCnvTool > m_xAODCnvTool;
};
#include "PileUpTools/PileUpMergeSvc.icc"
#endif /* PILEUPTOOLS_PILEUPMERGESVC_H */

