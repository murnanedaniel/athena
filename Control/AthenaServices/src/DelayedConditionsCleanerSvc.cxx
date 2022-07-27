/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/
/**
 * @file AthenaServices/src/DelayedConditionsCleanerSvc.cxx
 * @author scott snyder <snyder@bnl.gov>
 * @date May, 2018
 * @brief Clean conditions containers after a delay.
 */


#include "DelayedConditionsCleanerSvc.h"
#include "AthenaBaseComps/AthProperties.h"
#include "AthenaKernel/CondCont.h"
#include "AthenaKernel/IRCUSvc.h"
#include "CxxUtils/StrFormat.h"
#include "GaudiKernel/EventContext.h"
#include "GaudiKernel/ServiceHandle.h"
#include <algorithm>
#include <unordered_set>


// We wanted to try to allow cleaning to run as asynchronous tasks,
// but there are some race conditions that appear to be difficult
// to resolve.  These stem from the fact that conditions are only
// written from CondInputLoader.  If we clean after that, and don't
// have the full set of IOV keys for all executing events, then we can
// clean an item which the slot won't be able to recover.
// Leave the code commented-out for now while we think about this further.
#define USE_ASYNC_TASK 0


#if USE_ASYNC_TASK
#include "tbb/task.h"
#endif


namespace Athena {


class DelayedConditionsCleanerSvcProps
  : public AthProperties<DelayedConditionsCleanerSvc>
{
public:
  DelayedConditionsCleanerSvcProps (DelayedConditionsCleanerSvc* parent)
    : AthProperties<DelayedConditionsCleanerSvc>(parent) {}

  /// Property: Number of previous events for which to remember IOV history.
  Gaudi::Property<size_t> m_ringSize
    { parent(), "RingSize", 100,
      "Number of previous events for which to remember IOV history." };

  /// Property: Number of events after adding a conditions object
  ///           we try to clean its container.
  Gaudi::Property<size_t> m_cleanDelay
    { parent(), "CleanDelay", 100,
      "Number of events after adding a conditions object we try to clean its container." };

  /// Property: Maximum number of events to consolodate together when cleaning.
  Gaudi::Property<size_t> m_lookAhead
    { parent(), "LookAhead", 10,
      "Maximum number of events to consolodate together when cleaning." };


  /// Property: If true, run cleaning asynchronously in an MT job.
#if USE_ASYNC_TASK
  Gaudi::Property<bool> m_async
    { parent(), "Async", false,
      "If true, run cleaning asynchronously in an MT job." };
#else
  bool m_async = false;
#endif

  /// Property: RCU Service.
  ServiceHandle<Athena::IRCUSvc> m_rcu
    { parent(), "RCUSvc", "Athena::RCUSvc",
      "The RCU service." };
};


#if USE_ASYNC_TASK
/**
 * @brief TBB task for running cleaning asynchronously.
 */
class DelayedConditionsCleanerTask
  : public tbb::task
{
public:
  /**
   * @brief A TBB task to run cleaning asynchronously.
   * @param cleaner The parent service instance.
   * @param cis List of conditions containers to clean.
   * @param keyType Run+LBN or timestamp keys?
   * @param keys Set of IOV keys for recent events.
   */
  DelayedConditionsCleanerTask (DelayedConditionsCleanerSvc& cleaner,
                                std::vector<DelayedConditionsCleanerSvc::CondContInfo*>&& cis,
				DelayedConditionsCleanerSvc::twoKeys_t&& keys);

  /**
   * @brief Run asynchonous cleaning.
   */
  tbb::task* execute() override;


private:
  /// The parent service instance.
  DelayedConditionsCleanerSvc& m_cleaner;

  /// List of conditions containers to clean.
  std::vector<DelayedConditionsCleanerSvc::CondContInfo*> m_cis;

  /// Set of IOV keys for recent events.
  DelayedConditionsCleanerSvc::twoKeys_t m_keys;
};


/**
 * @brief A TBB task to run cleaning asynchronously.
 * @param cleaner The parent service instance.
 * @param cis List of conditions containers to clean.
 * @param keys Set of IOV keys for recent events.
 */
DelayedConditionsCleanerTask::DelayedConditionsCleanerTask
  (DelayedConditionsCleanerSvc& cleaner,
   std::vector<DelayedConditionsCleanerSvc::CondContInfo*>&& cis,
   DelayedConditionsCleanerSvc::twoKeys_t&& keys)
  : m_cleaner (cleaner),
    m_cis (cis),
    m_keys (keys)
{
}


/**
 * @brief Run asynchonous cleaning.
 */
tbb::task* DelayedConditionsCleanerTask::execute()
{
  // Do the cleaning.
  m_cleaner.cleanContainers (std::move (m_cis), std::move (m_keys));

  // This task is terminating.
  --m_cleaner.m_cleanTasks;
  return nullptr;
}
#endif // USE_ASYNC_TASK


/**
 * @brief Standard Gaudi constructor.
 * @param name Service name.
 * @param svc Service locator.
 */
DelayedConditionsCleanerSvc::DelayedConditionsCleanerSvc (const std::string& name,
                                                          ISvcLocator* svc)
  : base_class (name, svc),
    m_props (std::make_unique<DelayedConditionsCleanerSvcProps> (this))
{
}


/**
 * @brief Standard Gaudi initialize method.
 */
StatusCode DelayedConditionsCleanerSvc::initialize()
{
  // Set the ring buffer sizes.
  m_runlbn.reset (m_props->m_ringSize);
  m_timestamp.reset (m_props->m_ringSize);

  ATH_CHECK( m_props->m_rcu.retrieve() );
  size_t nslots = m_props->m_rcu->getNumSlots();
  m_slotLBN.resize (nslots);
  m_slotTimestamp.resize (nslots);

  return StatusCode::SUCCESS;
}


/**
 * @brief Called at the start of each event.
 * @param ctx The current event context.
 * @param allowAsync If true, then cleaning may be run
 *                   in an asynchronous TBB task.
 */
StatusCode
DelayedConditionsCleanerSvc::event (const EventContext& ctx, bool allowAsync)
{
  // Push the IOV key for the current event into the ring buffers.
  // Also save in the per-slot arrays.
  key_type key_lbn = CondContBase::keyFromRunLBN (ctx.eventID());
  key_type key_ts  = CondContBase::keyFromTimestamp (ctx.eventID());
  m_runlbn.push (key_lbn);
  m_timestamp.push (key_ts);
  EventContext::ContextID_t slot = ctx.slot();
  if (slot != EventContext::INVALID_CONTEXT_ID) {
    m_slotLBN[slot] = key_lbn;
    m_slotTimestamp[slot] = key_ts;
  }

  // Return now if an asynchronous cleaning task is still running ---
  // we don't want to start a new one yet.  We'll check pending work
  // on the next call.
  if (m_cleanTasks > 0) {
    return StatusCode::SUCCESS;
  }

  // Collect conditions containers in need of cleaning.
  std::vector<CondContInfo*> ci_vec;
  {
    lock_t lock (m_workMutex);
    // Is it time to clean the container at the top of the work queue?
    if (!m_work.empty() && m_work.top().m_evt <= ctx.evt()) {
      ++m_nEvents;
      size_t sz = m_work.size();
      m_queueSum += sz;
      m_maxQueue = std::max (m_maxQueue, sz);

      // Yes.  Put it on the correct list.  Also look ahead in the queue
      // a bit; if there are other containers that we want to clean soon,
      // go ahead and do them now.
      do {
        CondContInfo* ci = m_work.top().m_ci;
        switch (ci->m_cc.keyType()) {
        case KeyType::SINGLE:
          break;
        case KeyType::RUNLBN:
        case KeyType::MIXED:
        case KeyType::TIMESTAMP:
          ci_vec.push_back (ci);
          break;
        default:
          std::abort();
        }
        m_work.pop();
        ++m_workRemoved;
      } while (!m_work.empty() && m_work.top().m_evt <= ctx.evt() + m_props->m_lookAhead);
    }
  }

  // Clean the containers.
  if (!ci_vec.empty()) {
    scheduleClean (std::move (ci_vec), getKeys(m_runlbn,m_timestamp),
                   allowAsync);
  }
  return StatusCode::SUCCESS;
}


/**
 * @brief Called after a conditions object has been added.
 * @param ctx The current event context.
 * @param cc The container to which the object was added.
 */
StatusCode DelayedConditionsCleanerSvc::condObjAdded (const EventContext& ctx,
                                                      CondContBase& cc)
{
  // Add this container to the priority queue.
  lock_t lock (m_workMutex);
  CCInfoMap_t::iterator it = m_ccinfo.find (&cc);
  if (it == m_ccinfo.end()) {
    it = m_ccinfo.emplace (&cc, CondContInfo (cc)).first;
  }

  EventContext::ContextEvt_t evt = ctx.evt();
  m_work.emplace (evt + m_props->m_cleanDelay, it->second);
  return StatusCode::SUCCESS;
}


/**
 * @brief Print some statistics about the garbage collection.
 *        Would generally be called in finalize(), but broken out
 *        as a separate interface for testing/debugging purposes.
 */
StatusCode DelayedConditionsCleanerSvc::printStats() const
{
  // Suppress output if we didn't actually do anything.
  if (m_nEvents == 0) {
    return StatusCode::SUCCESS;
  }

  ATH_MSG_INFO( "Conditions container statistics" );
  ATH_MSG_INFO( CxxUtils::strformat ("  Work q: Max size: %zu  (%zu queries) ",
                                     m_maxQueue, m_nEvents) );
  size_t den = std::max (m_nEvents, 1lu);
  ATH_MSG_INFO( CxxUtils::strformat ("          Avg size: %.2f / Avg removed: %.2f",
                                     static_cast<float>(m_queueSum)/den,
                                     static_cast<float>(m_workRemoved)/den) );

  std::vector<const CondContInfo*> infos;
  for (const auto& p : m_ccinfo) {
    infos.push_back (&p.second);
  }
  std::sort (infos.begin(), infos.end(),
             [](const CondContInfo* a, const CondContInfo* b)
             { return a->m_cc.id().key() < b->m_cc.id().key(); });

  for (const CondContInfo* ci : infos) {
    ATH_MSG_INFO( CxxUtils::strformat ("  %-20s nInserts %6zu maxSize %3zu",
                                       ci->m_cc.id().key().c_str(),
                                       ci->m_cc.nInserts(),
                                       ci->m_cc.maxSize()) );
    den = std::max (ci->m_nClean, 1lu);
    ATH_MSG_INFO( CxxUtils::strformat ("    nClean %zu avgRemoved %.2f 0/1/2+ %zu/%zu/%zu",
                                       ci->m_nClean,
                                       static_cast<float> (ci->m_nRemoved) / den,
                                       ci->m_removed0,
                                       ci->m_removed1,
                                       ci->m_removed2plus) );
  }
                                     
  return StatusCode::SUCCESS;
}


/**
 * @brief Clear the internal state of the service.
 * Only for testing.  Don't call if any other thread may be touching the service.
 */
StatusCode DelayedConditionsCleanerSvc::reset()
{
  m_runlbn.reset (m_props->m_ringSize);
  m_timestamp.reset (m_props->m_ringSize);

  std::fill (m_slotLBN.begin(), m_slotLBN.end(), 0);
  std::fill (m_slotTimestamp.begin(), m_slotTimestamp.end(), 0);

  m_ccinfo.clear();
  std::priority_queue<QueueItem> tmp;
  m_work.swap (tmp);

  m_nEvents = 0;
  m_queueSum = 0;
  m_workRemoved = 0;
  m_maxQueue = 0;
  m_cleanTasks = 0;

  return StatusCode::SUCCESS;
}


DelayedConditionsCleanerSvc::twoKeys_t 
DelayedConditionsCleanerSvc::getKeys(const Ring& runLBRing, const Ring& TSRing) const {
  
  // Get a copy of the contents of the ring buffer holding runLumi and time-stamp keys
  std::vector<key_type> runLBKeys=runLBRing.getKeysDedup();
  std::vector<key_type> TSKeys=TSRing.getKeysDedup();

  // Add in the keys for the currently-executing slots.
  // These are very likely to already be in the ring, but that's
  // not absolutely guaranteed.
  // FIXME: This probably does another memory allocation, due to
  // growing the buffer.  Would be nice to avoid that.
  runLBKeys.insert (runLBKeys.end(), m_slotLBN.begin(), m_slotLBN.end());
  TSKeys.insert(TSKeys.end(), m_slotTimestamp.begin(), m_slotTimestamp.end());

  twoKeys_t result{runLBKeys, TSKeys};

  /// Sort the key array and remove duplicates.
  /// We expect that the key array is probably `almost' sorted.
  /// std::sort, at least in the gcc implementation, is designed
  /// to perform well in such cases.
  for ( auto& keys : result ) {
    std::sort (keys.begin(), keys.end());
    auto end = std::unique (keys.begin(), keys.end());
    keys.resize (end - keys.begin());
  }


  return result;
}


/**
 * @brief Do cleaning for a set of containers.
 * @param cis Set of containers to clean.
 * @param ring Ring buffer with recent IOV keys.
 * @param slotKeys Vector of current keys for all slots.
 * @param allowAsync Can this task run asynchronously?
 *
 * This will either run cleaning directly, or submit it as a TBB task.
 */
void
DelayedConditionsCleanerSvc::scheduleClean (std::vector<CondContInfo*>&& cis,
					    twoKeys_t&& twoKeys,
                                            bool allowAsync)
{
  // Remove any duplicates from the list of containers.
  std::sort (cis.begin(), cis.end());
  auto pos = std::unique (cis.begin(), cis.end());
  cis.resize (pos - cis.begin());

  if (allowAsync && m_props->m_async) {
#if USE_ASYNC_TASK
    // Queue cleaning as a TBB task.
    // Count that we have another executing task.
    ++m_cleanTasks;

    // Create the TBB task and queue it.
    // TBB will delete the task object after it completes.
    tbb::task* t = new (tbb::task::allocate_root())
      DelayedConditionsCleanerTask (*this, std::move (cis),
                                    std::move (twoKeys));
    tbb::task::enqueue (*t);
#endif
  }
  else
  {
    // Call cleaning directly.
    cleanContainers (std::move (cis), std::move (twoKeys));
  }
}


/**
 * @brief Clean a set of containers.
 * @param cis Set of containers to clean.
 * @param keys Set of IOV keys for recent events.
 */
void
DelayedConditionsCleanerSvc::cleanContainers (std::vector<CondContInfo*>&& cis,
					      twoKeys_t&& twoKeys)
{
  // FIXME: Some conditions objects have pointers to parts of other
  // conditions objects, which violates the lifetime guarantees of
  // conditions containers (a pointer you get from a conditions container
  // is guaranteed to be valid until the end of the current event, but
  // not past that).  In some cases, this can be dealt with by replacing
  // the pointers with CondLink, but that is sometimes rather inconvenient.
  // Try to work around this for now by ensuring that when we delete an object,
  // we also try to clean the containers of other objects that may depend
  // on it.
  std::vector<CondContInfo*> toclean = std::move (cis);
  std::unordered_set<CondContInfo*> cleaned (toclean.begin(), toclean.end());
  while (!toclean.empty()) {
    std::vector<CondContInfo*> newclean;
    for (CondContInfo* ci : toclean) {
      if (cleanContainer (ci, twoKeys)) {
        lock_t lock (m_workMutex);
        for (CondContBase* dep : ci->m_cc.getDeps()) {
          CCInfoMap_t::iterator it = m_ccinfo.find (dep);
          // If we don't find it, then dep must have no conditions objects.
          if (it != m_ccinfo.end()) {
            CondContInfo* ci_dep = &it->second;
            if (cleaned.insert (ci_dep).second) {
              newclean.push_back (ci_dep);
            }
          }
        }
      }
    }
    toclean = std::move (newclean);
  }
}


/**
 * @brief Clean a single container.
 * @param ci The container to clean.
 * @param keyType Run+LBN or timestamp keys?
 * @param keys Set of IOV keys for recent events.
 *
 * Returns true if anything was removed from the container,
 */
bool
DelayedConditionsCleanerSvc::cleanContainer (CondContInfo* ci,
                                             const  twoKeys_t& twoKeys) const
{
  size_t n = ci->m_cc.trim (twoKeys[0],twoKeys[1]);

  ++ci->m_nClean;
  ci->m_nRemoved += n;
  switch (n) {
  case 0:
    ++ci->m_removed0;
    break;
  case 1:
    ++ci->m_removed1;
    break;
  default:
    ++ci->m_removed2plus;
  }

  return n > 0;
}
  
/**
 * @brief Standard destructor.
 */

DelayedConditionsCleanerSvc::~DelayedConditionsCleanerSvc() {}

/**
 * @brief Standard Gaudi finalize method.
 */
StatusCode DelayedConditionsCleanerSvc::finalize()
{
  ATH_CHECK( printStats() );
  return StatusCode::SUCCESS;
}


} // namespace Athena
