/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// $Id$
/**
 * @file AthContainers/src/AuxStoreInternal.cxx
 * @author scott snyder <snyder@bnl.gov>
 * @date Apr, 2013
 * @brief An auxiliary data store that holds data internally.
 */


#include <iostream>
#include <sstream>

#include "AthContainers/AuxStoreInternal.h"
#include "AthContainers/AuxTypeRegistry.h"
#include "AthContainers/exceptions.h"
#include "AthContainers/tools/foreach.h"
#include "AthContainers/tools/error.h"


namespace SG {


/**
 * @brief Constructor.
 * @param standalone If true, then write this in standalone mode.
 */
AuxStoreInternal::AuxStoreInternal (bool standalone /*= false*/)
  : m_standalone (standalone),
    m_tick (0),
    m_locked (false)
{
}


/**
 * @brief Destructor.
 *
 * All contained data will be deleted.
 */
AuxStoreInternal::~AuxStoreInternal()
{
  ATHCONTAINERS_FOREACH (IAuxTypeVector* p, m_vecs)
    delete p;
}


/**
 * @brief Copy constructor.
 */
AuxStoreInternal::AuxStoreInternal (const AuxStoreInternal& other)
  : m_standalone (other.m_standalone),
    m_isDecoration (other.m_isDecoration),
    m_auxids (other.m_auxids),
    m_tick (1),
    m_locked (other.m_locked)
{
  size_t size = other.m_vecs.size();
  m_vecs.resize (size);
  for (size_t i = 0; i < size; i++) {
    if (other.m_vecs[i])
      m_vecs[i] = other.m_vecs[i]->clone();
  }
}


/**
 * @brief Return the standalone flag.
 */
bool AuxStoreInternal::standalone() const
{
  return m_standalone;
}


/**
 * @brief Return the data vector for one aux data item
 * @param auxid The identifier of the desired aux data item.
 *
 * Each aux data item is stored as a vector, with one entry
 * per entry in the owning container.  This returns a pointer
 * to the start of the vector.
 *
 * This should return 0 if the item doesn't exist.
 */
const void* AuxStoreInternal::getData (auxid_t auxid) const
{
  guard_t guard (m_mutex);
  if (auxid >= m_vecs.size() || !m_vecs[auxid]) {
    // With the new behavior of SG::AuxElement::Accessor::isAvailable,
    // we shouldn't print an error message here. Asking the store whether
    // it has an element using this function is not necessarily an
    // error condition by now. In any case, the DataVector code will
    // complain itself in case of an error.
    return 0;
  }
  return m_vecs[auxid]->toPtr();
}


/**
 * @brief Return the data vector for one aux data item
 * @param auxid The identifier of the desired aux data item.
 * @param size The current size of the container (in case the data item
 *             does not already exist).
 * @param capacity The current capacity of the container (in case
 *                 the data item does not already exist).
 *
 * Each aux data item is stored as a vector, with one entry
 * per entry in the owning container.  This returns a pointer
 * to the start of the vector.
 *
 * If the data item does not exist, it should be created and initialized
 * to default values.  @c size and @c capacity give the size for the
 * new aux data item vector.
 */
void* AuxStoreInternal::getData (auxid_t auxid, size_t size, size_t capacity)
{
  return getDataInternal (auxid, size, capacity, false);
}


/**
 * @brief Return the data vector for one aux data decoration item.
 * @param auxid The identifier of the desired aux data item.
 * @param size The current size of the container (in case the data item
 *             does not already exist).
 * @param capacity The current capacity of the container (in case
 *                 the data item does not already exist).
 *
 * Each aux data item is stored as a vector, with one entry
 * per entry in the owning container.  This returns a pointer
 * to the start of the vector.
 *
 * If the data item does not exist, it then it will be created and initialized
 * with default values.  If the container is locked, then the new
 * item will be marked as a decoration.  @c size and @c capacity give
 * the size for the new aux data item vector.
 *
 * If the data item already exists, then we return it if either the
 * container is not locked or the item is marked as a decoration.
 * Otherwise we throw an exception.
 */
void*
AuxStoreInternal::getDecoration (auxid_t auxid, size_t size, size_t capacity)
{
  guard_t guard (m_mutex);
  if (m_vecs.size() <= auxid) {
    m_vecs.resize (auxid+1);
    m_isDecoration.resize (auxid+1);
  }
  if (m_vecs[auxid] == 0) {
    m_vecs[auxid] = AuxTypeRegistry::instance().makeVector (auxid, size, capacity);
    addAuxID (auxid);
    if (m_locked)
      m_isDecoration[auxid] = true;
  }
  if (m_locked && !m_isDecoration[auxid])
    throw ExcStoreLocked (auxid);
  return m_vecs[auxid]->toPtr();
}


/**
 * @brief Change the size of all aux data vectors.
 * @param sz The new size.
 *
 * This should be called when the size of the container changes.
 * This should resize the vectors for all aux data items.
 *
 * If the size of the container grows, the new elements should
 * be default-initialized; if it shrinks, destructors should
 * be run as appropriate.
 */
void AuxStoreInternal::resize (size_t sz)
{
  guard_t guard (m_mutex);
  if (m_locked)
    throw ExcStoreLocked ("resize");
  ATHCONTAINERS_FOREACH (IAuxTypeVector* v, m_vecs) {
    if (v)
      v->resize (sz);
  }
}


/**
 * @brief Change the capacity of all aux data vectors.
 * @param sz The new capacity.
 *
 * This should be called when the capacity of the container changes
 * (by @c reserve).  This should change the capacity for the vectors
 * for all aux data items.
 */
void AuxStoreInternal::reserve (size_t sz)
{
  guard_t guard (m_mutex);
  if (m_locked)
    throw ExcStoreLocked ("reserve");
  ATHCONTAINERS_FOREACH (IAuxTypeVector* v, m_vecs) {
    if (v)
      v->reserve (sz);
  }
}


/**
 * @brief Shift the elements of the container.
 * @param pos The starting index for the shift.
 * @param offs The (signed) amount of the shift.
 *
 * This operation shifts the elements in the vectors for all
 * aux data items, to implement an insertion or deletion.
 * @c offs may be either positive or negative.
 *
 * If @c offs is positive, then the container is growing.
 * The container size should be increased by @c offs,
 * the element at @c pos moved to @c pos + @c offs,
 * and similarly for following elements.
 * The elements between @c pos and @c pos + @c offs should
 * be default-initialized.
 *
 * If @c offs is negative, then the container is shrinking.
 * The element at @c pos should be moved to @c pos + @c offs,
 * and similarly for following elements.
 * The container should then be shrunk by @c -offs elements
 * (running destructors as appropriate).
 */
void AuxStoreInternal::shift (size_t pos, ptrdiff_t offs)
{
  guard_t guard (m_mutex);
  if (m_locked)
    throw ExcStoreLocked ("shift");
  ATHCONTAINERS_FOREACH (IAuxTypeVector* v, m_vecs) {
    if (v)
      v->shift (pos, offs);
  }
}


/**
 * @brief Return a set of identifiers for existing data items
 *        in this store.
 *
 *        This should include identifiers for all items,
 *        const and non-const.
 */
const SG::auxid_set_t&
AuxStoreInternal::getAuxIDs() const
{
  guard_t guard (m_mutex);
  if (m_tsAuxids.get() == 0) {
    m_tsAuxids.reset (new TSAuxidSet (m_tick, m_auxids));
  }
  else if (m_tsAuxids->m_tick != m_tick) {
    m_tsAuxids->m_set = m_auxids;  // May need to optimize this!
    m_tsAuxids->m_tick = m_tick;
  }
  return m_tsAuxids->m_set;
}


/**
 * @brief Return a set of identifiers for writable data items
 *        in this store.
 *
 *        This should include only non-const identifiers.
 */
const SG::auxid_set_t&
AuxStoreInternal::getWritableAuxIDs() const
{
  return getAuxIDs();
}


/**
 * @brief Return a pointer to the data to be stored for one aux data item.
 * @param auxid The identifier of the desired aux data item.
 * @param quiet If true, then don't print an error on failure.
 *
 * This will usually be a pointer to a @c std::vector; however, it may
 * be something different for a standalone object.
 *
 * Returns 0 and reports an error if the requested aux data item
 * does not exist.
 */
const void* AuxStoreInternal::getIODataInternal (auxid_t auxid, bool quiet) const
{
  guard_t guard (m_mutex);
  if (auxid >= m_vecs.size() || !m_vecs[auxid]) {
    if (!quiet) {
      std::ostringstream ss;
      ss  << "Requested variable "
          << SG::AuxTypeRegistry::instance().getName (auxid)
          << " (" << auxid << ") doesn't exist";
      ATHCONTAINERS_ERROR("AuxStoreInternal::getIODataInternal", ss.str());
    }
    return 0;
  }

  if (m_standalone)
    return m_vecs[auxid]->toPtr();
  return m_vecs[auxid]->toVector();
}


/**
 * @brief Return a pointer to the data to be stored for one aux data item.
 * @param auxid The identifier of the desired aux data item.
 *
 * This will usually be a pointer to a @c std::vector; however, it may
 * be something different for a standalone object.
 *
 * Returns 0 and reports an error if the requested aux data item
 * does not exist.
 */
const void* AuxStoreInternal::getIOData (auxid_t auxid) const
{
  return getIODataInternal (auxid, false);
}


/**
 * @brief Return the type of the data to be stored for one aux data item.
 * @param auxid The identifier of the desired aux data item.
 *
 * For an aux data item of type @c T, this will usually be
 * @c std::vector<T>.  For standalone objects, however, it will
 * usually be @c T; and @c std::vector<char> will be used instead
 * of @c std::vector<bool>.
 *
 * Returns 0 if the requested aux data item does not exist.
 */
const std::type_info* AuxStoreInternal::getIOType (auxid_t auxid) const
{
  if (m_standalone)
    return SG::AuxTypeRegistry::instance().getType (auxid);
  return SG::AuxTypeRegistry::instance().getVecType (auxid);
}


/**
 * @brief Get the list of all variables that need to be handled.
 */
const SG::auxid_set_t&
AuxStoreInternal::getDynamicAuxIDs() const
{
  return getAuxIDs();
}


/**
 * @brief Lock the container.
 *
 * After this, only decorations can be changed/modified.
 * If the container is already locked, this is a no-op.
 */
void AuxStoreInternal::lock()
{
  guard_t guard (m_mutex);
  m_locked = true;
}


/**
 * @brief Clear all decorations.
 *
 * Erase all decorations from the store, restoring the state to when
 * @c lock was called.  Be sure to clear the cache of the referencing
 * container!
 */
void AuxStoreInternal::clearDecorations()
{
  guard_t guard (m_mutex);
  for (auxid_t id = 0; id < m_vecs.size(); id++) {
    if (m_isDecoration[id]) {
      m_isDecoration[id] = false;
      delete m_vecs[id];
      m_vecs[id] = 0;
      m_auxids.erase (id);
      ++m_tick;
    }
  }
}


/**
 * @brief Return the number of elements in the store.
 *
 * May return 0 for a store with no aux data.
 */
size_t AuxStoreInternal::size() const
{
  guard_t guard (m_mutex);
  for (SG::auxid_t id : m_auxids) {
    if (id < m_vecs.size() && m_vecs[id] && m_vecs[id]->size() > 0)
      return m_vecs[id]->size();
  }
  return 0;
}


/**
 * @brief Add a new auxid to the set of those being managed by this store.
 * @param auxid The auxid to add.
 */
void AuxStoreInternal::addAuxID (auxid_t auxid)
{
  m_auxids.insert (auxid);
  ++m_tick;
}


/**
 * @brief Return the data vector for one aux data item
 * @param auxid The identifier of the desired aux data item.
 * @param size The current size of the container (in case the data item
 *             does not already exist).
 * @param capacity The current capacity of the container (in case
 *                 the data item does not already exist).
 * @param no_lock_check If true, then skip the test for a locked container.
 *
 * Each aux data item is stored as a vector, with one entry
 * per entry in the owning container.  This returns a pointer
 * to the start of the vector.
 *
 * If the data item does not exist, it should be created and initialized
 * to default values.  @c size and @c capacity give the size for the
 * new aux data item vector.
 */
void* AuxStoreInternal::getDataInternal (auxid_t auxid,
                                         size_t size,
                                         size_t capacity,
                                         bool no_lock_check)
{
  guard_t guard (m_mutex);
  if (m_vecs.size() <= auxid) {
    m_vecs.resize (auxid+1);
    m_isDecoration.resize (auxid+1);
  }
  if (m_vecs[auxid] == 0) {
    if (m_locked && !no_lock_check)
      throw ExcStoreLocked (auxid);
    m_vecs[auxid] = AuxTypeRegistry::instance().makeVector (auxid, size, capacity);
    addAuxID (auxid);
  }
  return m_vecs[auxid]->toPtr();
}


} // namespace SG
