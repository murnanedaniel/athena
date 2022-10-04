//Dear emacs, this is -*- c++ -*-

/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/


#ifndef LARRAWCONDITIONS_LARCONDITIONSCONTAINERDB_H
#define LARRAWCONDITIONS_LARCONDITIONSCONTAINERDB_H

/** LArConditionsContainerDB template class
    Author: Walter Lampl
    Initial Version: April, 22nd 2005
    
    Container class to store any kind of data object 
    that is organized by FEB and FEB-channel number. 
	
    The data object stored in this container has to have
    a sensible default constructor. A default instance is
    returned by Get-function if the requested channel is 
    not known. 
    Furthermore, the data object must have a operator= 
    because this operator is used by the set function.

    About the iterator: 
    The iterator-class holds internally a iterator over the 
    FEB map and a iterator over the channels inside this FEB 
    as well as an interator pointing to the last feb (NOT 
    the end()!). begin() points to the first channel of the 
    first FEB, end() to the 129th channel of the last FEB. 
    The operator++ jumps form a last channel of FEB n to
    the first channel of FEB n+1 except we are already at
    the last FEB. The operator-- jumps from the first 
    channel of FEB n to the last one of FEB n-1. It will
    try to do this even if FEB n is the first one and 
    thus cause a segfault. 

**/

#include "LArRawConditions/LArConditionsSubset.h"
#include "LArIdentifier/LArOnlineID_Base.h"
#include <map>
#include <unordered_map>
#include <vector>

template<class T> 
class LArConditionsContainerDB 
{
public: 

    typedef          LArConditionsSubsetTraits<T>        Traits;
    typedef typename Traits::FebId                       FebId;
    typedef typename Traits::ChannelVector               ChannelVector;
    typedef typename Traits::ConstChannelVector          ConstChannelVector;
    typedef typename Traits::ChannelVectorPointer        ChannelVectorPointer;
    typedef typename std::map<FebId, ChannelVectorPointer >   ConditionsMap;
    typedef typename std::unordered_map<FebId, ChannelVectorPointer >   ConditionsHashMap;
    typedef typename ChannelVector::const_iterator      ConstChannelIt;
    typedef typename ChannelVector::iterator            ChannelIt;
    typedef typename ConditionsMap::const_iterator      ConstConditionsMapIterator;
    typedef typename ConditionsMap::iterator            ConditionsMapIterator;
    typedef typename ConditionsMap::size_type           size_type;
    typedef std::vector<FebId>                          FebIdVec;
    typedef typename Traits::Reference                  Reference;
    typedef typename Traits::ConstReference             ConstReference;
    typedef typename Traits::Pointer                    Pointer;
    typedef typename Traits::ConstPointer               ConstPointer;

    /// Declaration of const iterator
    template <class MAP_ITERATOR, class CHAN_ITERATOR, class POINTER, class REFERENCE>
    class iteratorT
    {
    public:
	iteratorT();
	iteratorT(MAP_ITERATOR febit, 
                  CHAN_ITERATOR chanit, 
                  MAP_ITERATOR febendit,
                  const LArOnlineID_Base* onlineHelper);
	iteratorT(MAP_ITERATOR febit, 
                  CHAN_ITERATOR chanit, 
                  MAP_ITERATOR febendit,
                  const LArOnlineID_Base* onlineHelper,
                  const FebIdVec& febIds);
	REFERENCE operator*() const;
        POINTER   operator->() const;
	iteratorT& operator++();
	iteratorT  operator++(int);
	bool           operator!=(const iteratorT& a) const;
	bool           operator==(const iteratorT& a) const;
	FebId          getFebId() const;
	int            getChannel() const;
	HWIdentifier   febId() const;
	HWIdentifier   channelId() const;

    private:
	CHAN_ITERATOR m_channelIt;
	MAP_ITERATOR m_febIt;
	// m_lastFebit is initialized with m_febMap.end()...
	//...and decremented to point to the last filled FEB.
	MAP_ITERATOR m_lastFebIt;
	const LArOnlineID_Base*         m_onlineHelper;
	FebIdVec                   m_febIds;
	unsigned int               m_febIdsIndex;

	
    };

    typedef iteratorT<ConditionsMapIterator,
                      ChannelIt,
                      Pointer,
                      Reference> iterator;

    typedef iteratorT<ConstConditionsMapIterator,
                      ConstChannelIt,
                      ConstPointer,
                      ConstReference> const_iterator;

    /// Constructor with gain
    LArConditionsContainerDB(unsigned int gain = 0);
    
    /// Setter
    void                         set(const FebId id, const int channel, const T& payload); 

    /// Getter
    ConstReference               get(const FebId id, const int channel) const;

    /// non-const get function (might extend FEB map)
    Reference                    getNonConst(const FebId id, const int channel);
  
    /// Do conditions for FebId? Used by writers to check whether add
    /// is needed. 
    bool                         exist(FebId id) const;

    /// Iterator over all channels 
    const_iterator               begin(const LArOnlineID_Base* onlineHelper) const;
    const_iterator               end  (const LArOnlineID_Base* onlineHelper) const;

    iterator                     begin(const LArOnlineID_Base* onlineHelper);
    iterator                     end  (const LArOnlineID_Base* onlineHelper);

    /// Iterator over all channels of selected FEBs
    const_iterator               begin(const LArOnlineID_Base* onlineHelper,
				       const FebIdVec& febIds) const;

    iterator                     begin(const LArOnlineID_Base* onlineHelper,
				       const FebIdVec& febIds);

    /// Size of map
    size_type                    size() const;

    /// access to gain
    int                          getGain() const  ; 

    /// set gain
    void                         setGain(int g) ; 

    ///  const ConditionsMap& getConditionsMap()const { return m_febMap ;} ; 
    void                         add(const LArConditionsContainerDB<T>* p); 

    ///  add pointer to vector of (febid/channel vector)
    void                         add(FebId id, ChannelVectorPointer channelVec);

    ///  erase element for this FebId
    void                         erase(FebId id);

protected:

    // Store the FEB -> ChannelVector mapping twice: once as an ordered map
    // (this will be used for iteration) and once as an unordered_map
    // (this will be used for fast lookup).
    ConditionsMap                m_febMap;
    ConditionsHashMap            m_febHashMap;

    unsigned int                 m_gain; 

private:
    //Dummy iterators to mark the beginning and the end of an empty container instance
    const iterator m_dummyIt;
    const const_iterator m_dummyConstIt;
};



template<class T> 
inline
LArConditionsContainerDB<T>::LArConditionsContainerDB(unsigned int gain) 
    :
    m_gain(gain)
{}


template<class T> 
inline
typename LArConditionsContainerDB<T>::const_iterator 
LArConditionsContainerDB<T>::begin(const LArOnlineID_Base* onlineHelper) const 
{
    if (m_febMap.size()==0) { //no elements yet, begin() and end() have to be identical
	return m_dummyConstIt;
    }
    else {
	return const_iterator(m_febMap.begin(),
			      m_febMap.begin()->second->begin(),
			      m_febMap.end(),
			      onlineHelper);    }
}

template<class T> 
inline
typename LArConditionsContainerDB<T>::iterator 
LArConditionsContainerDB<T>::begin(const LArOnlineID_Base* onlineHelper)
{
    if (m_febMap.size()==0) { //no elements yet, begin() and end() have to be identical
	return m_dummyIt;
    }
    else {
	return iterator(m_febMap.begin(),
                        m_febMap.begin()->second->begin(),
                        m_febMap.end(),
                        onlineHelper);    }
}

template<class T> 
inline
typename LArConditionsContainerDB<T>::const_iterator 
LArConditionsContainerDB<T>::begin(const LArOnlineID_Base* onlineHelper,
				   const FebIdVec& febIds) const 
{
    if (m_febMap.size()==0) { //no elements yet, begin() and end() have to be identical
	return m_dummyConstIt;
    }
    else {
      if(febIds.size()==0) return end(onlineHelper);

      ConstConditionsMapIterator it2 = m_febMap.end() ;
      for (FebId id : febIds) {
	it2 = m_febMap.find(id);
	if( it2 != m_febMap.end()) break; 
      }
      if( it2 == m_febMap.end()) return end(onlineHelper);

      return const_iterator( it2,
                             it2->second->begin(),
                             m_febMap.end(),
                             onlineHelper,
                             febIds);
    }
}


template<class T> 
inline
typename LArConditionsContainerDB<T>::iterator 
LArConditionsContainerDB<T>::begin(const LArOnlineID_Base* onlineHelper,
				   const FebIdVec& febIds)
{
    if (m_febMap.size()==0) { //no elements yet, begin() and end() have to be identical
	return m_dummyIt;
    }
    else {
      if(febIds.size()==0) return end(onlineHelper);

      ConditionsMapIterator it2 = m_febMap.end() ;
      for (FebId id : febIds) {
	it2 = m_febMap.find(id);
	if( it2 != m_febMap.end()) break; 
      }
      if( it2 == m_febMap.end()) return end(onlineHelper);

      return iterator( it2,
                       it2->second->begin(),
                       m_febMap.end(),
                       onlineHelper,
                       febIds);
    }
}



template<class T> 
inline
typename LArConditionsContainerDB<T>::const_iterator
LArConditionsContainerDB<T>::end(const LArOnlineID_Base* onlineHelper) const
{
    if (m_febMap.size()==0) {
	return m_dummyConstIt;
    }
    else {
	ConstConditionsMapIterator feb_end_it=m_febMap.end();
	ConstConditionsMapIterator last_feb_it=feb_end_it;
	last_feb_it--;
	return const_iterator(last_feb_it,
			      last_feb_it->second->end(),
			      feb_end_it,
			      onlineHelper);
    }
}


template<class T> 
inline
typename LArConditionsContainerDB<T>::iterator
LArConditionsContainerDB<T>::end(const LArOnlineID_Base* onlineHelper)
{
    if (m_febMap.size()==0) {
	return m_dummyIt;
    }
    else {
	ConditionsMapIterator feb_end_it=m_febMap.end();
	ConditionsMapIterator last_feb_it=feb_end_it;
	last_feb_it--;
	return iterator(last_feb_it,
                        last_feb_it->second->end(),
                        feb_end_it,
                        onlineHelper);
    }
}



#define ITERATORT iteratorT<MAP_ITERATOR, CHAN_ITERATOR, POINTER, REFERENCE>

template<class T> 
template <class MAP_ITERATOR, class CHAN_ITERATOR, class POINTER, class REFERENCE>
inline
LArConditionsContainerDB<T>::ITERATORT::iteratorT() : 
	m_channelIt(),
	m_febIt(),
	m_onlineHelper(0),
	m_febIdsIndex(0)
{}

template<class T> 
template <class MAP_ITERATOR, class CHAN_ITERATOR, class POINTER, class REFERENCE>
inline
LArConditionsContainerDB<T>::ITERATORT::iteratorT
  (MAP_ITERATOR febit, 
   CHAN_ITERATOR chanit, 
   MAP_ITERATOR febendit,
   const LArOnlineID_Base* onlineHelper) : 
	m_channelIt(chanit),
	m_febIt(febit),
	m_lastFebIt(febendit),
	m_onlineHelper(onlineHelper),
	m_febIdsIndex(0)
{ m_lastFebIt--; }

template<class T> 
template <class MAP_ITERATOR, class CHAN_ITERATOR, class POINTER, class REFERENCE>
inline
LArConditionsContainerDB<T>::ITERATORT::iteratorT
  (MAP_ITERATOR febit, 
   CHAN_ITERATOR chanit, 
   MAP_ITERATOR febendit,
   const LArOnlineID_Base* onlineHelper,
   const FebIdVec& febIds)
	:
	m_channelIt(chanit),
	m_febIt(febit),
	m_lastFebIt(febendit),
	m_onlineHelper(onlineHelper),
	m_febIds(febIds.size()),
	m_febIdsIndex(0)
{ 
    // Save last feb it for id match below
    MAP_ITERATOR lastFeb = m_lastFebIt;

    // decrement iterator to point to the last valid feb
    m_lastFebIt--; 

    // If there is a non-zero vector of feb ids, sort them, save them
    // into local vector and move feb iterator to the first valid
    // febId
    if (febIds.size() && m_febIt != lastFeb) {
        typename FebIdVec::const_iterator it    = febIds.begin();
	typename FebIdVec::const_iterator last  = febIds.end();
	std::set<FebId> febIdsSet;
	for (; it != last; ++it) febIdsSet.insert(*it);
	typename std::set<FebId>::const_iterator itSet   = febIdsSet.begin();
	typename std::set<FebId>::const_iterator lastSet = febIdsSet.end();
	for (unsigned int i = 0; itSet != lastSet; ++itSet, ++i) m_febIds[i] = (*itSet);

	// Move set index to point where id is <= to the first febId
	m_febIdsIndex = std::lower_bound(m_febIds.begin(), m_febIds.end(), m_febIt->first)
	    - m_febIds.begin(); 

	// We now require an id match, or one of the two iterators hits the end

	// If febId is not the one found, iterate until a match is found
	while (m_febIt != lastFeb && m_febIdsIndex < m_febIds.size() && 
	       m_febIds[m_febIdsIndex] != m_febIt->first) {
	    // If not equal, then we need to move iterator forward
	    // for the smaller of the two febids
	    if (m_febIt->first < m_febIds[m_febIdsIndex]) {
		++m_febIt; // Move to next feb in container
	    }
	    else {
		++m_febIdsIndex; // Move to next feb in selection
	    }

	}

	if (m_febIt == lastFeb || m_febIdsIndex == m_febIds.size()) {
	    // no match - set iterators to the end
	    m_channelIt = m_lastFebIt->second->end();
	    m_febIt     = m_lastFebIt;
	    return;
	}
	if (m_febIds[m_febIdsIndex] != m_febIt->first) {
	    std::cout << "LArConditionsContainerDB<T>::iteratorT constructor - ERROR: NO MATCH FOUND" 
		      << std::endl;
	}
    }
}


template<class T> 
template <class MAP_ITERATOR, class CHAN_ITERATOR, class POINTER, class REFERENCE>
inline
REFERENCE
LArConditionsContainerDB<T>::ITERATORT::operator*() const 
{return *m_channelIt; }

template<class T> 
template <class MAP_ITERATOR, class CHAN_ITERATOR, class POINTER, class REFERENCE>
inline
POINTER
LArConditionsContainerDB<T>::ITERATORT::operator->() const 
{
  auto& ref = *m_channelIt;
  return &ref;
}

template<class T> 
template <class MAP_ITERATOR, class CHAN_ITERATOR, class POINTER, class REFERENCE>
inline
typename LArConditionsContainerDB<T>::template ITERATORT&
LArConditionsContainerDB<T>::ITERATORT::operator++()
{ 

    // Increment to next channel within the currect feb vector
    ++m_channelIt;
    if (m_channelIt == m_febIt->second->end() && m_febIt != m_lastFebIt) { 

	// Have reached last channel of a FEB vector but not yet in LAST FEB vector
	
	// If there is a non-zero feb id vector, then we require a match
	// of feb ids with those in the vector
	if (m_febIds.size()) {
	    // Get map end - it is 'one past' the saved last m_lastFebIt
            MAP_ITERATOR lastFeb = m_lastFebIt;
	    ++lastFeb;

	    // Move the feb vector index forward
	    ++m_febIdsIndex;
	    
	    // Now if not at end, iterate until the next match is found
	    while (m_febIt != lastFeb && m_febIdsIndex < m_febIds.size() && 
		   m_febIds[m_febIdsIndex] != m_febIt->first) {
		// If not equal, then we need to move iterator forward
		// for the smaller of the two febids
		if (m_febIt->first < m_febIds[m_febIdsIndex]) {
		    ++m_febIt; // Move to next feb in container
		}
		else {
		    ++m_febIdsIndex; // Move to next feb in selection
		}

	    }
	    if (m_febIt == lastFeb || m_febIdsIndex == m_febIds.size()) {
		// no match - set iterators to the end
		m_channelIt = m_lastFebIt->second->end();
		m_febIt     = m_lastFebIt;
		return (*this);
	    }
	    
	    // Print out warning if no match found
	    if (m_febIds[m_febIdsIndex] != m_febIt->first) {
		std::cout << "LArConditionsContainerDB<T>::iteratorT constructor - ERROR: NO MATCH FOUND" 
			  << std::endl;
	    }
	}
	else {
	    // Simple iteration
	    ++m_febIt;
	}
	m_channelIt = m_febIt->second->begin();
    }
    return *this;
}

template<class T> 
template <class MAP_ITERATOR, class CHAN_ITERATOR, class POINTER, class REFERENCE>
inline
typename LArConditionsContainerDB<T>::template ITERATORT
LArConditionsContainerDB<T>::ITERATORT::operator++(int) 
{
  iteratorT tmp = *this;
  operator++();
  return tmp;
}

template<class T> 
template <class MAP_ITERATOR, class CHAN_ITERATOR, class POINTER, class REFERENCE>
inline
bool
LArConditionsContainerDB<T>::ITERATORT::operator==
  (const iteratorT& a) const
{ return (m_channelIt==a.m_channelIt) && (m_febIt==a.m_febIt); }

template<class T> 
template <class MAP_ITERATOR, class CHAN_ITERATOR, class POINTER, class REFERENCE>
inline
bool
LArConditionsContainerDB<T>::ITERATORT::operator!=
  (const iteratorT& a) const
{ return !(*this == a); }

template<class T> 
template <class MAP_ITERATOR, class CHAN_ITERATOR, class POINTER, class REFERENCE>
inline
typename LArConditionsContainerDB<T>::FebId 
LArConditionsContainerDB<T>::ITERATORT::getFebId() const {
  return m_febIt->first;
}

template<class T>
template <class MAP_ITERATOR, class CHAN_ITERATOR, class POINTER, class REFERENCE>
inline
int
LArConditionsContainerDB<T>::ITERATORT::getChannel() const {
  return (m_channelIt - m_febIt->second->begin());
}

template<class T>
template <class MAP_ITERATOR, class CHAN_ITERATOR, class POINTER, class REFERENCE>
inline
HWIdentifier
LArConditionsContainerDB<T>::ITERATORT::febId() const
{
    return (HWIdentifier(m_febIt->first));
}

template<class T>
template <class MAP_ITERATOR, class CHAN_ITERATOR, class POINTER, class REFERENCE>
inline
HWIdentifier
LArConditionsContainerDB<T>::ITERATORT::channelId() const
{
    // Construct channel id
    if (m_onlineHelper) {
	int chan = m_channelIt - m_febIt->second->begin();
	return (m_onlineHelper->channel_Id(febId(), chan));
    }
    HWIdentifier result;
    return (result); // return invalid one if helper is not available
}

#undef ITERATORT



template<class T> 
inline
void
LArConditionsContainerDB<T>::set(const FebId febId, const int channel, const T& payload)
{
  typename ConditionsHashMap::iterator it = m_febHashMap.find (febId);
  if (it == m_febHashMap.end()) {
    // Should never get here
    std::cout << "LArConditionsContainerDB<T>::set ERROR could not find FEB ID " 
              << febId << std::endl;
    return;
  }
    
  // Set payload
  (*(it->second))[channel] = payload;
}


template<class T> 
inline
typename LArConditionsContainerDB<T>::ConstReference
LArConditionsContainerDB<T>::get(const FebId febId, const int channel) const	
{
    // The get/set methods used to cache the iterator in m_febMap where the
    // last access was found.  That, however, is not at all thread-safe
    // (and was in fact observed to lead to incorrect results in MT jobs).
    // We need to get rid of the caching, but we'll also introduce
    // m_febHashMap to try to reduce the lookup overhead.
    // If this turns out to be a hot spot, then the iterator could
    // be returned to the caller via a cookie.
    typename ConditionsHashMap::const_iterator it = m_febHashMap.find (febId);
    if (it == m_febHashMap.end()) {
      // Unknown FEB
      return Traits::empty();
    }

    const ConstChannelVector& vec = *it->second;
    // First check the size - channel vec may be empty
    if (channel < static_cast<int>(vec.size())) {
	return vec[channel];
    }
    else {
        return Traits::empty();
    }
}

template<class T> 
inline
typename LArConditionsContainerDB<T>::Reference
LArConditionsContainerDB<T>::getNonConst(const FebId febId, const int channel)
{
    typename ConditionsHashMap::iterator it = m_febHashMap.find (febId);
    if (it == m_febHashMap.end()) {
      throw std::runtime_error ("Unknown FEB");
    }

    ChannelVector& vec = *it->second;
    // Note: channel vec should always be full for non-const access
    return vec[channel];
}

template<class T> 
inline
bool
LArConditionsContainerDB<T>::exist(FebId id) const
{
  return m_febHashMap.find (id) != m_febHashMap.end();
}

template<class T> 
inline
void  LArConditionsContainerDB<T>::setGain(int g)
{ 
  m_gain = g;
  return;
} 

template<class T> 
inline
typename LArConditionsContainerDB<T>::size_type      
LArConditionsContainerDB<T>::size() const
{
    return (m_febMap.size());
}

template<class T> 
inline
int LArConditionsContainerDB<T>::getGain() const
{ 
  return m_gain;
} 


template<class T> 
inline
void LArConditionsContainerDB<T>::add(const LArConditionsContainerDB<T>* p)
{ 
  // copy the data to this container. FEB by FEB. 
  //  if FEB id exist, the old data will be overwritten. 

  for (const auto& pp : p->m_febMap) {
    m_febMap[pp.first] = pp.second;
    m_febHashMap[pp.first] = pp.second;
  }
} 

template<class T> 
inline
void 
LArConditionsContainerDB<T>::add(FebId id, ChannelVectorPointer channelVec)
{
    m_febMap[id]   = channelVec;
    m_febHashMap[id]   = channelVec;
}

template<class T> 
inline
void 
LArConditionsContainerDB<T>::erase(FebId id)
{
    m_febMap.erase(id);
    m_febHashMap.erase(id);
}


#endif
