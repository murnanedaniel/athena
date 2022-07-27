/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef ViewAlgs_IDCCacheCreatorBase_h
#define ViewAlgs_IDCCacheCreatorBase_h

#include "AthenaBaseComps/AthReentrantAlgorithm.h"
#include <atomic>
#include "EventContainers/IdentifiableValueCache.h"
#include "EventContainers/IdentifiableCache.h"

class IDCCacheCreatorBase : public AthReentrantAlgorithm {
 public:
  /// Constructor
  IDCCacheCreatorBase(const std::string &name,ISvcLocator *pSvcLocator);
  /// Destructor
  virtual ~IDCCacheCreatorBase()=default;
protected:
  template<bool checkKey = true, typename T>
  StatusCode createContainer(const SG::WriteHandleKey<T>& , long unsigned int , const EventContext& ) const;
  template<bool checkKey = true, typename T, typename X>
  StatusCode createValueContainer(const SG::WriteHandleKey<T>& , long unsigned int , const EventContext&, const X& defaultValue ) const;
  mutable std::atomic_flag m_disableWarningCheck ATLAS_THREAD_SAFE = ATOMIC_FLAG_INIT;
  bool isInsideView(const EventContext&) const;
  StatusCode checkInsideViewOnce(const EventContext&) const;
};


template<bool checkKey, typename T>
StatusCode IDCCacheCreatorBase::createContainer(const SG::WriteHandleKey<T>& containerKey, long unsigned int size, const EventContext& ctx) const{
    static_assert(std::is_base_of<EventContainers::IdentifiableCacheBase, T>::value, "Expects a IdentifiableCache Class" );
    if constexpr (checkKey){
      if(containerKey.key().empty()){
        ATH_MSG_DEBUG( "Creation of container "<< containerKey.key() << " is disabled (no name specified)");
        return StatusCode::SUCCESS;
      }
    }
    SG::WriteHandle<T> ContainerCacheKey(containerKey, ctx);
    ATH_CHECK( ContainerCacheKey.recordNonConst ( std::make_unique<T>(IdentifierHash(size), nullptr) ));
    ATH_MSG_DEBUG( "Container "<< containerKey.key() << " created to hold " << size );
    return StatusCode::SUCCESS;
}

template<bool checkKey, typename T, typename X>
StatusCode IDCCacheCreatorBase::createValueContainer(const SG::WriteHandleKey<T>& containerKey, long unsigned int size, const EventContext& ctx, const X& defaultValue) const{
    static_assert(std::is_base_of<IdentifiableValueCache<X>, T>::value, "Expects a IdentifiableValueCache Class" );
    if constexpr (checkKey){
      if(containerKey.key().empty()){
        ATH_MSG_DEBUG( "Creation of container "<< containerKey.key() << " is disabled (no name specified)");
        return StatusCode::SUCCESS;
      }
    }
    SG::WriteHandle<T> ContainerCacheKey(containerKey, ctx);
    ATH_CHECK( ContainerCacheKey.recordNonConst ( std::make_unique<T>(size, defaultValue) ));
    ATH_MSG_DEBUG( "ValueContainer "<< containerKey.key() << " created to hold " << size );
    return StatusCode::SUCCESS;
}

#endif
