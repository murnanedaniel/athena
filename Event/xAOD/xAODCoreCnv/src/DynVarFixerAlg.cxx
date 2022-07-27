/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

// Local include(s):
#include "DynVarFixerAlg.h"

// EDM include(s):
#include "AthContainers/AuxTypeRegistry.h"
#include "StoreGate/ReadHandle.h"

// ROOT include(s):
#include <TClass.h>
#include <TVirtualCollectionProxy.h>
#include <TGenCollectionProxy.h>

// System include(s).
#include <algorithm>
#include <cassert>
#include <string>

namespace xAODMaker {

   StatusCode DynVarFixerAlg::initialize() {

      // Greet the user.
      ATH_MSG_INFO( "Initialising" );
      ATH_MSG_DEBUG( "Will be fixing dynamic variables in container(s):" );
      ATH_MSG_DEBUG( "  " << m_constKeys );

      // Initialise the read handle keys.
      ATH_CHECK( m_constKeys.initialize() );
      std::vector< std::string > keys( m_constKeys.size() );
      std::transform( m_constKeys.begin(), m_constKeys.end(), keys.begin(),
                      []( const SG::ReadHandleKey< SG::IConstAuxStore >& key ) {
                         return key.key(); } );
      ATH_CHECK( m_ioKeys.assign( keys ) );
      ATH_CHECK( m_ioKeys.initialize() );

      // Reset the cache variable.
      m_dicts.clear();

      // Return gracefully.
      return StatusCode::SUCCESS;
   }

   StatusCode DynVarFixerAlg::execute() {

      // Construct the SG::IConstAuxStore and SG::IAuxStoreIO read handles.
      auto constHandles = m_constKeys.makeHandles();
      auto ioHandles = m_ioKeys.makeHandles();
      assert( constHandles.size() == ioHandles.size() );

      // Loop over the specified containers:
      for( std::size_t i = 0; i < constHandles.size(); ++i ) {

         // Access the appropriate handles.
         SG::ReadHandle< SG::IConstAuxStore >& constHandle =
            constHandles.at( i );
         SG::ReadHandle< SG::IAuxStoreIO >& ioHandle =
            ioHandles.at( i );

         // If the store is empty, stop here.
         if( constHandle->size() == 0 ) {
            continue;
         }

         // Remember the store's name and size.
         const std::string& cname = constHandle.key();
         const size_t size = constHandle->size();
         ATH_MSG_VERBOSE( "Size of container \"" << cname << "\": " << size );

         // Loop over all dynamic variables of the container:
         for( SG::auxid_t auxid : ioHandle->getDynamicAuxIDs() ) {

            // Access the registry:
            static const auto& reg = SG::AuxTypeRegistry::instance();

            // Tell the user what's happening.
            ATH_MSG_VERBOSE( "Checking variable \"" << cname
                             << reg.getName( auxid ) << "\"" );

            // Access the std::vector variable:
            const void* vecPtr = ioHandle->getIOData( auxid );

            // Get a dictionary for the variable:
            ::TClass* cl = getClass( auxid );
            if( ! cl ) {
               continue;
            }

            // Get the collection proxy for the class:
            ::TVirtualCollectionProxy* proxy = cl->GetCollectionProxy();
            if( ! proxy ) {
               ATH_MSG_FATAL( "Couldn't get collection proxy for variable \""
                              << reg.getName( auxid ) << "\"" );
               return StatusCode::FAILURE;
            }

            // Attach the proxy to the auxiliary variable:
            ::TVirtualCollectionProxy::TPushPop helper( proxy,
                                                        ( void* ) vecPtr );

            // Check if the variable is of the right size. If it is, don't
            // bother the variable any longer.
            if( proxy->Size() == size ) {
               continue;
            }
            ATH_MSG_VERBOSE( "Size of variable \"" << cname
                             << reg.getName( auxid ) << "\": "
                             << proxy->Size() );

            // If it's not, then try to resize it using the TGenCollectionProxy
            // interface:
            ::TGenCollectionProxy* genProxy =
                  dynamic_cast< ::TGenCollectionProxy* >( proxy );
            if( ! genProxy ) {
               ATH_MSG_WARNING( "Variable \"" << reg.getName( auxid )
                                << "\" doesn't have a TGenCollectionProxy. "
                                << "Variable not fixed!" );
               continue;
            }

            // Resize the variable to the right size:
            static const ::Bool_t FORCE_DELETE = kFALSE;
            genProxy->Resize( size, FORCE_DELETE );
            ATH_MSG_VERBOSE( "Size of variable \"" << cname
                             << reg.getName( auxid ) << "\" after resize: "
                             << proxy->Size() );
            ATH_MSG_DEBUG( "Fixed variable: \"" << cname << reg.getName( auxid )
                           << "\"" );
         }
      }

      // Return gracefully.
      return StatusCode::SUCCESS;
   }

   ::TClass* DynVarFixerAlg::getClass( SG::auxid_t auxid ) {

      // Check if we already checked, and there's no dictionary available for
      // this type:
      if( ( m_noDict.size() > auxid ) && m_noDict[ auxid ] ) {
         return nullptr;
      }

      // Check if the cache already has this variable:
      if( ( m_dicts.size() > auxid ) && m_dicts[ auxid ] ) {
         return m_dicts[ auxid ];
      }

      // Make sure that the cache variables are large enough:
      if( m_dicts.size() <= auxid ) {
         m_dicts.resize( auxid + 1, nullptr );
      }
      if( m_noDict.size() <= auxid ) {
         m_noDict.resize( auxid + 1, false );
      }

      // If not, then look for a dictionary for this type:
      static const auto& reg = SG::AuxTypeRegistry::instance();
      static const ::Bool_t LOAD = kTRUE;
      static const ::Bool_t SILENT = kTRUE;
      ::TClass* cl = ::TClass::GetClass( reg.getVecTypeName( auxid ).c_str(),
                                         LOAD, SILENT );

      // Cache the result:
      if( cl ) {
         m_dicts[ auxid ] = cl;
      } else {
         ATH_MSG_INFO( "No dictionary available for variable \""
                       << reg.getName( auxid ) << "\"" );
         m_noDict[ auxid ] = true;
      }

      // Return what we found:
      return cl;
   }

} // namespace xAODMaker
