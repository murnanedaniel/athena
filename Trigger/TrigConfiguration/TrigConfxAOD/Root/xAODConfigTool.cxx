/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

// System include(s):
#include <stdexcept>

// Trigger include(s):
#include "TrigConfL1Data/BunchGroup.h"
#include "TrigConfHLTData/HLTSequenceList.h"

// Local include(s):
#include "TrigConfxAOD/xAODConfigTool.h"
#include "TrigConfxAOD/tools/prepareTriggerMenu.h"
#include "TrigConfxAOD/tools/xAODKeysMatch.h"

// Boost includes
#define BOOST_BIND_GLOBAL_PLACEHOLDERS // Needed to silence Boost pragma message
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace TrigConf {

   struct xAODConfigTool::Impl
   {
      /// Note: The run 3 JSON derives objects below are used to
      /// additionally populate these objects for interface backwards compatibility
      /// The "translated" LVL1 configuration object
      CTPConfig m_ctpConfig;
      /// The "translated" HLT configuration object
      HLTChainList m_chainList;
      /// The "translated" HLT configuration object
      HLTSequenceList m_sequenceList;
      /// The "translated" bunch group set object
      BunchGroupSet m_bgSet;

      /// The JSON decoded Run3 HLT menu
      HLTMenu m_currentHlt;
      /// The JSON decoded Run3 HLT monitoring
      HLTMonitoring m_currentHltmonitoring;
      /// The JSON decoded Run3 L1 menu
      L1Menu m_currentL1;
      /// The JSON decoded current HLT prescales set
      HLTPrescalesSet m_currentHltps;
      /// The JSON decoded current L1 prescales set
      L1PrescalesSet m_currentL1ps;
      /// The JSON decoded current Bunchgroup configuration
      L1BunchGroupSet m_currentBg;
   };

   xAODConfigTool::xAODConfigTool( const std::string& name )
      : asg::AsgMetadataTool( name ),
        m_tmc( nullptr ),
        m_hltJson( nullptr), m_hltmonitoringJson(nullptr), m_l1Json( nullptr ),
        m_hltpsJson( nullptr ), m_l1psJson( nullptr ), m_bgJson( nullptr ),
        m_menu( nullptr ),
        m_currentHltJson( nullptr ), m_currentHltmonitoringJson( nullptr ), m_currentL1Json( nullptr ),
        m_currentHltpsJson( nullptr ), m_currentL1psJson( nullptr ), m_currentBgJson( nullptr ),
        m_triggerMenuContainerAvailable(false),
        m_menuJSONContainerAvailable(false),
        m_impl (std::make_unique<Impl>()) {

      declareProperty( "EventObjectName", m_eventName = "TrigConfKeys" );
      declareProperty( "BGKeysObjectName", m_bgkeysName = "BunchConfKey" );
      declareProperty( "MetaObjectName", m_metaName_run2 = "TriggerMenu" );

      declareProperty( "JSONMetaObjectNameHLT", m_metaNameJSON_hlt = "TriggerMenuJson_HLT" );
      declareProperty( "JSONMetaObjectNameHLTMonitoring", m_metaNameJSON_hltmonitoring = "TriggerMenuJson_HLTMonitoring" );
      declareProperty( "JSONMetaObjectNameL1", m_metaNameJSON_l1 = "TriggerMenuJson_L1" );
      declareProperty( "JSONMetaObjectNameHLTPS", m_metaNameJSON_hltps = "TriggerMenuJson_HLTPS" );
      declareProperty( "JSONMetaObjectNameL1PS", m_metaNameJSON_l1ps = "TriggerMenuJson_L1PS" );
      declareProperty( "JSONMetaObjectNameBunchgroup", m_metaNameJSON_bg = "TriggerMenuJson_BG" );
   }

   // Defined here rather than in the header because it requires the knowledge of Impl defined above
   xAODConfigTool::~xAODConfigTool() = default;

   StatusCode xAODConfigTool::initialize() {

      // Greet the user:
      ATH_MSG_INFO( "Initialising..." );
      ATH_MSG_DEBUG( "EventObjectName = " << m_eventName );
      ATH_MSG_DEBUG( "-- Run 2 AOD Configuration Settings");
      ATH_MSG_DEBUG( "MetaObjectName  = " << m_metaName_run2 );
      ATH_MSG_DEBUG( "-- Run 3 AOD Configuration Settings");
      ATH_MSG_DEBUG( "JSONMetaObjectNameHLT = " << m_metaNameJSON_hlt );
      ATH_MSG_DEBUG( "JSONMetaObjectNameHLTMonitoring = " << m_metaNameJSON_hltmonitoring );
      ATH_MSG_DEBUG( "JSONMetaObjectNameL1 = " << m_metaNameJSON_l1 );
      ATH_MSG_DEBUG( "JSONMetaObjectNameHLTPS = " << m_metaNameJSON_hltps );
      ATH_MSG_DEBUG( "JSONMetaObjectNameL1PS = " << m_metaNameJSON_l1ps );
      ATH_MSG_DEBUG( "JSONMetaObjectNameBunchgroup = " << m_metaNameJSON_bg );

      // Reset the pointers:
      m_tmc = nullptr;
      m_menu = nullptr;
      //
      m_hltJson = nullptr;
      m_hltmonitoringJson = nullptr;
      m_l1Json = nullptr;
      m_hltpsJson = nullptr;
      m_l1psJson = nullptr;
      m_bgJson = nullptr;
      m_currentHltJson = nullptr;
      m_currentHltmonitoringJson = nullptr;
      m_currentL1Json = nullptr;
      m_currentHltpsJson = nullptr;
      m_currentL1psJson = nullptr;
      m_currentBgJson = nullptr;

      // Return gracefully:
      return StatusCode::SUCCESS;
   }

   const CTPConfig* xAODConfigTool::ctpConfig() const {
      // Check if the object is well prepared:
      if( m_impl->m_ctpConfig.menu().size() == 0 ) {
         ATH_MSG_FATAL( "Trigger menu not loaded" );
         throw std::runtime_error( "Tool not initialised correctly" );
      }

      // Return the pointer:
      return &m_impl->m_ctpConfig;
   }

   const BunchGroupSet* xAODConfigTool::bunchGroupSet() const {
      // Check if the object is well prepared:
      if( m_impl->m_bgSet.bunchGroups().size() == 0 ) {
         ATH_MSG_FATAL( "Trigger menu not loaded" );
         throw std::runtime_error( "Tool not initialised correctly" );
      }

      // Return the pointer:
      return &m_impl->m_bgSet;
   }

   uint32_t xAODConfigTool::lvl1PrescaleKey() const {
      if ( m_menuJSONContainerAvailable ) {

         return m_currentL1psJson->key();

      } else {

         // Check if the object is well prepared:
         if( ! m_menu ) {
            ATH_MSG_FATAL( "Trigger menu not loaded" );
            throw std::runtime_error( "Tool not initialised correctly" );
         }

         // Return the key from the metadata object:
         return m_menu->l1psk();

      }
   }

   const HLTChainList& xAODConfigTool::chains() const {
      // Check if the object is well prepared:
      if( m_impl->m_chainList.size() == 0 ) {
         ATH_MSG_FATAL( "Trigger menu not loaded" );
         throw std::runtime_error( "Tool not initialised correctly" );
      }

      // Return the object:
      return m_impl->m_chainList;
   }

   const HLTSequenceList& xAODConfigTool::sequences() const {
      // Check if the object is well prepared:
      if( m_impl->m_sequenceList.size() == 0 ) {
         ATH_MSG_FATAL( "Trigger menu not loaded" );
         throw std::runtime_error( "Tool not initialised correctly" );
      }

      // Return the object:
      return m_impl->m_sequenceList;
   }

   uint32_t xAODConfigTool::masterKey() const {
      if (m_menuJSONContainerAvailable) {

         return m_currentHltJson->key();

      } else {

         // Check if the object is well prepared:
         if( ! m_menu ) {
            ATH_MSG_FATAL( "Trigger menu not loaded" );
            throw std::runtime_error( "Tool not initialised correctly" );
         }

         // Return the key from the metadata object:
         return m_menu->smk();

      }
   }

   uint32_t xAODConfigTool::hltPrescaleKey() const {
      if (m_menuJSONContainerAvailable) {

         return m_currentHltpsJson->key();

      } else {

         // Check if the object is well prepared:
         if( ! m_menu ) {
            ATH_MSG_FATAL( "Trigger menu not loaded" );
            throw std::runtime_error( "Tool not initialised correctly" );
         }

         // Return the key from the metadata object:
         return m_menu->hltpsk();

      }
   }

   const HLTMenu& xAODConfigTool::hltMenu(const EventContext&) const {
      if( ! m_menuJSONContainerAvailable ) {
         ATH_MSG_FATAL( "Run 3 format Trigger menu not loaded" );
         throw std::runtime_error( "Tool not initialised correctly" );
      }
      return m_impl->m_currentHlt;
   }

   const HLTMonitoring& xAODConfigTool::hltMonitoring(const EventContext&) const {
      if( ! m_menuJSONContainerAvailable ) {
         ATH_MSG_FATAL( "Run 3 format Trigger menu not loaded" );
         throw std::runtime_error( "Tool not initialised correctly" );
      }
      return m_impl->m_currentHltmonitoring;
   }

   const L1Menu& xAODConfigTool::l1Menu(const EventContext&) const {
      if( ! m_menuJSONContainerAvailable ) {
         ATH_MSG_FATAL( "Run 3 format Trigger menu not loaded" );
         throw std::runtime_error( "Tool not initialised correctly" );
      }
      return m_impl->m_currentL1;
   }

   const HLTPrescalesSet& xAODConfigTool::hltPrescalesSet(const EventContext&) const {
      if( ! m_menuJSONContainerAvailable ) {
         ATH_MSG_FATAL( "Run 3 format Trigger menu not loaded" );
         throw std::runtime_error( "Tool not initialised correctly" );
      }
      return m_impl->m_currentHltps;
   }

   const L1PrescalesSet& xAODConfigTool::l1PrescalesSet(const EventContext&) const {
      if( ! m_menuJSONContainerAvailable ) {
         ATH_MSG_FATAL( "Run 3 format Trigger menu not loaded" );
         throw std::runtime_error( "Tool not initialised correctly" );
      }
      return m_impl->m_currentL1ps;
   }

   const L1BunchGroupSet& xAODConfigTool::l1BunchGroupSet(const EventContext&) const {
      if( ! m_menuJSONContainerAvailable ) {
         ATH_MSG_FATAL( "Run 3 format Trigger menu not loaded" );
         throw std::runtime_error( "Tool not initialised correctly" );
      }
      return m_impl->m_currentBg;
   }

   StatusCode xAODConfigTool::beginInputFile() {

      // Tell the user what's happening:
      ATH_MSG_DEBUG( "Loading the trigger menu from a new input file" );

      // Try to read the R2 metadata object:
      m_tmc = nullptr;
      m_triggerMenuContainerAvailable = true;
      if( !inputMetaStore()->contains<xAOD::TriggerMenuContainer>(m_metaName_run2)
          or inputMetaStore()->retrieve( m_tmc, m_metaName_run2 ).isFailure() )
      {
         m_triggerMenuContainerAvailable = false;
      }

      // Try to read the R3 metadata object:
      m_hltJson = nullptr;
      m_hltmonitoringJson = nullptr;
      m_l1Json = nullptr;
      m_hltpsJson = nullptr;
      m_l1psJson = nullptr;
      m_bgJson = nullptr;
      m_menuJSONContainerAvailable = true;
      if( !inputMetaStore()->contains<xAOD::TriggerMenuJsonContainer>(m_metaNameJSON_hlt)
          or inputMetaStore()->retrieve( m_hltJson, m_metaNameJSON_hlt ).isFailure() )
      {
         m_menuJSONContainerAvailable = false;
      }
      if( !inputMetaStore()->contains<xAOD::TriggerMenuJsonContainer>(m_metaNameJSON_hltmonitoring)
          or inputMetaStore()->retrieve( m_hltmonitoringJson, m_metaNameJSON_hltmonitoring ).isFailure() )
      {
         // m_menuJSONContainerAvailable = false;
         // Currently planning on only storing these data in MC. Hence this has to be an optional input.
      }
      if( !inputMetaStore()->contains<xAOD::TriggerMenuJsonContainer>(m_metaNameJSON_l1)
          or inputMetaStore()->retrieve( m_l1Json, m_metaNameJSON_l1 ).isFailure() )
      {
         m_menuJSONContainerAvailable = false;
      }
      if( !inputMetaStore()->contains<xAOD::TriggerMenuJsonContainer>(m_metaNameJSON_hltps)
          or inputMetaStore()->retrieve( m_hltpsJson, m_metaNameJSON_hltps ).isFailure() )
      {
         m_menuJSONContainerAvailable = false;
      }
      if( !inputMetaStore()->contains<xAOD::TriggerMenuJsonContainer>(m_metaNameJSON_l1ps)
          or inputMetaStore()->retrieve( m_l1psJson, m_metaNameJSON_l1ps ).isFailure() )
      {
         m_menuJSONContainerAvailable = false;
      }
      if( !inputMetaStore()->contains<xAOD::TriggerMenuJsonContainer>(m_metaNameJSON_bg)
          or inputMetaStore()->retrieve( m_bgJson, m_metaNameJSON_bg ).isFailure() ) {
         // m_menuJSONContainerAvailable = false;
         // This was not written up to the end of 2021, we have to make it optional for at least a while
      }


      if( !m_triggerMenuContainerAvailable && !m_menuJSONContainerAvailable ) {
         ATH_MSG_WARNING( "No trigger configurations are available on "
                          "the input" );
         return StatusCode::SUCCESS;
      }

      // Prefer run three format, if available
      if (m_menuJSONContainerAvailable) {

         // A little sanity check:
         if( m_hltJson->size() == 0 ) {
            // This can happen when we encounter empty input files. In which
            // case we should not bail, but continue, and only bail if by the
            // start of an event, we still don't see any configurations.
            ATH_MSG_WARNING( "No trigger configurations are available on "
                             "the input" );
            return StatusCode::SUCCESS;
         }

         // Load the first elements by default
         ATH_CHECK( m_hltJson->size() > 0 );
         //ATH_CHECK( m_hltmonitoringJson->size() > 0 ); // Optional - re-enable this check in the future when all used AODs contains this
         ATH_CHECK( m_l1Json->size() > 0 );
         ATH_CHECK( m_hltpsJson->size() > 0 );
         ATH_CHECK( m_l1psJson->size() > 0 );
         if (m_bgJson) { // Considered optional for now
            ATH_CHECK( m_bgJson->size() > 0 );
         }

         m_currentHltJson = m_hltJson->at( 0 );
         if (m_hltmonitoringJson && m_hltmonitoringJson->size()) { // Optional
            m_currentHltmonitoringJson = m_hltmonitoringJson->at( 0 );
         }
         m_currentL1Json = m_l1Json->at( 0 );
         m_currentHltpsJson = m_hltpsJson->at( 0 );
         m_currentL1psJson = m_l1psJson->at( 0 );
         if (m_bgJson) {  // Considered optional for now
            m_currentBgJson = m_bgJson->at( 0 );
         }

         ATH_CHECK( loadPtrees() );

         // Loading the payload doesn't additionally let the object know about its own key. We can load this in now too.
         m_impl->m_currentHlt.setSMK( m_currentHltJson->key() );
         if (m_currentHltmonitoringJson) { // Optional
            m_impl->m_currentHltmonitoring.setSMK( m_currentHltmonitoringJson->key() );
         }
         m_impl->m_currentL1.setSMK( m_currentL1Json->key() );
         m_impl->m_currentHltps.setPSK( m_currentHltpsJson->key() );
         m_impl->m_currentL1ps.setPSK( m_currentL1psJson->key() );
         if (m_bgJson) { // Optional (for now)
            m_impl->m_currentBg.setBGSK( m_currentBgJson->key() );
         }

         ATH_CHECK( prepareTriggerMenu( m_impl->m_currentHlt,
                                        m_impl->m_currentL1,
                                        m_impl->m_currentHltps,
                                        m_impl->m_currentL1ps,
                                        m_impl->m_currentBg,
                                        m_impl->m_ctpConfig,
                                        m_impl->m_chainList,
                                        m_impl->m_sequenceList,
                                        m_impl->m_bgSet, msg() ) );

         return StatusCode::SUCCESS;

      } else if (m_triggerMenuContainerAvailable) {

         // A little sanity check:
         if( m_tmc->size() == 0 ) {
            // This can happen when we encounter empty input files. In which
            // case we should not bail, but continue, and only bail if by the
            // start of an event, we still don't see any configurations.
            ATH_MSG_WARNING( "No trigger configurations are available on "
                             "the input" );
            return StatusCode::SUCCESS;
         }

         m_menu = m_tmc->at( 0 );
         // Cache the menu's configuration:
         ATH_CHECK( prepareTriggerMenu( m_menu, m_impl->m_ctpConfig,
                                        m_impl->m_chainList,
                                        m_impl->m_sequenceList,
                                        m_impl->m_bgSet, msg() ) );

         return StatusCode::SUCCESS;
      }

      ATH_MSG_ERROR( "Both m_menuJSONContainerAvailable and m_triggerMenuContainerAvailable are false");
      return StatusCode::FAILURE; // Should never get here as checked that one or the other is true above
   }

   StatusCode xAODConfigTool::beginEvent() {

      // It may be that the input file opening event is missed in standalone
      // mode at the very beginning of the application. (Depending on the
      // creation order of the objects.) So make sure that we have the
      // configuration objects at hand...
      if( !m_triggerMenuContainerAvailable and !m_menuJSONContainerAvailable ) {
         ATH_CHECK( beginInputFile() );
      }

      // Read the current event's trigger keys:
      const xAOD::TrigConfKeys* keys = nullptr;
      ATH_CHECK( evtStore()->retrieve( keys, m_eventName ) );

      const xAOD::BunchConfKey* bgKey = nullptr; // The BG key is currently optional, only files written from late 2021 will contain this 
      if (evtStore()->contains<xAOD::BunchConfKey>(m_bgkeysName)) {
         ATH_CHECK(evtStore()->retrieve(bgKey, m_bgkeysName));
      }

      // Prefer Run 3 data if available
      if (m_menuJSONContainerAvailable) {
         return beginEvent_Run3(keys, bgKey);
      } else if (m_triggerMenuContainerAvailable) {
         return beginEvent_Run2(keys);
      }

      ATH_MSG_ERROR( "Both m_menuJSONContainerAvailable and m_triggerMenuContainerAvailable are false");
      return StatusCode::FAILURE;

   }

   StatusCode xAODConfigTool::beginEvent_Run2(const xAOD::TrigConfKeys* keys) {
      // Check if we have the correct menu already:
      if( m_menu && xAODKeysMatch( keys, m_menu ) ) {
         return StatusCode::SUCCESS;
      }

      // If not, let's look for the correct configuration:
      xAOD::TriggerMenuContainer::const_iterator menu_itr = m_tmc->begin();
      xAOD::TriggerMenuContainer::const_iterator menu_end = m_tmc->end();
      for( ; menu_itr != menu_end; ++menu_itr ) {
         // Check if this is the menu we're looking for:
         if( ! xAODKeysMatch( keys, *menu_itr ) ) continue;
         // Remember it's pointer:
         m_menu = *menu_itr;
         // Cache the menu's configuration:
         ATH_CHECK( prepareTriggerMenu( m_menu, m_impl->m_ctpConfig,
                                        m_impl->m_chainList,
                                        m_impl->m_sequenceList,
                                        m_impl->m_bgSet, msg() ) );
         // We're done:
         return StatusCode::SUCCESS;
      }

      // Apparently we didn't find the correct menu!
      ATH_MSG_ERROR( "Couldn't find configuration for current event (SMK:"
                     << keys->smk() << ", L1PSK:" << keys->l1psk()
                     << ", HLTPSK:" << keys->hltpsk() << ")" );
      return StatusCode::FAILURE;
   }

   StatusCode xAODConfigTool::beginEvent_Run3(const xAOD::TrigConfKeys* keys, const xAOD::BunchConfKey* bgKey) {
      if (keys==nullptr) {
         ATH_MSG_ERROR("nullptr TrigConfKeys");
         return StatusCode::FAILURE;
      }

      // Check if we have the correct menu already:
      bool validConfig = true;
      if (m_currentHltJson->key() != keys->smk()) {
         validConfig = false;
      }
      // m_currentHltminitoringJson is checked by the m_currentHltJson check
      if (m_currentL1Json->key() != keys->smk()) {
         validConfig = false;
      }
      if (m_currentHltpsJson->key() != keys->hltpsk()) {
         validConfig = false;
      }
      if (m_currentL1psJson->key() != keys->l1psk()) {
         validConfig = false;
      }
      if (m_bgJson && m_currentBgJson && bgKey && m_currentBgJson->key() != static_cast<unsigned int>(bgKey->id())) {
          validConfig = false;
      }

      if (validConfig) {
         return StatusCode::SUCCESS;
      }

      // If not, load correct JSON menus from their respective containers, matching against the event's keys ...
      ATH_CHECK( loadJsonByKey("HLT Menu", m_hltJson, keys->smk(), m_currentHltJson) );
      if (m_hltmonitoringJson && m_hltmonitoringJson->size()) {
         ATH_CHECK( loadJsonByKey("HLT Monitoring", m_hltmonitoringJson, keys->smk(), m_currentHltmonitoringJson) );
      }
      ATH_CHECK( loadJsonByKey("L1 Menu", m_l1Json, keys->smk(), m_currentL1Json) );
      ATH_CHECK( loadJsonByKey("HLT Prescales", m_hltpsJson, keys->hltpsk(), m_currentHltpsJson) );
      ATH_CHECK( loadJsonByKey("L1 Prescales", m_l1psJson, keys->l1psk(), m_currentL1psJson) );
      if (m_bgJson && bgKey) {
         ATH_CHECK( loadJsonByKey("Bunchgroups", m_bgJson, bgKey->id(), m_currentBgJson) );
      }

      // ... and from these serialised JSON strings, populate the ptree data structures...
      ATH_CHECK( loadPtrees() ); 

      // ... and set the keys of the objects in the objects (not part of the JSON pyaload) ...
      m_impl->m_currentHlt.setSMK( m_currentHltJson->key() );
      if (m_currentHltmonitoringJson) {
         m_impl->m_currentHltmonitoring.setSMK( m_currentHltmonitoringJson->key() );
      }
      m_impl->m_currentL1.setSMK( m_currentL1Json->key() );
      m_impl->m_currentHltps.setPSK( m_currentHltpsJson->key() );
      m_impl->m_currentL1ps.setPSK( m_currentL1psJson->key() );
      if (m_bgJson && bgKey) {
         m_impl->m_currentBg.setBGSK( m_currentBgJson->key() );
      }

      // R3 interfaces now active

      // ... and finally populate the legacy interface from the ptree data
      ATH_CHECK( prepareTriggerMenu( m_impl->m_currentHlt,
                                     m_impl->m_currentL1,
                                     m_impl->m_currentHltps,
                                     m_impl->m_currentL1ps,
                                     m_impl->m_currentBg,
                                     m_impl->m_ctpConfig,
                                     m_impl->m_chainList,
                                     m_impl->m_sequenceList,
                                     m_impl->m_bgSet, msg() ) ); // R2 interfaces now active

      return StatusCode::SUCCESS;
   }

   StatusCode xAODConfigTool::loadJsonByKey(const std::string& humanName,
      const xAOD::TriggerMenuJsonContainer* metaContainer,
      const uint32_t keyToCheck,
      const xAOD::TriggerMenuJson*& ptrToSet)
   {
      xAOD::TriggerMenuJsonContainer::const_iterator menu_itr = metaContainer->begin();
      xAOD::TriggerMenuJsonContainer::const_iterator menu_end = metaContainer->end();
      for( ; menu_itr != menu_end; ++menu_itr ) {
         // Check if this is the menu we're looking for:
         if( keyToCheck != (*menu_itr)->key() ) continue;
         ptrToSet = *menu_itr;
         return StatusCode::SUCCESS;
      }

      ATH_MSG_FATAL("Couldn't find configuration for current event"
         << ", Requested key=" << keyToCheck
         << ", Requested menu=" << humanName);
      return StatusCode::FAILURE;
   }

   StatusCode xAODConfigTool::loadPtrees() {
      ATH_CHECK( loadPtree("HLT Menu", m_currentHltJson, m_impl->m_currentHlt) );
      if (m_currentHltmonitoringJson) {
         ATH_CHECK( loadPtree("HLT Monitoring", m_currentHltmonitoringJson, m_impl->m_currentHltmonitoring) );
      }
      ATH_CHECK( loadPtree("L1 Menu", m_currentL1Json, m_impl->m_currentL1) );
      ATH_CHECK( loadPtree("HLT Prescales", m_currentHltpsJson, m_impl->m_currentHltps) );
      ATH_CHECK( loadPtree("L1 Prescales", m_currentL1psJson, m_impl->m_currentL1ps) );
      if (m_bgJson) {
         ATH_CHECK( loadPtree("Bunchgroups", m_currentBgJson, m_impl->m_currentBg) );
      }
      return StatusCode::SUCCESS;
   }

   StatusCode xAODConfigTool::loadPtree(const std::string& humanName,
      const xAOD::TriggerMenuJson* menu,
      DataStructure& dataStructure)
   {
      std::stringstream rawData;
      rawData << menu->payload();
      dataStructure.clear();
      try {
         boost::property_tree::ptree pt;
         boost::property_tree::read_json(rawData, pt);
         dataStructure.setData(std::move(pt));
      } catch (const boost::property_tree::json_parser_error& e) {
         ATH_MSG_FATAL("Unable to decode a JSON trigger menu metadata payload for " << humanName << " with key " << menu->key());
         ATH_MSG_FATAL(e.what());
         return StatusCode::FAILURE;
      }
      return StatusCode::SUCCESS;
   }

} // namespace TrigConf
