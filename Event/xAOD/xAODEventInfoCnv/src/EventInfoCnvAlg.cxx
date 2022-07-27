/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/


// System include(s):
#include <memory>

// Gaudi/Athena include(s):
#include "AthenaKernel/errorcheck.h"
#include "GaudiKernel/ConcurrencyFlags.h"

// EDM include(s):
#include "EventInfo/EventInfo.h"
#include "EventInfo/PileUpEventInfo.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODEventInfo/EventAuxInfo.h"
#include "xAODEventInfo/EventInfoContainer.h"
#include "xAODEventInfo/EventInfoAuxContainer.h"

// Local include(s):
#include "EventInfoCnvAlg.h"

namespace xAODMaker {

   EventInfoCnvAlg::EventInfoCnvAlg( const std::string& name,
                                     ISvcLocator* svcLoc )
      : AthReentrantAlgorithm( name, svcLoc ),
        m_cnvTool( "xAODMaker::EventInfoCnvTool/EventInfoCnvTool", this ) {

      declareProperty( "AODKey", m_aodKey = "" );
      declareProperty( "xAODKey", m_xaodKey = "EventInfo" );
      declareProperty( "PileupKey", m_pileupKey = "" );
      declareProperty( "CnvTool", m_cnvTool );
   }

   StatusCode EventInfoCnvAlg::initialize() {

      // Greet the user:
      ATH_MSG_INFO( "Initializing " << name() );
      ATH_MSG_DEBUG( " AOD Key: " << m_aodKey );
      ATH_MSG_DEBUG( "xAOD Key: " << m_xaodKey.key() );

      // Retrieve the converter tool:
      CHECK( m_cnvTool.retrieve() );

      if (m_pileupKey.key().empty()) {
        m_pileupKey = "Pileup" + m_xaodKey.key();
      }

      CHECK( m_aodKey.initialize(SG::AllowEmpty) );
      CHECK( m_xaodKey.initialize() );

      /// FIXME: Should do this only if we have a PileUpEventInfo.
      CHECK( m_pileupKey.initialize(!m_pileupKey.key().empty()) );

      // Return gracefully:
      return StatusCode::SUCCESS;
   }

   StatusCode EventInfoCnvAlg::execute (const EventContext& ctx) const {

      // Check if anything needs to be done:
      // FIXME: Job configuration should be fixed so we don't need this test.
      if( evtStore()->contains< xAOD::EventInfo >( m_xaodKey.key() ) ) {
        ATH_MSG_WARNING( "xAOD::EventInfo with key \"" << m_xaodKey.key()
                          << "\" is already in StoreGate; "
                          << "EventInfoCnvAlg should not be scheduled.");
         return StatusCode::SUCCESS;
      }

      // Retrieve the AOD object:
      // FIXME: Use a ReadHandle.
      const EventInfo* aod = nullptr;
      if( m_aodKey.empty() ) {
         // If key has not been set, do a keyless retrieve instead.
         // This is not standard behavior, but is for compatibility
         // with existing configurations.
         CHECK( evtStore()->retrieve( aod ) );
      } else {
         SG::ReadHandle<EventInfo> ei (m_aodKey, ctx);
         aod = ei.cptr();
      }

      // Create the xAOD object(s):
      auto ei = std::make_unique<xAOD::EventInfo>();
      auto ei_aux = std::make_unique<xAOD::EventAuxInfo>();
      ei->setStore (ei_aux.get());

      // Check if this is a PileUpEventInfo object:
      const PileUpEventInfo* paod =
         dynamic_cast< const PileUpEventInfo* >( aod );
      if( paod ) {
        // Create an EventInfoContainer for the pileup events:
        auto puei = std::make_unique<xAOD::EventInfoContainer>();
        auto puei_aux = std::make_unique<xAOD::EventInfoAuxContainer>();
        puei->setStore (puei_aux.get());

        // Sub-events for the main EventInfo object:
        std::vector< xAOD::EventInfo::SubEvent > subEvents;

        // A convenience type declaration:
        typedef ElementLink< xAOD::EventInfoContainer > EiLink;

        // Create xAOD::EventInfo objects for each pileup EventInfo object:
        PileUpEventInfo::SubEvent::const_iterator pu_itr = paod->beginSubEvt();
        PileUpEventInfo::SubEvent::const_iterator pu_end = paod->endSubEvt();
        for( ; pu_itr != pu_end; ++pu_itr ) {
          // Create a new xAOD object:
          xAOD::EventInfo* ei = new xAOD::EventInfo();
          puei->push_back( ei );
          // Fill it with information:
          CHECK( m_cnvTool->convert( pu_itr->pSubEvt, ei, true, false, ctx ) );
          // And now add a sub-event to the temporary list:
          xAOD::EventInfo::PileUpType type = xAOD::EventInfo::Unknown;
          switch (pu_itr->type()) {
          case PileUpTimeEventIndex::Signal:
            type = xAOD::EventInfo::Signal;
            break;
          case PileUpTimeEventIndex::MinimumBias:
            type = xAOD::EventInfo::MinimumBias;
            break;
          case PileUpTimeEventIndex::Cavern:
            type = xAOD::EventInfo::Cavern;
            break;
          case PileUpTimeEventIndex::HaloGas:
            type = xAOD::EventInfo::HaloGas;
            break;
          case PileUpTimeEventIndex::HighPtMinimumBias:
            type = xAOD::EventInfo::HighPtMinimumBias;
            break;
          case PileUpTimeEventIndex::ZeroBias:
            type = xAOD::EventInfo::ZeroBias;
            break;
          default:
            break;
          }
          subEvents.emplace_back( pu_itr->time(),
                                  pu_itr->index(),
                                  type,
                                  EiLink( m_pileupKey.key(),
                                          puei->size() -
                                          1 ) );
        }

        // And now update the main EventInfo object with the sub-events:
        ei->setSubEvents( subEvents );

        // Record PU objects.
        SG::WriteHandle<xAOD::EventInfoContainer> puei_h (m_pileupKey, ctx);
        CHECK( puei_h.record (std::move(puei), std::move(puei_aux)) );
      }

      // Do the translation:
      CHECK( m_cnvTool->convert( aod, ei.get(), false, true, ctx ) );

      // Record EI objects.
      SG::WriteHandle<xAOD::EventInfo> ei_h (m_xaodKey, ctx);
      CHECK( ei_h.record (std::move(ei), std::move (ei_aux)) );

      // Return gracefully:
      return StatusCode::SUCCESS;
   }

} // namespace xAODMaker
