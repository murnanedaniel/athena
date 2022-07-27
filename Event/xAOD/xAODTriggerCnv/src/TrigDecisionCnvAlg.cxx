/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/


// Gaudi/Athena include(s):
#include "AthenaKernel/errorcheck.h"

// Local include(s):
#include "TrigDecisionCnvAlg.h"

namespace xAODMaker {

   TrigDecisionCnvAlg::TrigDecisionCnvAlg( const std::string& name,
                                           ISvcLocator* svcLoc )
      : AthReentrantAlgorithm( name, svcLoc ) {
   }


   TrigDecisionCnvAlg::~TrigDecisionCnvAlg() {
   }


   StatusCode TrigDecisionCnvAlg::initialize() {

       // Greet the user:
      ATH_MSG_INFO( "Initializing" );
      ATH_MSG_DEBUG( " AOD Key: " << m_aodKey );
      ATH_MSG_DEBUG( "xAOD Key: " << m_xaodKey );

      CHECK( m_eventInfoKey.initialize(SG::AllowEmpty) );
      CHECK( m_aodKey.initialize() );
      CHECK( m_xaodKey.initialize() );

      // Retrieve the converter tool:
      CHECK( m_cnvTool.retrieve() );

      // Return gracefully:
      return StatusCode::SUCCESS;
   }

   StatusCode TrigDecisionCnvAlg::execute(const EventContext& ctx) const {

      // Retrieve the AOD object:
      SG::ReadHandle<TrigDec::TrigDecision> aodRH{m_aodKey, ctx};
      CHECK( aodRH.isValid() );
      const TrigDec::TrigDecision* aod = aodRH.cptr();

      // trigger Info
      const TriggerInfo* trigInfo = nullptr;
      if (not m_eventInfoKey.empty()) {
	  SG::ReadHandle<EventInfo> eiRH{m_eventInfoKey, ctx};
	  trigInfo = eiRH.isValid() ? eiRH.cptr()->trigger_info() : nullptr;
      }

      // Create the xAOD object and its auxiliary store:
      std::unique_ptr<xAOD::TrigDecisionAuxInfo> aux = std::make_unique<xAOD::TrigDecisionAuxInfo>();
      std::unique_ptr<xAOD::TrigDecision> xaod = std::make_unique<xAOD::TrigDecision>();
      xaod->setStore( aux.get() );

      // Fill the xAOD object:
      CHECK( m_cnvTool->convert( aod, xaod.get(), trigInfo ) );

      // Record the xAOD objects:
      SG::WriteHandle<xAOD::TrigDecision> writeHandle{m_xaodKey, ctx};
      CHECK( writeHandle.record(std::move(xaod), std::move(aux)));

      // Return gracefully:
      return StatusCode::SUCCESS;
   }

} // namespace xAODMaker
