/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// $Id: CMMCPHitsCnvAlg.cxx 576052 2013-12-18 09:13:50Z morrisj $

// Gaudi/Athena include(s):
#include "AthenaKernel/errorcheck.h"

// EDM include(s):
#include "TrigT1CaloEvent/CMMCPHitsCollection.h"
#include "xAODTrigL1Calo/CMMCPHitsContainer.h"
#include "xAODTrigL1Calo/CMMCPHitsAuxContainer.h"

// Local include(s):
#include "CMMCPHitsCnvAlg.h"

namespace xAODMaker {

   CMMCPHitsCnvAlg::CMMCPHitsCnvAlg( const std::string& name,
                                     ISvcLocator* svcLoc )
      : AthAlgorithm( name, svcLoc ),
        m_cnvTool( "xAODMaker::CMMCPHitsCnvTool/CMMCPHitsCnvTool", this ) {

      declareProperty( "ESDKey" , m_esdKey  = "CMMCPHits" );
      declareProperty( "xAODKey", m_xaodKey = "CMMCPHits" );
      declareProperty( "CnvTool", m_cnvTool );
   }

   StatusCode CMMCPHitsCnvAlg::initialize() {

      // Greet the user:
      ATH_MSG_INFO( "Initializing - Package version: " << PACKAGE_VERSION );
      ATH_MSG_DEBUG( " ESD Key: " << m_esdKey );
      ATH_MSG_DEBUG( "xAOD Key: " << m_xaodKey );

      // Retrieve the converter tool:
      CHECK( m_cnvTool.retrieve() );

      // Return gracefully:
      return StatusCode::SUCCESS;
   }

   StatusCode CMMCPHitsCnvAlg::execute() {
      
      // Retrieve the ESD container:
      const CMMCPHitsCollection* esd = 0;
      CHECK( evtStore()->retrieve( esd, m_esdKey ) );
      
      // Create the xAOD container and its auxiliary store:
      xAOD::CMMCPHitsAuxContainer* aux = new xAOD::CMMCPHitsAuxContainer();
      xAOD::CMMCPHitsContainer* xaod = new xAOD::CMMCPHitsContainer();
      xaod->setStore( aux );
      
      // Fill the xAOD container:
      CHECK( m_cnvTool->convert( esd, xaod ) );     
      
      // Record the xAOD containers:
      CHECK( evtStore()->record( aux, m_xaodKey + "Aux." ) );
      CHECK( evtStore()->record( xaod, m_xaodKey ) );

      // Return gracefully:
      return StatusCode::SUCCESS;     
   }

} // namespace xAODMaker
