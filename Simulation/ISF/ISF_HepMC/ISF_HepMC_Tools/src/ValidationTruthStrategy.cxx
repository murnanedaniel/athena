/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
// ValidationTruthStrategy.cxx, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

// class header include
#include "ValidationTruthStrategy.h"

// ISF includes
#include "ISF_Event/ITruthIncident.h"
#include "ISF_Event/ISFParticle.h"

/** Constructor **/
ISF::ValidationTruthStrategy::ValidationTruthStrategy(const std::string& t, const std::string& n, const IInterface* p) :
  base_class(t,n,p),
  m_minParentP2(0.)
{
    // parent particle minimum momentum
    declareProperty("ParentMinP"        , m_minParentP2     );
    declareProperty("Regions"                   , m_regionListProperty );
}

/** Destructor **/
ISF::ValidationTruthStrategy::~ValidationTruthStrategy()
{
}

// Athena algtool's Hooks
StatusCode  ISF::ValidationTruthStrategy::initialize()
{
    ATH_MSG_VERBOSE("Initializing ...");

    // (*) setup parent particle cuts
    // -> compute p^2 for fast comparison
    m_minParentP2 *= m_minParentP2;

    for(auto region : m_regionListProperty.value()) {
      if(region < AtlasDetDescr::fFirstAtlasRegion || region >= AtlasDetDescr::fNumAtlasRegions) {
        ATH_MSG_ERROR("Unknown Region (" << region << ") specified. Please check your configuration.");
        return StatusCode::FAILURE;
      }
    }

    return StatusCode::SUCCESS;
}

StatusCode  ISF::ValidationTruthStrategy::finalize()
{
    ATH_MSG_VERBOSE("Finalizing ...");
    return StatusCode::SUCCESS;
}

bool ISF::ValidationTruthStrategy::pass( ITruthIncident& ti) const {

  // parent particle check
  bool pass =  ( ti.parentP2() >= m_minParentP2 );

  return pass;
}

bool ISF::ValidationTruthStrategy::appliesToRegion(unsigned short geoID) const
{
  return std::find( m_regionListProperty.begin(),
                    m_regionListProperty.end(),
                    geoID ) != m_regionListProperty.end();
}
