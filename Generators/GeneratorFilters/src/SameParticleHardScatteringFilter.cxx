/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "GeneratorFilters/SameParticleHardScatteringFilter.h"

SameParticleHardScatteringFilter::SameParticleHardScatteringFilter(const std::string& name, ISvcLocator* pSvcLocator)
  : GenFilter(name,pSvcLocator)
{
  declareProperty("PDGParent",m_PDGParent);
  declareProperty("PDGChild", m_PDGChild);
}

StatusCode SameParticleHardScatteringFilter::filterInitialize() {
 if (m_PDGParent.size() == 0) ATH_MSG_ERROR("PDGParent[] not set ");
 if (m_PDGChild.size() == 0) ATH_MSG_ERROR("PDGChild[] not set ");
 for (int i=0; i < int(m_PDGParent.size()); i++) ATH_MSG_DEBUG("PDGParent["<<i<<"] = " << m_PDGParent[i]);
 for (int i=0; i < int(m_PDGChild.size()); i++) ATH_MSG_DEBUG("PDGChild["<<i<<"] = " << m_PDGChild[i]);
 return StatusCode::SUCCESS;
}

StatusCode SameParticleHardScatteringFilter::filterEvent() {
  ATH_MSG_DEBUG(" SameParticleHardScattering filtering for: Parent --> " << m_PDGParent[0]
		<< " and parent " << -m_PDGParent[0]
                << ", Child --> " << m_PDGChild[0]);
  int N_Parent[2];
  N_Parent[0] = 0;
  N_Parent[1] = 0;

  for (const HepMC::GenEvent* genEvt : *events_const()) {
      for (auto pitr: *genEvt) 
	{
	  int id = pitr->pdg_id();
	  if (std::abs(id) != m_PDGChild[0]) continue; // searching for only b-quarks
	  
	  // a pointer to the production vertex
	  HepMC::ConstGenVertexPtr productionVtx = pitr->production_vertex();
	  
	  // Verify if we got a valid pointer and retrieve the number of parents
	  if (!productionVtx) continue;
	  // Incoming particle range check
#ifdef HEPMC3
	  if (productionVtx->particles_in().size() < 2) continue; //  we are looking for excited tau-leptons produced in b-quark b-antiquark scattering
	  for (auto thisParent:  productionVtx->particles_in()) {
#else
	  if (productionVtx->particles_in_size() < 2) continue; //  we are looking for excited tau-leptons produced in b-quark b-antiquark scattering
	  HepMC::GenVertex::particles_in_const_iterator firstParentIt = productionVtx->particles_in_const_begin();
	  HepMC::GenVertex::particles_in_const_iterator endParentIt = productionVtx->particles_in_const_end();
	  for (HepMC::GenVertex::particles_in_const_iterator thisParentIt = firstParentIt ; thisParentIt != endParentIt; ++thisParentIt) {
		auto thisParent= *thisParentIt;
#endif
	    ATH_MSG_DEBUG(" SelectBQuarkScattering Filter: parent ==> " <<thisParent->pdg_id() << " child ===> "  << pitr->pdg_id());
	    if ( thisParent->pdg_id()    == m_PDGParent[0] )
	      {
		N_Parent[0]++;
	      }
	    if ( thisParent->pdg_id()   == -m_PDGParent[0] )
	      {
		N_Parent[1]++;
	      }
	  }
	}
    }
  setFilterPassed(N_Parent[0] >= 1 && N_Parent[1] >= 1);
  return StatusCode::SUCCESS;
}

