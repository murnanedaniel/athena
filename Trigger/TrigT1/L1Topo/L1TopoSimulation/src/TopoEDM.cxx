// Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration

#include "./TopoEDM.h"

#include "xAODTrigger/L1TopoSimResults.h"

using namespace LVL1;

TopoEDM::TopoEDM(const std::string& type, const std::string& name, 
                                       const IInterface* parent) :
   AthAlgTool(type, name, parent)
{
   declareInterface<LVL1::IReadTopoEDM>( this );
}

TopoEDM::~TopoEDM()
{}

StatusCode
TopoEDM::initialize() {
  ATH_CHECK( m_legacyL1topoKey.initialize() );
  ATH_CHECK( m_l1topoKey.initialize() );
  
  return StatusCode::SUCCESS;
}


StatusCode
TopoEDM::Read(bool isLegacy) const {
  const xAOD::L1TopoSimResults* l1topo_dec = 0;
  
  SG::ReadHandle<xAOD::L1TopoSimResultsContainer> cont(isLegacy? m_legacyL1topoKey : m_l1topoKey);
  if(!cont.isValid()){
    ATH_MSG_FATAL("Could not retrieve EDM Container");
    return StatusCode::FAILURE;
  }

  ATH_MSG_DEBUG("----got container: " << cont.key());

  for(const xAOD::L1TopoSimResults* it : * cont){
    l1topo_dec = it;
    ATH_MSG_DEBUG( "Reading L1Topo EDM:: BoardName: " << l1topo_dec->boardName() << " Clock: " << l1topo_dec->clock() << " Word32: " << l1topo_dec->word32() << " Word64: " << l1topo_dec->word64() << " WordOpt: " << l1topo_dec->wordOptical() );
  }

  return StatusCode::SUCCESS;
}
