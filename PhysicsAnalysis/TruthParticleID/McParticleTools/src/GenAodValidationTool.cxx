/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

/////////////////////////////////////////////////////////////////// 
// GenAodValidationTool.cxx 
// Implementation file for class GenAodValidationTool
// Author: S.Binet<binet@cern.ch>
/////////////////////////////////////////////////////////////////// 


// STL includes
#include <sstream>

// FrameWork includes
#include "Gaudi/Property.h"

// StoreGate
#include "StoreGate/StoreGateSvc.h"

// CLHEP/HepMC includes
#include "GeneratorObjects/McEventCollection.h"
#include "AtlasHepMC/GenEvent.h"
#include "AtlasHepMC/GenVertex.h"
#include "AtlasHepMC/GenParticle.h"
#include "AtlasHepMC/Flow.h"
#include "AtlasHepMC/Polarization.h"
// McParticleKernel includes
#include "McParticleKernel/IIOHepMcTool.h"

// McParticleTools includes
#include "GenAodValidationTool.h"

/////////////////////////////////////////////////////////////////// 
/// Public methods: 
/////////////////////////////////////////////////////////////////// 

/// Constructors
////////////////
GenAodValidationTool::GenAodValidationTool( const std::string& type, 
					    const std::string& name, 
					    const IInterface* parent ) : 
  TruthParticleValidationBaseTool( type, name, parent ),
  m_outFile( nullptr )
{
  //
  // Property declaration
  // 
  
  declareProperty( "RefMcEvents",   
		   m_refMcEventsName = "TruthEvent",
		   "StoreGate location of the reference McEventCollection" );

  declareProperty( "CheckMcEvents", 
		   m_checkMcEventsName = "GEN_AOD",
		   "Location of the 'to-be-validated' McEventCollection" );

  declareProperty( "HardScatteringVtxOutputFile",
		   m_hardVtxOutFileName = "hepmc.hard.vtx.txt",
		   "Name of the output text file which will contain the " 
		   "hard-scattering vertices" );

  declareProperty( "RefHepMcWriterTool",
		   m_refMcEventWriter = HepMcTool_t( "HepMcWriterTool", 
						     this ),
		   "Tool to write the reference HepMC::GenEvent into a dedicated file" );

  declareProperty( "CheckHepMcWriterTool",
		   m_checkMcEventWriter = HepMcTool_t( "HepMcWriterTool", 
						       this ),
		   "Tool to write the \"to-be-checked\" HepMC::GenEvent "
		   "into a dedicated file" );

}

/// Destructor
///////////////
GenAodValidationTool::~GenAodValidationTool()
{ 
  ATH_MSG_DEBUG("Calling destructor");

  if ( m_outFile ) {
    m_outFile->close();
    delete m_outFile;
    m_outFile = nullptr;
  }

}


/////////////////////////////////////////////////////////////////// 
/// Non-const methods: 
/////////////////////////////////////////////////////////////////// 

StatusCode GenAodValidationTool::initializeTool()
{
  ATH_MSG_INFO("Initializing " << name() << "...");

  if ( m_outFile ) {
    delete m_outFile;
    m_outFile = nullptr;
  }
  m_outFile = new std::ofstream( m_hardVtxOutFileName.value().c_str(),
				 std::ios::trunc );
  if ( !m_outFile || !m_outFile->is_open() || m_outFile->bad() ) {
    ATH_MSG_ERROR("Could not open output Hard-Scattering vertices file at : " 
		  << m_hardVtxOutFileName.value());
    if ( m_outFile ) {
      delete m_outFile;
      m_outFile = nullptr;
    }
    return StatusCode::FAILURE;
  }

  // configure hard-scattering vertices
  m_ppFilter.setDecayPattern( "-1|-2|-3|-4|-5|-6|21 + 1|2|3|4|5|6|21 -> " );
  
  // filter which will be used to not select shower vertices while
  // looking for hard-scattering vertices
  m_showerFilter.setDecayPattern( "-> 91|92|94" );


  // retrieve and configure HepMC writer tools
  if ( setupHepMcWriterTools().isFailure() ) {
    ATH_MSG_ERROR("Could not configure the HepMC::GenEvent writer tools !!");
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

StatusCode GenAodValidationTool::finalizeTool()
{
  ATH_MSG_INFO("Finalizing " << name() << "...");

  return StatusCode::SUCCESS;
}

StatusCode GenAodValidationTool::executeTool()
{
  // retrieve reference McEventCollection
  const McEventCollection * refMcEvents = nullptr;
  if ( evtStore()->retrieve( refMcEvents, m_refMcEventsName ).isFailure() ||
       nullptr == refMcEvents ) {
    ATH_MSG_ERROR("Could not retrieve reference McEventCollection at ["
		  << m_refMcEventsName << "] !!");
    return StatusCode::FAILURE;
  } else {
    ATH_MSG_VERBOSE("Successfully retrieved reference McEventCollection at ["
		    << m_refMcEventsName << "]");

    if ( !m_refMcEventWriter || 
	 m_refMcEventWriter->execute().isFailure() ) {
      ATH_MSG_WARNING("Failed to write the reference GenEvent !!");
    }
  }

  // retrieve checking McEventCollection
  const McEventCollection * checkMcEvents = nullptr;
  if ( evtStore()->retrieve(checkMcEvents, m_checkMcEventsName).isFailure() ||
       nullptr == checkMcEvents ) {
    ATH_MSG_ERROR("Could not retrieve 'to-be-validated' McEventCollection at ["
		  << m_checkMcEventsName << "] !!");
    return StatusCode::FAILURE;
  } else {
    ATH_MSG_VERBOSE
      ("Successfully retrieved 'to-be-validated' McEventCollection at ["
       << m_checkMcEventsName << "]");

    if ( !m_checkMcEventWriter ||
	 m_checkMcEventWriter->execute().isFailure() ) {
      ATH_MSG_WARNING("Failed to write the reference GenEvent !!");
    }
  }

  return executeTool( refMcEvents, checkMcEvents );
}

StatusCode 
GenAodValidationTool::executeTool( const McEventCollection* refMcEvents,
				   const McEventCollection* checkMcEvents )
{
  if ( nullptr == refMcEvents ) {
    ATH_MSG_ERROR("NULL pointer to reference McEventCollection !!");
    return StatusCode::FAILURE;
  }

  if ( nullptr == checkMcEvents ) {
    ATH_MSG_ERROR("NULL pointer to the 'to-be-validated' McEventCollection !!");
    return StatusCode::FAILURE;
  }

  if ( refMcEvents->size() != checkMcEvents->size() ) {
    ATH_MSG_ERROR("McEventCollections to be compared don't have the same size !"
		  << endmsg
		  << "Ref:     " << refMcEvents->size() << endmsg
		  << "Current: " << checkMcEvents->size());
  }

  const unsigned int nMaxEvts = std::min( refMcEvents->size(),
					  checkMcEvents->size() );
  for ( unsigned int iEvt = 0; iEvt != nMaxEvts; ++iEvt ) {

    if ( executeTool( (*refMcEvents)[iEvt], 
		      (*checkMcEvents)[iEvt] ).isFailure() ) {
      ATH_MSG_WARNING("Could not VALIDATE this event ["
		      << (*refMcEvents)[iEvt]->event_number()
		      << "] !!");
    }
  }

  return StatusCode::SUCCESS;
}

StatusCode 
GenAodValidationTool::executeTool( const HepMC::GenEvent* refMcEvts,
				   const HepMC::GenEvent* checkMcEvts )
{
  if ( nullptr == refMcEvts ) {
    ATH_MSG_ERROR("NULL pointer to reference HepMC::GenEvent !!");
    return StatusCode::FAILURE;
  }

  if ( nullptr == checkMcEvts ) {
    ATH_MSG_ERROR("NULL pointer to the 'to-be-validated' HepMC::GenEvent !!");
    return StatusCode::FAILURE;
  }

  const unsigned int evtNbr = refMcEvts->event_number();
  (*m_outFile) << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
	       << std::endl
	       << "Event: " << evtNbr
	       << std::endl;

#ifdef HEPMC3   
  // loop over reference vertices
  for ( const auto&  vtx: refMcEvts->vertices()) {
    if ( m_ppFilter.isAccepted(vtx) &&
	 !m_showerFilter.isAccepted(vtx) ) {
      HepMC::Print::line(*m_outFile,vtx);
      HepMC::ConstGenVertexPtr checkVtx = HepMC::barcode_to_vertex(checkMcEvts,HepMC::barcode(vtx));
      if ( !checkVtx ) {
	ATH_MSG_WARNING
	  ("Output GenEvent is missing the selected HardScattering Vtx !!"
	   << " (" << vtx << ")");
      } else {
	(*m_outFile) << "---------" << std::endl;
	 HepMC::Print::line(*m_outFile,checkVtx);
	if ( !compareVtx( vtx, checkVtx ) ) {
	  ATH_MSG_WARNING("Selected HardScattering vertices are NOT the same !!"
			  << " at Event [" << evtNbr << "]"
			  << " refVtx = " << vtx);
	}
      }
    }
  }
  (*m_outFile) << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
  // loop over slimmed HepMC::GenEvent and check that vertices
  // are comparable
  // loop over reference vertices
  for (const auto& checkVtx: checkMcEvts->vertices()) {
    HepMC::ConstGenVertexPtr refVtx = HepMC::barcode_to_vertex(refMcEvts,HepMC::barcode(checkVtx));
    if (!refVtx) {
      ATH_MSG_WARNING("In Event [" << evtNbr
		      << "]: got null ref-vertex ( " 
		      << checkVtx << ")");
      continue;
    }
    if ( !compareVtx( refVtx, checkVtx ) ) {
      ATH_MSG_WARNING("In Event [" << evtNbr
		      << "]: vertices are not the SAME (" 
		      << refVtx << ")");
      std::stringstream refVtxStr;
      HepMC::Print::line(refVtxStr,refVtx);
      std::stringstream checkVtxStr;
      HepMC::Print::line(checkVtxStr,checkVtx);
      msg(MSG::WARNING) << std::endl
			<< "######### Ref vertex:" << std::endl
			<< refVtxStr.str()
			<< std::endl
			<< "######### Check vertex:" << std::endl
			<< checkVtxStr.str()
			<< endmsg;
    }
  }  
#else
  // loop over reference vertices
  for ( HepMC::GenEvent::vertex_const_iterator vtx = refMcEvts->vertices_begin();
	vtx != refMcEvts->vertices_end(); 
	++vtx ) {
    if ( m_ppFilter.isAccepted(*vtx) &&
	 !m_showerFilter.isAccepted(*vtx) ) {
      (*vtx)->print(*m_outFile);
      const HepMC::GenVertex* checkVtx = checkMcEvts->barcode_to_vertex((*vtx)->barcode());
      if ( 0 == checkVtx ) {
	ATH_MSG_WARNING
	  ("Output GenEvent is missing the selected HardScattering Vtx !!"
	   << " (" << (*vtx)->barcode() << ")");
      } else {
	(*m_outFile) << "---------" << std::endl;
	checkVtx->print(*m_outFile);
	if ( !compareVtx( *vtx, checkVtx ) ) {
	  ATH_MSG_WARNING("Selected HardScattering vertices are NOT the same !!"
			  << " at Event [" << evtNbr << "]"
			  << " refVtx = " << (*vtx)->barcode());
	}
      }
    }
  }

  (*m_outFile) << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;

  // loop over slimmed HepMC::GenEvent and check that vertices
  // are comparable
  // loop over reference vertices
  for ( HepMC::GenEvent::vertex_const_iterator vtx = checkMcEvts->vertices_begin();
	vtx != checkMcEvts->vertices_end(); 
	++vtx ) {
    const HepMC::GenVertex* checkVtx = *vtx;
    const HepMC::GenVertex* refVtx = refMcEvts->barcode_to_vertex(checkVtx->barcode());
    if (0 == refVtx) {
      ATH_MSG_WARNING("In Event [" << evtNbr
		      << "]: got null ref-vertex (barcode: " 
		      << checkVtx->barcode() << ")");
      continue;
    }
    if ( !compareVtx( refVtx, checkVtx ) ) {
      ATH_MSG_WARNING("In Event [" << evtNbr
		      << "]: vertices are not the SAME (" 
		      << refVtx->barcode() << ")");
      std::stringstream refVtxStr;
      refVtx->print(refVtxStr);
      std::stringstream checkVtxStr;
      checkVtx->print(checkVtxStr);
      msg(MSG::WARNING) << std::endl
			<< "######### Ref vertex:" << std::endl
			<< refVtxStr.str()
			<< std::endl
			<< "######### Check vertex:" << std::endl
			<< checkVtxStr.str()
			<< endmsg;
    }
  }  
#endif

  return StatusCode::SUCCESS;
}

bool GenAodValidationTool::compareVtx( const HepMC::ConstGenVertexPtr& vtx1, const HepMC::ConstGenVertexPtr& vtx2 ) const
{
  if ( !vtx1 || !vtx2 ) {
    ATH_MSG_ERROR("One of vertices is a NULL pointer !!" << endmsg
		  << " vtx1: " << vtx1 << endmsg
		  << " vtx2: " << vtx2);
    return false;
  }

#ifdef HEPMC3 
  auto inVtx1 = vtx1->particles_in();
  auto inVtx2 = vtx2->particles_in();
  if (inVtx1.size() !=  inVtx2.size()) {
    ATH_MSG_ERROR("Not the same number of branches !!" << endmsg  << " in:  " << inVtx1.size()  << "\t" << inVtx2.size());
    return false;
  }
  std::sort(inVtx1.begin(), inVtx1.end(), [](auto& a, auto& b) -> bool{a->momentum().pz() > b->momentum().pz(); });
  std::sort(inVtx2.begin(), inVtx2.end(), [](auto& a, auto& b) -> bool{a->momentum().pz() > b->momentum().pz(); });
  for (size_t i=0; i<inVtx1.size(); ++i){
     if ( compareParts( inVtx1.at(i), inVtx2.at(i) ) ) continue;
     ATH_MSG_ERROR("In-going particles are NOT matching !!");
     return false;
  }
  auto outVtx1 = vtx1->particles_out();
  auto outVtx2 = vtx2->particles_out();
  if (outVtx1.size() !=  outVtx2.size()) {
    ATH_MSG_ERROR("Not the same number of branches !!" << endmsg  << " out:  " << outVtx1.size()  << "\t" << outVtx2.size());
    return false;
  }
  std::sort(outVtx1.begin(), outVtx1.end(), [](auto& a, auto& b) -> bool{a->momentum().pz() > b->momentum().pz(); });
  std::sort(outVtx2.begin(), outVtx2.end(), [](auto& a, auto& b) -> bool{a->momentum().pz() > b->momentum().pz(); });
  for (size_t i=0; i<outVtx1.size(); ++i){
     if ( compareParts( outVtx1.at(i), outVtx2.at(i) ) ) continue;
     ATH_MSG_ERROR("Out-going particles are NOT matchiutg !!");
     return false;
  }
#else
  const int inVtx1 = vtx1->particles_in_size();
  const int inVtx2 = vtx2->particles_in_size();

  const int outVtx1 = vtx1->particles_out_size();
  const int outVtx2 = vtx2->particles_out_size();
  
  if (  inVtx1 !=  inVtx2 || outVtx1 != outVtx2 ) {
    ATH_MSG_ERROR("Not the same number of branches !!" << endmsg
		  << " in:  " << inVtx1  << "\t" << inVtx2  << endmsg
		  << " out: " << outVtx1 << "\t" << outVtx2);
    return false;
  }

  for ( HepMC::GenVertex::particles_in_const_iterator inPart1 = vtx1->particles_in_const_begin();
	inPart1 != vtx1->particles_in_const_end();
	++inPart1 ) {
    bool inParticleOK = false;
    for ( HepMC::GenVertex::particles_in_const_iterator inPart2 = vtx2->particles_in_const_begin();
	inPart2 != vtx2->particles_in_const_end();
	++inPart2 ) {
      if ( compareParts( *inPart1, *inPart2 ) ) {
	inParticleOK = true;
	break;
      }
    } //> end loop over in-part2

    if ( !inParticleOK ) {
      ATH_MSG_ERROR("In-going particles are NOT matching !!");
      return false;
    }

  }//> end loop over in-part1


  for ( HepMC::GenVertex::particles_out_const_iterator outPart1 = vtx1->particles_out_const_begin();
	outPart1 != vtx1->particles_out_const_end();
	++outPart1 ) {
    bool outParticleOK = false;
    for ( HepMC::GenVertex::particles_out_const_iterator outPart2 = vtx2->particles_out_const_begin();
	outPart2 != vtx2->particles_out_const_end();
	++outPart2 ) {
      if ( compareParts( *outPart1, *outPart2 ) ) {
	outParticleOK = true;
	break;
      }
    } //> end loop over out-part2

    if ( !outParticleOK ) {
      ATH_MSG_ERROR("Out-going particles are NOT matching !!");
      return false;
    }

  }//> end loop over out-part1
#endif

  return true;
}

bool 
GenAodValidationTool::compareParts( const HepMC::ConstGenParticlePtr& p1, const HepMC::ConstGenParticlePtr& p2 ) const
{
  if ( !p1 || !p2 ) {
    ATH_MSG_ERROR("One of particles is a NULL pointer !!" << endmsg
		  << " p1: " << p1 << endmsg
		  << " p2: " << p2);
    return false;
  }

  if ( HepMC::barcode(p1) != HepMC::barcode(p2) ) {
    return false;
  }

  const HepMC::FourVector hlv1   = p1->momentum();
  const int id1                  = p1->pdg_id();
  const int status1              = p1->status();
  auto flow1                     = HepMC::flow(p1);
  auto pol1                      = HepMC::polarization(p1);

  const HepMC::FourVector hlv2   = p2->momentum();
  const int id2                  = p2->pdg_id();
  const int status2              = p2->status();
  auto flow2                     = HepMC::flow(p2);
  auto pol2                      = HepMC::polarization(p2);

  const bool isOK = ( hlv1    == hlv2    &&
		      id1     == id2     && 
		      status1 == status2 &&
		      flow1   == flow2   &&
		      pol1    == pol2    );


  if ( !isOK ) {
    ATH_MSG_WARNING 
      ("=============================" << endmsg
       << "\thlv:    " << std::boolalpha << (hlv1    == hlv2   ) << endmsg
       << "\tid:     " << std::boolalpha << (id1     == id2    ) << endmsg
       << "\tstatus: " << std::boolalpha << (status1 == status2) << endmsg
       << "\tflow:   " << std::boolalpha << (flow1   == flow2  ) << endmsg
       << "\tpolar:  " << std::boolalpha << (pol1    == pol2   ) << endmsg
       << "=============================");
  }

  return isOK;
}

StatusCode GenAodValidationTool::setupHepMcWriterTools()
{
  if ( !m_refMcEventWriter.retrieve().isSuccess() ) {
    ATH_MSG_ERROR("Creation of algTool ["
		  << m_refMcEventWriter.type() << "/" 
		  << "] FAILED !");
    return StatusCode::FAILURE;
  }

  if ( !m_checkMcEventWriter.retrieve().isSuccess() ) {
    ATH_MSG_ERROR("Creation of algTool ["
		  << m_checkMcEventWriter.type() << "/" 
		  << "] FAILED !");
    return StatusCode::FAILURE;
  }

  // now we configure the tools
  StringProperty refProp( "McEvents", m_refMcEventsName.value() );
  if ( m_refMcEventWriter->setProperty( refProp ).isFailure() ) {
    ATH_MSG_ERROR("Could not set property [" << refProp.name() 
		  << "] for tool [" << m_refMcEventWriter.type() << "/" 
		  << "] !");
    return StatusCode::FAILURE;
  }

  StringProperty checkProp( "McEvents", m_checkMcEventsName.value() );
  if ( m_checkMcEventWriter->setProperty( checkProp ).isFailure() ) {
    ATH_MSG_ERROR("Could not set property [" << checkProp.name()
		  << "] for tool [" << m_checkMcEventWriter.type() << "/" 
		  << "] !");
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}
