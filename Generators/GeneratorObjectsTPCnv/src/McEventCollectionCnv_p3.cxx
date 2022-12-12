///////////////////////// -*- C++ -*- /////////////////////////////

/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

// McEventCollectionCnv_p3.cxx
// Implementation file for class McEventCollectionCnv_p3
// Author: S.Binet<binet@cern.ch>
///////////////////////////////////////////////////////////////////


// STL includes
#include <utility>
#include <cmath>

// GeneratorObjectsTPCnv includes
#include "GeneratorObjectsTPCnv/McEventCollectionCnv_p3.h"
#include "HepMcDataPool.h"


///////////////////////////////////////////////////////////////////
// Constructors
///////////////////////////////////////////////////////////////////


McEventCollectionCnv_p3::McEventCollectionCnv_p3() :
  Base_t( )
{}

McEventCollectionCnv_p3::McEventCollectionCnv_p3( const McEventCollectionCnv_p3& rhs ) :
  Base_t( rhs )
{}

McEventCollectionCnv_p3&
McEventCollectionCnv_p3::operator=( const McEventCollectionCnv_p3& rhs )
{
  if ( this != &rhs ) {
    Base_t::operator=( rhs );
  }
  return *this;
}

///////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////

McEventCollectionCnv_p3::~McEventCollectionCnv_p3()
{
}


void McEventCollectionCnv_p3::persToTrans( const McEventCollection_p3* persObj,
                                           McEventCollection* transObj,
                                           MsgStream& msg )
{
  msg << MSG::DEBUG << "Loading McEventCollection from persistent state..."
      << endmsg;

  // elements are managed by DataPool
  transObj->clear(SG::VIEW_ELEMENTS);
  HepMC::DataPool datapools;
  const unsigned int nVertices = persObj->m_genVertices.size();
  if ( datapools.vtx.capacity() - datapools.vtx.allocated() < nVertices ) {
    datapools.vtx.reserve( datapools.vtx.allocated() + nVertices );
  }
  const unsigned int nParts = persObj->m_genParticles.size();
  if ( datapools.part.capacity() - datapools.part.allocated() < nParts ) {
    datapools.part.reserve( datapools.part.allocated() + nParts );
  }
  const unsigned int nEvts = persObj->m_genEvents.size();
  if ( datapools.evt.capacity() - datapools.evt.allocated() < nEvts ) {
    datapools.evt.reserve( datapools.evt.allocated() + nEvts );
  }

  transObj->reserve( nEvts );
  for ( std::vector<GenEvent_p3>::const_iterator
          itr = persObj->m_genEvents.begin(),
          itrEnd = persObj->m_genEvents.end();
        itr != itrEnd;
        ++itr ) {
    const GenEvent_p3& persEvt = *itr;

    HepMC::GenEvent * genEvt        = datapools.getGenEvent();
#ifdef HEPMC3
    genEvt->add_attribute ("barcodes", std::make_shared<HepMC::GenEventBarcodes>());
    genEvt->add_attribute("signal_process_id",std::make_shared<HepMC3::IntAttribute>(persEvt.m_signalProcessId));
    genEvt->set_event_number(persEvt.m_eventNbr);
    genEvt->add_attribute("event_scale",std::make_shared<HepMC3::DoubleAttribute>(persEvt.m_eventScale));
    genEvt->add_attribute("alphaQCD",std::make_shared<HepMC3::DoubleAttribute>(persEvt.m_alphaQCD));
    genEvt->add_attribute("alphaQED",std::make_shared<HepMC3::DoubleAttribute>(persEvt.m_alphaQED));
    genEvt->weights()= persEvt.m_weights;
    genEvt->add_attribute("random_states",std::make_shared<HepMC3::VectorLongIntAttribute>(persEvt.m_randomStates));
    transObj->push_back( genEvt );

    ParticlesMap_t partToEndVtx( (persEvt.m_particlesEnd- persEvt.m_particlesBegin)/2 );

    // create the vertices
    const unsigned int endVtx = persEvt.m_verticesEnd;
    for ( unsigned int iVtx= persEvt.m_verticesBegin; iVtx != endVtx; ++iVtx ) {
      createGenVertex( *persObj, persObj->m_genVertices[iVtx],partToEndVtx, datapools, genEvt );
    } 
#else
    genEvt->m_signal_process_id     = persEvt.m_signalProcessId;
    genEvt->m_event_number          = persEvt.m_eventNbr;
    genEvt->m_event_scale           = persEvt.m_eventScale;
    genEvt->m_alphaQCD              = persEvt.m_alphaQCD;
    genEvt->m_alphaQED              = persEvt.m_alphaQED;
    genEvt->m_signal_process_vertex = 0;
    genEvt->m_weights               = persEvt.m_weights;
    genEvt->m_random_states         = persEvt.m_randomStates;
    genEvt->m_vertex_barcodes.clear();
    genEvt->m_particle_barcodes.clear();
    genEvt->m_pdf_info = 0;         //> not available at that time...

    transObj->push_back( genEvt );

    // create a temporary map associating the barcode of an end-vtx to its
    // particle.
    // As not all particles are stable (d'oh!) we take 50% of the number of
    // particles as an initial size of the hash-map (to prevent re-hash)
    ParticlesMap_t partToEndVtx( (persEvt.m_particlesEnd-
                                  persEvt.m_particlesBegin)/2 );

    // create the vertices
    const unsigned int endVtx = persEvt.m_verticesEnd;
    for ( unsigned int iVtx= persEvt.m_verticesBegin; iVtx != endVtx; ++iVtx ) {
      genEvt->add_vertex( createGenVertex( *persObj,
                                           persObj->m_genVertices[iVtx],
                                           partToEndVtx,
                                           datapools ) );
    } //> end loop over vertices
#endif

    // set the signal process vertex
    const int sigProcVtx = persEvt.m_signalProcessVtx;
    if ( sigProcVtx != 0 ) {
      HepMC::set_signal_process_vertex(genEvt, HepMC::barcode_to_vertex(genEvt, sigProcVtx ) );
    }

    // connect particles to their end vertices
    for ( ParticlesMap_t::iterator
            p = partToEndVtx.begin(),
            endItr = partToEndVtx.end();
          p != endItr;
          ++p ) {
      auto decayVtx=HepMC::barcode_to_vertex(genEvt, p->second );
      if ( decayVtx ) {
        decayVtx->add_particle_in( p->first );
      } else {
        msg << MSG::ERROR
            << "GenParticle points to null end vertex !!"
            << endmsg;
      }
    }
  } //> end loop over m_genEvents

  msg << MSG::DEBUG << "Loaded McEventCollection from persistent state [OK]"
      << endmsg;

  return;
}

void McEventCollectionCnv_p3::transToPers( const McEventCollection* /*transObj*/,
                                           McEventCollection_p3* /*persObj*/,
                                           MsgStream& msg )
{
  msg << MSG::ERROR
      << "This transient-to-persistent converter method has been RETIRED !!"
      << endmsg
      << "You are not supposed to end-up here ! Go away !"
      << endmsg;

  throw std::runtime_error( "Retired McEventCollectionCnv_p3::transToPers() !!" );
  return;
}


HepMC::GenVertexPtr
McEventCollectionCnv_p3::createGenVertex( const McEventCollection_p3& persEvt,
                                          const GenVertex_p3& persVtx,
                                          ParticlesMap_t& partToEndVtx,
                                          HepMC::DataPool& datapools, HepMC::GenEvent* parent ) const
{
  HepMC::GenVertexPtr vtx = datapools.getGenVertex();
  if (parent) parent->add_vertex(vtx);

#ifdef HEPMC3
  vtx->set_position( HepMC::FourVector(persVtx.m_x,persVtx.m_y,persVtx.m_z,persVtx.m_t) );
  vtx->set_status(persVtx.m_id);
  // cast from std::vector<float> to std::vector<double>
  std::vector<double> weights( persVtx.m_weights.begin(), persVtx.m_weights.end() );
  vtx->add_attribute("weights",std::make_shared<HepMC3::VectorDoubleAttribute>(weights));
  HepMC::suggest_barcode(vtx,persVtx.m_barcode);
  // handle the in-going (orphans) particles
  //Is this needed for HEPMC3?
  const unsigned int nPartsIn = persVtx.m_particlesIn.size();
  for ( unsigned int i = 0; i != nPartsIn; ++i ) {
    createGenParticle( persEvt.m_genParticles[persVtx.m_particlesIn[i]], partToEndVtx, datapools, vtx, false );
  }
  // now handle the out-going particles
  const unsigned int nPartsOut = persVtx.m_particlesOut.size();
  for ( unsigned int i = 0; i != nPartsOut; ++i ) {
     createGenParticle( persEvt.m_genParticles[persVtx.m_particlesOut[i]], partToEndVtx, datapools, vtx );
  }
#else
  vtx->m_position.setX( persVtx.m_x );
  vtx->m_position.setY( persVtx.m_y );
  vtx->m_position.setZ( persVtx.m_z );
  vtx->m_position.setT( persVtx.m_t );
  vtx->m_particles_in.clear();
  vtx->m_particles_out.clear();
  vtx->m_id      = persVtx.m_id;
  vtx->m_weights.m_weights.reserve( persVtx.m_weights.size() );
  vtx->m_weights.m_weights.assign ( persVtx.m_weights.begin(),
                                    persVtx.m_weights.end() );
  vtx->m_event   = 0;
  vtx->m_barcode = persVtx.m_barcode;

  // handle the in-going (orphans) particles
  const unsigned int nPartsIn = persVtx.m_particlesIn.size();
  for ( unsigned int i = 0; i != nPartsIn; ++i ) {
    createGenParticle( persEvt.m_genParticles[persVtx.m_particlesIn[i]],
                       partToEndVtx,
                       datapools );
  }

  // now handle the out-going particles
  const unsigned int nPartsOut = persVtx.m_particlesOut.size();
  for ( unsigned int i = 0; i != nPartsOut; ++i ) {
    vtx->add_particle_out( createGenParticle( persEvt.m_genParticles[persVtx.m_particlesOut[i]],
                                              partToEndVtx,
                                              datapools ) );
  }
#endif

  return vtx;
}

HepMC::GenParticlePtr
McEventCollectionCnv_p3::createGenParticle( const GenParticle_p3& persPart,
                                            ParticlesMap_t& partToEndVtx,
                                            HepMC::DataPool& datapools, HepMC::GenVertexPtr parent, bool add_to_output ) const
{
  HepMC::GenParticlePtr p    = datapools.getGenParticle();
  if (parent) add_to_output?parent->add_particle_out(p):parent->add_particle_in(p);

#ifdef HEPMC3
  p->set_pdg_id(persPart.m_pdgId);
  p->set_status(persPart.m_status);
  // Note: do the E calculation in extended (long double) precision.
  // That happens implicitly on x86 with optimization on; saying it
  // explicitly ensures that we get the same results with and without
  // optimization.  (If this is a performance issue for platforms
  // other than x86, one could change to double for those platforms.)
  double temp_e=0.0;
  if ( 0 == persPart.m_recoMethod ) {
    temp_e = std::sqrt( (long double)(persPart.m_px)*persPart.m_px +
                          (long double)(persPart.m_py)*persPart.m_py +
                          (long double)(persPart.m_pz)*persPart.m_pz +
                          (long double)(persPart.m_m) *persPart.m_m );
  } else {
    const int signM2 = ( persPart.m_m >= 0. ? 1 : -1 );
    const double persPart_ene =
      std::sqrt( std::abs((long double)(persPart.m_px)*persPart.m_px +
                (long double)(persPart.m_py)*persPart.m_py +
                (long double)(persPart.m_pz)*persPart.m_pz +
                signM2* (long double)(persPart.m_m)* persPart.m_m));
    const int signEne = ( persPart.m_recoMethod == 1 ? 1 : -1 );
    temp_e=signEne * persPart_ene;
  }
  p->set_momentum(HepMC::FourVector(persPart.m_px,persPart.m_py,persPart.m_pz,temp_e));
  // setup flow
  // fillin' the flow
  std::vector<int> flows;
  const unsigned int nFlow = persPart.m_flow.size();
  for ( unsigned int iFlow= 0; iFlow != nFlow; ++iFlow ) {
  flows.push_back(persPart.m_flow[iFlow].second );
  }
  //We construct it here as vector w/o gaps.
  p->add_attribute("flows", std::make_shared<HepMC3::VectorIntAttribute>(flows));
  HepMC::suggest_barcode(p,persPart.m_barcode);
#else
  p->m_pdg_id              = persPart.m_pdgId;
  p->m_status              = persPart.m_status;
  p->m_polarization.m_theta= static_cast<double>(persPart.m_thetaPolarization);
  p->m_polarization.m_phi  = static_cast<double>(persPart.m_phiPolarization  );
  p->m_production_vertex   = 0;
  p->m_end_vertex          = 0;
  p->m_barcode             = persPart.m_barcode;

  // Note: do the E calculation in extended (long double) precision.
  // That happens implicitly on x86 with optimization on; saying it
  // explicitly ensures that we get the same results with and without
  // optimization.  (If this is a performance issue for platforms
  // other than x86, one could change to double for those platforms.)
  if ( 0 == persPart.m_recoMethod ) {
    p->m_momentum.setPx( persPart.m_px );
    p->m_momentum.setPy( persPart.m_py );
    p->m_momentum.setPz( persPart.m_pz );
    double temp_e = std::sqrt( (long double)(persPart.m_px)*persPart.m_px +
                          (long double)(persPart.m_py)*persPart.m_py +
                          (long double)(persPart.m_pz)*persPart.m_pz +
                          (long double)(persPart.m_m) *persPart.m_m );
    p->m_momentum.setE( temp_e );
  } else {
    const int signM2 = ( persPart.m_m >= 0. ? 1 : -1 );
    const double persPart_ene =
      std::sqrt( std::abs((long double)(persPart.m_px)*persPart.m_px +
                (long double)(persPart.m_py)*persPart.m_py +
                (long double)(persPart.m_pz)*persPart.m_pz +
                signM2* (long double)(persPart.m_m)* persPart.m_m));
    const int signEne = ( persPart.m_recoMethod == 1 ? 1 : -1 );
    p->m_momentum.set( persPart.m_px,
                       persPart.m_py,
                       persPart.m_pz,
                       signEne * persPart_ene );
  }

  // setup flow
  const unsigned int nFlow = persPart.m_flow.size();
  p->m_flow.clear();
  for ( unsigned int iFlow= 0; iFlow != nFlow; ++iFlow ) {
    p->m_flow.set_icode( persPart.m_flow[iFlow].first,
                         persPart.m_flow[iFlow].second );
  }
#endif

  if ( persPart.m_endVtx != 0 ) {
    partToEndVtx[p] = persPart.m_endVtx;
  }

  return p;
}
