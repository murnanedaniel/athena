/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
// EntryLayerTool.cxx, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

// class header include
#include "ISF_Geant4CommonTools/EntryLayerTool.h"

// ISF includes
#include "ISF_Event/ISFParticle.h"
#include "ISF_Interfaces/IGeoIDSvc.h"

// ISF interfaces
#include "ISF_Interfaces/IParticleFilter.h"

/** Constructor **/
ISF::EntryLayerTool::EntryLayerTool(const std::string& t, const std::string& n, const IInterface* p) : 
  AthAlgTool(t,n,p),
  m_incidentSvc("IncidentSvc",n),
  m_geoIDSvc("GeoIDSvc",n),
  m_geoIDSvcQuick(0),
  m_particleFilterHandle(),
  m_particleFilter(),
  m_numParticleFilters(0),
  m_collection(),
  m_SGName(),
  m_volumeName()
{
  declareInterface<ISF::IEntryLayerTool>(this);

  // geometry identification service
  declareProperty( "GeoIDSvc",
                   m_geoIDSvc,
                   "AthenaService used to indentify sub-detector by (x,y,z) coordintes.");

   // particle filters
  declareProperty( "ParticleFilters",
                   m_particleFilterHandle,
                   "ISF Particle filters, defining whether a particle will be stored or not.");

  // storegate collection names
  declareProperty( "CaloEntryLayer",
                   m_SGName[ISF::fAtlasCaloEntry] = "CaloEntryLayer",
                   "Name of CaloEntryLayer collection on Storegate");
  declareProperty( "MuonEntryLayer",
                   m_SGName[ISF::fAtlasMuonEntry] = "MuonEntryLayer",
                   "Name of MuonEntryLayer collection on Storegate");
  declareProperty( "MuonExitLayer",
                   m_SGName[ISF::fAtlasMuonExit]  = "MuonExitLayer",
                   "Name of MuonExitLayer collection on Storegate");

  // volumeNames for TrackRecords
  declareProperty( "CaloEntryVolumeString",
                   m_volumeName[ISF::fAtlasCaloEntry] = "IDET::IDET",
                   "VolumeName in TrackRecords in CaloEntryLayer");
  declareProperty( "MuonEntryVolumeString",
                   m_volumeName[ISF::fAtlasMuonEntry] = "CALO::CALO",
                   "VolumeName in TrackRecords in MuonEntryLayer");
  declareProperty( "MuonExitVolumeString",
                   m_volumeName[ISF::fAtlasMuonExit]  = "MUONQ02::MUONQ02",
                   "VolumeName in TrackRecords in MuonExitLayer");
}


/** Destructor **/
ISF::EntryLayerTool::~EntryLayerTool()
{
}


/** Athena algtool's Hooks */
StatusCode  ISF::EntryLayerTool::initialize()
{
  ATH_MSG_INFO("initialize() ...");

  // incident service from Athena/Gaudi framework
  if (m_incidentSvc.retrieve().isFailure()) {
    ATH_MSG_FATAL("Could not retrieve IncidentService '" << m_incidentSvc << "'. Exiting.");
    return StatusCode::FAILURE;
  }
  // register to the incident service: BeginEvent for TrackRecordCollection
  m_incidentSvc->addListener( this, IncidentType::BeginEvent);


  // retrieve the GeoIDSvc
  if ( m_geoIDSvc.retrieve().isFailure() ){
    ATH_MSG_FATAL( "Could not retrieve GeoID service. Abort.");
    return StatusCode::FAILURE;
  }
  // store a quick-access pointer to the c++ class directly
  m_geoIDSvcQuick = &(*m_geoIDSvc);


  // retrieve the ParticleFilters
  if ( m_particleFilterHandle.retrieve().isFailure() ){
    ATH_MSG_FATAL( "Could not retrieve ParticleFilters. Abort.");
    return StatusCode::FAILURE;
  }
  // store a quick-access pointer to the c++ classes directly
  m_numParticleFilters = m_particleFilterHandle.size();
  m_particleFilter = new ISF::IParticleFilter*[ m_numParticleFilters ];
  for ( size_t curFilter = 0; curFilter<m_numParticleFilters; curFilter++) {
      // convert ToolHandle to standard c++ class pointer
      m_particleFilter[curFilter] = &(*m_particleFilterHandle[curFilter]);
  }

  ATH_MSG_INFO("initialize() successful");
  return StatusCode::SUCCESS;
}


/** Athena algtool's Hooks */
StatusCode  ISF::EntryLayerTool::finalize()
{
  ATH_MSG_INFO("finalize() successful");
  return StatusCode::SUCCESS;
}


/** framework handle */
void ISF::EntryLayerTool::handle( const Incident& inc ) {

  // check the incident type
  if ( inc.type() == IncidentType::BeginEvent ) {
    // initialize storegate collections at BeginEvent incident
    for ( int entryLayer=ISF::fFirstAtlasEntryLayer;
          entryLayer<ISF::fNumAtlasEntryLayers;
          entryLayer++) {
      m_collection[entryLayer] = setupSGCollection( m_SGName[entryLayer] );
    }
  }

  return;
}


/** Add the given particle to the corresponding Entry/Exit layer if applicable */
ISF::EntryLayer ISF::EntryLayerTool::registerParticle(const ISF::ISFParticle& particle, ISF::EntryLayer layerHit)
{
  // (1.) check whether the particle actually passes all the filters
  //      -> rather fast usually
  {
    bool pass = true;
    for ( size_t curFilter=0; pass && (curFilter<m_numParticleFilters); curFilter++) {
      // check current filter
      pass = m_particleFilter[curFilter]->passFilter( particle);
    }
    // return if not passed
    if (!pass)
      return ISF::fUnsetEntryLayer;
  }

  // at this point all filters are passed

  if (layerHit == ISF::fUnsetEntryLayer) {

    // (2.) check whether the particle lies on any entry surface
    //      -> this goes second because computation intensive routines
    //         are used for this

    const Amg::Vector3D &pos = particle.position();

    // check if particle on ID and/or Calo surface
    bool onIDSurface   = (m_geoIDSvcQuick->inside( pos, AtlasDetDescr::fAtlasID)   == ISF::fSurface);
    bool onCaloSurface = (m_geoIDSvcQuick->inside( pos, AtlasDetDescr::fAtlasCalo) == ISF::fSurface);

    // on CaloEntry layer ?
    if ( onIDSurface && onCaloSurface ) {
      layerHit = ISF::fAtlasCaloEntry;
    }

    // no surface hit yet -> test MS volume surface hit
    else {
      // check if particle on MS surface
      bool onMSSurface   = (m_geoIDSvcQuick->inside( pos, AtlasDetDescr::fAtlasMS) == ISF::fSurface);

      // on MuonEntry layer ?
      if (onCaloSurface && onMSSurface) {
        layerHit = ISF::fAtlasMuonEntry;
      }
      // on MuonExit layer ?
      else if (onMSSurface) {
        layerHit = ISF::fAtlasMuonExit;
      }
    }
  }

  // particle is on a boundary surface -> add it to TrackRecordCollection
  if ( layerHit != ISF::fUnsetEntryLayer) {
    ATH_MSG_VERBOSE( "Particle >>" << particle << "<< hit boundary surface, "
                     "adding it to '" << m_SGName[layerHit] << "' TrackRecord collection");
    TrackRecord tRec = convert(particle, m_volumeName[layerHit]);
    m_collection[layerHit]->push_back(tRec);
  }

  return layerHit;
}


/** used to setup a TrackRecordCollection on storegate */
TrackRecordCollection *ISF::EntryLayerTool::setupSGCollection(const std::string &collectionName)
{
  TrackRecordCollection *collection = 0;

  // check if storegate already contains the collection
  // (a)     if yes ... try to retrieve it
  if ( evtStore()->contains<TrackRecordCollection>(collectionName) ){
      if ( (evtStore()->retrieve( collection, collectionName)).isFailure() )
          ATH_MSG_ERROR( "[ --- ] Unable to retrieve TrackRecordCollection " << collectionName);
  // (b)     if no ... try to create it     
  } else {
      collection = new TrackRecordCollection( collectionName);
      if ( (evtStore()->record( collection, collectionName, true)).isFailure() ) {
          ATH_MSG_ERROR( "[ --- ] Unable to record SiHitCollection " << collectionName);
          delete collection;
          collection=0;
      }
  }

  return collection;
}


/** convert ISFParticle to TrackRecord */
TrackRecord& ISF::EntryLayerTool::convert( const ISF::ISFParticle &p,
                                               std::string &volumeName )
{
  // definitions for quick access:
  double mass            = p.mass();
  CLHEP::Hep3Vector      mom( p.momentum().x(), p.momentum().y(), p.momentum().z() );// not optimal, but required by TrackRecord
  CLHEP::Hep3Vector      pos( p.position().x(), p.position().y(), p.position().z() );// not optimal, but required by TrackRecord
  double energy          = sqrt(mass*mass + mom.mag2());
  double time            = 0.; //!< @TODO: fix time

  // finally create the TrackRecord
  TrackRecord tRec( p.pdgCode(),
                    energy,
                    mom,
                    pos,
                    time,
                    p.barcode(),
                    volumeName );
  return tRec;
}

