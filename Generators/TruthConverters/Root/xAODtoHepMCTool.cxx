/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/


#include "TruthConverters/xAODtoHepMCTool.h"
#ifndef XAOD_STANDALONE
#include "AthAnalysisBaseComps/AthAnalysisHelper.h"
#endif

xAODtoHepMCTool::xAODtoHepMCTool( const std::string& name ) 
  : asg::AsgTool( name ),
    m_momFac(1.), //<-- MeV/GeV conversion factor
    m_lenFac(1.), //<-- length conversion factor
    m_signalOnly(true),
    m_maxCount(0)
{
}


StatusCode xAODtoHepMCTool::initialize() {
  ATH_MSG_INFO ("Initializing xAODtoHepMCTool "<< name() << "...");
  ATH_MSG_INFO ("SignalOnly         = " << m_signalOnly);
  m_evtCount   = 0;
  m_badSuggest = 0;
  m_noProdVtx  = 0;
  m_badBeams   = 0;
  return StatusCode::SUCCESS;
}



StatusCode xAODtoHepMCTool :: finalize ()
{
  ATH_MSG_INFO ("==============================================================");
  ATH_MSG_INFO ("==========    xAOD -> HepMC Tool :: Run Summary     ==========");
  ATH_MSG_INFO ("==============================================================");
  if( m_badSuggest ){
    ATH_MSG_INFO ("Number of suggest_barcode failures = " << m_badSuggest);
    ATH_MSG_INFO (" ... check input particle/vertex barcodes.");
  } else {
    ATH_MSG_INFO ("No suggest_barcode failures");
  }
  if ( m_noProdVtx ) {
    ATH_MSG_INFO ("Number of events that have a missing production vertex = " << m_noProdVtx);
  }else{
    ATH_MSG_INFO ("No missing production vertices");
  }
  if ( m_badBeams ){
    ATH_MSG_INFO ("Number events with improperly defined beamparticles = " << m_badBeams);
    ATH_MSG_INFO ("Inconsistencies with the beam particles storage in the event record!");
    ATH_MSG_INFO ("... check xAOD::TruthEvent record ...");
  }else{
    ATH_MSG_INFO ("No events with undefined beams.");
  }
  ATH_MSG_INFO ("==============================================================");
  ATH_MSG_INFO ("===================    End Run Summary     ===================");
  ATH_MSG_INFO ("==============================================================");
  return StatusCode::SUCCESS;
}

std::vector<HepMC::GenEvent> xAODtoHepMCTool :: getHepMCEvents(const xAOD::TruthEventContainer* xTruthEventContainer, 
                                                               const unsigned long long evtNum) {
  ++m_evtCount;
  bool doPrint = m_evtCount < m_maxCount;
  // CREATE HEPMC EVENT COLLECTION
  std::vector<HepMC::GenEvent> mcEventCollection;
  // Loop over xAOD truth container
  // Signal event is always first, followed by pileup events
  ATH_MSG_DEBUG("Start event loop");
  for (const xAOD::TruthEvent* xAODEvent : *xTruthEventContainer) {
    if( doPrint ) ATH_MSG_DEBUG("XXX Printing xAOD Event");
    if( doPrint ) printxAODEvent( xAODEvent,evtNum );
    // Create GenEvent for each xAOD truth event
    ATH_MSG_DEBUG("Create new GenEvent");
    const HepMC::GenEvent& hepmcEvent = createHepMCEvent(xAODEvent,evtNum);
    // Insert into McEventCollection
    mcEventCollection.push_back(hepmcEvent);    
    if( doPrint ) ATH_MSG_DEBUG("XXX Printing HepMC Event");
    if( doPrint ) hepmcEvent.print();
    // Quit if signal only
    if( m_signalOnly ) break;
  }  
  return mcEventCollection;
}

std::vector<HepMC::GenEvent> xAODtoHepMCTool :: getHepMCEvents(const xAOD::TruthEventContainer* xTruthEventContainer, 
                                                               const xAOD::EventInfo* eventInfo) {
  return getHepMCEvents(xTruthEventContainer, eventInfo->eventNumber());
}


const HepMC::GenEvent xAODtoHepMCTool::createHepMCEvent(const xAOD::TruthEvent* xEvt, const unsigned long long evtNum) {
  
  ///EVENT LEVEL
  HepMC::GenEvent genEvt;

  genEvt.set_event_number(evtNum);
  ATH_MSG_DEBUG("Start createHepMCEvent for event " <<evtNum);
  
  //Beam particles
  std::pair<const xAOD::TruthParticle*,const xAOD::TruthParticle*> beamParticles = xEvt->beamParticles();
  

  xAOD::TruthEvent::PdfInfo pdfInfo = xEvt->pdfInfo();
  HepMC::PdfInfo pdfInfo_hepmc;
  pdfInfo_hepmc.set_id1( pdfInfo.pdgId1 );
  pdfInfo_hepmc.set_id2( pdfInfo.pdgId2 );
  pdfInfo_hepmc.set_pdf_id1( pdfInfo.pdfId1 );
  pdfInfo_hepmc.set_pdf_id2( pdfInfo.pdfId2 );
  pdfInfo_hepmc.set_x1( (double)pdfInfo.x1 );
  pdfInfo_hepmc.set_x2( (double)pdfInfo.x2 );
  pdfInfo_hepmc.set_pdf1( (double)pdfInfo.xf1 );
  pdfInfo_hepmc.set_pdf2( (double)pdfInfo.xf2 );
  //pdfInfo_hepmc.set_scalePDF( m_momFac * (double)pdfInfo.Q );
  genEvt.set_pdf_info( pdfInfo_hepmc );

  //Cross section <- Can't see this read in Reader
  // Need to look in xTruthEvent
  //HepMC::GenCrossSection* xsec = new HepMC::GenCrossSection();
  //xsec->set_cross_section((double)xEvt->crossSection());
  //xsec->set_cross_section_error((double)xEvt->crossSectionError());
  //genEvt.set_cross_section( *xsec ); // this is currently always -1 ... something wrong
  
  //Weights
  #ifndef XAOD_STANDALONE
  const std::vector<float> weights = xEvt->weights();
  std::map<std::string, int> weightNameMap;
  if (AthAnalysisHelper::retrieveMetadata("/Generation/Parameters","HepMCWeightNames", weightNameMap).isFailure()) {
    ATH_MSG_DEBUG("Couldn't find meta-data for weight names.");
  }
  if (weightNameMap.size()) {
    HepMC::WeightContainer& wc = genEvt.weights();
    wc.clear();
    auto itEnd = weightNameMap.end();
    for (int idx = 0; idx < int(weights.size()); ++idx) {
      for (auto it = weightNameMap.begin(); it != itEnd; ++it) {
	if (it->second == idx) {
	  wc[ it->first ] = weights[idx];
	  break;
	}
      }
    }
  }
  else {
    for ( std::vector<float>::const_iterator wgt = weights.begin(); wgt != weights.end(); ++wgt ) {
      genEvt.weights().push_back(*wgt);
    }
  }
  #endif

  // PARTICLES AND VERTICES  
  // Map of existing vertices - needed for the tree linking
  std::map<const xAOD::TruthVertex*,HepMC::GenVertex*> vertexMap;
  std::pair<HepMC::GenParticle*,HepMC::GenParticle*> HepMC_beamParticles;

  // Loop over all of the particles in the event, call particle builder
  // Call suggest_barcode only after insertion!
  //int i=0;
  for (auto tlink : xEvt->truthParticleLinks()) {
    if (!tlink.isValid()) continue;
    const xAOD::TruthParticle* xPart = *tlink;
    //xAOD::TruthParticle* xPart = *tlink;
    // sanity check
    if (xPart == nullptr) {
      ATH_MSG_WARNING("xAOD TruthParticle is equal to NULL. This should not happen!");
      continue;
    }

    if( !xPart->hasProdVtx() && !xPart->hasDecayVtx() ){
      ATH_MSG_WARNING("xAOD particle with no vertices, bc = "<<xPart->barcode());
      continue;
    }

    // skip particles with barcode > 200000 --> Geant4 secondaries 
    if ( xPart->barcode() > 200000 ) continue; 

    // Create GenParticle
    HepMC::GenParticle* hepmcParticle( createHepMCParticle(xPart) );
    int bcpart = xPart->barcode();
    //i+=1;
    //std::cout<< "i: "<< i <<std::endl;
    //std::cout << "xPart pdgId: " << xPart->pdgId() << std::endl;
    //std::cout << "xPart status 1: " << xPart->status() << std::endl;
    //std::cout << "xPart energy: " << xPart->e() << std::endl;
    //std::cout<< "barcode: " << bcpart <<std::endl;
    //std::cout<< "Energy" << xPart-><<std::endl;
    //std::cout<< "beamParticles.first: "<< beamParticles.first << std::endl;
    //std::cout<< "beamParticles.second: "<< beamParticles.second << std::endl;
    if( xPart == beamParticles.first  )  HepMC_beamParticles.first  = hepmcParticle;
    if( xPart == beamParticles.second )  HepMC_beamParticles.second = hepmcParticle;
    //std::cout << "hepmcParticle status: "<< hepmcParticle->status() << std::endl;
    //if(hepmcParticle->status()!=4 && hepmcParticle->pdg_id()!=2212){
    //  std::cout<< "INCORRECT BEAM PARTICLE "<< std::endl;
    //}
    //if(xPart->status()!=4 && (i==1 || i==2)){
    //  xAOD::TruthParticle* xPart_initial;
    //  xPart_initial->setPdgId(4);
      //xPart->setPdgId(4);
    //}
    /*
    if (xPart->status()==21 && (i==1 || i==2) ){
      j+=1;
      std::cout<<"j: "<< j << std::endl;
      std::cout <<"WRONG INCOMING BEAM PARTICLE!! " << std::endl;
      continue;
      }*/
    //if (xPart->status() == 4){ 
    //  xAOD::TruthParticle* xPart_four;
    //  xPart_four->setStatus(2);
    //  std::cout << "xPart status 2: " << xPart_four->status() << std::endl;
    //}

    // status 10902 should be treated just as status 2
    if ( hepmcParticle->status() == 10902 ) hepmcParticle->set_status(2);

    // Get the production and decay vertices
    if( xPart->hasProdVtx() ) {
      const xAOD::TruthVertex* xAODProdVtx = xPart->prodVtx();
      // skip production vertices with barcode > 200000 --> Geant4 secondaries 
      if ( std::abs(xAODProdVtx->barcode()) > 200000 ) continue; 
      bool prodVtxSeenBefore(false); // is this new?
      HepMC::GenVertex* hepmcProdVtx = vertexHelper(xAODProdVtx,vertexMap,prodVtxSeenBefore);
      // Set the decay/production links
      hepmcProdVtx->add_particle_out(hepmcParticle);
      // Insert into Event
      if (!prodVtxSeenBefore){ 
        genEvt.add_vertex(hepmcProdVtx);
        if( !hepmcProdVtx->suggest_barcode(xAODProdVtx->barcode()) ){
          ATH_MSG_WARNING("suggest_barcode failed for vertex "<<xAODProdVtx->barcode());
          ++m_badSuggest;
        }
      }
      if( !hepmcParticle->suggest_barcode(bcpart) ){
        ATH_MSG_DEBUG("suggest_barcode failed for particle " <<bcpart);
        ++m_badSuggest;
      }
      bcpart = 0;
    } else {
      ATH_MSG_DEBUG("No production vertex found for particle "<<xPart->barcode());
    }
    
    if( xPart->hasDecayVtx() ){
      const xAOD::TruthVertex* xAODDecayVtx = xPart->decayVtx();
      // skip decay vertices with barcode > 200000 --> Geant4 secondaries 
      if ( fabs(xAODDecayVtx->barcode()) > 200000 ) continue; 
      bool decayVtxSeenBefore(false); // is this new?
      HepMC::GenVertex* hepmcDecayVtx = vertexHelper(xAODDecayVtx,vertexMap,decayVtxSeenBefore);
      // Set the decay/production links
      hepmcDecayVtx->add_particle_in(hepmcParticle);
      // Insert into Event
      if (!decayVtxSeenBefore){ 
        genEvt.add_vertex(hepmcDecayVtx);
        if( !hepmcDecayVtx->suggest_barcode(xAODDecayVtx->barcode()) ){
          ATH_MSG_WARNING("suggest_barcode failed for vertex " << xAODDecayVtx->barcode());
          ++m_badSuggest;
        }
      }
      if( bcpart != 0 ){
        if( !hepmcParticle->suggest_barcode(bcpart) ){
          ATH_MSG_DEBUG("suggest_barcode failed for particle " <<bcpart);
          ++m_badSuggest;
        }
        bcpart = 0;
      }
    }

  } // end of particle loop

  genEvt.set_beam_particles(HepMC_beamParticles);
  //std::cout<< "HepMC_beamParticles first: "<<HepMC_beamParticles.first<<std::endl;
  //std::cout<< "HepMC_beamParticles second: "<<HepMC_beamParticles.second<<std::endl;
  const HepMC::GenEvent constGenEvt(genEvt);
  ATH_MSG_DEBUG("Returning const GenEvent");
  return constGenEvt;
  
}

// Helper to check whether a vertex exists or not using a map; 
// calls createHepMCVertex if not
HepMC::GenVertex* xAODtoHepMCTool::vertexHelper(const xAOD::TruthVertex* xaodVertex,
						std::map<const xAOD::TruthVertex*,HepMC::GenVertex*> &vertexMap,
						bool &seenBefore) {
  
  HepMC::GenVertex* hepmcVertex;
  std::map<const xAOD::TruthVertex*,HepMC::GenVertex*>::iterator vMapItr;
  vMapItr=vertexMap.find(xaodVertex);
  // Vertex seen before?
  if (vMapItr!=vertexMap.end()) {
    // YES: use the HepMC::Vertex already in the map
    hepmcVertex = (*vMapItr).second;
    seenBefore = true;
  } else {
    // NO: create a new HepMC::Vertex
    hepmcVertex = createHepMCVertex(xaodVertex);
    vertexMap[xaodVertex] = hepmcVertex;
    seenBefore = false;
  }
  return hepmcVertex;
  
}

// Create the HepMC GenParticle
// Call suggest_barcode after insertion!
HepMC::GenParticle* xAODtoHepMCTool::createHepMCParticle(const xAOD::TruthParticle* particle) {
  ATH_MSG_VERBOSE("Creating GenParticle for barcode " <<particle->barcode());
  const HepMC::FourVector fourVec( m_momFac * particle->px(), m_momFac * particle->py(), m_momFac * particle->pz(), m_momFac * particle->e() );
  HepMC::GenParticle* hepmcParticle = new HepMC::GenParticle(fourVec, particle->pdgId(), particle->status());
  //std::cout<< "particle pdgID: " << particle->pdgId() << std::endl;
  //std::cout<< "particle mass: " << particle->m() << std::endl;
  //std::cout<< "particle status: " << particle->status() <<"particle mass: "<< particle->m()<< std::endl;

  hepmcParticle->set_generated_mass( m_momFac * particle->m());
  return hepmcParticle;
}

// Create the HepMC GenVertex
// Call suggest_barcode after insertion!
HepMC::GenVertex* xAODtoHepMCTool::createHepMCVertex(const xAOD::TruthVertex* vertex) {
  ATH_MSG_VERBOSE("Creating GenVertex for barcode " <<vertex->barcode());
  HepMC::FourVector prod_pos( m_lenFac * vertex->x(), m_lenFac * vertex->y(),m_lenFac * vertex->z(), m_lenFac * vertex->t() );
  HepMC::GenVertex* genVertex = new HepMC::GenVertex(prod_pos);
  return genVertex;
}

// Print xAODTruth Event. The printout is particle oriented, unlike the
// HepMC particle/vertex printout. Geant and pileup particles with
// barcode>100000 are omitted.
void xAODtoHepMCTool::printxAODEvent(const xAOD::TruthEvent* event, const unsigned long long evtNum) const{
  
  std::vector<int> bcPars;
  std::vector<int> bcKids;

  //long long int evtNum = eventInfo->eventNumber();

  std::cout <<"======================================================================================" <<std::endl;
  std::cout <<"xAODTruth Event " << evtNum <<std::endl;
  std::cout <<"   Barcode      PDG Id  Status   px(GeV)   py(GeV)   pz(GeV)    E(GeV)   Parent: Decay" <<std::endl;
  std::cout <<"   -----------------------------------------------------------------------------------" <<std::endl;

  int nPart = event->nTruthParticles();
  for(int i=0; i<nPart; ++i){
    const xAOD::TruthParticle* part = event->truthParticle(i);
    if (part==nullptr) continue;
    int bc = part->barcode();
    if( bc > 100000 ) continue;
    int id = part->pdgId();
    if ( id != 25 ) continue;
    int stat = part->status();
    float px = part->px()/1000.;
    float py = part->py()/1000.;
    float pz = part->pz()/1000.;
    float e = part->e()/1000.;
    bcPars.clear();
    bcKids.clear();

    if( part->hasProdVtx() ){
      const xAOD::TruthVertex* pvtx = part->prodVtx();
      if( pvtx ) bcPars.push_back(pvtx->barcode());
    }

    if( part->hasDecayVtx() ){
      const xAOD::TruthVertex* dvtx = part->decayVtx();
      if( dvtx ) bcKids.push_back(dvtx->barcode());
    }

    std::cout <<std::setw(10)<<bc <<std::setw(12)<<id
              <<std::setw(8)<<stat
              <<std::setprecision(2)<<std::fixed
              <<std::setw(10)<<px <<std::setw(10)<<py
              <<std::setw(10)<<pz <<std::setw(10)<<e <<"   ";
    std::cout <<"P: ";
    for(unsigned int k=0; k<bcPars.size(); ++k){
      std::cout <<bcPars[k] <<" ";
    }
    std::cout <<"  D: ";
    for(unsigned int k=0; k<bcKids.size(); ++k){
      std::cout <<bcKids[k] <<" ";
    }
    std::cout <<std::endl;
  }
  std::cout <<"======================================================================================" <<std::endl;
}

