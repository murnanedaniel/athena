//*****************************************************************************
//
// Generators/EvtGen_i/EvtInclusiveDecay.h
//
// $Id: EvtInclusiveDecay.cxx,v 1.6 2007-03-01 23:23:44 binet Exp $
//
// EvtInclusiveDecay is a TopAlg that takes HepMC events from StoreGate and
// generates particle decays using EvtGen. Depending on job options either all or
// only a subset of the particles which have decays defined in the EvtGen
// decay files will be handed to EvtGen. Both undecayed particles and particles
// with an existing decay tree can be handled (in the latter case,
// EvtInclusiveDecay will remove the existing decay tree).
//
// Written in February 2006 by Juerg Beringer, based in part on the existing
// EvtDecay module.
//
//*****************************************************************************

#include "EvtGen_i/EvtInclusiveDecay.h"

#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGen_i/EvtGenExternal/EvtExternalGenList.hh"

#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "EvtGen/EvtGen.hh"
#include "EvtGenBase/EvtRandomEngine.hh"
#include "EvtGenBase/EvtDecayTable.hh"

#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/IPartPropSvc.h"
#include "AthenaKernel/IAtRndmGenSvc.h"
#include "StoreGate/StoreGateSvc.h"
#include "StoreGate/DataHandle.h"
#include "GeneratorObjects/McEventCollection.h"
#include "HepPDT/ParticleData.hh"
#include "HepPID/ParticleName.hh"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Vector/LorentzVector.h"

#include <stdlib.h>
#include <sstream>
#include <map>


EvtInclusiveDecay::EvtInclusiveDecay(const std::string& name, ISvcLocator* pSvcLocator):Algorithm( name, pSvcLocator ) {
  m_evtAtRndmGen = 0;
  m_myEvtGen = 0;

  // Basic EvtGen configuration: decay and particle definition files, random number stream
  declareProperty("pdtFile", m_pdtFile = "inclusive.pdt");
  declareProperty("decayFile", m_decayFile = "2014inclusive.dec");
  declareProperty("userDecayFile", m_userDecayFile = "");
  declareProperty("randomStreamName", m_randomStreamName = "EVTGEN");
  declareProperty("inputKeyName", m_inputKeyName = "GEN_EVENT");
  declareProperty("outputKeyName",m_outputKeyName = "GEN_EVENT_EVTGEN");
  declareProperty("readExisting",m_readExisting=false);

  // Selection of particles to be decayed
  declareProperty("prohibitFinalStateDecay", m_prohibitFinalStateDecay=false);
  declareProperty("prohibitReDecay", m_prohibitReDecay=false);
  declareProperty("prohibitUnDecay", m_prohibitUnDecay=true);
  declareProperty("prohibitRemoveSelfDecay", m_prohibitRemoveSelfDecay=false);
  declareProperty("blackList",m_blackList);
  declareProperty("allowAllKnownDecays", m_allowAllKnownDecays=true);
  declareProperty("allowDefaultBDecays", m_allowDefaultBDecays=true);
  declareProperty("whiteList",m_whiteList);

  // Level of output
  declareProperty("printHepMCBeforeEvtGen", m_printHepMCBeforeEvtGen=false);
  declareProperty("printHepMCAfterEvtGen", m_printHepMCAfterEvtGen=false);
  declareProperty("printHepMCHighlighted", m_printHepMCHighlighted=true);
  declareProperty("printHepMCHighLightTopLevelDecays", m_printHepMCHighLightTopLevelDecays=true);

  // Optional checks
  declareProperty("checkDecayTree", m_checkDecayTree=false);
  declareProperty("checkDecayChannels", m_checkDecayChannels=false);
  
  // Repeated decays
  declareProperty("maxNRepeatedDecays", m_maxNRepeatedDecays=1);
  
  // User selection
  declareProperty("applyUserSelection", m_applyUserSelection=false);
  declareProperty("userSelRequireOppositeSignedMu", m_userSelRequireOppositeSignedMu=true);
  declareProperty("userSelMu1MinPt", m_userSelMu1MinPt=0.);
  declareProperty("userSelMu2MinPt", m_userSelMu2MinPt=0.);
  declareProperty("userSelMu1MaxEta", m_userSelMu1MaxEta=102.5);
  declareProperty("userSelMu2MaxEta", m_userSelMu2MaxEta=102.5);
  declareProperty("userSelMinDimuMass", m_userSelMinDimuMass=0.);
  declareProperty("userSelMaxDimuMass", m_userSelMaxDimuMass=-1.); // set to negative to not apply cut

  m_atRndmGenSvc = 0;
  m_mcEvtColl = 0;
  m_sgSvc = 0;

  // We have decided to blacklist Tau decays because we are not sure whether the polarization
  // would be properly passed to EvtGen
  m_blackList.push_back(15);
}



EvtInclusiveDecay::~EvtInclusiveDecay() {
  delete m_myEvtGen;
  delete m_evtAtRndmGen;
}



StatusCode EvtInclusiveDecay::initialize() {
  StatusCode sc;
  MsgStream log(messageService(), name());
  log << MSG::INFO << "EvtInclusiveDecay initialize" << endreq;
  log << MSG::INFO << "Particle properties definition file = " << m_pdtFile << endreq;
  log << MSG::INFO << "Main decay file                     = " << m_decayFile << endreq;
  log << MSG::INFO << "User decay file                     = " << m_userDecayFile << endreq;
  log << MSG::INFO << "Max number of repeated decays       = " << m_maxNRepeatedDecays << endreq;
  log << MSG::INFO << "EvtInclusiveDecay selection parameters:" << endreq;
  log << MSG::INFO << "* prohibitFinalStateDecay = " << m_prohibitFinalStateDecay << endreq;
  log << MSG::INFO << "* prohibitReDecay         = " << m_prohibitReDecay << endreq;
  log << MSG::INFO << "* prohibitUnDecay         = " << m_prohibitUnDecay << endreq;
  log << MSG::INFO << "* prohibitRemoveSelfDecay = " << m_prohibitRemoveSelfDecay << endreq;
  log << MSG::INFO << "* allowAllKnownDecays     = " << m_allowAllKnownDecays << endreq;
  log << MSG::INFO << "* allowDefaultBDecays     = " << m_allowDefaultBDecays << endreq;
  log << MSG::INFO << "User selection parameters:" << endreq;
  log << MSG::INFO << "* applyUserSelection             = " << m_applyUserSelection << endreq;
  log << MSG::INFO << "* userSelRequireOppositeSignedMu = " << m_userSelRequireOppositeSignedMu << endreq;
  log << MSG::INFO << "* userSelMu1MinPt                = " << m_userSelMu1MinPt << endreq;
  log << MSG::INFO << "* userSelMu2MinPt                = " << m_userSelMu2MinPt << endreq;
  log << MSG::INFO << "* userSelMu1MaxEta               = " << m_userSelMu1MaxEta << endreq;
  log << MSG::INFO << "* userSelMu2MaxEta               = " << m_userSelMu2MaxEta << endreq;
  log << MSG::INFO << "* userSelMinDimuMass             = " << m_userSelMinDimuMass << endreq;
  log << MSG::INFO << "* userSelMaxDimuMass             = " << m_userSelMaxDimuMass << endreq;

  // Initialize and print blackList
  m_blackListSet.insert(m_blackList.begin(),m_blackList.end());
  log << MSG::INFO << "* blackList               = ";
  for (std::set<int>::iterator i = m_blackListSet.begin(); i!=m_blackListSet.end(); ++i)
    log << (*i) << " ";
  log << endreq;

  // Initialize and print whiteList
  m_whiteListSet.insert(m_whiteList.begin(),m_whiteList.end());
  log << MSG::INFO << "* whiteList               = ";
  for (std::set<int>::iterator i = m_whiteListSet.begin(); i!=m_whiteListSet.end(); ++i)
    log << (*i) << " ";
  log << endreq;

  // Obtain random number generator for EvtGen
  static const bool CREATEIFNOTTHERE(true);
  sc = service("AtRndmGenSvc", m_atRndmGenSvc, CREATEIFNOTTHERE);
  if (!sc.isSuccess() || 0==m_atRndmGenSvc) {
    log << MSG::ERROR << "Failed to initialize random number service" << endreq;
    return sc;
  }
  m_evtAtRndmGen = new EvtInclusiveAtRndmGen(m_atRndmGenSvc,m_randomStreamName);

  // Create an instance of EvtGen and read particle properties and decay files
  EvtExternalGenList genList(true,xmlpath(),"gamma");
  EvtAbsRadCorr* radCorrEngine = genList.getPhotosModel();
  std::list<EvtDecayBase*> extraModels = genList.getListOfModels();

   // Create the EvtGen generator object
  // EvtGen myGenerator("decayFile.dec", "evt.pdl", randomEnginePointer,
  //                    radCorrEngine, &extraModels);


  m_myEvtGen = new EvtGen( m_decayFile.c_str(), m_pdtFile.c_str(), m_evtAtRndmGen, radCorrEngine, &extraModels);
  if(!m_userDecayFile.empty())
    m_myEvtGen->readUDecay(m_userDecayFile.c_str());

  // Get a handle to StoreGate
  sc = service("StoreGateSvc", m_sgSvc);
  if (sc.isFailure()) {
    log << MSG::ERROR << "Could not find StoreGateSvc" << endreq;
    return sc;
  } 
  
  m_nRepeatedDecays = 0;

  return StatusCode::SUCCESS;
}



StatusCode EvtInclusiveDecay::execute() {
  MsgStream log(messageService(), name());
  log << MSG::DEBUG << "EvtInclusiveDecay executing" << endreq;

  std::string   key = m_inputKeyName;

  // retrieve event from Transient Store (Storegate)
  if(m_readExisting) {
    const McEventCollection* oldmcEvtColl;
    if(m_sgSvc->retrieve(oldmcEvtColl, key).isFailure()) {
      log << MSG::ERROR << "Could not retrieve const McEventCollection" << endreq;
      return StatusCode::FAILURE;
    }
 
    m_mcEvtColl = new McEventCollection;
    std::size_t evtIdx = 0;
    std::size_t pushCount = 0;
    for ( McEventCollection::const_iterator iEvt = oldmcEvtColl->begin();
	iEvt != oldmcEvtColl->end();
	++iEvt, ++evtIdx ) {
        m_mcEvtColl->push_back( new HepMC::GenEvent(**iEvt) );
        ++pushCount;
    } //> end loop over HepMC::GenEvent

    // If we are using the same key for the old and new collections, we have
    // to delete the old one
    if(m_outputKeyName==key) {
      if ( m_sgSvc->remove(oldmcEvtColl).isFailure() ) {
        log<< "Could not remove old HepMC collection" << endreq;
        return StatusCode::FAILURE;
      }
    }
    /*
  if (m_deleteOldCollection) {
    if ( evtStore()->remove(m_mcEvtColl).isFailure() ) {
      ATH_MSG_ERROR ( "Could not remove old HepMC collection" );
      return StatusCode::FAILURE;
    }
  

  if ( evtStore()->record(m_newMcEvtColl, m_newCollectionKey).isFailure() ) {
    ATH_MSG_ERROR ( "Could not add new HepMC collection with key " << m_newCollectionKey );
    return StatusCode::FAILURE;
  }
    */
  }
  else  if ( m_sgSvc->retrieve(m_mcEvtColl, key).isFailure() ) {
    log << MSG::ERROR << "Could not retrieve McEventCollection" << endreq;
    return StatusCode::FAILURE;
  }

  McEventCollection::iterator mcItr;
  for( mcItr = m_mcEvtColl->begin(); mcItr != m_mcEvtColl->end(); mcItr++ )   {
    HepMC::GenEvent* hepMC = *mcItr;

    // Search HepMC record for particles to be decayed by EvtGen
    // NOTE: In order to ensure repeatability, we use a std::set of barcodes to obtain
    //       an ordered list of particles to be decayed by EvtGen.
    std::set<HepMC::GenVertex*> visited;
    std::set<int> toBeDecayed;
    for (HepMC::GenEvent::particle_iterator itp = hepMC->particles_begin(); itp != hepMC->particles_end(); ++itp) {
      HepMC::GenParticle* p = *itp;
      if ( (!p->production_vertex()) ||
	   (p->production_vertex()->particles_in_size() == 0) ) {
	StatusCode sc = traverseDecayTree(p,false,visited,toBeDecayed);
	if (sc.isFailure())
	  return StatusCode::FAILURE;
      }
    }

    // Print HepMC in tree format if desired (before doing anything)
    if (m_printHepMCBeforeEvtGen) {
      log << MSG::DEBUG << "Printing HepMC record at " << hepMC << " BEFORE running EvtGen:" << endreq;
      if (m_printHepMCHighLightTopLevelDecays)
	printHepMC(hepMC,&toBeDecayed);
      else
	printHepMC(hepMC);
    }

    // Decay selected particles
    bool eventPassesCuts(false);
    int loopCounter(0);
    while( !eventPassesCuts && loopCounter < m_maxNRepeatedDecays ) {
      for (std::set<int>::iterator itb = toBeDecayed.begin(); itb!=toBeDecayed.end(); ++itb) {
        HepMC::GenParticle* p = hepMC->barcode_to_particle(*itb);
        if (p==0) {
          log << MSG::ERROR << "Overlapping decay tree encountered for barcode " << *itb << endreq;
          return StatusCode::FAILURE;
        }
        decayParticle(hepMC,p);
      }
      
      if(m_applyUserSelection)
        eventPassesCuts = passesUserSelection(hepMC);
      else
        eventPassesCuts = true;
      
      m_nRepeatedDecays++;
      loopCounter++;
    } 

    // Store the number of decay attempts in event weights std::map, only if repeated decays enabled
    if(m_maxNRepeatedDecays > 1)
      hepMC->weights()["nEvtGenDecayAttempts"] = loopCounter;
    //hepMC->weights().write();

    // Print HepMC in tree format if desired (after finishing all EvtGen decays)
    if (m_printHepMCAfterEvtGen) {
      log << MSG::DEBUG << "Printing HepMC record at " << hepMC << " AFTER running EvtGen:" << endreq;
      if (m_printHepMCHighLightTopLevelDecays)
	printHepMC(hepMC,&toBeDecayed);
      else
	printHepMC(hepMC);
    }
  }
  if(m_readExisting) {
    if ( m_sgSvc->record( m_mcEvtColl,m_outputKeyName).isFailure() ) {
      log << "Could not add new HepMC collection with key " << m_outputKeyName << endreq;
    return StatusCode::FAILURE;
    }
  }
  return StatusCode::SUCCESS;
}



StatusCode EvtInclusiveDecay::finalize() {
  MsgStream log(messageService(), name());
  if (m_checkDecayChannels) {
    log << MSG::INFO << "The following particles were checked and didn't have any decay channels:" << endreq;
    if (log.level() <= MSG::INFO) {
      std::cout << std::endl;
      std::cout << " Particle code    Name from HepPDT        # Occurences" << std::endl;
      std::cout << "------------------------------------------------------"  << std::endl;
      for (std::map<int,long>::iterator p = m_noDecayChannels.begin(); p!=m_noDecayChannels.end(); ++p) {
	int id = p->first;
	int count = p->second;
	std::cout << std::setw(14) << id
		  << std::setw(20) << HepPID::particleName(id)
		  << std::setw(20) << count
		  << std::endl;
      }
      std::cout << std::endl;
    }
  }
  log << MSG::INFO << "Total number of repeated decays: " << m_nRepeatedDecays << endreq;
  log << MSG::INFO << "EvtInclusiveDecay finalized" << endreq;
  return StatusCode::SUCCESS;
}



//
// Recursively traverse the decay tree of a particle p, looking for particles to be
// decayed by EvtGen. Care is taken to traverse each particle only once.
// Note that since in each decay tree only the top-most particle will be decayed
// by EvtGen (with its decay tree being deleted beforehand), we cannot use HepMC's
// "descendant" iterator.
//
StatusCode EvtInclusiveDecay::traverseDecayTree(HepMC::GenParticle* p,
						bool isToBeRemoved,
						std::set<HepMC::GenVertex*>& visited,
						std::set<int>& toBeDecayed) {
  MsgStream log(messageService(), name());
  log << MSG::VERBOSE << "Inspecting: " << pdgName(p) << "   barcode: " << p->barcode() << endreq;
  if (!isToBeRemoved) {
    if (isToBeDecayed(p,true)) {
      toBeDecayed.insert(p->barcode());
      isToBeRemoved = true;
      log << MSG::VERBOSE << "Selected particle for decay: " << pdgName(p) << " (barcode " << p->barcode() << ")" << endreq;
      // In principle we could stop the recursion here. However, to prevent
      // pathological cases in certain decay trees (in particular from Herwig),
      // we continue in order to mark all descendants of this particle
      // as visited. Thus none of these descendants can be flagged for further
      // decay, even if it has several mothers.
    }
  }
  HepMC::GenVertex* v = p->end_vertex();
  if (v) {
    if (visited.insert(v).second) {
      if ( isToBeRemoved && (v->particles_in_size()>1) && m_checkDecayTree ) {
	// This is normal for Herwig but should not occur for Pythia
	log << MSG::WARNING << "Found particle to be decayed with vertex with >1 incoming mother particles in decay tree" << endreq;
	if (log.level() <= MSG::WARNING) {
	  std::cout << std::endl;
	  p->print();
	  v->print();
	  std::cout << std::endl;
	  //return StatusCode::FAILURE;
	}
      }
      for (HepMC::GenVertex::particle_iterator itp = v->particles_begin(HepMC::children);
	                                       itp != v->particles_end(HepMC::children);
                                               ++itp) {
	StatusCode sc = traverseDecayTree(*itp,isToBeRemoved,visited,toBeDecayed);
	if (sc.isFailure())
	  return StatusCode::FAILURE;
      }
    }
  }
  return StatusCode::SUCCESS;
}



//
// Remove an existing decay tree
//
void EvtInclusiveDecay::removeDecayTree(HepMC::GenEvent* hepMC, HepMC::GenParticle* p) {
  HepMC::GenVertex* v = p->end_vertex();
  if (v) {
    MsgStream log(messageService(), name());
    std::set<int> vtxBarCodesToDelete;
    vtxBarCodesToDelete.insert(v->barcode());
    for (HepMC::GenVertex::vertex_iterator itv = v->vertices_begin(HepMC::descendants);
                                           itv != v->vertices_end(HepMC::descendants);
                                           ++itv)
      vtxBarCodesToDelete.insert((*itv)->barcode());
    for (std::set<int>::iterator itb = vtxBarCodesToDelete.begin(); itb != vtxBarCodesToDelete.end(); ++itb) {
      HepMC::GenVertex* vdel = hepMC->barcode_to_vertex(*itb);
      hepMC->remove_vertex(vdel);
      delete vdel;
    }
    p->set_status(1);   // For now, flag particle as undecayed (stable)
    log << MSG::DEBUG << "Removed existing " << pdgName(p) << " (barcode " << p->barcode() << ")" 
	<< " decay tree with " << vtxBarCodesToDelete.size() << " vertices" << endreq;
  }
}



//
// Decay a particle with EvtGen. Any existing decay tree will be removed.
//
// The following status codes are used: 
//
// status == 1     - undecayed particle (also for particles that are not supposed to decay)
// status == 999   - particle decayed by EvtGen
// status == 899   - particle was supposed to be decayed by EvtGen, but found no decay channel
//
// Note that if a particle with an existing decay tree but no defined decay channels
// in EvtGen is passed to this routine, the net effect will be to "undecay" this
// particle. The particle will be flagged by status code 899, but MAY BE LOST further
// downstream in the simulation chain. The default job parameters are for routine
// isToBeDecayed() to never enable decays of such particles by EvtGen.
//
void EvtInclusiveDecay::decayParticle(HepMC::GenEvent* hepMC, HepMC::GenParticle* part) {
  MsgStream log(messageService(), name());
  log << MSG::DEBUG << "Decaying particle " << pdgName(part) << " (barcode " << part->barcode() << ")" << endreq;
  if (log.level() <= MSG::VERBOSE) part->print();

  // Remove existing decay tree, if any, and flag particle as being decayed by EvtGen
  removeDecayTree(hepMC,part);
  part->set_status(899);

  // Create EvtGen version of part and have EvtGen decay it.
  // Since EvtGen uses GeV, convert particles momentum from MeV to GeV.
  int id = part->pdg_id();
  EvtId evtId=EvtPDL::evtIdFromStdHep(id);
  double en =(part->momentum()).e()/1000.;
  double px=(part->momentum()).px()/1000.;
  double py=(part->momentum()).py()/1000.;
  double pz=(part->momentum()).pz()/1000.;
  EvtVector4R evtP(en,px,py,pz);
  EvtParticle* evtPart = EvtParticleFactory::particleFactory(evtId,evtP);
  m_myEvtGen->generateDecay(evtPart);
  if (log.level() <= MSG::VERBOSE) evtPart->printTree();

  // Add new decay tree to hepMC, converting back from GeV to MeV.
  addEvtGenDecayTree(hepMC, part, evtPart, 1000.);
  if(evtPart->getNDaug() !=0) part->set_status(999);  
  evtPart->deleteTree();
}



void EvtInclusiveDecay::addEvtGenDecayTree(HepMC::GenEvent* hepMC, HepMC::GenParticle* part,
				  EvtParticle* evtPart, double momentumScaleFactor) {  
  if(evtPart->getNDaug()!=0) {  

    // Add decay vertex
    double ct=(evtPart->getDaug(0)->get4Pos()).get(0);
    double x=(evtPart->getDaug(0)->get4Pos()).get(1);
    double y=(evtPart->getDaug(0)->get4Pos()).get(2);
    double z=(evtPart->getDaug(0)->get4Pos()).get(3);
    HepMC::GenVertex* end_vtx = new HepMC::GenVertex(CLHEP::HepLorentzVector(x,y,z,ct));
    hepMC->add_vertex(end_vtx);
    end_vtx->add_particle_in(part);

    // Add decay daughter with their own decay trees
    for(uint it=0; it<evtPart->getNDaug(); it++) {
      double e=(evtPart->getDaug(it)->getP4Lab()).get(0) * momentumScaleFactor;
      double px=(evtPart->getDaug(it)->getP4Lab()).get(1) * momentumScaleFactor;
      double py=(evtPart->getDaug(it)->getP4Lab()).get(2) * momentumScaleFactor;
      double pz=(evtPart->getDaug(it)->getP4Lab()).get(3) * momentumScaleFactor;
      int id=EvtPDL::getStdHep(evtPart->getDaug(it)->getId());
      int status=1;
      if(evtPart->getDaug(it)->getNDaug() != 0) status=999;
      HepMC::GenParticle* daughter = new HepMC::GenParticle(CLHEP::HepLorentzVector(px,py,pz,e),id,status);
      end_vtx->add_particle_out(daughter);
      addEvtGenDecayTree(hepMC, daughter, evtPart->getDaug(it), momentumScaleFactor);
    }
  }
}



//
// isToBeDecayed returns true if we want the particle p to be decayed by
// EvtGen based on the job options selected by the user.
// The parameter doCrossChecks is used to prevent double-counting for cross-checks
// if isToBeDecayed is called more than once for the same particle.
//
bool EvtInclusiveDecay::isToBeDecayed(const HepMC::GenParticle* p, bool doCrossChecks) {
  MsgStream log(messageService(), name());
  int id = p->pdg_id();
  int stat = p->status();
  int nDaughters = 0;
  HepMC::GenVertex* v = p->end_vertex();
  if (v) nDaughters = v->particles_out_size();

  // Ignore documentation lines
  if (stat == 3) return false;

  // Particularly for Herwig, try to ignore particles that really should
  // be flagged as documentation lines
  double m2 = p->momentum().m2();
  if (m2 < -1.0E-3) {
    log << MSG::DEBUG << "Ignoring particle " << pdgName(p) << " with m^2 = " << m2 << endreq;
    return false;
  }

  // Check whether EvtGen has any decay channels defined for this particle
  EvtId evtId = EvtPDL::evtIdFromStdHep(id);
  // std::cout << "EVTID: " << evtId.getId() << " alias " << evtId.getAlias() << std::endl;
  int nModes = 0;
  if (evtId.getId()>=0)
    //    nModes = EvtDecayTable::getNMode(evtId.getAlias());
    nModes = EvtDecayTable::getInstance()->getNMode(evtId.getAlias());
  if (doCrossChecks) {
    log << MSG::VERBOSE << "Checking particle " << pdgName(p)
	<< " (status = " << stat
	<<") -- " << nModes << " decay modes found" << endreq;
    if (m_checkDecayChannels && nModes==0) {
      std::map<int,long>::iterator pos = m_noDecayChannels.find(id);
      if (pos != m_noDecayChannels.end())
	(pos->second)++;
      else
	m_noDecayChannels[id] = 1;
    }
  }

  // Check prohibit* settings
  if (m_prohibitFinalStateDecay && stat==1) return false;
  if (m_prohibitReDecay && nDaughters>0) return false;
  if (m_prohibitUnDecay && nModes==0) return false;
  if (m_prohibitRemoveSelfDecay && nDaughters>0) {
    // For now, check only children - this should be sufficient and checking all
    // descendants would be very expensive.
    for (HepMC::GenVertex::particle_iterator itd = v->particles_begin(HepMC::children);
	                                     itd != v->particles_end(HepMC::children);
                                             ++itd) {
      if (abs((*itd)->pdg_id()) == abs(id)) return false;
    }
  }

  // Check blackList
  if (m_blackListSet.count(abs(id))>0) return false;

  // Check allow* settings
  if (m_allowAllKnownDecays && nModes>0) return true;
  if (m_allowDefaultBDecays && isDefaultB(id)) return true;

  // Check whiteList
  if (m_whiteListSet.count(abs(id))>0) return true;

  return false;   // Default is NOT to decay through EvtGen
}



//
// The following mimicks the particle selection implemented in EvtDecay.
//
bool EvtInclusiveDecay::isDefaultB(const int pId) const {
  int id = abs(pId);
  if ( id == 511   || 
       id == 521   ||
       id == 531   ||
       id == 541   || 
       id == 5122  ||      
       id == 5132  ||                         
       id == 5232  || 
       id == 5112  ||        
       id == 5212  ||        
       id == 5222 )
    return true;
  else
    return false;
}

//
// Function to apply the user selection after repeated decay
// Now the selection is based on di-muon kinematics only
// TODO: to be replaced by something more configurable
//
bool EvtInclusiveDecay::passesUserSelection(HepMC::GenEvent* hepMC) {
//   MsgStream log(messageService(), name());
  bool passed(false);
  std::vector<HepMC::GenParticle*> *muons = new std::vector<HepMC::GenParticle*>;
  for (HepMC::GenEvent::particle_iterator itp = hepMC->particles_begin(); itp != hepMC->particles_end(); ++itp) {
    HepMC::GenParticle* p = *itp;
    if( abs(p->pdg_id()) == 13 )
      muons->push_back(p);
  }
//   log << MSG::INFO << "Found " << muons->size() << " muons" << endreq;
  
  for (std::vector<HepMC::GenParticle*>::iterator muItr1 = muons->begin(); muItr1 != muons->end(); ++muItr1) {
    for (std::vector<HepMC::GenParticle*>::iterator muItr2 = muItr1+1; muItr2 != muons->end(); ++muItr2) {
      if( m_userSelRequireOppositeSignedMu && (*muItr1)->pdg_id() * (*muItr2)->pdg_id() > 0)
        continue;
//       log << MSG::INFO << "Opp sign ok" << endreq;
      if( !( (*muItr1)->momentum().perp() > m_userSelMu1MinPt && fabs((*muItr1)->momentum().pseudoRapidity()) < m_userSelMu1MaxEta && 
             (*muItr2)->momentum().perp() > m_userSelMu2MinPt && fabs((*muItr2)->momentum().pseudoRapidity()) < m_userSelMu2MaxEta ) &&
          !( (*muItr2)->momentum().perp() > m_userSelMu1MinPt && fabs((*muItr2)->momentum().pseudoRapidity()) < m_userSelMu1MaxEta && 
             (*muItr1)->momentum().perp() > m_userSelMu2MinPt && fabs((*muItr1)->momentum().pseudoRapidity()) < m_userSelMu2MaxEta ) )
        continue;
//       log << MSG::INFO << "Pt ok" << endreq;
      double dimuMass = invMass((*muItr1),(*muItr2));
      if( !( dimuMass > m_userSelMinDimuMass && (dimuMass < m_userSelMaxDimuMass || m_userSelMaxDimuMass < 0.) ) )
        continue;
//       log << MSG::INFO << "Mass ok" << endreq;
      passed = true;
    }
  }
  
  delete muons;
  
  return passed;
}

double EvtInclusiveDecay::invMass(const HepMC::GenParticle* p1, const HepMC::GenParticle* p2) {
  double p1Px = p1->momentum().px();
  double p1Py = p1->momentum().py();
  double p1Pz = p1->momentum().pz();
  double p1E  = p1->momentum().e();
  double p2Px = p2->momentum().px();
  double p2Py = p2->momentum().py();
  double p2Pz = p2->momentum().pz();
  double p2E  = p2->momentum().e();
  double dimuE = p2E + p1E;
  double dimuPx = p2Px + p1Px;
  double dimuPy = p2Py + p1Py;
  double dimuPz = p2Pz + p1Pz;
  double invMass = sqrt(dimuE*dimuE - dimuPx*dimuPx - dimuPy*dimuPy - dimuPz*dimuPz);
  
  return invMass;
}

//
// Utility functions to print a HepMC event record in a tree-like format, using
// colors to denote the status of particles and to indicate which particles
// are selected by the job options to be decayed by EvtGen.
//
void EvtInclusiveDecay::printHepMC(HepMC::GenEvent* hepMC, std::set<int>* barcodeList) {
  std::set<HepMC::GenVertex*> visited;
  unsigned int nParticlesFound = 0;
  unsigned int nTreesFound = 0;
  for (HepMC::GenEvent::particle_iterator itp = hepMC->particles_begin(); itp != hepMC->particles_end(); ++itp) {
    HepMC::GenParticle* p = *itp;
    if ( (!p->production_vertex()) ||
	 (p->production_vertex()->particles_in_size() == 0) ) {
      nTreesFound++;
      std::cout << "\n    Found new partial decay tree:\n" << std::endl;
      unsigned int nParticlesVisited = printTree(p,visited,1,barcodeList);
      std::cout << "\n    " << nParticlesVisited << " particles in this subtree" << std::endl;
      nParticlesFound += nParticlesVisited;
    }
  }
  std::cout << "\n    Total of " << nParticlesFound << " particles found in "
	    << nTreesFound << " decay subtrees in HepMC event record\n" << std::endl;
}

unsigned int EvtInclusiveDecay::printTree(HepMC::GenParticle* p,
				 std::set<HepMC::GenVertex*>& visited, int level, std::set<int>* barcodeList) {
  unsigned int nParticlesVisited = 1;
  for (int i=0; i<level; i++) std::cout << "    ";
  std::cout << pdgName(p,m_printHepMCHighlighted,barcodeList);
  HepMC::GenVertex* v = p->end_vertex();
  if (v) {
    if (v->particles_in_size() > 1)
      std::cout << " [interaction: " << v->particles_in_size() << " particles, barcode " << v->barcode() << "]    -->   ";
    else
      std::cout << "   -->   ";
    if (visited.insert(v).second) {
      for (HepMC::GenVertex::particle_iterator itp = v->particles_begin(HepMC::children);
	                                       itp != v->particles_end(HepMC::children);
	                                       ++itp) {
	std::cout << pdgName(*itp,m_printHepMCHighlighted,barcodeList) << "   ";
      }
      std::cout << std::endl;
      for (HepMC::GenVertex::particle_iterator itp = v->particles_begin(HepMC::children);
	                                       itp != v->particles_end(HepMC::children);
                                               ++itp) {
	if ((*itp)->end_vertex())
	  nParticlesVisited += printTree(*itp, visited, level+1, barcodeList);
	else
	  nParticlesVisited++;
      }
    } else
      std:: cout << "see above" << std::endl;
  } else
    std::cout << "   no decay vertex\n" << std::endl;
  return nParticlesVisited;
}

std::string EvtInclusiveDecay::pdgName(const HepMC::GenParticle* p, bool statusHighlighting, std::set<int>* barcodeList) {
  // Note: HepPDT doesn't seem to know anti-particle names
  // const HepPDT::ParticleData* pData = m_pdt->particle(HepPDT::ParticleID(abs(p->pdg_id())));
  std::ostringstream buf;
  if (statusHighlighting) {
    if ( ((barcodeList!=0) && (barcodeList->find(p->barcode()) != barcodeList->end())) ||
         ((barcodeList==0) && isToBeDecayed(p,false)) )
      buf << "\033[7m";   // reverse
    if (p->status() != 1) {
      if (p->status() == 999)
	buf << "\033[33m";   // yellow
      else
	buf << "\033[31m";   // red
    }
  }
  buf << p->pdg_id();
  //if (pData)
  //  buf << "/" << pData->name();
  buf << "/" << HepPID::particleName(p->pdg_id());
  if (statusHighlighting) {
    buf << "\033[0m";   // revert color attributes
  }
  return buf.str();
}



//
// Interface between Athena random number service and EvtGen's EvtRandomEngine class
//
EvtInclusiveAtRndmGen::EvtInclusiveAtRndmGen(IAtRndmGenSvc* atRndmGenSvc, std::string streamName)
  : m_atRndmGenSvc(atRndmGenSvc)
  , m_streamName(streamName)
{}

EvtInclusiveAtRndmGen::~EvtInclusiveAtRndmGen() {
}

double EvtInclusiveAtRndmGen::random() {
  CLHEP::HepRandomEngine* engine = m_atRndmGenSvc->GetEngine(m_streamName);
  return CLHEP::RandFlat::shoot(engine);
}

////////////////////////////////////////////////////////////////////////////////
std::string EvtInclusiveDecay::xmlpath(){

  char *cmtpath = getenv("CMTPATH");
  char *cmtconfig = getenv("CMTCONFIG");

  std::string foundpath = "";

  if(cmtpath != 0 && cmtconfig != 0){

    std::vector<std::string> cmtpaths;
    boost::split(cmtpaths, cmtpath, boost::is_any_of(std::string(":")));

    std::string installPath = "/InstallArea/" + std::string(cmtconfig) + "/share/Pythia8/xmldoc";

    for(std::vector<std::string>::const_iterator path = cmtpaths.begin();
        path != cmtpaths.end() && foundpath == ""; ++path){
      std::string testPath = *path + installPath;
      std::ifstream testFile(testPath.c_str());
      if(testFile.good()) foundpath = testPath;
      testFile.close();
    }

  }

  return foundpath;
}

