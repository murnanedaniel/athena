/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
// SimKernel.cxx, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

// ISF_Algs includes
#include "SimKernel.h"
// ISF_Event includes
#include "ISF_Event/ISFParticle.h"
// ISF_Interfaces includes
#include "ISF_Interfaces/IParticleBroker.h"
#include "ISF_Interfaces/ITruthSvc.h"
#include "ISF_Interfaces/ISimHitSvc.h"
#include "ISF_Interfaces/ISimulationSvc.h"
#include "ISF_Interfaces/IMonitoringTool.h"
#include "ISF_Interfaces/IEventFilterTool.h"
// FrameWork includes
#include "GaudiKernel/Property.h"
// Boost
#include <boost/lexical_cast.hpp>
// ATLAS cxx utils
#include "CxxUtils/make_unique.h"
#include "PmbCxxUtils/CustomBenchmark.h"
// ROOT includes
#include "TTree.h"
// DetectorDescription
#include "AtlasDetDescr/AtlasRegionHelper.h"
// McEventCollection
#include "GeneratorObjects/McEventCollection.h"

///////////////////////////////////////////////////////////////////
// Public methods:
///////////////////////////////////////////////////////////////////

// Constructors
////////////////
ISF::SimKernel::SimKernel( const std::string& name, ISvcLocator* pSvcLocator ) :
  ::AthAlgorithm( name, pSvcLocator ),
  m_inputHardScatterEvgen(),
  m_inputPileupEvgen(),
  m_outputHardScatterTruth(),
  m_outputPileupTruth(),
  m_inputConverter("",name),
  m_validationOutput(false),
  m_thistSvc("THistSvc",name),
  m_validationStream("ISFSimKernel"),
  m_t_simParticles(0),
  m_val_x(0.0),
  m_val_y(0.0),
  m_val_z(0.0),
  m_val_p(0.0),
  m_val_px(0.0),
  m_val_py(0.0),
  m_val_pz(0.0),
  m_val_meta(0.0),
  m_val_peta(0.0),
  m_val_pdg(0),
  m_val_simID(0),
  m_val_geoID(0),
  m_val_sc(0),
  m_particleBroker("ISF_ParticleBroker", name),
  m_truthRecordSvc("ISF_TruthRecordSvc", name),
  m_simHitSvc("ISF_SimHitSvc", name),
  m_doMemMon(true),
  m_memMon("MemMonitoringTool"),
  m_memUsageEvts(1000),
  m_simSvcs(ISF::fMaxNumAtlasSimIDs),
  m_simSvcNames(ISF::fMaxNumAtlasSimIDs),
  m_numSimSvcs(ISF::fFirstAtlasSimID), // ==1 since UndefinedSimID is always there
  m_numISFEvents(0),
  m_screenOutputPrefix("isf >> "),
  m_screenEmptyPrefix(""),
  m_doCPUMon(true),
  //m_benchPDGCode(0), TODO: implement this if feasible
  //m_benchGeoID(0), TODO: implement this if feasible
  m_benchSimID(0),
  m_numParticles(0),
  m_maxParticleVectorSize(10240)
{
    declareProperty("InputHardScatterCollection",
                    m_inputHardScatterEvgen,
                    "Input Hard Scatter EVGEN collection.");
    declareProperty("InputPileupCollection",
                     m_inputPileupEvgen,
                    "Input Pileup EVGEN collection.");
    declareProperty("OutputHardScatterTruthCollection",
                    m_outputHardScatterTruth,
                    "Output Hard Scatter Truth collection.");
    declareProperty("OutputPileupTruthCollection",
                     m_outputPileupTruth,
                    "Output Pileup Truth collection.");
    declareProperty("InputConverter",
                    m_inputConverter,
                    "Input McEventCollection->ISFParticleContainer conversion service.");

    // validation output section
    declareProperty( "ValidationOutput",
                     m_validationOutput = false,
                     "If turned on, write out a ROOT tree.");
    declareProperty("ValidationStreamName",
                     m_validationStream = "ISFSimKernel",
                     "Name of the output stream" );
    declareProperty("THistService",
                     m_thistSvc,
                     "The THistSvc" );

    // the general services and tools needed
    declareProperty("ParticleBroker"             , m_particleBroker                  );
    declareProperty("TruthRecordService"         , m_truthRecordSvc                  );
    declareProperty("SimHitService"              , m_simHitSvc                       );
    declareProperty("DoCPUMonitoring"            , m_doCPUMon                        );
    declareProperty("DoMemoryMonitoring"         , m_doMemMon                        );
    declareProperty("MemoryMonitoringTool"       , m_memMon                          );
    declareProperty("SummarizeMemUsageEveryNEvts", m_memUsageEvts                    );
    // routing tool
    declareProperty("BeamPipeSimulationSelectors", m_simSelectors[AtlasDetDescr::fAtlasForward]  );
    declareProperty("IDSimulationSelectors"      , m_simSelectors[AtlasDetDescr::fAtlasID]       );
    declareProperty("CaloSimulationSelectors"    , m_simSelectors[AtlasDetDescr::fAtlasCalo]     );
    declareProperty("MSSimulationSelectors"      , m_simSelectors[AtlasDetDescr::fAtlasMS]       );
    declareProperty("CavernSimulationSelectors"  , m_simSelectors[AtlasDetDescr::fAtlasCavern]   );
    // event filter
    declareProperty("EventFilterTools"           , m_eventFilters                      );
    // tuning parameters
    declareProperty("MaximumParticleVectorSize"  , m_maxParticleVectorSize             );
    // refine the screen output for debugging
    declareProperty("ScreenOutputPrefix"      , m_screenOutputPrefix  );
}

// Destructor
///////////////
ISF::SimKernel::~SimKernel()
{}

// Athena Algorithm's Hooks
////////////////////////////
StatusCode ISF::SimKernel::initialize()
{

  // Screen output, part 1
  for (size_t prl = 0; prl < m_screenOutputPrefix.size(); ++prl) m_screenEmptyPrefix += " ";
  ATH_MSG_VERBOSE ( m_screenOutputPrefix << "--------------------------------------------------------");
  ATH_MSG_INFO( m_screenOutputPrefix << "Initializing the ISF KERNEL " );

  // setup memory monitoring Tool
  if ( m_doMemMon) {
    // memory monitoring tool given -> do memory monitoring
    if ( m_memMon.retrieve().isFailure() ){
        ATH_MSG_FATAL( m_screenOutputPrefix <<  "Could not retrieve MemoryMonitoring Service. Abort.");
        return StatusCode::FAILURE;
    } else
        ATH_MSG_INFO( m_screenEmptyPrefix <<  "- MemoryMonitoring  : " << m_memMon.typeAndName() );
    // record current memory usage
    m_memMon->recordCurrent("at beginning of SimKernel initialize()");
  }

  // setup for validation mode
  if ( m_validationOutput) {

    // retrieve the histogram service
    if ( m_thistSvc.retrieve().isSuccess() ) {
      // Create the prefix of histogram names for the THistSvc
      std::string prefix = "/" + m_validationStream + "/";
      std::string treeName="ISF Simulated Particles";
      m_t_simParticles = new TTree( "particles", treeName.c_str() );
      m_t_simParticles->Branch("x"         , &m_val_x    , "x/F"         );
      m_t_simParticles->Branch("y"         , &m_val_y    , "y/F"         );
      m_t_simParticles->Branch("z"         , &m_val_z    , "z/F"         );
      m_t_simParticles->Branch("p"         , &m_val_p    , "p/F"         );
      m_t_simParticles->Branch("pdg"       , &m_val_pdg  , "pdg/I"       );
      m_t_simParticles->Branch("statuscode", &m_val_sc   , "statuscode/S");
      m_t_simParticles->Branch("px"        , &m_val_px   , "px/F"        );
      m_t_simParticles->Branch("py"        , &m_val_py   , "py/F"        );
      m_t_simParticles->Branch("pz"        , &m_val_pz   , "pz/F"        );
      m_t_simParticles->Branch("meta"      , &m_val_meta , "meta/F"     );
      m_t_simParticles->Branch("peta"      , &m_val_peta , "peta/F"     );
      m_t_simParticles->Branch("simID"     , &m_val_simID, "simID/I"   );
      //m_t_simParticles->Branch("geoID"     , &m_val_geoID, "geoID/I"   );

      // register the Tree to the THistSvc and return it's StatusCode
      ATH_CHECK(m_thistSvc->regTree( prefix+treeName, m_t_simParticles) );
      ATH_MSG_INFO( m_screenOutputPrefix <<"Validation mode creating ISF Simulated Particles tree");

    }

    // error when trying to retrieve the THistSvc
    else {
      // -> turn off validation output
      ATH_MSG_ERROR( m_screenOutputPrefix << "Validation mode turned on but unable to retrieve THistService. Will not write out ROOT histograms/Trees.");
      m_validationOutput = false;
    }

  } // end if (m_validationOutput)
  else {
    //ATH_MSG_ERROR( m_screenOutputPrefix <<"m_validationOutput didn't work");
  }

  // setup CPU Benchmarks
  if (m_doCPUMon) {
    //if (!m_benchPDGCode)
    //  m_benchPDGCode = new PMonUtils::CustomBenchmark(ISF::fMaxBenchmarkPDGCode);
    //if (!m_benchGeoID)
    //  m_benchGeoID   = new PMonUtils::CustomBenchmark(AtlasDetDescr::fNumAtlasRegions     );
    if (!m_benchSimID)
      m_benchSimID   = new PMonUtils::CustomBenchmark(ISF::fMaxNumAtlasSimIDs  );
  }

  // retrieve the stack service
  if ( m_particleBroker.retrieve().isFailure() ){
      ATH_MSG_FATAL( m_screenOutputPrefix <<  "Could not retrieve ParticleBroker Service. Abort.");
      return StatusCode::FAILURE;
  } else
      ATH_MSG_INFO( m_screenEmptyPrefix <<  "- ParticleBroker   : " << m_particleBroker.typeAndName() );

  // the truth service
  if ( m_truthRecordSvc.retrieve().isFailure() ){
      ATH_MSG_FATAL( m_screenOutputPrefix <<  "Could not retrieve ParticleStack Service. Abort.");
      return StatusCode::FAILURE;
  } else
      ATH_MSG_INFO( m_screenEmptyPrefix <<  "- TruthRecordSvc   : " << m_truthRecordSvc.typeAndName() );

  // and the simhit service
  if ( m_simHitSvc.retrieve().isFailure() ){
      ATH_MSG_FATAL( m_screenOutputPrefix <<  "Could not retrieve SimHit Service. Abort.");
      return StatusCode::FAILURE;
  } else
      ATH_MSG_INFO( m_screenEmptyPrefix <<  "- SimHitSvc   : " << m_simHitSvc.typeAndName() );

  // initialize all SimulationServices
  //
  for ( short geoID=AtlasDetDescr::fFirstAtlasRegion; geoID<AtlasDetDescr::fNumAtlasRegions ; ++geoID) {
    if ( initSimSvcs(m_simSelectors[geoID]).isFailure())
    return StatusCode::FAILURE;
  }

  // initialize all the EventFilterTools
  if( m_eventFilters.retrieve().isFailure() ) {
      ATH_MSG_FATAL( m_screenOutputPrefix <<  "Failed to retrieve Event Filters. Abort.");
      return StatusCode::FAILURE;
  } else {
    ATH_MSG_INFO( m_screenOutputPrefix << "The following Event Filters are defined:");
    ATH_MSG_INFO( m_screenOutputPrefix << m_eventFilters);
  }

  // free unused space
  m_simSvcs.resize( m_numSimSvcs);
  m_simSvcNames.resize( m_numSimSvcs);
  // some screen output
  ATH_MSG_INFO ( m_screenOutputPrefix << "The following SimulationSvc are registered to ISF:");
  for (SimSvcID id=ISF::fFirstAtlasSimID; id<m_numSimSvcs; id++)
    ATH_MSG_INFO ( m_screenOutputPrefix << " ID: " << id << "\t  Name: '" << m_simSvcNames[id]
                   << "'");

  // setup the simulation selectors
  //
  for ( short geoID=AtlasDetDescr::fFirstAtlasRegion; geoID<AtlasDetDescr::fNumAtlasRegions ; ++geoID) {
    if ( m_particleBroker->registerSimSelector( m_simSelectors[geoID], (AtlasDetDescr::AtlasRegion)geoID).isFailure()) {
      ATH_MSG_ERROR( m_screenOutputPrefix << "Unable to register SimulationSelectors for GeoID="
                     << AtlasDetDescr::AtlasRegionHelper::getName(geoID));
      return StatusCode::FAILURE;
    }
  }
  // screen output
  ATH_MSG_INFO( m_screenOutputPrefix << "The following routing chains are defined:");
  for ( short geoID = 0; geoID<AtlasDetDescr::fNumAtlasRegions ; ++geoID) {
    ATH_MSG_INFO( m_screenOutputPrefix << AtlasDetDescr::AtlasRegionHelper::getName(geoID)
                  << " (GeoID=" << geoID << "): \t" << m_simSelectors[geoID]);
  }

  // record current memory usage
  if (m_doMemMon) m_memMon->recordCurrent("at end of ISF SimKernel initialize()");

  // intialziation successful
  return StatusCode::SUCCESS;
}


StatusCode ISF::SimKernel::finalize()
{
  ATH_MSG_INFO ( m_screenOutputPrefix << "Finalizing ...");

  // record current memory usage
  if (m_doMemMon) m_memMon->recordCurrent("at beginning of ISF SimKernel finalize()");

  // statistics: number of particles handled
  ATH_MSG_INFO(" Number of particles handled by the ISF SimKernel: " << m_numParticles );

  ATH_MSG_INFO(" ========================= ISF Timing Stats  =========================");
  // Benchmarking: by SimID
  if (m_benchSimID) {
    ATH_MSG_INFO("Breakdown of simulation loop by SimulatorID:");
    //TODO: for (unsigned simID=0;simID<m_benchSimID->size();++simID) {
    for (unsigned simID=0;simID<ISF::fMaxNumAtlasSimIDs;++simID) {
      uint64_t count;
      double time_ms;
      m_benchSimID->getData(simID, count, time_ms);
      if (count>0)
       ATH_MSG_INFO("  "<<std::setprecision(4)
                        <<m_simSvcNames[simID]<<" (id="<<simID<<")"
                        <<"\t\tn="<<count<<"\t\tt=" <<time_ms<<" ms\t\tt/n="<<time_ms/count<<" ms"
                        <<std::setprecision(-1) );
    }

    delete m_benchSimID;
    m_benchSimID=0;
  }

  //TODO: implement this if feasible
  // Benchmarking: by GeoID
  //if (m_benchGeoID) {
  //  ATH_MSG_INFO("Breakdown of simulation loop by GeoID:");
  //  //TODO: for (unsigned geoID=0;geoID<m_benchGeoID->size();++simID) {
  //  for (unsigned geoID=AtlasDetDescr::fFirstAtlasRegion;geoID<AtlasDetDescr::fNumAtlasRegions;++geoID) {
  //    uint64_t count;
  //    double time_ms;
  //    m_benchGeoID->getData(geoID, count, time_ms);
  //    if (count>0)
  //     ATH_MSG_INFO("  "<<std::setprecision(4)
  //                      <<AtlasDetDescr::AtlasRegionHelper::getName(geoID)<<" (id="<<geoID<<")"
  //                      <<"\t\tn="<<count<<"\t\tt=" <<time_ms<<" ms\t\tt/n="<<time_ms/count<<" ms"
  //                      <<std::setprecision(-1) );
  //  }

  //  delete m_benchGeoID;
  //  m_benchGeoID=0;
  //}

  //TODO: implement this if feasible
  // Benchmarking: by PDGCode
  //if (m_benchPDGCode) {
  //  ATH_MSG_INFO("Breakdown of simulation loop by PDG Particle Code:");
  //  //TODO: for (unsigned geoID=0;geoID<m_benchGeoID->size();++simID) {
  //  for (int pdgCode=ISF::fUndefinedPDGCode;pdgCode<ISF::fMaxBenchmarkPDGCode;++pdgCode) {
  //    uint64_t count;
  //    double time_ms;
  //    m_benchPDGCode->getData(pdgCode, count, time_ms);
  //    if (count>0)
  //     ATH_MSG_INFO( std::setprecision(4)
  //                   << "  |PDGCode|="<<pdgCode<<", n="
  //                   <<count<<", t="<<time_ms<<" ms, t/n="<<time_ms/count<<" ms"
  //                   <<std::setprecision(-1) );
  //  }

  //  delete m_benchPDGCode;
  //  m_benchPDGCode=0;
  //}

  // call the memory monitoring tool to print some memory stats
  if (m_doMemMon) {
    ATH_MSG_INFO(" ====================== ISF Memory Usage Stats =======================");
    m_memMon->dumpSummary("end of ISF event");
  }

  ATH_MSG_INFO(" =====================================================================");

  return StatusCode::SUCCESS;
}


StatusCode ISF::SimKernel::initSimSvcs( SimSelectorToolArray &simSelectorTools)
{
  // (1.) retrieve all SimulationSelector tools in the array
  if ( simSelectorTools.retrieve().isFailure() ) {
      ATH_MSG_FATAL( m_screenOutputPrefix <<  "Could not retrieve SimulatorSelector Tool Array. Abort.");
      return StatusCode::FAILURE;
  }

  // (2.) loop over SimulationSelector tool array and retrieve simulators
  SimSelectorToolArray::iterator fSimSelectorIter    = simSelectorTools.begin();
  SimSelectorToolArray::iterator fSimSelectorIterEnd = simSelectorTools.end();
  for ( ; fSimSelectorIter != fSimSelectorIterEnd; ++fSimSelectorIter ) {

    // take the simulator from the current SimulationSelector
    ServiceHandle<ISimulationSvc> *curSimulator = (*fSimSelectorIter)->simulator();

    if ( (*curSimulator).retrieve().isFailure() ){
        ATH_MSG_FATAL( m_screenOutputPrefix << "Could not retrieve SimulatorSelector Tool. Abort.");
        return StatusCode::FAILURE;
    } else
        ATH_MSG_INFO( m_screenEmptyPrefix <<  "- SimulationSelector   : " << fSimSelectorIter->typeAndName() );

    // hand over particle broker to simulator
    if ( (*curSimulator)->setParticleBroker( &*m_particleBroker).isFailure() ){
      ATH_MSG_FATAL( m_screenOutputPrefix << "Unable to register ParticleService to SimulationService "
                                          << *curSimulator );
      return StatusCode::FAILURE;
    }

    // get the unique ID assigned to the Simulation Service
    SimSvcID curID = (*curSimulator)->simSvcID();
    // if no ID assigned yet -> new simulator
    if ( curID == ISF::fUndefinedSimID) {
      // assign a new new ID to the simulator
      (*curSimulator)->assignSimSvcID( m_numSimSvcs);
      // register current simulator to the simulatorArray
      m_simSvcs[m_numSimSvcs]     = (&**curSimulator);
      m_simSvcNames[m_numSimSvcs] = (*curSimulator)->simSvcDescriptor();
      ATH_MSG_DEBUG( m_screenOutputPrefix << "Assigned SimSvcID=" << m_numSimSvcs
                    << " to simulator '" << m_simSvcNames[m_numSimSvcs]
                    << "'");
      // increment the total number of simulators registered (=used as IDs)
      ++m_numSimSvcs;
    }

  } // loop over simulation Selectors


  return StatusCode::SUCCESS;
}


StatusCode ISF::SimKernel::execute()
{

  ATH_MSG_DEBUG ( m_screenOutputPrefix << "Executing ...");

  // dump and record current memory stats
  if ( m_doMemMon && (m_numISFEvents==0) ) {
    m_memMon->dumpCurrent( "before 1st event", false );
    m_memMon->recordCurrent("before 1st event");
  }

  // read and convert input
  //  a. hard-scatter
  ISFParticleContainer simParticles{}; // particles for ISF simulation
  ATH_CHECK( prepareInput(m_inputHardScatterEvgen, m_outputHardScatterTruth, simParticles) );
  //  b. pileup
  if (!m_inputPileupEvgen.key().empty()) {
    bool isPileup = true;
    ATH_CHECK( prepareInput(m_inputPileupEvgen, m_outputPileupTruth, simParticles, isPileup) );
  }

  // -----------------------------------------------------------------------------------------------
  // Step 1: Initialize the particle stack and the TruthManager, ABORT if failure
  StatusCode sc = m_particleBroker->initializeEvent( std::move(simParticles) );
  if ( sc.isFailure() ){
      ATH_MSG_FATAL( m_screenOutputPrefix <<  "Failed to initialize Particle Broker. Abort.");
      return StatusCode::FAILURE;
  }
  if ( (m_truthRecordSvc->initializeTruthCollection()).isFailure() ){
      ATH_MSG_FATAL( m_screenOutputPrefix <<  "Failed to initialize TruthService. Abort.");
      return StatusCode::FAILURE;
  }
  // -----------------------------------------------------------------------------------------------


  // -----------------------------------------------------------------------------------------------
  // Step 2: Initialize the Event
  {
    if (m_simHitSvc->initializeEvent().isFailure() ){
      ATH_MSG_ERROR( m_screenOutputPrefix << "Event initialize failed for "<< m_simHitSvc );
      return StatusCode::FAILURE;
    } else
      ATH_MSG_DEBUG( m_screenOutputPrefix << "Event initialize done for "<< m_simHitSvc );

    std::vector<ISimulationSvc*>::iterator fSimSvcIter     = m_simSvcs.begin();
    std::vector<ISimulationSvc*>::iterator fSimSvcIterEnd  = m_simSvcs.end();
    for ( ; fSimSvcIter != fSimSvcIterEnd; ++fSimSvcIter ){
      ISimulationSvc *curSimSvc = (*fSimSvcIter);
      // if simulation with current flavour is registered
      //  -> setupEvent
      if ( curSimSvc){
        if( curSimSvc->setupEvent().isFailure() ) {
          ATH_MSG_WARNING( m_screenOutputPrefix <<  "Event setup failed for "
                           << curSimSvc->simSvcDescriptor() );
        } else {
          ATH_MSG_DEBUG  ( m_screenOutputPrefix <<  "Event setup done for "
                           << curSimSvc->simSvcDescriptor() );
        }
      }
    }
  }
  // -----------------------------------------------------------------------------------------------



  // -----------------------------------------------------------------------------------------------
  // Step 3: ISimulation KERNEL : loop over particle stack, until empty
  ATH_MSG_DEBUG( m_screenOutputPrefix << "Starting simulation loop, initial particle stack size: " << m_particleBroker->numParticles());
  while ( true ) {
    // get next vector of particles for simulation
    const ISF::ConstISFParticleVector &particles = m_particleBroker->popVector(m_maxParticleVectorSize);
    int numParticles = particles.size();
    // particle vector empty -> end simulation
    if (numParticles==0) break;

    // for job statistics
    m_numParticles += numParticles;

    // retrieve the particle destination simulator (and geoID)
    const ISFParticle *firstP = particles.front();
    ISF::SimSvcID simID = firstP->nextSimID();
    //AtlasDetDescr::AtlasRegion geoID = firstP->nextGeoID();

    ATH_MSG_DEBUG  ( m_screenOutputPrefix <<  "Took " << numParticles << " particles from queue (remaining: " << m_particleBroker->numParticles() << ")" );
    ATH_MSG_VERBOSE( m_screenOutputPrefix <<  " -> All particles will be sent to '" << m_simSvcNames[simID] << "' simulator (SimSvcID=" << simID << ")"  );

    // block defines scope for Benchmarks
    //  -> benchmarks will be stared/stopped automatically via the CustomBenchmarkGuard
    //     constructor and destructor, respectively
    StatusCode simSC;
    {
      // setup sim svc benchmarks
      //PMonUtils::CustomBenhmarkGuard benchPDG  ( m_benchPDGCode, pdgCode );
      //PMonUtils::CustomBenchmarkGuard benchGeoID( m_benchGeoID  , geoID , numParticles );
      PMonUtils::CustomBenchmarkGuard benchSimID( m_benchSimID  , simID , numParticles );

      // ===> simulate particle
      simSC = m_simSvcs[simID]->simulateVector( particles);
      if ( simSC.isFailure())
        ATH_MSG_WARNING( m_screenOutputPrefix << "Simulation of particles failed in Simulator: " << m_simSvcNames[simID]);
    }

    // validation output
    if (m_validationOutput) {
      // loop over all particles in the 'particles' vector
      ISF::ConstISFParticleVector::const_iterator partIt    = particles.begin();
      ISF::ConstISFParticleVector::const_iterator partItEnd = particles.end();
      for ( ; partIt != partItEnd; partIt++) {
        const ISFParticle *curPart = (*partIt);
        m_val_x     = curPart->position().x();
        m_val_y     = curPart->position().y();
        m_val_z     = curPart->position().z();
        m_val_p     = curPart->momentum().mag2();
        m_val_px    = curPart->momentum().x();
        m_val_py    = curPart->momentum().y();
        m_val_pz    = curPart->momentum().z();
        m_val_meta  = curPart->momentum().eta();
        m_val_peta  = curPart->position().eta();
        m_val_pdg   = curPart->pdgCode();
        m_val_sc    = simSC;
        m_val_simID = simID;
        //m_val_geoID = geoID;

        m_t_simParticles->Fill();
      }
    }
  }
  ATH_MSG_VERBOSE( m_screenOutputPrefix << "Ending simulation loop, no more particles in the stack");
  // -----------------------------------------------------------------------------------------------



  // Step 4: Finalize the Event
  //  -> stack service
  if ( m_particleBroker->finalizeEvent().isFailure()) {
    ATH_MSG_WARNING( m_screenOutputPrefix <<  "ParticleBroker returned with an error in event finalization.");
  }
  //  -> simulator services
  {
    std::vector<ISimulationSvc*>::iterator fSimSvcIter     = m_simSvcs.begin();
    std::vector<ISimulationSvc*>::iterator fSimSvcIterEnd  = m_simSvcs.end();
    for ( ; fSimSvcIter != fSimSvcIterEnd; ++fSimSvcIter ){
      ISimulationSvc *curSimSvc = (*fSimSvcIter);
      // if simulation with current flavour is registered
      //  -> releaseEvent()
      if ( curSimSvc){
        if( curSimSvc->releaseEvent().isFailure() ) {
          ATH_MSG_WARNING( m_screenOutputPrefix <<  "Event setup failed for "
                           << curSimSvc->simSvcDescriptor() );
        } else {
          ATH_MSG_DEBUG  ( m_screenOutputPrefix <<  "Event setup done for "
                           << curSimSvc->simSvcDescriptor() );
        }
      }
    } // -> loop over SimSvcs
  }

  if ( m_truthRecordSvc->releaseEvent().isFailure() ){
      ATH_MSG_FATAL( m_screenOutputPrefix <<  "Event finalize failed for TruthService. Abort.");
      return StatusCode::FAILURE;
  }
  if ( m_simHitSvc->releaseEvent().isFailure() ){
      ATH_MSG_ERROR( m_screenOutputPrefix << "Event finalize failed for "<< m_simHitSvc );
      return StatusCode::FAILURE;
  } else
      ATH_MSG_DEBUG( m_screenOutputPrefix << "Event finalize done for "<< m_simHitSvc );
  // -----------------------------------------------------------------------------------------------

  // Step 5: Check Any Filters
  ToolHandleArray<IEventFilterTool>::iterator eventFilter(m_eventFilters.begin());
  const ToolHandleArray<IEventFilterTool>::iterator endOfEventFilters(m_eventFilters.end());
  while (eventFilter != endOfEventFilters) {
    if (!((**eventFilter).eventPassesFilter())) {
      setFilterPassed(false);
      ATH_MSG_INFO("This event failed the " << (**eventFilter).name() << " Filter. Therefore it will not be recorded.");
      break;
    }
    ++eventFilter;
  }


  // -----------------------------------------------------------------------------------------------

  // dump current memory monitoring information
  if (m_doMemMon) {
    std::string evtStr = boost::lexical_cast<std::string>( m_numISFEvents );
    std::string descr("after event " + evtStr);
    m_memMon->dumpCurrent( descr.c_str(), true);

    // ISF internal event counting
    m_numISFEvents++;

    // memory monitoring records for the final summary
    if ( !(m_numISFEvents%m_memUsageEvts) ) m_memMon->recordCurrent( descr.c_str() );
    else if ( m_numISFEvents==1)   m_memMon->recordCurrent("after   1st event");
    else if ( m_numISFEvents==2)   m_memMon->recordCurrent("after   2nd event");
    else if ( m_numISFEvents==10)  m_memMon->recordCurrent("after  10th event");
    else if ( m_numISFEvents==100) m_memMon->recordCurrent("after 100th event");
  }

  return StatusCode::SUCCESS;
}


/** Convert input generator particles to ISFParticles and copy input
    generator truth collection into output simulation truth collection */
StatusCode ISF::SimKernel::prepareInput(SG::ReadHandle<McEventCollection>& inputTruth,
                                        SG::WriteHandle<McEventCollection>& outputTruth,
                                        ISFParticleContainer& simParticles,
                                        bool isPileup) const {

  if (!inputTruth.isValid()) {
    ATH_MSG_FATAL("Unable to read input GenEvent collection '" << inputTruth.key() << "'");
    return StatusCode::FAILURE;
  }

  // create copy
  outputTruth = CxxUtils::make_unique<McEventCollection>(*inputTruth);

  ATH_CHECK( m_inputConverter->convert(*outputTruth, simParticles, isPileup) );

  return StatusCode::SUCCESS;
}
