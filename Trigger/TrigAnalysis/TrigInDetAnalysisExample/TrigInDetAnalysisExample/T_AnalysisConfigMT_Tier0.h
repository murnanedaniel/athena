/* emacs: this is -*- c++ -*- */
/**
 **     @file    T_AnalysisConfigMT_Tier0.h
 **
 **     @brief   baseclass template so that we can use in different contexts 
 **              in different ways in the monitoring
 ** 
 **       NB: this will be configured to run *either* a standard 
 **       analysis, or a "purity" analysis. If a purity analysis, 
 **       the trigger tracks become the reference (with all the 
 **       selection) and the offline or truth the "test" tracks
 **       This would be a simple switch if the reference tracks
 **       were in the RoI, but as they are not we need to move the 
 **       RoI filtering to the test filter and *not* the reference 
 **       filter grrrrrr  
 **
 **     @author  mark sutton
 **     @date    Tue 16 May 2017 09:28:55 CEST 
 **
 **     Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
 **/

#ifndef TrigInDetAnalysisExample_T_AnalysisConfigMT_Tier0_H
#define TrigInDetAnalysisExample_T_AnalysisConfigMT_Tier0_H

#include "TrigInDetAnalysis/TIDAEvent.h"
#include "TrigInDetAnalysis/TIDAVertex.h"
#include "TrigInDetAnalysis/TrackSelector.h"     
#include "TrigInDetAnalysisUtils/T_AnalysisConfig.h"
#include "TrigInDetAnalysisUtils/TagNProbe.h"

#include "TrigInDetAnalysisExample/Analysis_Tier0.h"
#include "TrigInDetAnalysisExample/AnalysisR3_Tier0.h"
#include "TrigInDetAnalysisExample/VtxAnalysis.h"
#include "TrigInDetAnalysisExample/ChainString.h"
#include "TrigInDetAnalysisExample/TIDATools.h"

#include "TTree.h"
#include "TFile.h"

#include "GaudiKernel/ToolHandle.h"
#include "AthenaMonitoringKernel/GenericMonitoringTool.h"
 

// McParticleEvent includes
#include "McParticleEvent/TruthParticleContainer.h"

#include "GeneratorObjects/McEventCollection.h"
#include "AtlasHepMC/GenEvent.h"
#include "AtlasHepMC/GenVertex.h"
#include "AtlasHepMC/GenParticle.h"

#include "xAODEventInfo/EventInfo.h"



#include "TrigInDetAnalysis/TIDDirectory.h"

#include "TrigInDetAnalysis/Filter_AcceptAll.h"

#include "TrigInDetAnalysisUtils/TIDARoiDescriptorBuilder.h"
#include "TrigInDetAnalysisUtils/Filter_etaPT.h"
#include "TrigInDetAnalysisUtils/Filter_RoiSelector.h"
#include "TrigInDetAnalysisUtils/Associator_BestMatch.h"
#include "TrigInDetAnalysisUtils/Filters.h"

                            

#include "VxVertex/VxContainer.h"

#include "muonEvent/MuonContainer.h"

#include "egammaEvent/ElectronContainer.h"

#include "tauEvent/TauJetContainer.h"


#include "TrigSteeringEvent/HLTResult.h"
#include "TrigDecisionTool/ExpertMethods.h"

// #include "TrigSteeringEvent/TrigRoiDescriptorCollection.h"

#include "xAODTracking/TrackParticle.h"
#include "xAODTracking/TrackParticleContainer.h"

#include "TrigCompositeUtils/TrigCompositeUtils.h"





template<typename T, typename A=Analysis_Tier0>
class T_AnalysisConfigMT_Tier0 : public T_AnalysisConfig<T> {

public:

  // Full constructor: test/reference/selection
  // - analysisInstanceName: the name of the analysis chain being created
  // - xxxChainName: the name of the chain to be used as test/reference/selection; must be "StoreGate" in case of direct access to SG containers
  // - xxxType: the type of tracks to be retrieved from the test/reference/selection chain or container
  // - xxxKey:  the key for tracks to be retrieved from the test/reference/selection chain or container
  // - all standard operations are performed in loops over 0=test 1=reference 2=selection
  T_AnalysisConfigMT_Tier0(const std::string& analysisInstanceName,
			 const std::string& testChainName,      const std::string& testType,      const std::string& testKey,
			 const std::string& referenceChainName, const std::string& referenceType, const std::string& referenceKey,
			 TrackFilter*     testFilter,  TrackFilter*     referenceFilter, 
			 TrackAssociator* associator,
			 TrackAnalysis*   analysis,
			 TagNProbe*      TnP_tool = 0) :
    T_AnalysisConfig<T>( analysisInstanceName,
			 testChainName,      testType,      testKey,
			 referenceChainName, referenceType, referenceKey,
			 testFilter, referenceFilter,
			 associator,
			 analysis),
    m_useBeamCondSvc(false),
    m_doOffline(true),
    m_doMuons(false),
    m_doElectrons(false),
    m_doTaus(false),
    m_doBjets(false),
    m_hasTruthMap(false),
    m_NRois(0),
    m_NRefTracks(0),
    m_NTestTracks(0),
    m_runPurity(false),
    m_shifter(false),
    m_pTthreshold(0),
    m_first(true),
    m_containTracks(false),
    m_TnP_tool(TnP_tool),
    m_invmass(0),
    m_invmass_obj(0)
  {
    /// FIXME: the m_event should not be needed, we need to make this a local variable
    m_event = new TIDA::Event();
    m_chainNames.push_back(testChainName);

#if 0
    ChainString& chain = m_chainNames.back(); 

    std::cout << "\nT_AnalysisConfigMT_Tier0::name:                " << name() << "\t" << this << std::endl;
    std::cout <<  "T_AnalysisConfigMT_Tier0::chain specification: " << testChainName << " -> " << chain << "\t" << chain.raw() << std::endl;
    std::cout << "\tchain: " << chain.head()    << std::endl;
    std::cout << "\tkey:   " << chain.tail()    << std::endl;
    std::cout << "\troi:   " << chain.roi()     << std::endl;
    std::cout << "\tvtx:   " << chain.vtx()     << std::endl;
    std::cout << "\tte:    " << chain.element() << std::endl;

    std::cout << "\tpost:  " << chain.post()          << std::endl; 
    std::cout << "\tpt:    " << chain.postvalue("pt") << std::endl;

    std::cout << "\tcontainTracks: " << m_containTracks << std::endl; 
#endif
    
    m_testType = testType;
  }

  virtual ~T_AnalysisConfigMT_Tier0() {
    delete m_event;
    if ( m_TnP_tool != 0 ) delete m_TnP_tool ;
  }

  void setRunPurity( bool b ) { m_runPurity=b; }

  void setShifter( bool b )    { m_shifter=b; }

  void useBeamCondSvc( bool b ) { m_useBeamCondSvc = b; }

  void containTracks( bool b ) { m_containTracks = b; }

public:

  A* _analysis;

protected:

  using T_AnalysisConfig<T>::m_provider;
  using T_AnalysisConfig<T>::m_tdt;
  using T_AnalysisConfig<T>::m_mcTruth;

  using T_AnalysisConfig<T>::name;
  using T_AnalysisConfig<T>::m_analysis;

  using T_AnalysisConfig<T>::m_selectorTest;
  using T_AnalysisConfig<T>::m_selectorRef;
  using T_AnalysisConfig<T>::m_associator;
  using T_AnalysisConfig<T>::m_filters;
  
  
  virtual void loop() {
    
    bool TnP_flag = (m_TnP_tool != 0) ; // flag for tag and probe analysis
    
    if( m_provider->msg().level() <= MSG::VERBOSE) {
      m_provider->msg(MSG::VERBOSE) <<  "AnalysisConfigMT_Tier0::loop() for " << T_AnalysisConfig<T>::m_analysisInstanceName <<  endmsg;
    }

    // get (offline) beam position
    double xbeam = 0;
    double ybeam = 0;

    if ( m_first ) {      

      m_first = false;
      
      if ( m_provider->msg().level() <= MSG::VERBOSE ) {
	m_provider->msg(MSG::VERBOSE) << " using beam position\tx=" << xbeam << "\ty=" << ybeam << endmsg;
	
	std::vector<std::string> configuredChains  = (*(m_tdt))->getListOfTriggers("L2_.*, EF_.*, HLT_.*");
	
       	for ( unsigned i=0 ; i<configuredChains.size() ; i++ ) {
	  //	std::cout << "Configured chain " << configuredChains[i]  << std::endl;
	  m_provider->msg(MSG::VERBOSE)  << "Chain " << configuredChains[i]  << endmsg;
	}
      }
    
      //      std::cout << "\tloop() analyse chains " << m_chainNames.size() << std::endl;
	
      std::vector<ChainString>::iterator chainitr = m_chainNames.begin();
      
      std::vector<ChainString> chains;
      
      /// handle wildcard chain selection - but only the first time
      /// NB: also check all other chains as well - only set up an
      ///     analysis for configured chains
      while ( chainitr!=m_chainNames.end() ) {
	
        /// get chain
        ChainString& chainName = (*chainitr);

	if (m_provider->msg().level() <= MSG::VERBOSE) {
	  m_provider->msg(MSG::VERBOSE) << "process chain " << chainName << endmsg;
	}	

	if ( chainName.head() == "" ) { 
	  
	  std::string selectChain = chainName.raw(); 

	  chains.push_back( ChainString(selectChain) );
	    
	}
	else { 

	  /// get matching chains
	  std::vector<std::string> selectChains  = (*(m_tdt))->getListOfTriggers( chainName.head() );
	  	  
	  for ( unsigned iselected=0 ; iselected<selectChains.size() ; iselected++ ) {

	    selectChains[iselected] = chainName.subs( selectChains[iselected] ); 
	    
#if 0
	    std::cout << "sorting:: chain specification: " << chainName << "\traw:" << chainName.raw() << std::endl;
	    std::cout << "\tchain: " << chainName.head()    << std::endl;
	    std::cout << "\tkey:   " << chainName.tail()    << std::endl;
	    std::cout << "\tind:   " << chainName.extra()   << std::endl;
	    std::cout << "\troi:   " << chainName.roi()     << std::endl;
	    std::cout << "\tvtx:   " << chainName.vtx()     << std::endl;
	    std::cout << "\tte:    " << chainName.element() << std::endl;
#endif
	    
	    /// replace wildcard with actual matching chains ...
	    chains.push_back( ChainString(selectChains[iselected]) );
	    
	    if(m_provider->msg().level() <= MSG::VERBOSE) {
	      m_provider->msg(MSG::VERBOSE) << "Matching chain " << selectChains[iselected] << " (" << chainName.head() << ")" << endmsg;
	    }

	    //	    std::cout<< "\tMatching chain " << sc << " (" << chainName.head() << ")" << std::endl;
	    
	  }
	}

        chainitr++;
      }

      // m_chainNames.insert( m_chainNames.end(), chains.begin(), chains.end() );
      m_chainNames = chains;
      
      if(m_provider->msg().level() <= MSG::VERBOSE) {
	for ( size_t ic=m_chainNames.size() ; ic-- ;  ) m_provider->msg(MSG::VERBOSE) << "Analyse chain " << m_chainNames[ic] << endmsg;
      }      

    } /// end of first event setup
    
    //    std::cout << "\tloop() event analysis ..." << std::endl;

    /// all this should perhaps be class variables
    Filter_True filter;
    
    Filter_etaPT     filter_etaPT(5,200);
    Filter_Combined filter_truth( &filter_etaPT,   &filter_etaPT);
    
    /// will need to add a vertex filter at some point probably
    // Filter_Combined filterRef (&filter_offline, &filter_vertex);

    int iRefFilter  = 1;
    int iTestFilter = 0;

    if ( m_runPurity ) { 
      iRefFilter  = 0;
      iTestFilter = 1;
    }
    
    Filter_Combined filterRef(   m_filters[iRefFilter][0],  &filter );
    Filter_Combined filterTest(  m_filters[iTestFilter][0], &filter );


    TrigTrackSelector selectorTruth( &filter_truth );
    TrigTrackSelector selectorRef( &filterRef );
    m_selectorRef = &selectorRef;
    TrigTrackSelector selectorTest( &filterTest );
    m_selectorTest = &selectorTest;

    if ( xbeam!=0 || ybeam!=0 ) { 
      m_selectorRef->setBeamline(  xbeam, ybeam );
    }  

    /// now start everything going for this event properly ...

    // clear the ntuple TIDA::Event class
    m_event->clear();

    /// (obviously) get the event info

    const xAOD::EventInfo* pEventInfo;

    unsigned           run_number        = 0;
    unsigned long long event_number      = 0;
    unsigned           lumi_block        = 0;
    unsigned           bunch_crossing_id = 0;
    unsigned           time_stamp        = 0;
    double             mu_val            = 0;

    //    std::cout << "\tloop() get EventInfo" << std::endl;

    if ( m_provider->evtStore()->retrieve(pEventInfo).isFailure() ) {
      m_provider->msg(MSG::WARNING) << "Failed to get EventInfo " << endmsg;
    } else {

      run_number        = pEventInfo->runNumber();
      event_number      = pEventInfo->eventNumber();
      lumi_block        = pEventInfo->lumiBlock();
      time_stamp        = pEventInfo->timeStamp();
      bunch_crossing_id = pEventInfo->bcid();
      mu_val            = pEventInfo->averageInteractionsPerCrossing();
    }
  
    if(m_provider->msg().level() <= MSG::VERBOSE){
      m_provider->msg(MSG::VERBOSE) << "run "     << run_number
                                    << "\tevent " << event_number
                                    << "\tlb "    << lumi_block << endmsg;
    }
    
    //      m_provider->msg(MSG::INFO) << "run "     << run_number
    //				 << "\tevent " << event_number
    //				 << "\tlb "    << lumi_block << endmsg;

    //    std::cout << "\trun "     << run_number  << "\tevent " << event_number  << "\tlb "    << lumi_block << std::endl;


    // clear the ntuple TIDA::Event class
    m_event->clear();

    m_event->run_number(run_number);
    m_event->event_number(event_number);
    m_event->lumi_block(lumi_block);
    m_event->time_stamp(time_stamp);
    m_event->bunch_crossing_id(bunch_crossing_id);
    m_event->mu(mu_val);

    /// first check whether the chains have actually run, otherwise there's no point
    /// doing anything

    bool analyse = false;
  
    // Check HLTResult

    //    std::cout << "\tloop() loop over trigger chains to determine whether to process this event ..." << std::endl;  

    for ( unsigned ichain=0 ; ichain<m_chainNames.size() ; ichain++ ) {

      const std::string& chainname = m_chainNames[ichain].head();
 
      if ( chainname == "" ) analyse = true;
      else { 

	//Only for trigger chains
	if ( chainname.find("L2")  == std::string::npos &&
	     chainname.find("EF")  == std::string::npos &&
	     chainname.find("HLT") == std::string::npos ) continue;
	
	if ( m_provider->msg().level() <= MSG::DEBUG ) {
	  m_provider->msg(MSG::DEBUG) << "Chain "  << chainname
				      << "\tpass " << (*m_tdt)->isPassed(chainname)
				      << "\tpres " << (*m_tdt)->getPrescale(chainname) << endmsg;
	}
	
	//	std::cout << "\tChain "  << chainname << "\tpass " << (*m_tdt)->isPassed(chainname)
	//		  << "\tpres " << (*m_tdt)->getPrescale(chainname) << std::endl;
	
	if ( (*(m_tdt))->isPassed(chainname) || (*(m_tdt))->getPrescale(chainname) ) analyse = true;
	
      }
    }
    

    //    if ( (*m_tdt)->ExperimentalAndExpertMethods().isHLTTruncated() ) {
    //    m_provider->msg(MSG::WARNING) << "HLTResult truncated, skipping event" << endmsg;
    //   return;
    //  }
    
    if ( !this->m_keepAllEvents && !analyse ) {
      //     m_provider->msg(MSG::VERBOSE) << "No chains passed unprescaled - not processing this event" << endmsg;
      if(m_provider->msg().level() <= MSG::VERBOSE) {
        m_provider->msg(MSG::VERBOSE) << "No chains passed unprescaled - not processing this event" << endmsg;
      }
      //    std::cout << "\tNo chains passed unprescaled - not processing this event" << std::endl;
      return;
    }

    

    /// for Monte Carlo get the truth particles if requested to do so

    selectorTruth.clear();

    

    if(m_provider->msg().level() <= MSG::VERBOSE)
      m_provider->msg(MSG::VERBOSE) << "MC Truth flag " << m_mcTruth << endmsg;

    const TrigInDetTrackTruthMap* truthMap = 0;

    if ( m_mcTruth ) {
      if(m_provider->msg().level() <= MSG::VERBOSE ) m_provider->msg(MSG::VERBOSE) << "getting Truth" << endmsg;

      if ( m_provider->evtStore()->retrieve(truthMap, "TrigInDetTrackTruthMap").isFailure()) {
        if(m_provider->msg().level() <= MSG::VERBOSE)
          m_provider->msg(MSG::VERBOSE) << "TrigInDetTrackTruthMap not found" << endmsg;
        m_hasTruthMap = false;
      }
      else {
        if(m_provider->msg().level() <= MSG::VERBOSE)
          m_provider->msg(MSG::VERBOSE) << "TrigInDetTrackTruthMap found" << endmsg;
        m_hasTruthMap = true;
      }
    }
    

    /// get the offline vertices into our structure

    std::vector<TIDA::Vertex> vertices;
    std::vector<TIDA::Vertex> vertices_rec;

    std::vector<double> refbeamspot;
    std::vector<double> testbeamspot;

    /// fetch offline vertices ...

    m_provider->msg(MSG::VERBOSE) << "fetching AOD Primary vertex container" << endmsg;

    if ( !this->select( vertices, "PrimaryVertices" ) ) { 
      m_provider->msg(MSG::VERBOSE) << "could not retrieve vertex collection " "PrimaryVertices" << std::endl;
    }

    /// add the truth particles if needed

    if ( m_mcTruth ) {
      m_event->addChain( "Truth" );
      m_event->back().addRoi(TIDARoiDescriptor());
      m_event->back().back().addTracks(selectorTruth.tracks());
    }

    /// now add the vertices

    if ( m_doOffline ) {
      for ( unsigned i=0 ; i<vertices.size() ; i++ )  {
        if(m_provider->msg().level() <= MSG::VERBOSE)
          m_provider->msg(MSG::VERBOSE) << "vertex " << i << " " << vertices[i] << endmsg;
        m_event->addVertex(vertices[i]);
      }
    }

    /// now add the offline tracks and reco objects

    // int Noff = 0;
    std::vector<TIDA::Track*> offline_tracks;
    std::vector<TIDA::Track*> electron_tracks;
    std::vector<TIDA::Track*> muon_tracks;

    std::vector<TIDA::Track*> ref_tracks;
    std::vector<TIDA::Track*> test_tracks;

    offline_tracks.clear();
    electron_tracks.clear();
    muon_tracks.clear();

    ref_tracks.clear();
    test_tracks.clear();

    // offline track retrieval now done once for each chain rather than each roi
    if ( m_provider->msg().level() <= MSG::VERBOSE )
      m_provider->msg(MSG::VERBOSE) << "MC Truth flag " << m_mcTruth << endmsg;

    bool foundTruth = false;

    if ( !m_doOffline && m_mcTruth ) {

      filter_truth.setRoi( 0 ); // don't filter on RoI yet (or until needed)  

      selectorTruth.clear();

      if ( m_provider->msg().level() <= MSG::VERBOSE )
	m_provider->msg(MSG::VERBOSE) << "getting Truth" << endmsg;

      if ( m_provider->evtStore()->template contains<TruthParticleContainer>("INav4MomTruthEvent") ) {
	//ESD
	this->template selectTracks<TruthParticleContainer>( &selectorTruth, "INav4MomTruthEvent" );
	foundTruth = true;
      }
      else if ( m_provider->evtStore()->template contains<TruthParticleContainer>("SpclMC") ) {
	/// AOD
	this->template selectTracks<TruthParticleContainer>( &selectorTruth, "SpclMC");
	foundTruth = true;
      }
      else if ( m_provider->evtStore()->template contains<xAOD::TruthParticleContainer>("TruthParticles") ) {
	/// xAOD::TruthParticles
	this->template selectTracks<xAOD::TruthParticleContainer>( &selectorTruth, "TruthParticles");
	foundTruth = true;
      }
      else if ( m_provider->evtStore()->template contains<TruthParticleContainer>("") ) {
	/// anything else?
	this->template selectTracks<TruthParticleContainer>( &selectorTruth, "");
	foundTruth = true;
      }
      else
	if ( m_provider->msg().level() <= MSG::VERBOSE ) {
	  m_provider->msg(MSG::VERBOSE) << "Truth not found - none whatsoever!" << endmsg;
	}
    }
      
    if ( !m_doOffline && m_mcTruth && !foundTruth ) {
          
      if ( m_provider->msg().level() <= MSG::VERBOSE ) { 
	m_provider->msg(MSG::VERBOSE) << "getting Truth" << endmsg;
      }

      /// selectTracks<TruthParticleContainer>( &selectorTruth, "INav4MomTruthEvent" );

      const DataHandle<McEventCollection> mcevent;

      /// now as a check go through the GenEvent collection

      std::string keys[4] = { "GEN_AOD", "TruthEvent", "", "G4Truth" };

      std::string key = "";

      bool foundcollection = false;

      for ( int ik=0 ; ik<4 ; ik++ ) {
         
	if ( m_provider->msg().level() <= MSG::VERBOSE ) {
	  m_provider->msg(MSG::VERBOSE) << "Try McEventCollection: " << keys[ik] << endmsg;
	}

	if ( !m_provider->evtStore()->template contains<McEventCollection>(keys[ik]) ) {
	  if( m_provider->msg().level() <= MSG::VERBOSE )
	    m_provider->msg(MSG::VERBOSE) << "No McEventCollection: " << keys[ik] << endmsg;
	  continue;
	}

	if ( m_provider->msg().level() <= MSG::VERBOSE )
	  m_provider->msg(MSG::VERBOSE) << "evtStore()->retrieve( mcevent, " << keys[ik] << " )" << endmsg;

	if ( m_provider->evtStore()->template retrieve( mcevent, keys[ik] ).isFailure() ) {
	  if ( m_provider->msg().level() <= MSG::VERBOSE )
	    m_provider->msg(MSG::VERBOSE) << "Failed to get McEventCollection: " << keys[ik] << endmsg;
	}
	else {
	  /// found this key
	  key = keys[ik];
	  if(m_provider->msg().level() <= MSG::VERBOSE)
	    m_provider->msg(MSG::VERBOSE) << "Found McEventCollection: " << key << endmsg;
	  foundcollection = true;
	  break;
	}
      }

      /// not found any truth collection
      if ( !foundcollection ) {
	if(m_provider->msg().level() <= MSG::VERBOSE)
	  m_provider->msg(MSG::WARNING) << "No MC Truth Collections of any sort, whatsoever!!!" << endmsg;

	//    m_tree->Fill();
	//    return StatusCode::FAILURE;

	return;
      }

      if ( m_provider->msg().level() <= MSG::VERBOSE ) {
	m_provider->msg(MSG::VERBOSE) << "Found McEventCollection: " << key << "\tNevents " << mcevent->size() << endmsg;
      }

      McEventCollection::const_iterator evitr = mcevent->begin();
      McEventCollection::const_iterator evend = mcevent->end();

      unsigned ie = 0; /// count of "events" - or interactions
      unsigned ip = 0; /// count of particles

      unsigned ie_ip = 0; /// count of "events with some particles"

      while ( evitr!=evend ) {

	int _ip = 0; /// count of particles in this interaction

	int pid = HepMC::signal_process_id((*evitr));

	//The logic should be clarified here
	if ( pid!=0 ) { /// hooray! actually found a sensible event

	  for (auto pitr: *(*evitr) ) { 

	    selectorTruth.selectTrack( pitr );

	    ++_ip;
                
	  }

	}
	++ie;
	++evitr;

	if ( _ip>0 ) {
	  /// if there were some particles in this interaction ...
	  //      m_provider->msg(MSG::VERBOSE) << "Found " << ie << "\tpid " << pid << "\t with " << ip << " TruthParticles (GenParticles)" << endmsg;
	  ++ie_ip;
	  ip += _ip;
	}
      }

      if(m_provider->msg().level() <= MSG::VERBOSE){
	m_provider->msg(MSG::VERBOSE) << "Found " << ip << " TruthParticles (GenParticles) in " << ie_ip << " GenEvents out of " << ie << endmsg;
	m_provider->msg(MSG::VERBOSE) << "selected " << selectorTruth.size() << " TruthParticles (GenParticles)" << endmsg;
      }

      if(selectorTruth.size() > 0) foundTruth = true;

      if ( !(ip>0) ) {
	if (m_provider->msg().level() <= MSG::VERBOSE) m_provider->msg(MSG::WARNING) << "NO TRUTH PARTICLES - returning" << endmsg;
	return; /// need to be careful here, if not requiring truth *only* should not return
      }
	  
    }

    // m_provider->msg(MSG::VERBOSE) << " Offline tracks " << endmsg;
	
    if ( m_doOffline ) {
	  
      if ( m_provider->evtStore()->template contains<xAOD::TrackParticleContainer>("InDetTrackParticles") ) {
	this->template selectTracks<xAOD::TrackParticleContainer>( m_selectorRef, "InDetTrackParticles" );
	refbeamspot = this->template getBeamspot<xAOD::TrackParticleContainer>( "InDetTrackParticles" );
      }
      else if (m_provider->evtStore()->template contains<Rec::TrackParticleContainer>("TrackParticleCandidate") ) {
	this->template selectTracks<Rec::TrackParticleContainer>( m_selectorRef, "TrackParticleCandidate" );
      }
      else if ( m_provider->msg().level() <= MSG::WARNING ) {
	m_provider->msg(MSG::WARNING) << " Offline tracks not found " << endmsg;
      }

    }


    //    std::cout << "\tloop() loop over chains proper ..." << std::endl;

    /// now loop over all relevant chains to get the trigger tracks...
    for ( unsigned ichain=0 ; ichain<m_chainNames.size() ; ichain++ ) {

      /// create chains for ntpl

      // std::string& chainname = chains[ichain];
      const std::string& chainname = m_chainNames[ichain].head();
      const std::string&       key = m_chainNames[ichain].tail();
      const std::string&  vtx_name = m_chainNames[ichain].vtx();
      //Not used left just in case
      //const std::string&  roi_name = m_chainNames[ichain].roi();
      //const std::string&  te_name = m_chainNames[ichain].element();
      m_pTthreshold = 0;  /// why does this need to be a class variable ???

      if ( m_chainNames[ichain].postcount() ) { 
        std::string ptvalue = m_chainNames[ichain].postvalue("pt");
	if ( ptvalue!="" ) m_pTthreshold = std::stod(ptvalue);
      }


      unsigned _decisiontype = TrigDefs::Physics;
      unsigned  decisiontype;
     
      if ( this->requireDecision() ) _decisiontype =  TrigDefs::requireDecision;
      
      
      if ( m_chainNames[ichain].passed() ) decisiontype = _decisiontype;
      else                                 decisiontype = TrigDefs::alsoDeactivateTEs;

      /// useful debug information to be kept in for the time being 
      //      if ( decisiontype==TrigDefs::requireDecision ) std::cout << "\tSUTT TrigDefs::requireDecision " << decisiontype << std::endl;
      //      if ( decisiontype==TrigDefs::Physics )         std::cout << "\tSUTT TrigDefs::Physics "         << decisiontype << std::endl;

      if ( chainname!="" && m_provider->msg().level() <= MSG::VERBOSE ) {

        m_provider->msg(MSG::VERBOSE) << "status for chain " << chainname
                                      << "\tpass "           << (*m_tdt)->isPassed(chainname)
                                      << "\tprescale "       << (*m_tdt)->getPrescale(chainname) << endmsg;

        m_provider->msg(MSG::VERBOSE) << "fetching features for chain " << chainname << endmsg;
    
        m_provider->msg(MSG::VERBOSE) << chainname << "\tpassed: " << (*m_tdt)->isPassed( chainname ) << endmsg;
      }

      /// useful debug information to be kept in
      //      std::cout << "\tstatus for chain " << chainname
      //		<< "\tpass "           << (*m_tdt)->isPassed( chainname )
      //		<< "\tpassdt "         << (*m_tdt)->isPassed( chainname, decisiontype )
      //		<< "\tprescale "       << (*m_tdt)->getPrescale( chainname ) << std::endl;
	

      //   m_provider->msg(MSG::INFO) << chainname << "\tpassed: " << (*m_tdt)->isPassed( chainname ) << "\t" << m_chainNames[ichain] << "\trun " << run_number << "\tevent " << event_number << endmsg;


      if ( chainname!="" && !this->m_keepAllEvents && !(*m_tdt)->isPassed( chainname, decisiontype ) ) continue;


      /// Get chain combinations and loop on them
      /// - loop made on chain selected as the one steering RoI creation
      // Trig::FeatureContainer f = (*m_tdt)->features( chainname, TrigDefs::alsoDeactivateTEs);
     
      /// only use the TDT for extracting collections if this was a trigger analysis 
      /// for fullscan "offline" type analyses (ie fullscan FTK) do not use this

      // tag and probe analysis processes multiple chains passed in the tag and probe tool at once so loop over vector of chains
      std::vector<std::string> chainNames ;
            
      if ( !TnP_flag ) {
	chainNames.push_back(m_chainNames[ichain].raw()) ;
      }
      else {
	chainNames.push_back(m_TnP_tool->tag()) ;
	chainNames.push_back(m_TnP_tool->probe()) ;
      }

      // loop over new chainNames vector but doing the same stuff
      for ( unsigned i = 0 ; i < chainNames.size() ; ++i ) {

	ChainString chainConfig = chainNames[i] ;
	std::string chainName = chainConfig.head();
	
	m_event->addChain( chainNames[i] ); 
	
	TIDA::Chain& chain = m_event->back();
	
	if ( chainName == "" ) { 

	  /// do we still want the blind chain access for track collections ???

	  m_selectorTest->clear();
	
	  /// dummy full scan chain 

	  TIDARoiDescriptor* roiInfo = new TIDARoiDescriptor(true);
		
	  chain.addRoi( *roiInfo );
	
	  if ( m_provider->evtStore()->template contains<xAOD::TrackParticleContainer>(key) ) {
	    this->template selectTracks<xAOD::TrackParticleContainer>( m_selectorTest, key );
	    refbeamspot = this->template getBeamspot<xAOD::TrackParticleContainer>( key );
	  }

	  const std::vector<TIDA::Track*>& testtracks = m_selectorTest->tracks();

	  chain.back().addTracks(testtracks);
	
	  if ( vtx_name!="" ) { 
	  
	    /// MT Vertex access

	    m_provider->msg(MSG::VERBOSE) << "\tFetch xAOD::VertexContainer with key " << vtx_name << endmsg;

	    std::vector<TIDA::Vertex> tidavertices;

	    if ( this->select( tidavertices, vtx_name ) ) chain.back().addVertices( tidavertices );	 
	  }
	

	  if ( roiInfo ) delete roiInfo;

	}
	else {

	  /// new Roi based feature access
	  
	  //std::string roi_key = m_chainNames[ichain].roi();

	  std::string roi_key = chainConfig.roi();
	
	  unsigned feature_type =TrigDefs::lastFeatureOfType;

	  if ( roi_key!="" ) feature_type= TrigDefs::allFeaturesOfType;

	  /// new FeatureRequestDescriptor with leg access
	  
	  int leg = -1;
	  
	  if ( chainConfig.element()!="" ) { 
	    leg = std::atoi(chainConfig.element().c_str());
	  }
	  
	  std::vector< TrigCompositeUtils::LinkInfo<TrigRoiDescriptorCollection> > rois = 
	    (*m_tdt)->template features<TrigRoiDescriptorCollection>( Trig::FeatureRequestDescriptor( chainName,  
												      decisiontype,
												      roi_key, 
												      feature_type,
												      "roi", 
												      leg ) );
	  
	  //	const unsigned int featureCollectionMode = const std::string& navElementLinkKey = "roi") const;
	  
	  int iroi = 0; /// count of how many rois processed so far
	  
	  for ( const TrigCompositeUtils::LinkInfo<TrigRoiDescriptorCollection>& roi_info : rois ) {
	  
	    iroi++;

	    /// don't extract any additional rois if a superRoi is requested: 
	    /// In this case, the superRoi would be shared between the different 
	    /// chains 

	    if ( roi_key=="SuperRoi" && iroi>1 ) continue; 
	
	    //	  std::cout << "\troi: get link " << roi_key << " ..." << std::endl;

	    const ElementLink<TrigRoiDescriptorCollection> roi_link = roi_info.link;

	    /// check this is not a spurious TDT match
	    if ( roi_key!="" && roi_link.dataID()!=roi_key ) continue;

	    const TrigRoiDescriptor* const* roiptr = roi_link.cptr();

	    if ( roiptr == 0 ) { 
	      continue;
	    }  

	    //	  std::cout << "\troi: link deref ..." << *roiptr << std::endl;

	    if (m_provider->msg().level() <= MSG::VERBOSE) {
	      m_provider->msg(MSG::VERBOSE) << " RoI descriptor for seeded chain " << chainname << " " << **roiptr << endmsg;
	    }
	 
	    TIDARoiDescriptor* roiInfo = new TIDARoiDescriptor( TIDARoiDescriptorBuilder(**roiptr) );
	 
	    //	  if ( dbg ) std::cout << "\troi " << iroi << " " << *roiInfo << std::endl;

	    /// get the tracks 

	    m_selectorTest->clear();

	    if ( this->template selectTracks<xAOD::TrackParticleContainer>( m_selectorTest, roi_link,  key ) ) { } 

	    // beamspot stuff not needed for xAOD::TrackParticles

	    /// create analysis chain

	    chain.addRoi( *roiInfo );
	
	    /// get tracks 

	    const std::vector<TIDA::Track*>& testtracks = m_selectorTest->tracks();

	    chain.back().addTracks(testtracks);
	

	    /// now get the vertices 
	  
	    if ( vtx_name!="" ) { 

	      std::vector<TIDA::Vertex> tidavertices;

	      this->select( tidavertices, roi_link, vtx_name );

	      chain.back().addVertices( tidavertices );
	    
	    } /// retrieve online vertices
	  

#if 0
	    if ( dbg ) { 
	      std::cout << "\tTIDA analysis for chain: " << chainname << "\t key: " << key << "\t" << **roiptr << std::endl;
	      std::cout << "\tcollections: " << chain.back() << std::endl; 
	    }
#endif
	  
	    if ( roiInfo ) delete roiInfo;
	  
	  }


	} /// "offline" of "roi" type chains

      } // end of loop chainNames vector loop
	
	
      if ( m_provider->msg().level() <= MSG::VERBOSE ) {
	m_provider->msg(MSG::VERBOSE) << "event: " << *m_event << endmsg;
      }

    }

    // close previous loop over chains and open new one
    
    for ( unsigned ichain=0 ; ichain<m_event->size() ; ichain++ ) {
      
      TIDA::Chain& chain = (*m_event)[ichain];
      ChainString chainConfig(chain.name());        
      const std::string&  vtx_name = chainConfig.vtx();

      // skip tag chains to avoid performing standard analysis on them (done for tnp at the same time as probes)
      if ( TnP_flag && chainConfig.extra().find("_tag")!=std::string::npos ) continue ;
      
      std::vector<TIDA::Roi*> rois ;
      
      if (TnP_flag) {
	rois = m_TnP_tool->GetRois( m_event->chains(), m_selectorRef, &filterRef, m_invmass, m_invmass_obj );
	// needs to be done AFTER retrieving offline tracks as m_selectorRef passed as arguement, hence restructuring
      }
      else {
        rois.reserve( chain.size() );
        for ( size_t ir=0 ; ir<chain.size() ; ir++ ) {
	  rois.push_back( &(chain.rois()[ir]) );
	}
      }

      // now loop over the rois (again)

      for ( unsigned iroi=0 ; iroi<rois.size() ; iroi++ ) {

	if ( this->filterOnRoi() ) { 
	  filterRef.setRoi( &(rois.at(iroi)->roi() ) ); 
	  filterRef.containtracks( m_containTracks ); 
	}
	else filterRef.setRoi( 0 );
	
	test_tracks.clear();
	
	// this block is before the track retrieval in the original, is it working the same here?
	
	/// This is nonsense and needs restructuring - why is the truth and offline selection 
	/// done within this RoI loop? It means the complete offline and truth tracks will be 
	/// retrieved for every RoI ! really we should have the structure 
	///   
	///   - check_a_trigger_chain_has_passed
	///   - get_offline_or_truth_particles
	///   - loop_over_rois
	///     - get_trigger_tracks
	///     - filter_offline_or_truth_reference
	///     - match_tracks
	///     - call_analyis_routine
	///
	/// will leave as it is for the time being

        if ( m_provider->msg().level() <= MSG::VERBOSE )
          m_provider->msg(MSG::VERBOSE) << "MC Truth flag " << m_mcTruth << endmsg;

	if ( !m_doOffline && m_mcTruth ) {
	  if ( this->filterOnRoi() )  filter_truth.setRoi( &(rois.at(iroi)->roi() ) );
	  ref_tracks = m_selectorRef->tracks(&filter_truth);
	}
	else { // ie. if ( m_doOffline )
	  ref_tracks = m_selectorRef->tracks(&filterRef) ; 
	}
	
	if ( m_provider->msg().level() <= MSG::VERBOSE ) {
	  m_provider->msg(MSG::VERBOSE) << "ref tracks.size() " << m_selectorRef->tracks().size() << endmsg;
	  for ( int ii=m_selectorRef->tracks().size() ; ii-- ; ) {
	    m_provider->msg(MSG::VERBOSE) << "  ref track " << ii << " " << *m_selectorRef->tracks()[ii] << endmsg;
	  }
	}

        test_tracks.clear();


        for ( unsigned itrk=0 ; itrk<rois.at(iroi)->tracks().size() ; itrk++ ) {
          test_tracks.push_back(&(rois.at(iroi)->tracks().at(itrk)));
        }	
	
	//	std::cout << "sutt track multiplicities: offline " << offline_tracks.size() << "\ttest " << test_tracks.size() << std::endl;

        _analysis->setvertices( vertices.size() );  /// what is this for ???

	if ( refbeamspot.size()>0 )  _analysis->setBeamRef(   refbeamspot ); 
	if ( testbeamspot.size()>0 ) _analysis->setBeamTest( testbeamspot ); 

	/// if we want a purity, we need to swap round which tracks are the 
	/// reference tracks and which the test tracks

	if ( m_runPurity ) {

	  if ( this->getUseHighestPT() ) HighestPTOnly( test_tracks );

	  if ( m_pTthreshold>0 ) FilterPT( test_tracks, m_pTthreshold );

	  /// stats book keeping 
	  m_NRois++;
	  m_NRefTracks  += test_tracks.size();
	  m_NTestTracks += ref_tracks.size();

	  	  	  
	  /// match test and reference tracks
	  m_associator->match( test_tracks, ref_tracks );
	  
	  _analysis->execute( test_tracks, ref_tracks, m_associator );

	}
	else { 

	  /// filter on highest pt track only if required
	  if ( this->getUseHighestPT() ) HighestPTOnly( ref_tracks );

	  /// ignore all tracks belong the specific analysis pt threshold if set
	  
	  if ( m_pTthreshold>0 )         FilterPT( ref_tracks, m_pTthreshold );

	  /// stats book keeping 
	  m_NRois++;
	  m_NRefTracks  += ref_tracks.size();
	  m_NTestTracks += test_tracks.size();
	  
	  /// match test and reference tracks
	  m_associator->match( ref_tracks, test_tracks );
	
	  _analysis->setroi( &rois.at(iroi)->roi() );  
	  _analysis->execute( ref_tracks, test_tracks, m_associator );

	  if ( vtx_name!="" ) { 
	    /// get vertices for this roi - have to copy to a vector<Vertex*>
	    std::vector<TIDA::Vertex> vr = rois.at(iroi)->vertices();
	    std::vector<TIDA::Vertex*> vtx_rec;    
	    for ( unsigned iv=0 ; iv<vr.size() ; iv++ ) vtx_rec.push_back( &vr[iv] );

	    std::vector<TIDA::Vertex*> vtx;
	    if ( this->getVtxIndex()<0 ) { 
	      for ( unsigned iv=0 ; iv<vertices.size() ; iv++ ) vtx.push_back( &vertices[iv] );
	    }
	    else { 
	      if ( vertices.size()>unsigned(this->getVtxIndex()) ) vtx.push_back( &vertices[this->getVtxIndex()] );
	    }

	    _analysis->execute_vtx( vtx, vtx_rec, m_event );
	  }

	}
 
	if ( _analysis->debug() ) { 
	  m_provider->msg(MSG::INFO) << "Missing track for " << m_chainNames[ichain]  
				     << "\trun "             << run_number 
				     << "\tevent "           << event_number 
				     << "\tlb "              << lumi_block << endmsg;     
	}

      }

    }
    
    if ( m_provider->msg().level() <= MSG::VERBOSE ) {
      m_provider->msg(MSG::VERBOSE) << "\n\nEvent " << *m_event << endmsg;
    }
  }



  virtual void book() {
    
    if(m_provider->msg().level() <= MSG::VERBOSE)
      m_provider->msg(MSG::VERBOSE) << "AnalysisConfigMT_Tier0::book() " << name() << endmsg;

    // get the TriggerDecisionTool

    if( m_tdt->retrieve().isFailure() ) {
      if(m_provider->msg().level() <= MSG::ERROR)
        m_provider->msg(MSG::ERROR) << " Unable to retrieve the TrigDecisionTool: Please check job options file" << endmsg;
      // return StatusCode::FAILURE;
      return;
    }

    if(m_provider->msg().level() <= MSG::VERBOSE) {
      m_provider->msg(MSG::VERBOSE) << " Successfully retrived the TrigDecisionTool"  << endmsg;
    }


    /// get list of configured triggers
    if (m_provider->msg().level() <= MSG::VERBOSE) {
      std::vector<std::string> configuredChains  = (*(m_tdt))->getListOfTriggers("L2_.*, EF_.*, HLT_.*");

      m_provider->msg(MSG::VERBOSE)  << "Configured chains" << endmsg;
      for ( unsigned i=0 ; i<configuredChains.size() ; i++ ) {
        if( m_provider->msg().level() <= MSG::VERBOSE)
          m_provider->msg(MSG::VERBOSE)  << " Chain " << configuredChains[i]  << endmsg;
      }
    }


    // std::vector<std::string> chains;
    // chains.reserve( m_chainNames.size() );

    std::vector<ChainString>::iterator chainitr = m_chainNames.begin();

    std::vector<ChainString> chains;

    /// handle wildcard chain selection - but only the first time
    /// NB: also check all other chains as well - only set up an
    ///     analysis for configured chains
    while ( chainitr!=m_chainNames.end() ) {

      /// get chain
      ChainString& chainName = (*chainitr);

      if ( chainName.head() == "" ) {
	
	std::string selectChain = chainName.raw();

	/// replace wildcard with actual matching chains ...
	chains.push_back( selectChain );

      }
      else { 
	
	/// get matching chains
	std::vector<std::string> selectChains  = (*(m_tdt))->getListOfTriggers( chainName.head() );
	
	// std::cout << "selected chains " << selectChains.size() << std::endl;
	
	// if ( selectChains.size()==0 ) m_provider->msg(MSG::WARNING) << "No chains matched for  " << chainName << endmsg;
	
	for ( unsigned iselected=0 ; iselected<selectChains.size() ; iselected++ ) {
	  
	  selectChains[iselected] = chainName.subs( selectChains[iselected] ); 

	  /// replace wildcard with actual matching chains ...
	  chains.push_back( ChainString(selectChains[iselected]) );
	  
	  if(m_provider->msg().level() <= MSG::VERBOSE) {
	    m_provider->msg(MSG::VERBOSE) << "Matching chain " << selectChains[iselected] << " (" << chainName.head() << endmsg;
	  }
	  
	}
      }
	
      chainitr++;
    }

    m_chainNames = chains;

    for ( unsigned ic=0 ; ic<m_chainNames.size() ; ic++ ) {

      if ( ic>0 ) { 
	m_provider->msg(MSG::WARNING) << "more than one chain configured for this analysis - skipping " << m_chainNames[ic] << endmsg;
	continue;
      }

      m_provider->msg(MSG::VERBOSE) << "Analyse chain " << m_chainNames[ic] << endmsg;

      // m_provider->msg(MSG::VERBOSE)  << "--------------------------------------------------" << endmsg;
      
      std::string folder_name = "";
      
      if ( name()!="" )  folder_name = name(); 
      else               folder_name = "HLT/TRIDT/IDMon";  
      
      // unsigned decisiontype;
      // if ( m_chainNames.at(0).passed() ) decisiontype = TrigDefs::Physics;
      // else                               decisiontype = TrigDefs::alsoDeactivateTEs;
      
      /// folder_name.erase(0,3); // erase "L2_" or "EF_" so histograms all go in same chain folder - NO!!!! this means 
      /// they will over write, unless eg L2_, EF_, HLT_ etc is include in the histogram name 
      
      /// don't use test_type now ? 
      if( m_testType != "" ) folder_name = folder_name + "/" + m_testType;
      
      std::string mongroup;

      
      if ( name().find("Shifter")!=std::string::npos || m_shifter ) {
	/// shifter histograms - do not encode chain names
	if      ( m_chainNames.at(ic).tail().find("_FTF") != std::string::npos )              mongroup = folder_name + "/FTF";
	else if ( m_chainNames.at(ic).tail().find("_IDTrig") != std::string::npos || 
		  m_chainNames.at(ic).tail().find("_EFID") != std::string::npos )             mongroup = folder_name + "/EFID";
	else if ( m_chainNames.at(ic).tail().find("InDetTrigParticle") != std::string::npos ) mongroup = folder_name + "/EFID_RUN1";
	else if ( m_chainNames.at(ic).tail().find("_GSF")      != std::string::npos )         mongroup = folder_name + "/GSF";
	else                                                                                  mongroup = folder_name + "/Unknown";

	if ( m_chainNames.at(ic).vtx()!="" ) mongroup += "/" + m_chainNames.at(ic).vtx();

      }
      else { 
	/// these are the Expert / non-Shifter histograms - encode the full chain names

	if ( m_chainNames[ic].head() == "" ) mongroup = folder_name + "/Fullscan";
	else                                 mongroup = folder_name + "/" + m_chainNames[ic].head();
	std::string track_collection = ""; 

	if ( m_chainNames.at(ic).tail()!="" )  { 
	  track_collection =  "/" + m_chainNames.at(ic).tail();
	  if ( m_chainNames.at(ic).extra()!="" ) track_collection += "_" + m_chainNames.at(ic).extra();
	}

	if ( m_chainNames.at(ic).roi()!="" ) { 
	  if ( track_collection!="" ) track_collection += "_" + m_chainNames[ic].roi();
	  else                        track_collection = "/" + m_chainNames[ic].roi();
	}

	if ( m_chainNames.at(ic).vtx()!="" ) { 
	  if ( track_collection!="" ) track_collection += "_" + m_chainNames[ic].vtx();
	  else                        track_collection  = "/" + m_chainNames[ic].vtx();
	}

	/// add trigger element and roi descriptor names
	if ( m_chainNames.at(ic).element()!="" ) { 
	  if ( track_collection!="" ) track_collection += "_" + m_chainNames[ic].element();
	  else                        track_collection  = "/" + m_chainNames[ic].element();
	}
	
	if ( track_collection!="" )  mongroup += track_collection;

	if ( !m_chainNames.at(ic).passed() )      mongroup += "/DTE";

	//	std::cout << "\n SUTT chain " << m_chainNames.at(ic) << "\tvtx " << m_chainNames.at(ic).vtx() << "\tmongroup " << mongroup << std::endl;
	
      }

      //      std::cout << "SUTT chain " << "\tvtx " << m_chainNames.at(ic).vtx() << "\tmongroup " << mongroup << std::endl;
      
      m_provider->msg(MSG::VERBOSE) << " book mongroup " << mongroup << endmsg;
      
#     ifdef ManagedMonitorToolBase_Uses_API_201401
      m_provider->addMonGroup( new ManagedMonitorToolBase::MonGroup( m_provider, mongroup, ManagedMonitorToolBase::run ) );
#     else
      m_provider->addMonGroup( new ManagedMonitorToolBase::MonGroup( m_provider, mongroup, 
								     ManagedMonitorToolBase::shift, 
								     ManagedMonitorToolBase::run ) );
#   endif
      
      _analysis = dynamic_cast<A*>(m_analysis);
 
      if ( monTool() ) _analysis->set_monTool( monTool() );

      m_analysis->initialise();
      
      _analysis->setevent( m_event ); 

      
      std::map<std::string, TH1*>::const_iterator hitr = _analysis->THbegin();
      std::map<std::string, TH1*>::const_iterator hend = _analysis->THend();
      
      //    std::cout << "\tsutt adding to mongroup   " << mongroup << std::endl;

      // booking invariant mass histograms for tag and probe analysis
      if ( m_TnP_tool != 0 ) {
	m_invmass = new TH1F( "invmass", "invariant mass;mass [GeV]", 320, 0, 200 );                                                                                                      
	m_provider->addHistogram( m_invmass );                            

	m_invmass_obj = new TH1F( "invmass_obj", "invariant mass;mass [GeV]", 320, 0, 200 );                                                                                                      
	m_provider->addHistogram( m_invmass_obj );                            
	///}       
	//	m_TnP_tool->BookMinvHisto() ;
	//	m_provider->addHistogram( m_TnP_tool->GetMinvHisto(), mongroup ) ;
	//	m_provider->addHistogram( m_TnP_tool->GetMinvObjHisto(), mongroup ) ;
      }
      
      while ( hitr!=hend ) {
	//  std::cout << "\tsutt addHisto " << hitr->second->GetName() << std::endl;
	m_provider->addHistogram( hitr->second, mongroup );
	hitr++;
      }
      
      std::map<std::string, TProfile*>::const_iterator effitr = _analysis->TEffbegin();
      std::map<std::string, TProfile*>::const_iterator effend = _analysis->TEffend();
      
      while ( effitr!=effend ) {
	// std::cout << "\tsutt addProfile " << effitr->second->GetName() << std::endl;
	m_provider->addHistogram( effitr->second, mongroup );
	effitr++;
      }
      
      if(m_provider->msg().level() <= MSG::VERBOSE) {
	m_provider->msg(MSG::VERBOSE) << "AnalysisConfigMT_Tier0::book() done" << endmsg;
      }
    }

  }



  virtual void finalize() {

    if(m_provider->msg().level() <= MSG::VERBOSE){
      m_provider->msg(MSG::VERBOSE) << "AnalysisConfigMT_Tier0::finalise() " << m_provider->name() << endmsg;
    }
    
    m_analysis->finalise();

    // deleting instance of TnP_tool and setting pointer to null
    if ( m_TnP_tool != 0 ) {
      delete m_TnP_tool ;
      m_TnP_tool = 0 ;
    }
    
    m_provider->msg(MSG::INFO) << m_provider->name() << " " << m_chainNames[0] << "   \tNRois processed: " << m_NRois << "\tRef tracks: " << m_NRefTracks << "\tTestTracks: " << m_NTestTracks << endmsg;

    if(m_provider->msg().level() <= MSG::VERBOSE) {
      m_provider->msg(MSG::VERBOSE) << m_provider->name() << " finalised" << endmsg;
    }
  }

  void set_monTool( ToolHandle<GenericMonitoringTool>* m ) { m_monTool=m; }

  ToolHandle<GenericMonitoringTool>* monTool() { return m_monTool; }

protected:

  bool           m_useBeamCondSvc;

  TIDA::Event*  m_event;

  TFile*    mFile;
  TTree*    mTree;

  std::vector<ChainString>     m_chainNames;
  std::vector<A*> m_analyses;
  std::string m_testType;

  bool m_doOffline;
  bool m_doMuons;
  bool m_doElectrons;
  bool m_doTaus;
  bool m_doBjets;
  bool m_hasTruthMap;
  bool m_doTauThreeProng;
  bool m_tauEtCutOffline;

  std::string m_outputFileName;

  /// output stats
  int m_NRois;
  int m_NRefTracks;
  int m_NTestTracks;

  bool m_runPurity;

  bool m_shifter;

  double m_pTthreshold;

  bool   m_first;
  
  bool   m_containTracks;

  TagNProbe* m_TnP_tool; 

  TH1F*      m_invmass;
  TH1F*      m_invmass_obj;

  ToolHandle<GenericMonitoringTool>* m_monTool;

};



#endif  // TrigInDetAnalysisExample_T_AnalysisConfigMT_Tier0_H

