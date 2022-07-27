/**
 **     @file    AnalysisConfigMT_Ntuple.cxx
 **
 **     @author  mark sutton
 **
 **     Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
 **/

#include "TrigInDetAnalysisExample/AnalysisConfigMT_Ntuple.h"

#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthVertexContainer.h"

#include "xAODEventInfo/EventInfo.h"

#include "TrigDecisionTool/FeatureRequestDescriptor.h"

#include "TrigInDetAnalysis/Filter_AcceptAll.h"
#include "TrigInDetAnalysisUtils/Filter_etaPT.h"
#include "TrigInDetAnalysisUtils/Filter_RoiSelector.h"
#include "TrigInDetAnalysisUtils/Filters.h"
#include "TrigInDetAnalysisUtils/TIDAVertexBuilder.h"

#include "TrigInDetAnalysis/TIDDirectory.h"
#include "TrigInDetAnalysisUtils/TIDARoiDescriptorBuilder.h"


std::string date();


//function to find true taus
HepMC::ConstGenParticlePtr fromParent( int pdg_id, HepMC::ConstGenParticlePtr p, bool printout=false );
  


template<class T>
void remove_duplicates(std::vector<T>& vec) {
  std::sort(vec.begin(), vec.end());
  vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}




void AnalysisConfigMT_Ntuple::loop() {

        m_provider->msg(MSG::DEBUG) << "[91;1m" << "AnalysisConfigMT_Ntuple::loop() for " << m_analysisInstanceName 
				   << " compiled " << __DATE__ << " " << __TIME__ << "\t: " << date() << "[m" << endmsg;


	bool foundOffline = false;

	// get (offline) beam position
	double xbeam = 0;
	double ybeam = 0;
	double zbeam = 0;
	std::vector<double> beamline;


	// get (online) beam position
	double xbeam_online = 0;
	double ybeam_online = 0;
	double zbeam_online = 0;

	std::vector<double> beamline_online;

	//	m_provider->msg(MSG::INFO) << " offline beam position\tx=" << xbeam        << "\ty=" << ybeam        << "\tz=" << zbeam        << endmsg; 
	//	m_provider->msg(MSG::INFO) << " online  beam position\tx=" << xbeam_online << "\ty=" << ybeam_online << "\tz=" << zbeam_online << endmsg; 

	/// list the configured chains

	static std::set<std::string> configuredHLTChains;

	std::vector<ChainString> chainNames;

	if ( m_tida_first ) { 

		std::vector<std::string> configuredChains  = (*m_tdt)->getListOfTriggers("L2_.*, EF_.*, HLT_.*");

		if (m_provider->msg().level() <= MSG::VERBOSE) {
		  m_provider->msg(MSG::VERBOSE) << "[91;1m" << configuredChains.size() << " Configured Chains" << "[m" << endmsg;
		}

		for ( unsigned i=0 ; i<configuredChains.size() ; i++ ) { 
		  if (m_provider->msg().level() <= MSG::VERBOSE) {
		    m_provider->msg(MSG::VERBOSE) << "[91;1m" << "Chain " << configuredChains[i] << "   (ACN)[m" << endmsg;
		  }
		  configuredHLTChains.insert( configuredChains[i] );
		}

		m_tida_first = false;

		std::vector<ChainString>::iterator chainitr = m_chainNames.begin();

		/// handle wildcard chain selection - but only the first time
		while ( chainitr!=m_chainNames.end() ) {

		  /// get chain
		  ChainString& chainName = (*chainitr);

		  /// get matching chains
		  
		  std::vector<std::string> selectChains;
		  selectChains.clear();
		  if ( chainitr->head()=="" ) selectChains.push_back("");
		  else                        selectChains = (*m_tdt)->getListOfTriggers( chainName.head() );
		  
		  for ( unsigned iselected=0 ; iselected<selectChains.size() ; iselected++ ) {
 
		      selectChains[iselected] = chainName.subs( selectChains[iselected] );

		      /// replace wildcard with actual matching chains ...
		      chainNames.push_back( ChainString(selectChains[iselected]) );

		      m_provider->msg(MSG::INFO) << "[91;1m" << "Matching chain " << selectChains[iselected] << "[m" << endmsg;

		      /// if this has a cosmic chain, set the fiducial radius to be very large to 
		      /// allow the production vertex of the cosmic to be included
		      if ( selectChains[iselected].find("cosmic")!=std::string::npos ) m_fiducial_radius = 1e10; 
		     
		  }
		  
		  chainitr++;
		}
		
		m_chainNames = chainNames;

	}

	Filter_AcceptAll filter;
	/// FIXME: should really have hardcoded limits encoded as 
	///        const variables 
	Filter_etaPT       filterRef(5,500);   
	Filter_etaPT       filter_etaPT(3.0,m_ptmin);
	Filter_pdgIdpTeta  filter_pdgIdpTeta(m_TruthPdgId,3.0,m_ptmin); // |eta|<3, pt>1GeV 

	TrackFilter*        truthFilter = &filter_etaPT;

	TrigTrackSelector selectorTruth( truthFilter, m_fiducial_radius, m_TruthPdgId, m_parentTruthPdgId); 

	TrigTrackSelector selectorRef( &filter_etaPT ); 
	TrigTrackSelector selectorTest( &filter ); 

	TIDAVertexBuilder vertexBuilder;

	if ( xbeam!=0 || ybeam!=0 ) { 
	  selectorTruth.setBeamline( xbeam, ybeam, zbeam ); 
	  selectorRef.setBeamline( xbeam, ybeam, zbeam );
	}

	if ( xbeam_online!=0 || ybeam_online!=0 ) { 
	    selectorTest.setBeamline( xbeam_online, ybeam_online, zbeam_online ); 
	}

	selectorTruth.correctTracks( true );
	selectorRef.correctTracks( true );
	selectorTest.correctTracks( true );
	

	// clear the ntuple TIDA::Event class
	m_event->clear();

	const xAOD::EventInfo* pEventInfo = 0;

	unsigned           run_number         = 0;
	unsigned long long event_number       = 0;
	unsigned           lumi_block         = 0;
	unsigned           bunch_crossing_id  = 0;
	unsigned           time_stamp         = 0;
	double             mu_val             = 0;

	if ( m_provider->evtStore()->retrieve(pEventInfo).isFailure() ) {
		m_provider->msg(MSG::DEBUG) << "Failed to get EventInfo " << endmsg;
	} 
	else {
		run_number        = pEventInfo->runNumber();
		event_number      = pEventInfo->eventNumber();
		lumi_block        = pEventInfo->lumiBlock();
		time_stamp        = pEventInfo->timeStamp();
		bunch_crossing_id = pEventInfo->bcid();
		mu_val            = pEventInfo->averageInteractionsPerCrossing();
	}

	m_provider->msg(MSG::DEBUG) << "run "     << run_number 
				    << "\tevent " << event_number 
				    << "\tlb "    << lumi_block << endmsg;

	m_event->run_number(run_number);
	m_event->event_number(event_number);
	m_event->lumi_block(lumi_block);
	m_event->time_stamp(time_stamp);
	m_event->bunch_crossing_id(bunch_crossing_id);
	m_event->mu(mu_val);

	// build a chain group on the fly and use the reference
	const Trig::ChainGroup* L2chain=(*m_tdt)->getChainGroup("L2_.*");
	const Trig::ChainGroup* EFchain=(*m_tdt)->getChainGroup("EF_.*");
	const Trig::ChainGroup* HLTchain=(*m_tdt)->getChainGroup("HLT_.*");

	m_provider->msg(MSG::DEBUG) << "[91;1m" 
		<< "L2 pass " << L2chain->isPassed()  << "\t" 
		<< "EF pass " << EFchain->isPassed()  << "\t" 
		<< "HLT pass " << HLTchain->isPassed() << "[m"
		<< endmsg;
       
	/// check whether the chains have actually run, otherwise there's no point
	/// doing anything

	bool analyse = false;

	unsigned decisiontype = TrigDefs::Physics;
	
	/// bomb out if no chains passed and not told to keep all events  

	int passed_chains = 0;

	m_provider->msg(MSG::DEBUG) << "Checking " << m_chainNames.size() << " chains" << endmsg;
	
	if ( m_chainNames.empty() ) {
	  m_provider->msg(MSG::WARNING) << "No chains to check" << endmsg;
	  return;
	}

	for ( unsigned ichain=0 ; ichain<m_chainNames.size() ; ichain++ ) {
  
		std::string chainName = m_chainNames[ichain].head();

		// Only for trigger chains

		if (chainName.find("L2")  == std::string::npos && 
		    chainName.find("EF")  == std::string::npos && 
		    chainName.find("HLT") == std::string::npos ) continue;

		if ( configuredHLTChains.find(chainName)==configuredHLTChains.end() ) {
			m_provider->msg(MSG::WARNING) << "[91;1m" << "Chain " << chainName 
				<< " is not configured for this event" << "[m"
				<< endmsg;
			continue;
		}
		
		if ( m_chainNames[ichain].passed() ) decisiontype = TrigDefs::Physics;
		else                                 decisiontype = TrigDefs::includeFailedDecisions;
		
		std::string roistring = "";
		if ( m_chainNames[ichain].roi()!=""  ) roistring += "\troi " + m_chainNames[ichain].roi();  

		bool passPhysics = (*m_tdt)->isPassed(chainName); 

		m_provider->msg(MSG::DEBUG) << "Chain "  << chainName << "\troi " << roistring 
					   << "\tpres " << (*m_tdt)->getPrescale(chainName)
					   << ( passPhysics ? "[91;1m" : "" ) << "\tpass physics  " <<  passPhysics << ( passPhysics ? "[m" : "" ) 
					   << "\t: ( pass " << (*m_tdt)->isPassed(chainName, decisiontype ) << "\tdec type " << decisiontype << " ) " << endmsg;

		if ( (*m_tdt)->isPassed(chainName, decisiontype ) ||  !m_chainNames[ichain].passed() ) { 
		  analyse = true;
		  passed_chains++;
		}

	} /// finished loop over chains
 


	/// bomb out if no chains passed and not told to keep all events and found no 
	/// offline objects 
	if ( !analyse && !m_keepAllEvents && !foundOffline ) { 
	  m_provider->msg(MSG::DEBUG) << "No chains passed unprescaled - not processing this event: " << run_number << " " << event_number << " " << lumi_block << endmsg; 
	  return;
	}
	

	m_provider->msg(MSG::DEBUG) << "Chains passed " << passed_chains << endmsg;


	/// for Monte Carlo get the truth particles if requested to do so

	//  const TruthParticleContainer*  mcpartTES = 0;

	selectorTruth.clear();

	m_provider->msg(MSG::DEBUG) << "MC Truth flag " << m_mcTruth << endmsg; 
	const TrigInDetTrackTruthMap* truthMap = 0;

	if ( m_mcTruth) { 
		m_provider->msg(MSG::DEBUG) << "getting Truth" << endmsg; 
		if ( m_provider->evtStore()->retrieve(truthMap, "TrigInDetTrackTruthMap").isFailure()) {
			m_hasTruthMap = false;
		}
		else {
			m_hasTruthMap = true;
		}
		if (m_provider->evtStore()->contains<TruthParticleContainer>("INav4MomTruthEvent")) {
			//ESD
			selectTracks<TruthParticleContainer>( &selectorTruth, "INav4MomTruthEvent" );
		}
		else if (m_provider->evtStore()->contains<TruthParticleContainer>("SpclMC")) {
			/// AOD
			selectTracks<TruthParticleContainer>( &selectorTruth, "SpclMC");
		}
		else if (m_provider->evtStore()->contains<TruthParticleContainer>("")) {
			/// anything else?
			selectTracks<TruthParticleContainer>( &selectorTruth, "");
		}
		else if (m_provider->evtStore()->contains<xAOD::TruthParticleContainer>("TruthParticles")) {
			/// anything else?
		        selectTracks<xAOD::TruthParticleContainer>( &selectorTruth, "TruthParticles" );
		}
		else if (m_provider->evtStore()->contains<xAOD::TruthParticleContainer>("")) {
			/// anything else?
		        selectTracks<xAOD::TruthParticleContainer>( &selectorTruth, "" );
		}
		else { 
			m_provider->msg(MSG::WARNING) << "Truth not found - none whatsoever!" << endmsg; 
		}
	}


	// clear the ntuple TIDA::Event class
	m_event->clear();

	/// get (offline) reference tracks 

	/// get offline tracks

	m_provider->msg(MSG::DEBUG) << " Offline tracks " << endmsg;

	selectorRef.clear();

	if      (m_provider->evtStore()->contains<xAOD::TrackParticleContainer>("InDetTrackParticles")) {
	  selectTracks<xAOD::TrackParticleContainer>( &selectorRef, "InDetTrackParticles" );
	}
	else if (m_provider->evtStore()->contains<Rec::TrackParticleContainer>("TrackParticleCandidate")) {
	  selectTracks<Rec::TrackParticleContainer>( &selectorRef, "TrackParticleCandidate" );
        }
	else { 
	  m_provider->msg(MSG::WARNING) << " Offline tracks not found " << endmsg;
	}
	

	/// get the offline vertices into our structure
	for ( size_t iv=0; iv<m_vertexType.size(); iv++ ) {

		std::vector<TIDA::Vertex> vertices;
		
		std::string vertexType = "PrimaryVertices";
		std::string vertexChainname = "Vertex";
		if ( m_vertexType[iv]!="" ) { 
			vertexType = m_vertexType[iv];
			vertexChainname += ":" + vertexType;
		}

		m_provider->msg(MSG::VERBOSE) << "fetching offline AOD vertex container with key " << vertexType << endmsg;

		const xAOD::VertexContainer* xaodVtxCollection = 0;

		if ( m_provider->evtStore()->retrieve( xaodVtxCollection, vertexType ).isFailure()) {
		  if (m_provider->msg().level() <= MSG::WARNING) m_provider->msg(MSG::WARNING) << "xAOD vertex container not found with key " << vertexType <<  endmsg;
		}
		
		if ( xaodVtxCollection!=0 ) { 
		
			m_provider->msg(MSG::DEBUG) << "xAOD vertex container " << vertexType << " found with " << xaodVtxCollection->size() <<  " entries" << endmsg;
			
			// Vertex types in some secondary vertex collections are not properly set and are all 0, 
			// allow these vertices if primary vertices are not used
			if ( vertexType.find("SecVtx") != std::string::npos ) {
				vertices = vertexBuilder.select( xaodVtxCollection, &selectorRef.tracks(), true );
			}
			else {
				vertices = vertexBuilder.select( xaodVtxCollection, &selectorRef.tracks() );
			}
		}

		// now add the offline vertices
		if ( m_doOffline || m_doVertices ) { 
			m_event->addChain( vertexChainname );
			m_event->back().addRoi(TIDARoiDescriptor(true));
			m_event->back().back().addVertices( vertices );
		}	 
	}

	/// add offline Vertices to the Offline chain

	/// add the truth particles if needed
	
	if ( m_mcTruth ) {
	  m_event->addChain( "Truth" ); 
	  m_event->back().addRoi(TIDARoiDescriptor(true));
	  m_event->back().back().addTracks(selectorTruth.tracks());
	}

#if 0
	/// don't add them to the event - since now we store them in the Vertex chain ...
	for ( unsigned i=0 ; i<vertices.size() ; i++ )  { 
	  m_provider->msg(MSG::DEBUG) << "vertex " << i << " " << vertices[i] << endmsg;
	  m_event->addVertex(vertices[i]);
	}
#endif	

	/// offline object counters 		//   std::vector<TrackTrigObject> tidaVertexTracks;

	int Noff  = 0;
	int Nmu   = 0;
	int Nel   = 0;
        int Ntau  = 0;

	/// now add the offline tracks

	if ( m_doOffline ) { 
	  
	  m_event->addChain( "Offline" );
	  m_event->back().addRoi(TIDARoiDescriptor(true));
	  m_event->back().back().addTracks(selectorRef.tracks());
	
	  if ( selectorRef.getBeamX()!=0 || selectorRef.getBeamY()!=0 || selectorRef.getBeamZ()!=0 ) { 
	    std::vector<double> beamline_;
	    beamline_.push_back( selectorRef.getBeamX() );
	    beamline_.push_back( selectorRef.getBeamY() );
	    beamline_.push_back( selectorRef.getBeamZ() );
	    m_event->back().back().addUserData(beamline_);
	  }


	  Noff = selectorRef.tracks().size();
	  
	  m_provider->msg(MSG::DEBUG) << "ref tracks.size() " << selectorRef.tracks().size() << endmsg; 
	  for ( int ii=selectorRef.tracks().size() ; ii-- ; ) m_provider->msg(MSG::DEBUG) << "  ref track " << ii << " " << *selectorRef.tracks()[ii] << endmsg;  
	  
	}

	/// navigate through the requested storegate TEST chains
	for ( unsigned ichain=0 ; ichain<m_chainNames.size() ; ichain++ ) {  
	  
	  /// keep this printout here, but commented for usefull debug purposes ...
	  //	  m_provider->msg(MSG::INFO)<< "chain:\t" << m_chainNames[ichain] << endmsg;

	  /// get the chain, collection and TE names and track index
  	  std::string chainname      = m_chainNames[ichain].head();
	  std::string collectionname = m_chainNames[ichain].tail();
	  std::string vtx_name       = m_chainNames[ichain].vtx();


	  if ( chainname!="" )      continue;
	  if ( collectionname=="" ) continue;

	  chainname = collectionname;
	  if ( vtx_name!="" ) chainname += ":" + vtx_name; 
  
	  // useful debug information - leave this here

	  /// here we *only* want collections with no specified chain
	  /// name, then we look in storegate for the collections directly

	  selectorTest.clear();

	  bool found = false;
	  
	  std::string collection_test = collectionname;
	  size_t pos = collectionname.find("/");
	  if ( pos!=std::string::npos ) collection_test = collectionname.substr( pos+1, collectionname.size()-pos );

	  if (m_provider->evtStore()->contains<Rec::TrackParticleContainer>(collection_test)) {
	    found = selectTracks<Rec::TrackParticleContainer>( &selectorTest, collectionname );
	  }
	  else if (m_provider->evtStore()->contains<xAOD::TrackParticleContainer>(collection_test)) {
	    found = selectTracks<xAOD::TrackParticleContainer>( &selectorTest, collectionname );
	  }
	  else if (m_provider->evtStore()->contains<TrigInDetTrackCollection>(collection_test)) {
	    found = selectTracks<TrigInDetTrackCollection>( &selectorTest, collectionname );
	  }
	  else if (m_provider->evtStore()->contains<TrackCollection>(collection_test)) {
	    found = selectTracks<TrackCollection>( &selectorTest, collectionname );
	  }
	  else { 
	    m_provider->msg(MSG::WARNING) << "\tcollection " << collectionname << " not found" << endmsg;
	  }
	  

	  /// now retrieve any verttices for the analysis

	  std::vector<TIDA::Vertex> tidavertices;

	  m_provider->msg(MSG::DEBUG) << "\tFetch xAOD::VertexContainer with key " << vtx_name << endmsg;
	    
	  if ( vtx_name!="" ) { 
	        
	    m_provider->msg(MSG::DEBUG) << "\tFetch xAOD::VertexContainer with key " << vtx_name << endmsg;
	        
	    /// MT Vertex access
	        
	    const xAOD::VertexContainer* xaodVtxCollection = 0;
	    
	    if ( m_provider->evtStore()->retrieve( xaodVtxCollection, vtx_name ).isFailure() ) {
	      if (m_provider->msg().level() <= MSG::WARNING) m_provider->msg(MSG::WARNING) << "xAOD vertex container not found with key " << vtx_name <<  endmsg;
	    }
	    
	    if ( xaodVtxCollection!=0 ) { 
	            
	      m_provider->msg(MSG::DEBUG) << "\txAOD::VertexContainer found with size  " << xaodVtxCollection->size()
					  << "\t" << vtx_name << endmsg;

		  // Vertex types in some secondary vertex collections are not properly set and are all 0, 
	      // allow these vertices if primary vertices are not used
		  if ( vtx_name.find("SecVtx") != std::string::npos ) {
		    tidavertices = vertexBuilder.select( xaodVtxCollection, 0, true );
		  }
		  else {
		    tidavertices = vertexBuilder.select( xaodVtxCollection );
		  }		   
		}

	  }


	  if ( found ) { 
	    
	    m_event->addChain( chainname );
	    m_event->back().addRoi(TIDARoiDescriptor(true));
	    if ( vtx_name!="" ) m_event->back().back().addVertices( tidavertices );
	    m_event->back().back().addTracks(selectorTest.tracks());

	    if ( selectorTest.getBeamX()!=0 || selectorTest.getBeamY()!=0 || selectorTest.getBeamZ()!=0 ) { 
	      std::vector<double> beamline_;
	      beamline_.push_back( selectorTest.getBeamX() );
	      beamline_.push_back( selectorTest.getBeamY() );
	      beamline_.push_back( selectorTest.getBeamZ() );
	      m_event->back().back().addUserData(beamline_);
	    }
	    
	    int Ntest = selectorTest.tracks().size();
	    
	    m_provider->msg(MSG::DEBUG) << "collection " << collectionname << "\ttest tracks.size() " << Ntest << endmsg; 
	    for ( int ii=Ntest ; ii-- ; ) m_provider->msg(MSG::DEBUG) << "  test track " << ii << " " << *selectorTest.tracks()[ii] << endmsg;
	  }
	}	  


	std::string ElectronRef[7] =  { 
	  "", 
	  "TightCB", "MediumCB", "LooseCB",
	  "TightLH", "MediumLH", "LooseLH" };
 

	/// new electron selection 

	for ( size_t ielec=0 ; ielec<m_electronType.size() ; ielec++ ) {
	  /// hmm, if we stored the types as a map it would be more 
	  /// straightforward than having to stick all this in a loop

	  int itype = -1;
	  for ( int it=0 ; it<7 ; it++ ) if ( m_electronType[ielec]==ElectronRef[it] ) itype = it; 
	  if ( itype<0 ) continue;

	  std::vector<TrackTrigObject> elevec;
	  
	  int Nel_ = processElectrons( selectorRef, &elevec, itype, ( m_rawElectrons[ielec]=="raw" ? true : false ) );
	  
	  if ( Nel_ < 1 ) continue;
      
          Nel += Nel_;	

	  std::string echain = std::string("Electrons");
          if ( m_electronType[ielec]!="" )    echain += "_" + m_electronType[ielec];
	  if ( m_rawElectrons[ielec]=="raw" ) echain += "_raw";
	  
	  m_event->addChain( echain );
	  m_event->back().addRoi(TIDARoiDescriptor(true));
	  m_event->back().back().addTracks(selectorRef.tracks());
	  m_event->back().back().addObjects( elevec );

	  if ( selectorRef.getBeamX()!=0 || selectorRef.getBeamY()!=0 || selectorRef.getBeamZ()!=0 ) { 
	    std::vector<double> beamline_;
	    beamline_.push_back( selectorRef.getBeamX() );
	    beamline_.push_back( selectorRef.getBeamY() );
	    beamline_.push_back( selectorRef.getBeamZ() );
	    m_event->back().back().addUserData(beamline_);
	  }

	}
	
       

	std::string MuonRef[5] =  { "", "Tight", "Medium", "Loose", "VeryLoose" };

	/// get muons 
	for ( size_t imuon=0 ; imuon<m_muonType.size() ; imuon++ ) {
	  
	  m_provider->msg(MSG::DEBUG) << "fetching offline muons " << endmsg; 

          int muonType = -1;
          for ( int it=0 ; it<5 ; it++ ) if ( m_muonType[imuon] == MuonRef[it] ) muonType=it; 
          if ( muonType<0 ) continue; 

	  int Nmu_ = processMuons( selectorRef, muonType );

          if ( Nmu_ < 1 ) continue;

          Nmu += Nmu_;

	  m_provider->msg(MSG::DEBUG) << "found " << Nmu << " offline muons " << endmsg; 

          std::string mchain = "Muons";
          if ( m_muonType[imuon]!="" )  mchain += "_" + m_muonType[imuon];

	  m_event->addChain(mchain);
	  m_event->back().addRoi(TIDARoiDescriptor(true));
	  m_event->back().back().addTracks(selectorRef.tracks());

	  if ( selectorRef.getBeamX()!=0 || selectorRef.getBeamY()!=0 || selectorRef.getBeamZ()!=0 ) { 
	      std::vector<double> beamline_;
	      beamline_.push_back( selectorRef.getBeamX() );
	      beamline_.push_back( selectorRef.getBeamY() );
	      beamline_.push_back( selectorRef.getBeamZ() );
	      m_event->back().back().addUserData(beamline_);
	  }

	  m_provider->msg(MSG::DEBUG) << "ref muon tracks.size() " << selectorRef.tracks().size() << endmsg; 
	  for ( int ii=selectorRef.tracks().size() ; ii-- ; ) m_provider->msg(MSG::DEBUG) << "  ref muon track " << ii << " " << *selectorRef.tracks()[ii] << endmsg;  
	}
	



	/// get muons 
	if ( m_doMuonsSP ) { 
	  
	  m_provider->msg(MSG::DEBUG) << "fetching offline muons " << endmsg; 

          int muonType = 0;

	  Nmu += processMuons( selectorRef, muonType );

	  m_provider->msg(MSG::DEBUG) << "found " << Nmu << " offline muons " << endmsg; 

	  m_event->addChain("MuonsSP");
	  m_event->back().addRoi(TIDARoiDescriptor(true));
	  m_event->back().back().addTracks(selectorRef.tracks());

	  m_provider->msg(MSG::DEBUG) << "ref muon tracks.size() " << selectorRef.tracks().size() << endmsg; 
	  for ( int ii=selectorRef.tracks().size() ; ii-- ; ) m_provider->msg(MSG::DEBUG) << "  ref muon track " << ii << " " << *selectorRef.tracks()[ii] << endmsg;  
	}
	


	/// new tau selection 
	std::string TauRef[4] = { "", "Tight", "Medium", "Loose" };
	

	for ( size_t itau=0 ; itau<m_tauType.size() ; itau++ ) {
	  /// hmm, if we stored the types as a map it would be more 
	  /// straightforward than having to stick all this in a loop

	  int itype = -1;
	  for ( int it=0 ; it<4 ; it++ ) if ( m_tauType[itau]==TauRef[it] ) itype = it; 
	  if ( itype<0 ) continue;

	  /// use same threshold for 1 and 3 prong ??
	  int requireNtracks = 0;
	  if  ( m_tauProngs[itau]=="3Prong" ) requireNtracks = 3;	
	  if  ( m_tauProngs[itau]=="1Prong" ) requireNtracks = 1;

	  std::vector<TrackTrigObject> tauvec; 

	  int Ntau_ = processTaus( selectorRef, &tauvec, itype, requireNtracks, 20000 ); 

	  Ntau += Ntau_;

	  if ( Ntau_ > 0 ) { 
	    /// only add a tau collection if there are actually the 
	    /// relevant tausCH

	    std::string tchain = std::string("Taus");
	    if (   m_tauType[itau] != "" ) tchain += "_" + m_tauType[itau];
	    if ( m_tauProngs[itau] != "" ) tchain += "_" + m_tauProngs[itau];
	    
	    m_event->addChain( tchain );
	    m_event->back().addRoi(TIDARoiDescriptor(true));
	    m_event->back().back().addTracks(selectorRef.tracks());
	    m_event->back().back().addObjects( tauvec ) ; 

	    if ( selectorRef.getBeamX()!=0 || selectorRef.getBeamY()!=0 || selectorRef.getBeamZ()!=0 ) { 
	      std::vector<double> beamline_;
	      beamline_.push_back( selectorRef.getBeamX() );
	      beamline_.push_back( selectorRef.getBeamY() );
	      beamline_.push_back( selectorRef.getBeamZ() );
	      m_event->back().back().addUserData(beamline_);
	    }

	  }
	}
	
	if ( Nmu==0 && Noff==0 && Nel==0 && Ntau==0 ) m_provider->msg(MSG::DEBUG) << "No offline objects found " << endmsg;
	else foundOffline = true;


	// now loop over all relevant chains to get the trigger tracks...

	for ( unsigned ichain=0 ; ichain<m_chainNames.size() ; ichain++ ) {  

		// create chains for ntpl

		/// get the chain name
		const std::string&  chainName = m_chainNames[ichain].head();

		/// and the name of the collection (if any)    
		const std::string& collectionName = m_chainNames[ichain].tail();

		if( chainName.find("L2_")==std::string::npos && 
		    chainName.find("EF_")==std::string::npos && 
		    chainName.find("HLT_")==std::string::npos ) continue;

		if ( m_chainNames[ichain].passed() ) decisiontype = TrigDefs::Physics;
		else                                 decisiontype = TrigDefs::includeFailedDecisions;


		m_provider->msg(MSG::DEBUG) << "chain " << chainName 
					    << "\tprescale " << (*m_tdt)->getPrescale(chainName)
					    << "\tpass "     << (*m_tdt)->isPassed(chainName) << " physics " 
					    << "  (req dec " << (*m_tdt)->isPassed(chainName, decisiontype ) << " dec type " << decisiontype << ")"
					    << endmsg;
		
		/// now decide whether we want all the TEs for this chain, or just those 
		/// that are still active

		/// if the chain did not pass, skip this chain completely 
		if ( !(*m_tdt)->isPassed( chainName, decisiontype ) ) continue;

		/// new MT TDT feature access  
		

		std::string roi_key  = m_chainNames[ichain].roi();
		std::string vtx_name = m_chainNames[ichain].vtx();


#if 0
		/// this code needs to be here to be eventually replaced when 
		/// the TDT combination feature retrieval has been implemented
		/// at that point it can be replaced by the appropriate 
		/// code using the new TDT feature access


		if ( roi_name!="" ) { 

		  std::string roi_name_tmp = roi_name;
		  std::string roi_tename   = "";
		  
		  if ( roi_name.find("/")!=std::string::npos ) { 
		    roi_name_tmp = roi_name.substr( roi_name.find("/")+1, roi_name.size()-roi_name.find("/") );
		    roi_tename   = roi_name.substr( 0, roi_name.find("/") );
		  }
		  
		  roist = comb->get<TrigRoiDescriptor>( roi_name_tmp, decisiontype, roi_tename );
		  
		  if ( roist.size()>0 ) { 
		    for ( unsigned ir=0 ; ir<roist.size() ; ir++ ) m_provider->msg(MSG::DEBUG) << "\t\tRetrieved roi  " << roi_name << "\t" << *roist[ir].cptr() << endmsg; 
		  }
		  else { 
		    m_provider->msg(MSG::WARNING) << "\t\tRequested roi  " << roi_name << " not found" << endmsg; 
		  }
		  
		}
		else { 
		  roist = comb->get<TrigRoiDescriptor>("forID1"); 
		  if ( roist.empty() ) roist = comb->get<TrigRoiDescriptor>("forID"); 
		  if ( roist.empty() ) roist = comb->get<TrigRoiDescriptor>(""); 
		  if ( roist.empty() ) roist = comb->get<TrigRoiDescriptor>("initialRoI"); 
		}			  
#endif
		

		unsigned feature_type = TrigDefs::lastFeatureOfType;

		if ( roi_key!="" ) feature_type = TrigDefs::allFeaturesOfType;

		int leg = -1;

		if ( m_chainNames[ichain].element()!="" ) { 
		  leg = std::atoi(m_chainNames[ichain].element().c_str());
		}

		std::vector< TrigCompositeUtils::LinkInfo<TrigRoiDescriptorCollection> > rois = 
		  (*m_tdt)->template features<TrigRoiDescriptorCollection>( Trig::FeatureRequestDescriptor( chainName,  
													    decisiontype, 
													    roi_key, 
													    feature_type,
													    "roi", 
													    leg ) );
		int iroi = 0; /// count of how many rois processed so far

		/// if no rois for this chain then move along

		if ( rois.size()==0 ) continue;

		/// create the analysis chain

		m_event->addChain( m_chainNames[ichain] );

		TIDA::Chain& chain = m_event->back();

		/// I really, *really* hate range based loops ...
		for ( const TrigCompositeUtils::LinkInfo<TrigRoiDescriptorCollection>& roi_info : rois ) {
		  
		  iroi++;
		  
		  /// don't extract any additional rois if a superRoi is requested: 
		  /// In this case, the superRoi would be shared between the different 
		  /// chains 

		  if ( roi_key=="SuperRoi" && iroi>1 ) continue;
 
		  if ( roi_key.find("JetSuper")!=std::string::npos && iroi>1 ) continue; 
		  
		  const ElementLink<TrigRoiDescriptorCollection> roi_link = roi_info.link;

		  /// check this is not a spurious TDT match
		  if ( roi_key!="" && roi_link.dataID()!=roi_key ) continue;

		  const TrigRoiDescriptor* const* roiptr = roi_link.cptr();

		  if ( roiptr == 0 ) { 
		    //    std::cerr << "\treadback link is null  DAMMIT !!!" << std::endl;
		    continue;
		  }  

		  if (m_provider->msg().level() <= MSG::VERBOSE) {
		    m_provider->msg(MSG::VERBOSE) << " RoI descriptor for seeded chain " << chainName << " " << **roiptr << endmsg;
		  }

		  TIDARoiDescriptor* roi_tmp = new TIDARoiDescriptor( TIDARoiDescriptorBuilder(**roiptr) );		   

		  /// get the tracks 

		  /// useful diagnostic - leave in place ...
		  ///		  m_provider->msg(MSG::INFO) << "TIDARoi " << *roi_tmp << "\tcollectionName: " << collectionName << endmsg;
			      
		  /// this should *never* be the case, and we should only run this 
		  /// bit of code once the first time round the loop anyhow
			      
		  selectorTest.clear();


		  if ( chainName.find("HLT_")!=std::string::npos ) {
		    if ( selectTracks<xAOD::TrackParticleContainer>( &selectorTest, roi_link,  collectionName ) ); 
		    else {
		      if (m_provider->msg().level() <= MSG::DEBUG) {
			m_provider->msg(MSG::WARNING) << "\tNo track collection " << collectionName << " found"  << endmsg;  
		      }
		    }
		  }
		  
		  /// fetch vertices if available ...
		  
		  std::vector<TIDA::Vertex> tidavertices;	
		  
		  if ( vtx_name!="" ) { 
		    
		    m_provider->msg(MSG::DEBUG) << "\tFetch xAOD::VertexContainer for chain " << chainName << " with key " << vtx_name << endmsg;
		    
		    /// MT Vertex access
		    
		    std::pair< xAOD::VertexContainer::const_iterator, 
			       xAOD::VertexContainer::const_iterator > vtx_itrpair = getCollection<xAOD::VertexContainer>( roi_link, vtx_name );
		    
		    if ( vtx_itrpair.first == vtx_itrpair.second ) { 
		      if ( m_provider->msg().level() <= MSG::DEBUG ) {
			m_provider->msg(MSG::WARNING) << "\tNo xAOD::Vertex for chain " << chainName << " for key " << vtx_name << endmsg;
		      }
		    }
		    else {
		      
		      m_provider->msg(MSG::DEBUG) << "\txAOD::VertexContainer found with size  " << (vtx_itrpair.second - vtx_itrpair.first) 
						 << "\t" << vtx_name << endmsg;

			  // Vertex types in some secondary vertex collections are not properly set and are all 0, 
			  // allow these vertices if primary vertices are not used
			  if ( vtx_name.find("SecVtx") != std::string::npos ) {
				tidavertices = vertexBuilder.select( vtx_itrpair.first, vtx_itrpair.second, &selectorRef.tracks(), true );
			  }
			  else {
				tidavertices = vertexBuilder.select( vtx_itrpair.first, vtx_itrpair.second, &selectorTest.tracks() );
			  }		      
			}
		  }

#if 0 
		  //// not yet ready to get the jet yet - this can come once everything else is working 
		  // now get the jets if they are present
		  std::vector<TrackTrigObject> jets; 
		  if ( chainName.find("HLT_j")!=std::string::npos ) { 
		    if ( get_jets( comb, jets ) == 0 ) m_provider->msg(MSG::WARNING) << "\tjets could not be retrieved " << endmsg; 
		  }			  		    
#endif
		  
		  const std::vector<TIDA::Track*>& testTracks = selectorTest.tracks();
		  m_provider->msg(MSG::DEBUG) << "\ttest tracks.size() " << testTracks.size() << endmsg; 
		  for (unsigned int ii=0; ii < testTracks.size(); ii++) {
		    m_provider->msg(MSG::DEBUG) << "  test track " << ii << "for chain " << chainName + ":" + collectionName << " " << *testTracks[ii] << endmsg;  
		  }
		  
		  
		  // only add chain if there are any rois - also add beamline position for postprocessing
		  
		  
		  if ( roi_tmp == 0 ) { 
		    if ( testTracks.size()>0 ) m_provider->msg(MSG::WARNING) << "\ttest tracks.size() " << testTracks.size() << "found but no roi!!!" << endmsg;
		    roi_tmp = new TIDARoiDescriptor(true);
		  }
		  
		  chain.addRoi( *roi_tmp );

		  chain.back().addTracks(testTracks);
		  chain.back().addVertices(tidavertices);

#if 0
		  /// jets can't be added yet
		  if ( chainName.find("HLT_j")!=std::string::npos ) chain.back().addObjects( jets );
#endif
		  
		  if ( selectorTest.getBeamX()!=0 || selectorTest.getBeamY()!=0 || selectorTest.getBeamZ()!=0 ) { 
		    std::vector<double> beamline_;
		    beamline_.push_back( selectorTest.getBeamX() );
		    beamline_.push_back( selectorTest.getBeamY() );
		    beamline_.push_back( selectorTest.getBeamZ() );
		    chain.back().addUserData(beamline_);
		  }
	  
		  if ( roi_tmp ) delete roi_tmp;
		  roi_tmp = 0;

		}
		
	}

#if 0
	/// don't include this code yet - it is still being validated ...

	{ 
	  /// strip out the offline tracks not in any Roi ...

	  if ( filterOnRoi() || m_ptmin>0 ) { 
	    
	    TIDA::Chain* offline = 0;
	    
	    std::vector<std::string> chainnames = m_event->chainnames();
	    
	    /// get the offline chain
	    
	    for ( size_t ic=chainnames.size() ; ic-- ; ) {
	      if ( chainnames[ic] == "Offline" ) {
		offline = &(m_event->chains()[ic]);
		break;
	      }
	    }
	    
	    if ( offline ) { 
	      
	      std::vector<TIDA::Chain>& chains = m_event->chains();
	      std::vector<TIDA::Chain>::iterator citr = chains.begin();
	      
	      std::vector<std::pair<double,double> > philims;
	      
	      for ( ; citr!=chains.end() ; citr++ ) {
		if ( citr->name().find("HLT_")!=std::string::npos ) { 
		  for ( size_t ir=0 ; ir<citr->size() ; ir++ ) {
		    TIDARoiDescriptor& roi = citr->rois()[ir].roi();
		    if ( roi.composite() ) { 
		      for ( size_t isub=0 ; isub<roi.size() ; isub++ ) { 
			philims.push_back( std::pair<double,double>( roi[isub]->phiMinus(), roi[isub]->phiPlus() ) ); 
		      }
		    }
		    else philims.push_back( std::pair<double,double>( roi.phiMinus(), roi.phiPlus() ) ); 
		  }
		}
	      }
	      
	      remove_duplicates( philims );

	      for ( size_t iroi=0 ; iroi<offline->size() ; iroi++ ) {
		
		std::vector<TIDA::Track>& tracks = offline->rois()[iroi].tracks();
		
		/// may well put the reporting back in, so leaving this 
		/// this in place  
		//  size_t Noffline = tracks.size();

		for ( std::vector<TIDA::Track>::iterator it=tracks.begin() ; it<tracks.end() ; ) {
		  bool inc = true;
		  if ( m_ptmin>0 ) { 
		    if ( std::fabs(it->pT())<m_ptmin ) { inc=false; tracks.erase( it ); }
		  }
		  if ( inc && filterOnRoi() ) { 
		    bool remove_track = true;
		    for ( size_t isub=0 ; isub<philims.size() ; isub++ ) { 
		      
		      if ( philims[isub].first < philims[isub].second ) { 
			if ( it->phi()>=philims[isub].first && it->phi()<=philims[isub].second ) { 
			  remove_track = false; 
			  break;
			}
		      }
		      else  {
			if ( it->phi()>=philims[isub].first || it->phi()<=philims[isub].second ) { 
			  remove_track = false; 
			  break;
			}
		      }
		    }
		    if ( remove_track ) { inc=false; tracks.erase( it ); }
		  }
		  if ( inc ) it++;
		}
		
		/// may well put the reporting back in, so leaving this 
		/// this in place  		
		//  m_provider->msg(MSG::DEBUG) << "TIDA::Roi offline track reduction: " << Noffline << " -> " << tracks.size() << endmsg;
		
	      }
	     
	    }
	    	    
	  }
	}

#endif

	if ( mTree ) mTree->Fill();
	
}



