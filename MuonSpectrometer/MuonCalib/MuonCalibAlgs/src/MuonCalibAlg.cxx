/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#include "MuonCalibAlgs/MuonCalibAlg.h"

// gaudi
#include "xAODEventInfo/EventInfo.h"
#include "EventInfo/TagInfo.h"
#include "TrackRecord/TrackRecordCollection.h"
#include "TileEvent/TileContainer.h"
#include "TileIdentifier/TileTBID.h"
#include "LArRecEvent/LArCollisionTime.h"
#include "MuonPrepRawData/MuonPrepDataContainer.h"
#include "MuonTrigCoinData/TgcCoinDataContainer.h"
#include "MuonTrigCoinData/RpcCoinDataContainer.h"
#include "MuonDigitContainer/TgcDigit.h"
#include "MuonPrepRawData/CscStripPrepDataContainer.h"
#include "TrkEventPrimitives/ParamDefs.h"
#include "MuonReadoutGeometry/MdtReadoutElement.h"
#include "MuonSimData/MuonSimData.h"
#include "MuonSimData/MuonSimDataCollection.h"
#include "MuonSimData/CscSimDataCollection.h"
#include "GeneratorObjects/McEventCollection.h"
#include "AtlasHepMC/GenEvent.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "GeoPrimitives/GeoPrimitives.h"

// calibration 
#include "MuonCalibEvent/MuonCalibPatternCollection.h"
#include "MuonCalibNtuple/RootFileManager.h"
#include "MuonCalibIdentifier/MuonFixedId.h"
#include "MuonCalibStl/DeleteObject.h"
#include "MuonCalibEventBase/MuonCalibRawMdtHit.h"
#include "MuonCalibEventBase/MuonCalibRawRpcHit.h"
#include "MuonCalibEventBase/MuonCalibRawCscHit.h"
#include "MuonCalibEventBase/MuonCalibRawTgcHit.h"
#include "MuonCalibEventBase/MuonCalibRawTgcCoin.h"
#include "MuonCalibEventBase/MuonCalibRawRpcCoin.h"
#include "MuonCalibEventBase/MuonCalibRawRpcTriggerHit.h"
#include "MuonCalibEventBase/MuonCalibTruth.h"
#include "MuonCalibEventBase/MuonCalibMdtTruthHit.h"
#include "MuonCalibEventBase/MuonCalibRpcTruthHit.h"
#include "MuonCalibEventBase/MuonCalibTgcTruthHit.h"
#include "MuonCalibEventBase/MuonCalibCscTruthHit.h"
#include "MuonCalibEventBase/MuonCalibSegment.h"
#include "MuonCalibEventBase/MuonCalibEventInfo.h"
#include "MuonCalibEventBase/MuonCalibTriggerTimeInfo.h"
#include "MuonCalibEventBase/MuonCalibEvent.h"
#include "MuonCalibEventBase/MuonCalibTruthCollection.h"
#include "MuonCalibEventBase/MuonCalibRawHitCollection.h"
#include "MuonCalibEventBase/MuonCalibRawTriggerHitCollection.h"
#include "MuonCalibEventBase/MuonCalibPattern.h"

// trigger includes:
#include "TrigT1Result/CTP_RDO.h"
#include "TrigT1Result/CTP_Decoder.h"

#include "MuonRDO/RpcPad.h"
#include "MuonRDO/RpcPadContainer.h"
#include "MuonRDO/RpcSectorLogicContainer.h"

#include <algorithm>
#include <iostream>
#include "iomanip"
#include <map>

namespace MuonCalib { 

  MuonCalibAlg::MuonCalibAlg(const std::string& name, ISvcLocator* pSvcLocator) :
    AthAlgorithm(name, pSvcLocator),
    m_muonCalibTool("MuonCalib::MdtCalibTool/MdtCalibTool"),
    m_idToFixedIdTool("MuonCalib::IdToFixedIdTool"),
    m_stripFitter("CalibCscStripFitter/CalibCscStripFitter"),
    m_muonIdCutTool("MuonIdCutTool/MuonIdCutTool"),
    m_tileTBID(nullptr),
    m_ntupFileOpen(false)
  {
    // Contents in the ntuple: with or without truth?
    declareProperty("doTruth", m_doTruth = true);

    declareProperty("doTruthMuonOnly", m_doTruthMuonOnly = true);

    // pattern and mdtsegment location in storegate
    declareProperty("GlobalPatternLocation", m_globalPatternLocation="PatternsForCalibration");

    declareProperty("MuonCalibTool", m_muonCalibTool);

    declareProperty("EventTag", m_eventTag = "MoMu");

    declareProperty("DoRawData", m_doRawData = true);

    // Delete event at end of execute or at finalize, event needed if immediate calibration performed
    declareProperty("deleteEvent", m_doDeleteEvent = true);  

    declareProperty("addTriggerTag", m_addTriggerTag = true);

    declareProperty("MDT_TDC_cut",m_mdt_tdc_cut = true);

    // set by job options to turn on or off cuts on truth
    declareProperty("DoPrdSelect", m_doPrdSelect = false);

    //for Trigger time info
    declareProperty("MBTS_Threshold",  m_mbts_threshold = 40.0/222.0 );  // Value in pC

    // flags for ntuple production
    m_createRootFile = true;
    //declareProperty("CreateFile",m_createRootFile);
      
    declareProperty("NtupleName",m_ntupleName="CalibrationNtuple.root");

    declareProperty("doMDTs", m_doMDTs = true );
    
    declareProperty("doCSCs", m_doCSCs = true );
    declareProperty("CSCStripFitter", m_stripFitter);

    declareProperty("doRPCs", m_doRPCs = true );
    declareProperty("doTGCs", m_doTGCs = true );
    declareProperty("doRPCCoinData", m_doRPCCoin = true );
    declareProperty("doTGCCoinData", m_doTGCCoin = true );
    declareProperty("doRPCSectorLogic", m_doRpcSectorLogic = true );
    declareProperty("doRPCTriggerHits", m_doRPCTriggerHits = true );

    m_trigger_waning_printed=false;
    m_lar_waning_printed=false;
    m_mbts_warning_printed=false;
    m_eventNumber=0;
  }

  StatusCode MuonCalibAlg::initialize()
  {
 //-----------------------------------------------------
  
    ATH_MSG_INFO("Initialisation started     ");
    ATH_MSG_INFO("================================");
    ATH_MSG_INFO("= Proprieties are:");
    ATH_MSG_INFO("=   GlobalPatternLocation " << m_globalPatternLocation);
    ATH_MSG_INFO("=   MdtCalibTool          " << m_muonCalibTool);
    ATH_MSG_INFO("=   NtupleName            " << m_ntupleName   << " write " << m_createRootFile);
    if(m_doTruth) ATH_MSG_INFO("=  Truth information stored on ntuple");
    else ATH_MSG_INFO("=  Truth information NOT stored on ntuple");
    if(m_doRawData) ATH_MSG_INFO("=  Raw Hit information stored on ntuple");
    else ATH_MSG_INFO("=  Raw Hit information NOT stored on ntuple");
    ATH_MSG_INFO("================================");

    ATH_CHECK(m_idHelperSvc.retrieve());

    if (!m_stripFitter.empty()) {
      if (!m_idHelperSvc->hasCSC()) ATH_MSG_WARNING("The current layout does not have any CSC chamber but you gave a CSCStripFitter, ignoring it, but double-check configuration");
      else ATH_CHECK( m_stripFitter.retrieve() );
    }
    if (m_doCSCs && !m_idHelperSvc->hasCSC()) {
      ATH_MSG_WARNING("The current layout does not have any CSC chamber but you have set doCSCs to true, setting it to false, but double-check configuration");
      m_doCSCs = false;
    }

    m_eventNumber = 0;

    ATH_CHECK( m_idToFixedIdTool.retrieve() );

    ATH_CHECK( m_muonIdCutTool.retrieve() );

    // retrieve TileTBID helper 
    // (The MBTS was added to the Test Beam (TB) list.)
    if(detStore()->retrieve(m_tileTBID).isFailure()) {
      ATH_MSG_WARNING("Unable to retrieve TileTBID helper from DetectorStore");
    }
    
    // retrieve MuonDetectorManager from the conditions store
    ATH_CHECK(m_DetectorManagerKey.initialize());

    ATH_MSG_INFO("Initialization ended     ");
    return StatusCode::SUCCESS;
  }

  // Execute
  StatusCode MuonCalibAlg::execute(){
    if( m_createRootFile ) {
      if(!m_ntupFileOpen) {
        // create the root file: 
        // was originally in the constructor - from where it was not possible to 
        // configure the file name - please let me know if in some use-case 
        // having it here causes troubles - domizia.orestano@cern.ch
        RootFileManager::getInstance()->openFile(m_ntupleName);
        m_ntupFileOpen=true;
        ATH_CHECK(m_muonCalibTool.retrieve());
      }

      ATH_MSG_INFO("Write metadata");
      const TagInfo* tagInfo = nullptr;
      StatusCode sc = detStore()->retrieve(tagInfo);
      if (sc.isFailure() || !tagInfo) {
        ATH_MSG_FATAL("No TagInfo in DetectorStore");
        return StatusCode::FAILURE;
      }
      TagInfo::NameTagPairVec tags;
      tagInfo->getTags(tags);
      std::map<std::string, std::string> metadata;
      for (TagInfo::NameTagPairVec::const_iterator it=tags.begin(); it!=tags.end(); it++) metadata[it->first]=it->second;
      RootFileManager::getInstance()->WriteMetaData(metadata, m_ntupleName);
      m_createRootFile=false;
    }
	
    const MuonCalibEvent* event = retrieveEvent();
    m_muonCalibTool->handleEvent( event ) ;
    if (m_doDeleteEvent == true) delete event;
    else m_events.push_back(event);
    return StatusCode::SUCCESS;
  }

  StatusCode MuonCalibAlg::finalize(){
    ATH_MSG_INFO("Finalisation started     ");

    if(m_ntupFileOpen) {
      // perform analysis
      m_muonCalibTool->analyse();
      
      std::for_each(m_events.begin(),m_events.end(),DeleteObject());
      
      // close root file 
      if( m_createRootFile ) RootFileManager::getInstance()->closeFile();
    }

    return StatusCode::SUCCESS;
  }

  const MuonCalibEvent::MCPVec MuonCalibAlg::retrievePatterns() const 
  {
    MuonCalibEvent::MCPVec patternVec;
    //  retrieve MuonCalibPatternCollection
    if (!evtStore()->contains<MuonCalibPatternCollection>(m_globalPatternLocation)) {
      ATH_MSG_DEBUG("MuonCalibPatternCollection not contained in StoreGate at:"  << m_globalPatternLocation);
      return patternVec;
    }
    const MuonCalibPatternCollection* pats = 0;
    StatusCode sc = evtStore()->retrieve(pats,m_globalPatternLocation);
    if (sc.isFailure()) {
      ATH_MSG_WARNING("Could not find MuonCalibPatternCollection at " << m_globalPatternLocation);
    }else{
      ATH_MSG_DEBUG("retrieved MuonCalibPatternCollection "  << pats->size());
      patternVec.reserve(pats->size());
      
      MuonCalibPatternCollection::const_iterator pat_it     = pats->begin();
      MuonCalibPatternCollection::const_iterator pat_it_end = pats->end();
      for( ;pat_it!=pat_it_end; ++pat_it ){      
	MuonCalibPattern* pat = new MuonCalibPattern(*(*pat_it)) ;
	
	patternVec.push_back( pat );
      }    
    }
    
    return patternVec;
  }
  
  const MuonCalibTruthCollection* MuonCalibAlg::retrieveTruthCollection() const
  {  
    const TrackRecordCollection* truthCollection = nullptr;
    std::string location = "MuonEntryLayer";
    bool found = false;
    bool filled = false;
    if ( evtStore()->contains<TrackRecordCollection>(location) ) found = true ;
    if( found ) {
      if (evtStore()->retrieve(truthCollection,location).isFailure()) {
        ATH_MSG_WARNING("MuonCalibAlg::retrieveTruthCollection: could not retrieve TrackRecordCollection "<<location<<", trying MuonEntryRecord next");
      }
      if( truthCollection->size() == 0 ) filled = false ;
      else filled = true ;
    } else {
      ATH_MSG_DEBUG("location " << location << " discarded");
    }

    if(!found){
      location ="MuonEntryRecord";
      if ( evtStore()->contains<TrackRecordCollection>(location) ) found = true ;
      if( found ) {
        if (evtStore()->retrieve(truthCollection,location).isFailure()) {
          ATH_MSG_WARNING("MuonCalibAlg::retrieveTruthCollection: could not retrieve TrackRecordCollection "<<location);
        }
        if( truthCollection->size() == 0 ) filled = false ;
        else filled = true ;
      } else {
        ATH_MSG_DEBUG("location " << location << " discarded");
      }
    }

    if(found && filled){
      ATH_MSG_DEBUG("TracksRecordCollections collected at location: " << location);
    }
    if(!found){
      ATH_MSG_DEBUG("TracksRecordCollections nowhere to found in Storegate, nor filled");
    }
    if(found && !filled){
      ATH_MSG_DEBUG("TracksRecordCollections collected at location: " << location << " are empty ");
    }

    MuonCalibTruthCollection* MCtruthCollection = new MuonCalibTruthCollection();
  
    std::set < int > muonBarCode; 
    bool muonfound = false;  
    if(found && filled ){
      const McEventCollection* mcEventCollection = 0;
      std::string MClocation = "TruthEvent";
      if(evtStore()->contains<McEventCollection>(MClocation)) {
	StatusCode scmc = evtStore()->retrieve(mcEventCollection, MClocation);
	if (StatusCode::SUCCESS == scmc ){
	  ATH_MSG_DEBUG(" McEventCollection collection retrieved ");
	}

	ATH_MSG_DEBUG("retrieved TrackRecordCollection "  << truthCollection->size());
	//Convert TrackRecordCollection to MuonCalibTruthCollection
	
	TrackRecordConstIterator tr_it = truthCollection->begin();
	TrackRecordConstIterator tr_it_end = truthCollection->end();
	for(;tr_it!=tr_it_end; ++tr_it){
	  Amg::Vector3D  pos( (*tr_it).GetPosition().x(), (*tr_it).GetPosition().y(), (*tr_it).GetPosition().z() );
	  Amg::Vector3D mom( (*tr_it).GetMomentum().x(), (*tr_it).GetMomentum().y(), (*tr_it).GetMomentum().z() );
	  double kinEnergy( (*tr_it).GetEnergy() ) ;
	  int PDGCode( (*tr_it).GetPDGCode() ) ;
	  int barcode((*tr_it).GetBarCode() );	
	  double prec = std::sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2]);  
	  int newbarcode = 0; 	
	  if (std::abs(PDGCode) == 13 ) {
	    muonfound = true;  

	    // Barcode is not set correctly in TrackRecordCollection for earlier versions before June 2006
	    // Do matching with MC event collection

	    double dotprodbest=-1.;
           
	    if (StatusCode::SUCCESS == scmc ) {
	      const HepMC::GenEvent*    myGenEvent = *(mcEventCollection -> begin());
	      HepMC::GenEvent::particle_const_iterator particle = myGenEvent->particles_begin();
	      HepMC::GenEvent::particle_const_iterator particle_end = myGenEvent->particles_end();
	      for(; particle!=particle_end;++particle){
		if( (*particle)->pdg_id() ==  PDGCode ) {
		  Amg::Vector3D moms( (*particle)->momentum().x(), (*particle)->momentum().y(), (*particle)->momentum().z() );
		  double dotprod =  moms[0]*mom[0] + moms[1]*mom[1] + moms[2]*mom[2];
		  double psim = std::sqrt(moms[0]*moms[0]+moms[1]*moms[1]+moms[2]*moms[2]); 
		  psim =std::sqrt(psim);
		  dotprod = dotprod/(psim*prec);
		  if (dotprod > 0.5&&dotprod > dotprodbest) {
		    dotprodbest = dotprod;
		    newbarcode = (*particle)->barcode();
		  }
		}
	      }
	    }
	  }
	  // barcode from TrackRecord
	  MuonCalibTruth* truth = new MuonCalibTruth();
	  truth->setPosition(pos);
	  truth->setMomentum(mom);
	  truth->setKinEnergy(kinEnergy);
	  truth->setPDGCode(PDGCode); 
	  truth->setBarCode(barcode); 
	  //will destroyed in MuonCalibTruthCollection destructor
	  MCtruthCollection->addTruth(truth);
	  muonBarCode.insert(barcode);
	  // Store also newbarcode from MCEvent (negative value);
	  if (newbarcode>0 && barcode != newbarcode) {
	    truth = new MuonCalibTruth();
	    truth->setPosition(pos);
	    truth->setMomentum(mom);
	    truth->setKinEnergy(kinEnergy);
	    truth->setPDGCode(PDGCode); 
	    truth->setBarCode(-newbarcode); 
	    MCtruthCollection->addTruth(truth);
	    muonBarCode.insert(newbarcode);
	  }
	}
      } 
    }
    bool addmuonatIP = true; 
    if (!muonfound||addmuonatIP) { //retrieve truth track by other means: via McEventCollection...
      
      const McEventCollection* mcEventCollection = nullptr;
      std::string MClocation = "TruthEvent";
      if(evtStore()->contains<McEventCollection>(MClocation)) {
	StatusCode sc = evtStore()->retrieve(mcEventCollection, MClocation);

	if (StatusCode::SUCCESS == sc ) {
	  ATH_MSG_DEBUG(" McEventCollection collection retrieved ");
     
	  const HepMC::GenEvent*    myGenEvent = *(mcEventCollection -> begin());
	  HepMC::GenEvent::particle_const_iterator particle = myGenEvent->particles_begin();
	  HepMC::GenEvent::particle_const_iterator particle_end = myGenEvent->particles_end();
	  for(; particle!=particle_end;++particle){
	    if( std::abs( (*particle)->pdg_id() ) == 13) {
	      CLHEP::HepLorentzVector p((*particle)->momentum().px(),
				 (*particle)->momentum().py(),
				 (*particle)->momentum().pz(),
				 (*particle)->momentum().e());
	      Amg::Vector3D  pos(999.,999.,999.);
	      if ((*particle)->production_vertex()) { 
		pos[0] = (*particle)->production_vertex()->position().x();
		pos[1] = (*particle)->production_vertex()->position().y();
		pos[2] = (*particle)->production_vertex()->position().z();
	      }
	      Amg::Vector3D mom( (*particle)->momentum().x(), (*particle)->momentum().y(), (*particle)->momentum().z() );
	      double kinEnergy( p[3] ) ;
	      int PDGCode( (*particle)->pdg_id() ) ;
	    
	      MuonCalibTruth* truth = new MuonCalibTruth();
	    
	      truth->setPosition(pos);
	      truth->setMomentum(mom);
	      truth->setKinEnergy(kinEnergy);
	      truth->setPDGCode(PDGCode); 
	      int barcode = (*particle)->barcode();
	      truth->setBarCode(barcode);
	      MCtruthCollection->addTruth(truth);
	      muonBarCode.insert(barcode);
	    }
	  }
	} else {
	  ATH_MSG_WARNING(" McEventCollection collection NOT retrieved ");
	}
      }//retrieved other MCtruth
    }
  
    SG::ReadCondHandle<MuonGM::MuonDetectorManager> DetectorManagerHandle{m_DetectorManagerKey};
    const MuonGM::MuonDetectorManager* MuonDetMgr = DetectorManagerHandle.cptr(); 
    if(MuonDetMgr==nullptr){
      ATH_MSG_ERROR("Null pointer to the read MuonDetectorManager conditions object");
      return 0; 
    } 
  
    std::string sdoKey;
    const MuonSimDataCollection* sdoContainer = 0;
    if ( m_doMDTs ) {
      //Retrieving MDT truth hits from MDT_SDO container and store them in MuonCalibMdtTruthHits
      // and feed MuonCalibMdtTruthHits into MuonCalibTruhtCollection
      sdoKey= "MDT_SDO";
      bool addDummy = false;
      if(!evtStore()->contains<MuonSimDataCollection>(sdoKey)) {
	addDummy = true;
      }else{
	if (!evtStore()->retrieve(sdoContainer,sdoKey).isSuccess()) {
	  ATH_MSG_DEBUG("No MDT Sdo Container found");
	  addDummy = true;
	}
      }

      if( addDummy ){
	// Store one entry with zero's to recognize it offline 
	MuonCalibMdtTruthHit* mdtTruth = new MuonCalibMdtTruthHit();
	mdtTruth->setBarCode(0) ;
	mdtTruth->setDriftRadius(0.) ;
	mdtTruth->setPositionAlongTube(0.);
	MCtruthCollection->addTruth( mdtTruth );
      } else { 
	ATH_MSG_DEBUG("MDT Sdo Container found");
	MuonSimDataCollection::const_iterator mdt_it = sdoContainer->begin();
	MuonSimDataCollection::const_iterator mdt_it_end = sdoContainer->end();
	bool drop = false;
	int nmdtSDO = 0;
	int nsdo =0;
	for( ; mdt_it!=mdt_it_end; ++mdt_it ){
	  std::vector< MuonSimData::Deposit >::const_iterator deposits_it = (*mdt_it).second.getdeposits().begin() ;
	  std::vector< MuonSimData::Deposit >::const_iterator deposits_it_end = (*mdt_it).second.getdeposits().end() ;
	  for( ; deposits_it!= deposits_it_end; ++deposits_it) {
	    int barcode = (*deposits_it).first.barcode();
	    if (muonBarCode.count(barcode+10000)!=1&&muonBarCode.count(barcode)!=1&&muonBarCode.count(barcode%10000)!=1&&m_doTruthMuonOnly) continue;
	    if (nmdtSDO>999) drop = true;
	    if (!drop || (drop && (*deposits_it).first.barcode() != 0)) {
	      nmdtSDO++; 

	      const Identifier& id = (*mdt_it).first;
	      MuonFixedId mfi = m_idToFixedIdTool->idToFixedId( id );
	      if( m_doPrdSelect && m_muonIdCutTool->isCut(mfi) ) continue;

	      MuonCalibMdtTruthHit* mdtTruth = new MuonCalibMdtTruthHit();
	      mdtTruth->setBarCode( (*deposits_it).first.barcode() ) ;
	      mdtTruth->setIdentifier( mfi ) ;
	      nsdo++;
	      mdtTruth->setDriftRadius( (*deposits_it).second.firstEntry() ) ;
	      double posAlongAMDB = (*deposits_it).second.secondEntry();
	      Amg::Vector3D geoModelLocPos(0.,0.,posAlongAMDB);
	      const MuonGM::MdtReadoutElement* detEl = MuonDetMgr->getMdtReadoutElement(id);
	      if( detEl ) {
		Amg::Vector3D gpos = detEl->localToGlobalCoords(geoModelLocPos,id);
		Amg::Vector3D localAMDBPos = detEl->GlobalToAmdbLRSTransform()*gpos;
		posAlongAMDB = localAMDBPos.x();
	      }
	      mdtTruth->setPositionAlongTube( posAlongAMDB );
	      mdtTruth->identify();
	      MCtruthCollection->addTruth( mdtTruth );
	    }
	  }
	}
      }
    } // if ( m_doMDTs )

    if ( m_doRPCs ) {
      //Retrieving RPC truth hits from RPC_SDO container and store them in MuonCalibRpcTruthHits
      // and feed MuonCalibRpcTruthHits into MuonCalibTruthCollection
      sdoKey= "RPC_SDO";
      bool addDummy = false;
      if(!evtStore()->contains<MuonSimDataCollection>(sdoKey)) {
	addDummy = true;
      }else{
	if (!evtStore()->retrieve(sdoContainer,sdoKey).isSuccess()) {
	  addDummy = true;
	}
      }

      if(addDummy) {
	ATH_MSG_DEBUG("No RPC Sdo Container found");

	// Store one entry with zero's to recognize it offline 

	MuonCalibRpcTruthHit* rpcTruth = new MuonCalibRpcTruthHit();
	rpcTruth->setBarCode(0) ;
	rpcTruth->setTime(0.) ;
	MCtruthCollection->addTruth( rpcTruth );
      } else { 
	ATH_MSG_DEBUG("RPC Sdo Container found");
	MuonSimDataCollection::const_iterator rpc_it = sdoContainer->begin();
	MuonSimDataCollection::const_iterator rpc_it_end = sdoContainer->end();
	bool drop = false;
	int nrpcSDO = 0;
	for( ; rpc_it!=rpc_it_end; ++rpc_it ){
	  std::vector< MuonSimData::Deposit >::const_iterator deposits_it = (*rpc_it).second.getdeposits().begin() ;
	  std::vector< MuonSimData::Deposit >::const_iterator deposits_it_end = (*rpc_it).second.getdeposits().end() ;
	  for( ; deposits_it!= deposits_it_end; ++deposits_it) {
	    int barcode = (*deposits_it).first.barcode();
	    // Check that deposit is from truth muon 
	    if (muonBarCode.count(barcode)!=1&&muonBarCode.count(barcode%10000)!=1&&m_doTruthMuonOnly) continue;

	    if (nrpcSDO>999) drop = true;
	    if (!drop || (drop && (*deposits_it).first.barcode() != 0)) {
	      nrpcSDO++; 

	      MuonFixedId mfi = m_idToFixedIdTool->idToFixedId( (*rpc_it).first );
	      if( m_doPrdSelect && m_muonIdCutTool->isCut(mfi) ) continue;

	      MuonCalibRpcTruthHit* rpcTruth = new MuonCalibRpcTruthHit();
	      rpcTruth->setBarCode( (*deposits_it).first.barcode() ) ;
	      rpcTruth->setIdentifier( mfi ) ;
	      rpcTruth->setTime( (*deposits_it).second.firstEntry() ) ;
	      MCtruthCollection->addTruth( rpcTruth );
	    }
	  }
	}
      }
    } // if ( m_doRPCs )

    if ( m_doTGCs ) {
      //Retrieving TGC truth hits from TGC_SDO container and store them in MuonCalibTgcTruthHits
      // and feed MuonCalibTgcTruthHits into MuonCalibTruthCollection
      sdoKey= "TGC_SDO";
      bool addDummy = false;
      if(!evtStore()->contains<MuonSimDataCollection>(sdoKey)) {
	addDummy = true;
      }else{
	if (!evtStore()->retrieve(sdoContainer,sdoKey).isSuccess()) {
	  ATH_MSG_DEBUG("No TGC Sdo Container found");
	  addDummy = true;
	}
      }
      if (addDummy) {

	// Store one entry with zero's to recognize it offline 

	MuonCalibTgcTruthHit* tgcTruth = new MuonCalibTgcTruthHit();
	tgcTruth->setBarCode(0) ;
	tgcTruth->setTime(0.) ;
	MCtruthCollection->addTruth( tgcTruth );
      } else { 
	ATH_MSG_DEBUG("Sdo TGC Container found");
	MuonSimDataCollection::const_iterator tgc_it = sdoContainer->begin();
	MuonSimDataCollection::const_iterator tgc_it_end = sdoContainer->end();
	bool drop = false;
	int ntgcSDO = 0;
	for( ; tgc_it!=tgc_it_end; ++tgc_it ){
	  std::vector< MuonSimData::Deposit >::const_iterator deposits_it = (*tgc_it).second.getdeposits().begin() ;
	  std::vector< MuonSimData::Deposit >::const_iterator deposits_it_end = (*tgc_it).second.getdeposits().end() ;
	  for( ; deposits_it!= deposits_it_end; ++deposits_it) {
	    int barcode = (*deposits_it).first.barcode();
	    // Check that deposit is from truth muon 
	    if (muonBarCode.count(barcode+10000)!=1&&muonBarCode.count(barcode)!=1&&muonBarCode.count(barcode%10000)!=1&&m_doTruthMuonOnly) continue;
	    if (ntgcSDO>999) drop = true;
	    if (!drop || (drop && (*deposits_it).first.barcode() != 0)) {

	      MuonFixedId mfi =  m_idToFixedIdTool->idToFixedId( (*tgc_it).first );
	      if( m_doPrdSelect && m_muonIdCutTool->isCut(mfi) ) continue;

	      ntgcSDO++; 
	      MuonCalibTgcTruthHit* tgcTruth = new MuonCalibTgcTruthHit();
	      tgcTruth->setBarCode( (*deposits_it).first.barcode() ) ;
	      tgcTruth->setIdentifier( mfi ) ;
	      tgcTruth->setTime( (*deposits_it).second.firstEntry() ) ;
	      MCtruthCollection->addTruth( tgcTruth );
	    }
	  }
	}
      }
    } // if ( m_doTGCs )


    if ( m_doCSCs ) {
      //Retrieving CSC truth hits from CSC_SDO container and store them in MuonCalibCscTruthHits
      // and feed MuonCalibCscTruthHits into MuonCalibTruthCollection
      sdoKey= "CSC_SDO";
      if (evtStore()->contains<MuonSimDataCollection>(sdoKey))
	{
	  if (evtStore()->retrieve(sdoContainer,sdoKey).isSuccess()) {
	    ATH_MSG_DEBUG("CSC Sdo Container found");
	    MuonSimDataCollection::const_iterator csc_it = sdoContainer->begin();
	    MuonSimDataCollection::const_iterator csc_it_end = sdoContainer->end();
	    bool drop = false;
	    int ncscSDO = 0;
	    for( ; csc_it!=csc_it_end; ++csc_it ){
	      std::vector< MuonSimData::Deposit >::const_iterator deposits_it = (*csc_it).second.getdeposits().begin() ;
	      std::vector< MuonSimData::Deposit >::const_iterator deposits_it_end = (*csc_it).second.getdeposits().end() ;
	      for( ; deposits_it!= deposits_it_end; ++deposits_it) {
		int barcode = (*deposits_it).first.barcode();
		// Check that deposit is from truth muon 
		if (muonBarCode.count(barcode+10000)!=1&&muonBarCode.count(barcode)!=1&&muonBarCode.count(barcode%10000)!=1&&m_doTruthMuonOnly) continue;
		if (ncscSDO>999) drop = true;
		if (!drop || (drop && (*deposits_it).first.barcode() != 0)) {
		  ncscSDO++; 
		  MuonFixedId mfi = m_idToFixedIdTool->idToFixedId( (*csc_it).first );
		  if( m_doPrdSelect && m_muonIdCutTool->isCut(mfi) ) continue;

		  MuonCalibCscTruthHit* cscTruth = new MuonCalibCscTruthHit();
		  cscTruth->setBarCode( (*deposits_it).first.barcode() ) ;
		  cscTruth->setIdentifier( mfi ) ;
		  cscTruth->setTime( (*deposits_it).second.firstEntry() ) ;
		  MCtruthCollection->addTruth( cscTruth );
		}
	      }
	    }
	  }
	} else { 
	  ATH_MSG_DEBUG("No CSC MuonSimHit Sdo Container found");

	  const CscSimDataCollection* cscSdoContainer;
	  bool addDummy = false;
	  if(!evtStore()->contains<CscSimDataCollection>(sdoKey)) {
	    addDummy = true;
	  }else{
	    if (!evtStore()->retrieve(cscSdoContainer,sdoKey).isSuccess()) {
	      addDummy = true;
	    }
	  }
	  if (addDummy) {
	    ATH_MSG_DEBUG("No CSC Sdo Container found");

	    // Store one entry with zero's to recognize it offline 

	    MuonCalibCscTruthHit* cscTruth = new MuonCalibCscTruthHit();
	    cscTruth->setBarCode(0) ;
	    cscTruth->setTime(0.) ;
	    MCtruthCollection->addTruth( cscTruth );
	  } else { 
	    ATH_MSG_DEBUG("CSC Sdo Container found");
	    CscSimDataCollection::const_iterator csc_it = cscSdoContainer->begin();
	    CscSimDataCollection::const_iterator csc_it_end = cscSdoContainer->end();
	    bool drop = false;
	    int ncscSDO = 0;
	    for( ; csc_it!=csc_it_end; ++csc_it ){
	      std::vector< CscSimData::Deposit >::const_iterator deposits_it = (*csc_it).second.getdeposits().begin() ;
	      std::vector< CscSimData::Deposit >::const_iterator deposits_it_end = (*csc_it).second.getdeposits().end() ;
	      for( ; deposits_it!= deposits_it_end; ++deposits_it) {
		int barcode = (*deposits_it).first.barcode();
		// Check that deposit is from truth muon 
		if (muonBarCode.count(barcode)!=1&&muonBarCode.count(barcode%10000)!=1&&m_doTruthMuonOnly) continue;
		if (ncscSDO>999) drop = true;
		if (!drop || (drop && (*deposits_it).first.barcode() != 0)) {
		  ncscSDO++; 
		  Identifier id = (*csc_it).first; 
		  const MuonGM::CscReadoutElement* detEl = MuonDetMgr->getCscReadoutElement(id);
		  if( !detEl ){
		    ATH_MSG_WARNING("Found CSC Identifier which seems to have no readout element " 
			<< m_idHelperSvc->mdtIdHelper().print_to_string(id));
		    continue;
		  }
		  MuonFixedId mfi = m_idToFixedIdTool->idToFixedId(id);
		  if( m_doPrdSelect && m_muonIdCutTool->isCut(mfi) ) continue;

		  MuonCalibCscTruthHit* cscTruth = new MuonCalibCscTruthHit();
		  cscTruth->setBarCode( (*deposits_it).first.barcode() ) ;
		  cscTruth->setIdentifier( mfi ) ;
		  double ypos_local = 	 (*deposits_it).second.ypos() ;
		  double zpos_local = 	 (*deposits_it).second.zpos() ;
		  const Amg::Vector3D locpos3d = Amg::Vector3D(0.,ypos_local, zpos_local);
		  const Amg::Vector3D globalPos = (detEl->transform())*locpos3d;
		  double phipos = std::atan2(globalPos.y(),globalPos.x());
		  // use cscTime for phi, since time is no longer set for new SDO format!!
		  cscTruth->setTime( phipos ) ;
		  MCtruthCollection->addTruth( cscTruth );
		}
	      }
	    }
	  }
	}
    } // m_doCSCs

     ATH_MSG_DEBUG("retrieved MCtruthCollection size " << MCtruthCollection->numberOfTruth());
  
    return MCtruthCollection;
  }

  MuonCalibEventInfo MuonCalibAlg::retrieveEventInfo() const
    { 
      const xAOD::EventInfo* eventInfo=nullptr;
      MuonCalibEventInfo MCeventInfo;
      StatusCode sc = evtStore()->retrieve(eventInfo);
      if ( sc.isFailure() ) {
	ATH_MSG_ERROR("Could not find eventInfo ");
	return MCeventInfo;
      }else{
	ATH_MSG_DEBUG("retrieved eventInfo");
      }

      //Cast eventID into MuonCalibEventInfo class:
   
      MCeventInfo.setRunNumber( eventInfo->runNumber() ) ;
      MCeventInfo.setEventNumber( eventInfo->eventNumber() );
      MCeventInfo.setTimeStamp( eventInfo->timeStamp() );
      MCeventInfo.setLumiBlock( eventInfo->lumiBlock() );
      MCeventInfo.setBcId( eventInfo->bcid() );
      std::string eventTag=m_eventTag;
      MCeventInfo.setTag( eventTag );
      if ( m_addTriggerTag ) {
	const DataHandle< CTP_RDO > ctpRDO;
	if ( ! evtStore()->contains(ctpRDO, "CTP_RDO") )
                {
		if(!m_trigger_waning_printed)
		{
      ATH_MSG_INFO("No trigger info, not added to EventTag");
		m_trigger_waning_printed=true;
		}
		else
		{
		ATH_MSG_DEBUG("No trigger info, not added to EventTag");
		}
		
		return MCeventInfo;
		}
	if ( evtStore()->retrieve( ctpRDO, "CTP_RDO" ).isFailure() ) {
	if(!m_trigger_waning_printed)
		{
      ATH_MSG_INFO("CTP_RDO trigger info missing, not added to EventTag");
		m_trigger_waning_printed=true;
		}
	else
		{
		ATH_MSG_DEBUG("CTP_RDO trigger info missing, not added to EventTag");
		}
	return MCeventInfo;
	}
	CTP_Decoder ctp;
	ctp.setRDO(ctpRDO);
	uint16_t l1aPos = ctpRDO->getL1AcceptBunchPosition();
	if(l1aPos >= ctp.getBunchCrossings().size()) return MCeventInfo;
	const CTP_BC& bunch = ctp.getBunchCrossings().at(l1aPos);
	MCeventInfo.setNumberOfTriggerBits(bunch.getTAV().size() + bunch.getTAP().size());
	int nth_bit(0);
	for(unsigned int i=0; i<bunch.getTAV().size(); i++)
	  {
	    MCeventInfo.setTriggerBit(nth_bit, bunch.getTAV().test(i));
	    nth_bit++;
	  }
	for(unsigned int i=0; i<bunch.getTAP().size(); i++)
	  {
	    MCeventInfo.setTriggerBit(nth_bit, bunch.getTAP().test(i));
	    nth_bit++;
	  }
      }
      return MCeventInfo;
    }

  ///
  const MuonCalibTriggerTimeInfo* MuonCalibAlg::retrieveTriggerTimeInfo() const
  { 
    MuonCalibTriggerTimeInfo trigTimeInfo;
    ///////////////////////////////////////////////////////
    // MBTS word - copied from: 
    // http://alxr.usatlas.bnl.gov/lxr/source/atlas/TriggerAlgorithms/TrigT2MinBias/src/T2MbtsFex.cxx
    const DataHandle<TileCellContainer> tileCellCnt;
    uint32_t mbtsBits = 0;
    float timeA = 0.F;
    float timeC = 0.F;
    int countA = 0;
    int countC = 0;
    
    if(evtStore()->contains(tileCellCnt, "MBTSContainer"))
      {
      StatusCode sc=evtStore()->retrieve(tileCellCnt, "MBTSContainer");
      if (!sc.isFailure()) {
      
	uint32_t ibit, bit_pos = 0;
	float charge = 0;
	// Discriminate the signals
	TileCellContainer::const_iterator itr = tileCellCnt->begin();
	TileCellContainer::const_iterator itr_end = tileCellCnt->end();
	for(; itr != itr_end; ++itr) {
	  charge = (*itr)->energy();
	  ATH_MSG_DEBUG("Energy =" << charge << "pC");
	  if(charge > m_mbts_threshold) {
	    Identifier id=(*itr)->ID();
	    // cache type, module and channel
	    // MBTS Id type is  "side"  +/- 1
	    int type_id = m_tileTBID->type(id);
	    // MBTS Id module is "phi"  0-7
	    int module_id = m_tileTBID->module(id);
	    // MBTS Id channel is "eta"  0-1   zero is closer to beam pipe
	    int channel_id = m_tileTBID->channel(id);

	    if (type_id > 0) {
	      timeA +=  (*itr)->time();
	      countA++ ; 
	    } else {
	      timeC +=  (*itr)->time();
	      countC++;
	    }

	    // Catch errors
	    if( std::abs(type_id) != 1 ){
	      ATH_MSG_WARNING("MBTS identifier type is out of range");
	      continue;
	    }
	    if( channel_id < 0 || channel_id > 1 ){
	      ATH_MSG_WARNING("MBTS identifier channel is out of range");
	      continue;
	    }
	    if( module_id < 0 || module_id > 7 ){
	      ATH_MSG_WARNING("MBTS identifier module is out of range");
	      continue;
	    }      
	    bit_pos = 0; // The position of the bit
	    if(type_id == -1) {
	      bit_pos += 16;
	    }
	    bit_pos += channel_id*8;
	    bit_pos += module_id;
	    ibit = 1; // The mbts bit for this TileCell
	    ibit <<= bit_pos;
	    mbtsBits += ibit;
	  }
	}
      }}
     else
     	{
	if(!m_mbts_warning_printed)
		{
		ATH_MSG_INFO("No MBTS info in store gate");
		m_mbts_warning_printed=true;
		}
	}

    float timeDiff = -999.; // Or whatever default you want
    if (countA > 1 && countC > 1) timeDiff = (timeA/countA) - (timeC/countC);

    //LAr time diff
    float LArECtimeDiff=-999.;
//    LArCollisionTime* tps;
    const DataHandle< LArCollisionTime > tps;
    if ( evtStore()->contains(tps, "LArCollisionTime") )
    	{
	StatusCode sc = evtStore()->retrieve(tps,"LArCollisionTime");
	if (!sc.isFailure()) {
	const       int nMin=2;
	if (tps->ncellA() > nMin && tps->ncellC() > nMin)       LArECtimeDiff =   tps->timeA()-tps->timeC();
	}	
	}
    else
    	{
	if(!m_lar_waning_printed)
		{
		ATH_MSG_INFO("LArCollisionTime not in store gate.");
		m_lar_waning_printed=true;
		}
	}
    //Cast eventID into MuonCalibEventInfo class:
    trigTimeInfo.setdMbtsTime( timeDiff ) ;
    trigTimeInfo.setdLarTime( LArECtimeDiff ) ;
    
    MuonCalibTriggerTimeInfo* ttInfo = new MuonCalibTriggerTimeInfo(trigTimeInfo);

    return ttInfo;
  }
  ////
  
  const MuonCalibRawHitCollection* MuonCalibAlg::retrieveRawHits( const MuonCalibEvent::MCPVec &patterns ) const{    
    MuonCalibRawHitCollection* rawHits = new MuonCalibRawHitCollection();
    
    //First, construct maps of the hits on segment. They are ordered by MuonFixedId,
    //keeping track on how many times they're assigned to a segment.
    std::map<MuonFixedId,int> mdtMap, rpcMap, tgcMap, cscMap ;
    
    MuonCalibEvent::MCPVecCit pat_it     = patterns.begin();
    MuonCalibEvent::MCPVecCit pat_it_end = patterns.end(); 
    for( ;pat_it!=pat_it_end; ++pat_it ){
      MuonCalibPattern::MuonSegCit seg_it     = (*pat_it)->muonSegBegin();
      MuonCalibPattern::MuonSegCit seg_it_end = (*pat_it)->muonSegEnd();
      for( ;seg_it!=seg_it_end;++seg_it ){ 
	
	//Filling the mdtMap with the hits
	MuonCalibSegment::MdtHitCit mdt_it = (*seg_it)->mdtHOTBegin();
	MuonCalibSegment::MdtHitCit mdt_it_end = (*seg_it)->mdtHOTEnd();
	for( ;mdt_it!=mdt_it_end;++mdt_it){ 
	  std::map<MuonFixedId, int>::iterator position = mdtMap.find( (*mdt_it)->identify() );
	  if(position == mdtMap.end() ){
	    mdtMap[ (*mdt_it)->identify() ] = 1;
	  } else {
	    ++mdtMap[ (*mdt_it)->identify() ] ; 
	  }	    
	}
	
	//Filling the rpcMap with the hits
	MuonCalibSegment::RpcHitCit rpc_it = (*seg_it)->rpcHOTBegin();
	MuonCalibSegment::RpcHitCit rpc_it_end = (*seg_it)->rpcHOTEnd();
	for( ;rpc_it!=rpc_it_end;++rpc_it){ 
	  std::map<MuonFixedId, int>::iterator position = rpcMap.find( (*rpc_it)->identify() );
	  if(position == rpcMap.end() ){
	    rpcMap[ (*rpc_it)->identify() ] = 1;
	  } else {
	    ++rpcMap[ (*rpc_it)->identify() ] ; 
	  }	    
	}
	
	//Filling the tgcMap with the hits
	MuonCalibSegment::TgcHitCit tgc_it = (*seg_it)->tgcHOTBegin();
	MuonCalibSegment::TgcHitCit tgc_it_end = (*seg_it)->tgcHOTEnd();
	for( ;tgc_it!=tgc_it_end;++tgc_it){
	  std::map<MuonFixedId, int>::iterator position = tgcMap.find( (*tgc_it)->identify() );
	  if(position == tgcMap.end() ){
	    tgcMap[ (*tgc_it)->identify() ] = 1;
	  } else {
	    ++tgcMap[ (*tgc_it)->identify() ] ; 
	  }	    
	}
	
	//Filling the cscMap with the hits
	MuonCalibSegment::CscHitCit csc_it = (*seg_it)->cscHOTBegin();
	MuonCalibSegment::CscHitCit csc_it_end = (*seg_it)->cscHOTEnd();
	for( ;csc_it!=csc_it_end;++csc_it){ 
	  std::map<MuonFixedId, int>::iterator position = cscMap.find( (*csc_it)->identify() );
	  if(position == cscMap.end() ){
	    cscMap[ (*csc_it)->identify() ] = 1;
	  } else {
	    ++cscMap[ (*csc_it)->identify() ] ; 
	  }	    
	}
      }
    }
    
    //Now, get PRDs from storegate, cast them into MuonCalibRawHits and use the maps
    //constructed above to decide which 'occupancy' flag they hold.
    
    StatusCode sc;
    sc.ignore();
    if ( m_doMDTs ) {
      //MDT raw hits....
      
      std::vector<const Muon::MdtPrepDataCollection*> mdtCols;
      const Muon::MdtPrepDataContainer* mdtPrds = 0;      
      if( !evtStore()->contains<Muon::MdtPrepDataContainer>("MDT_DriftCircles")){
	ATH_MSG_DEBUG("MdtPrepDataContainer MDT_DriftCircles not contained in SG");
      }else if( evtStore()->retrieve(mdtPrds,"MDT_DriftCircles").isFailure() ){
	ATH_MSG_DEBUG("Cannot retrieve MdtPrepDataContainer MDT_DriftCircles");
      }else{
	
	ATH_MSG_DEBUG("Retrieved MdtPrepDataContainer ");
	
	mdtCols.reserve(mdtPrds->size()); // number of mdt chambers	
	Muon::MdtPrepDataContainer::const_iterator it = mdtPrds->begin();
	Muon::MdtPrepDataContainer::const_iterator it_end = mdtPrds->end();
	for( ; it!=it_end; ++it ) {
	  mdtCols.push_back( *it );
	}
      }
      
      std::vector< const Muon::MdtPrepDataCollection*>::const_iterator mdtCollection = mdtCols.begin();
      std::vector< const Muon::MdtPrepDataCollection*>::const_iterator mdtlastColl = mdtCols.end();
      for ( ; mdtCollection != mdtlastColl ; ++mdtCollection ) 
	{
	  Amg::Transform3D globalToStation;
	  
	  Muon::MdtPrepDataCollection::const_iterator mdt_it = (*mdtCollection)->begin();
	  Muon::MdtPrepDataCollection::const_iterator mdt_it_end = (*mdtCollection)->end();
	  bool amdbtransform = false;
	  for( ; mdt_it!=mdt_it_end; ++ mdt_it)
	    {
	      if (m_mdt_tdc_cut == true && (*mdt_it)->status()==Muon::MdtStatusMasked) continue;
	      if((*mdt_it)->localPosition()[Trk::locR]==0) continue;
	      
	      MuonCalibRawMdtHit* rawMdtHit = new MuonCalibRawMdtHit();     
	      MuonFixedId fID = m_idToFixedIdTool->idToFixedId( (*mdt_it)->identify() ) ;
	      
	      if ( amdbtransform == false ) { // shouldn't this be done all the time? since there are 2 detector elements (one for each ml) per collection!
		const MuonGM::MdtReadoutElement* detEl = (*mdt_it)->detectorElement() ;
		globalToStation = detEl->GlobalToAmdbLRSTransform();
		amdbtransform = true;
	      }
	      
	      Amg::Vector3D tubePos = (*mdt_it)->globalPosition();
	      Amg::Vector3D tubePosLoc = globalToStation*tubePos;
	      
	      rawMdtHit->setId( fID );
	      rawMdtHit->setLocalPosition( tubePosLoc );
	      rawMdtHit->setGlobalPosition( tubePos );
	      rawMdtHit->setAdc( (*mdt_it)->adc() );
	      rawMdtHit->setTdc( (*mdt_it)->tdc() );
	      rawMdtHit->setDriftTime( 0. ); 
	      rawMdtHit->setDriftRadius( (*mdt_it)->localPosition()[Trk::locR] );
	      rawMdtHit->setDriftRadiusError( 1./std::sqrt((*mdt_it)->localCovariance()(Trk::locR,Trk::locR)) );
	      
	      int occupancy = 0;
	      std::map<MuonFixedId, int>::const_iterator position = mdtMap.find( fID );
	      
	      if( position == mdtMap.end() )
		occupancy = -1; //ID of RawHit not found on any segment
	      else
		occupancy = position->second;             //ID of RawHit found n times on a segment in this event
	      rawMdtHit->setOccupancy(occupancy);
	      
	      rawHits->addMuonCalibRawHit( rawMdtHit );
	    }
	}
    } // if ( m_doMDTs )
    
    if ( m_doRPCs ) {
      //RPC raw hits...
      
      std::vector<const Muon::RpcPrepDataCollection*> rpcCols;
      const Muon::RpcPrepDataContainer* rpcPrds = 0;      
      if( !evtStore()->contains<Muon::RpcPrepDataContainer>("RPC_Measurements") ){
	ATH_MSG_DEBUG("RpcPrepDataContainer RPC_Measurements not contained in SG");
      }else if(evtStore()->retrieve(rpcPrds,"RPC_Measurements").isFailure()) {
	ATH_MSG_DEBUG("Cannot retrieve RpcPrepDataContainer RPC_Measurements");
      }else{
	ATH_MSG_DEBUG("Retrieved RpcPrepDataContainer ");
	
	rpcCols.reserve(rpcPrds->size()); // number of rpc chambers	
	Muon::RpcPrepDataContainer::const_iterator it = rpcPrds->begin();
	Muon::RpcPrepDataContainer::const_iterator it_end = rpcPrds->end();
	for( ; it!=it_end; ++it ) {
	  rpcCols.push_back( *it );
	}
      }
      
      std::vector< const Muon::RpcPrepDataCollection*>::const_iterator rpcCollection = rpcCols.begin();
      std::vector< const Muon::RpcPrepDataCollection*>::const_iterator rpclastColl = rpcCols.end();
      for ( ; rpcCollection != rpclastColl ; ++rpcCollection ) 
	{
	  Muon::RpcPrepDataCollection::const_iterator rpc_it = (*rpcCollection)->begin();
	  Muon::RpcPrepDataCollection::const_iterator rpc_it_end = (*rpcCollection)->end();
	  for( ; rpc_it!=rpc_it_end; ++ rpc_it){
	    
	    MuonCalibRawRpcHit* rawRpcHit = createRawRpcHit(**rpc_it);
	    
	    int occupancy = 0;
	    std::map<MuonFixedId, int>::const_iterator position = rpcMap.find( rawRpcHit->identify() );
	    if( position == rpcMap.end() ) occupancy = -1; //ID of RawHit not found on any segment
	    else occupancy = position->second;             //ID of RawHit found n times on a segment in this event
	    rawRpcHit->setOccupancy(occupancy);
	    
	    rawHits->addMuonCalibRawHit( rawRpcHit );
	  }
	}
    } // if ( m_doRPCs )

    if ( m_doCSCs ) {
      //CSC raw hits...
      
      std::vector<const Muon::CscStripPrepDataCollection*> cscCols;
      const Muon::CscStripPrepDataContainer* cscPrds = 0;      
      if( !evtStore()->contains<Muon::CscStripPrepDataContainer>("CSC_Measurements") ){
	ATH_MSG_DEBUG("CscPrepDataContainer CSC_Measurements not contained in SG");
      }else if(evtStore()->retrieve(cscPrds,"CSC_Measurements").isFailure()) {
	ATH_MSG_DEBUG("Cannot retrieve CscPrepDataContainer CSC_Measurements");
      }else{
	ATH_MSG_DEBUG("Retrieved CscPrepDataContainer ");
	
	cscCols.reserve(cscPrds->size());
	Muon::CscStripPrepDataContainer::const_iterator it = cscPrds->begin();
	Muon::CscStripPrepDataContainer::const_iterator it_end = cscPrds->end();
	for( ; it!=it_end; ++it ) {
	  cscCols.push_back( *it );
	}
      }
      
      std::vector< const Muon::CscStripPrepDataCollection*>::const_iterator cscCollection = cscCols.begin();
      std::vector< const Muon::CscStripPrepDataCollection*>::const_iterator csclastColl = cscCols.end();
      
      for ( ; cscCollection != csclastColl ; ++cscCollection ) 
	{
	  Muon::CscStripPrepDataCollection::const_iterator csc_it = (*cscCollection)->begin();
	  Muon::CscStripPrepDataCollection::const_iterator csc_it_end = (*cscCollection)->end();
	  for( ; csc_it!=csc_it_end; ++ csc_it){
	    
	    MuonCalibRawCscHit* rawCscHit = new MuonCalibRawCscHit(); 
	    MuonFixedId fID = m_idToFixedIdTool->idToFixedId( (*csc_it)->identify() ) ;    
	    rawCscHit->setGlobalPosition( (*csc_it)->globalPosition() );
	    rawCscHit->setId( fID );
        // time of first sample is not very useful. we need fitted time
        // invoke fitter first and then set the time
	    int measuresPhi    = m_idHelperSvc->cscIdHelper().measuresPhi((*csc_it)->identify());
	    int chamberLayer   = m_idHelperSvc->cscIdHelper().chamberLayer((*csc_it)->identify());
	    float stripWidth   = (*csc_it)->detectorElement()->cathodeReadoutPitch( chamberLayer, measuresPhi );
	    rawCscHit->setWidth( stripWidth );
        // invoke the strip fitter to fit the time samples (which is a vector of 4 elements)
        // fitter outputs the sample charge, time etc. which we can store here        
        ICscStripFitter::Result fitres; // fit result will be here
        fitres = m_stripFitter->fit(**csc_it);
	    rawCscHit->setCharge(fitres.charge);
        rawCscHit->setT(fitres.time);
	    
	    int occupancy = 0;
	    std::map<MuonFixedId, int>::const_iterator position = cscMap.find( fID );
	    if( position == cscMap.end() ) occupancy = -1; //ID of RawHit not found on any segment
	    else occupancy = position->second;             //ID of RawHit found n times on a segment in this event
	    rawCscHit->setOccupancy(occupancy);
	    
	    rawHits->addMuonCalibRawHit( rawCscHit );
	  }
	}
    } // if ( m_doCSCs )
    
    if ( m_doTGCs ) {
      //TGC raw hits...
      
      std::vector<const Muon::TgcPrepDataCollection*> tgcCols;
      const Muon::TgcPrepDataContainer* tgcPrds[3] = { 0, 0, 0 };
      
      for(int ibc=0; ibc<3; ibc++){
	int bcTag=ibc+1;
	std::ostringstream location;
	location << "TGC_Measurements" << (bcTag==TgcDigit::BC_PREVIOUS ? "PriorBC" : "")
		 << (bcTag==TgcDigit::BC_NEXT ? "NextBC" : "");	
	
	if( !evtStore()->contains<Muon::TgcPrepDataContainer>(location.str()) ){
	  ATH_MSG_DEBUG("Cannot retrieve TgcPrepDataContainer " << location.str());
	}else if(evtStore()->retrieve(tgcPrds[ibc],location.str()).isFailure()) {
	  ATH_MSG_DEBUG("Cannot retrieve TgcPrepDataContainer "  << location.str());
	}else{
	  ATH_MSG_DEBUG("Retrieved TgcPrepDataContainer "  << location.str());
	  
	  tgcCols.clear();
	  tgcCols.reserve(tgcPrds[ibc]->size()); // number of tgc chambers
	  Muon::TgcPrepDataContainer::const_iterator it = tgcPrds[ibc]->begin();
	  Muon::TgcPrepDataContainer::const_iterator it_end = tgcPrds[ibc]->end();
	  for( ; it!=it_end; ++it ) {
	    tgcCols.push_back( *it );
	  }
	}
	
	std::vector< const Muon::TgcPrepDataCollection*>::const_iterator tgcCollection = tgcCols.begin();
	std::vector< const Muon::TgcPrepDataCollection*>::const_iterator tgclastColl = tgcCols.end();
	for ( ; tgcCollection != tgclastColl ; ++tgcCollection ) 
	  {
	    Muon::TgcPrepDataCollection::const_iterator tgc_it = (*tgcCollection)->begin();
	    Muon::TgcPrepDataCollection::const_iterator tgc_it_end = (*tgcCollection)->end();
	    for( ; tgc_it!=tgc_it_end; ++ tgc_it){
	      MuonCalibRawTgcHit* rawTgcHit = new MuonCalibRawTgcHit();      
	      MuonFixedId fID = m_idToFixedIdTool->idToFixedId( (*tgc_it)->identify() ) ; 
	      bool measuresPhi = (bool) m_idHelperSvc->tgcIdHelper().isStrip((*tgc_it)->identify());
	      int gasGap = m_idHelperSvc->tgcIdHelper().gasGap((*tgc_it)->identify());
	      int channel = m_idHelperSvc->tgcIdHelper().channel((*tgc_it)->identify());
	      
	      const MuonGM::TgcReadoutElement* detEl = (*tgc_it)->detectorElement();
	      double width       = -999.;
	      double short_width = -999.;
	      double long_width  = -999.;
	      double length      = -999.;
	      if(measuresPhi){
		width = detEl->stripWidth(gasGap,channel);
		short_width = detEl->stripShortWidth(gasGap, channel);
		long_width = detEl->stripLongWidth(gasGap, channel);
		length = detEl->stripLength(gasGap, channel);
	      }else{   //measures Eta: gangs, not strips
		width = detEl->gangLength(gasGap,channel);
		short_width = detEl->gangShortWidth(gasGap, channel);
		long_width = detEl->gangLongWidth(gasGap, channel);
		length = detEl->gangThickness(gasGap, channel);
	      }
	      
	      rawTgcHit->setId( fID );
	      rawTgcHit->setGlobalPosition( (*tgc_it)->globalPosition() );
	      rawTgcHit->setStation( m_idHelperSvc->tgcIdHelper().stationName((*tgc_it)->identify()) );
	      rawTgcHit->setEta( m_idHelperSvc->tgcIdHelper().stationEta((*tgc_it)->identify()) );
	      rawTgcHit->setPhi( m_idHelperSvc->tgcIdHelper().stationPhi((*tgc_it)->identify()) );
	      rawTgcHit->setGasGap( gasGap );
	      rawTgcHit->setIsStrip( static_cast<int>(measuresPhi) ); 
	      rawTgcHit->setChannel( channel );
	      rawTgcHit->setBcTag( bcTag );
	      
	      rawTgcHit->setWidth( width );
	      rawTgcHit->setShortWidth( short_width );
	      rawTgcHit->setLongWidth( long_width );
	      rawTgcHit->setLength( length );
	      
	      int occupancy = 0;
	      std::map<MuonFixedId, int>::const_iterator position = tgcMap.find( fID );
	      if( position == tgcMap.end() ) occupancy = -1; //ID of RawHit not found on any segment
	      else occupancy = position->second;          //ID of RawHit found n times on a segment in this event
	      rawTgcHit->setOccupancy(occupancy);//Will be filled after loop over hits assigned to segment;
	      
	      rawHits->addMuonCalibRawHit( rawTgcHit );
	      
	    }
	  }
      }
    } // if ( m_doTGCs )
    
    // TGC raw coin...
    if( m_doTGCCoin ){
      std::vector<const Muon::TgcCoinDataCollection*> tgcCoinCols;
      const Muon::TgcCoinDataContainer* tgcCoinPrds[3] = { 0, 0, 0 };
      
      for(int ibc=0; ibc<3; ibc++){
	int bcTag=ibc+1;
	std::ostringstream location;
	location << "TrigT1CoinDataCollection" << (bcTag==TgcDigit::BC_PREVIOUS ? "PriorBC" : "")
		 << (bcTag==TgcDigit::BC_NEXT ? "NextBC" : "");	
	if( !evtStore()->contains<Muon::TgcCoinDataContainer>(location.str()) ){
	  ATH_MSG_DEBUG("TgcCoinDataContainer not contained in SG at" << location.str());
	}else if(evtStore()->retrieve(tgcCoinPrds[ibc],location.str()).isFailure()) {
	  ATH_MSG_DEBUG("Cannot retrieve TgcCoinDataContainer" << location.str());
	}else{
	  ATH_MSG_DEBUG("Retrieved TgcCoinDataContainer " << location.str());
	  
	  tgcCoinCols.clear();
	  tgcCoinCols.reserve(tgcCoinPrds[ibc]->size()); // number of tgcCoin chambers
	  Muon::TgcCoinDataContainer::const_iterator it = tgcCoinPrds[ibc]->begin();
	  Muon::TgcCoinDataContainer::const_iterator it_end = tgcCoinPrds[ibc]->end();
	  for( ; it!=it_end; ++it ) {
	    tgcCoinCols.push_back( *it );
	  }
	}
	
	std::vector< const Muon::TgcCoinDataCollection*>::const_iterator tgcCoinCollection = tgcCoinCols.begin();
	std::vector< const Muon::TgcCoinDataCollection*>::const_iterator tgcCoinlastColl = tgcCoinCols.end();
	
	for ( ; tgcCoinCollection != tgcCoinlastColl ; ++tgcCoinCollection ) 
	  {
	    Muon::TgcCoinDataCollection::const_iterator tgcCoin_it = (*tgcCoinCollection)->begin();
	    Muon::TgcCoinDataCollection::const_iterator tgcCoin_it_end = (*tgcCoinCollection)->end();
	    for( ; tgcCoin_it!=tgcCoin_it_end; ++ tgcCoin_it){
	      
	      if((*tgcCoin_it)->type()==Muon::TgcCoinData::TYPE_UNKNOWN || 
		 (*tgcCoin_it)->type()==Muon::TgcCoinData::TYPE_TRACKLET_EIFI) continue;
	      
	      MuonCalibRawTgcCoin* rawTgcCoin = new MuonCalibRawTgcCoin();      
	      MuonFixedId fID = m_idToFixedIdTool->idToFixedId( (*tgcCoin_it)->identify() ) ; 
	      rawTgcCoin->setId( fID );
	      
	      //*** calculation of sector ***//
	      int sector=-1;
	      if((*tgcCoin_it)->isForward()){//forward
		sector=(*tgcCoin_it)->phi()/2+1;//(phi,sector)= (1,1), (2,2), (3,2), (4,3)
	      }else{//endcap
		sector=((*tgcCoin_it)->phi()+1)/4+1;//(phi,sector)= (1,1), (2,1), (3,2), (4,2)
	      }
	      if(sector>12)sector=1;
	      //*** calculation of sector end ***//
	      
	      if((*tgcCoin_it)->type()==Muon::TgcCoinData::TYPE_TRACKLET){
		rawTgcCoin->setGlobalPositionIn((*tgcCoin_it)->globalposIn());
		rawTgcCoin->setGlobalPositionOut((*tgcCoin_it)->globalposOut());
		rawTgcCoin->setType(0);
		rawTgcCoin->setEta(m_idHelperSvc->tgcIdHelper().stationEta((*tgcCoin_it)->identify()));
		rawTgcCoin->setPhi((*tgcCoin_it)->phi());
		rawTgcCoin->setSector(sector);
		rawTgcCoin->setIsForward(static_cast<int>((*tgcCoin_it)->isForward()));
		rawTgcCoin->setIsStrip(static_cast<int>((*tgcCoin_it)->isStrip()));
		rawTgcCoin->setTrackletId((*tgcCoin_it)->trackletId());
		rawTgcCoin->setTrackletIdStrip(0);
		rawTgcCoin->setBcTag(bcTag);
		rawTgcCoin->setWidthIn((*tgcCoin_it)->widthIn());
		rawTgcCoin->setWidthOut((*tgcCoin_it)->widthOut());
		rawTgcCoin->setWidthR(0);
		rawTgcCoin->setWidthPhi(0);
		rawTgcCoin->setDelta((*tgcCoin_it)->delta());
		rawTgcCoin->setRoi(0);
		rawTgcCoin->setPt(0);
		rawTgcCoin->setSub((*tgcCoin_it)->sub());
		rawTgcCoin->setIsPositiveDeltaR(0);
	      } else if((*tgcCoin_it)->type()==Muon::TgcCoinData::TYPE_HIPT){
		rawTgcCoin->setGlobalPositionIn((*tgcCoin_it)->globalposIn());
		rawTgcCoin->setGlobalPositionOut((*tgcCoin_it)->globalposOut());
		rawTgcCoin->setType(1);
		rawTgcCoin->setEta(m_idHelperSvc->tgcIdHelper().stationEta((*tgcCoin_it)->identify()));
		rawTgcCoin->setPhi((*tgcCoin_it)->phi());
		rawTgcCoin->setSector(sector);
		rawTgcCoin->setIsForward(static_cast<int>((*tgcCoin_it)->isForward()));
		rawTgcCoin->setIsStrip(static_cast<int>((*tgcCoin_it)->isStrip()));
		rawTgcCoin->setTrackletId((*tgcCoin_it)->trackletId());
		rawTgcCoin->setTrackletIdStrip(0);
		rawTgcCoin->setBcTag(bcTag);
		rawTgcCoin->setWidthIn((*tgcCoin_it)->widthIn());
		rawTgcCoin->setWidthOut((*tgcCoin_it)->widthOut());
		rawTgcCoin->setWidthR(0);
		rawTgcCoin->setWidthPhi(0);
		rawTgcCoin->setDelta((*tgcCoin_it)->delta());
		rawTgcCoin->setRoi(0);
		rawTgcCoin->setPt(0);
		rawTgcCoin->setSub((*tgcCoin_it)->sub());
		rawTgcCoin->setIsPositiveDeltaR(0);
	      } else if((*tgcCoin_it)->type()==Muon::TgcCoinData::TYPE_SL){
		Amg::Vector3D tmp(0.0,0.0,0.0);
		rawTgcCoin->setGlobalPositionIn(tmp);
		rawTgcCoin->setGlobalPositionOut((*tgcCoin_it)->globalposOut());
		rawTgcCoin->setType(2);
		rawTgcCoin->setEta(m_idHelperSvc->tgcIdHelper().stationEta((*tgcCoin_it)->identify()));
		rawTgcCoin->setPhi((*tgcCoin_it)->phi());
		rawTgcCoin->setSector(sector);
		rawTgcCoin->setIsForward(static_cast<int>((*tgcCoin_it)->isForward()));
		rawTgcCoin->setIsStrip(static_cast<int>((*tgcCoin_it)->isStrip()));
		rawTgcCoin->setTrackletId((*tgcCoin_it)->trackletId());
		rawTgcCoin->setTrackletIdStrip((*tgcCoin_it)->trackletIdStrip());
		rawTgcCoin->setBcTag(bcTag);
		rawTgcCoin->setWidthIn(0);
		rawTgcCoin->setWidthOut(0);
		double w1 = 1./std::sqrt((*tgcCoin_it)->errMat()(Trk::loc1,Trk::loc2));
		double w2 = 1./std::sqrt((*tgcCoin_it)->errMat()(Trk::loc2,Trk::loc2));
		rawTgcCoin->setWidthR(w1);
		rawTgcCoin->setWidthPhi(w2);
		rawTgcCoin->setDelta(0);
		rawTgcCoin->setRoi((*tgcCoin_it)->roi());
		rawTgcCoin->setPt((*tgcCoin_it)->pt());	      	      
		rawTgcCoin->setSub(0);
		rawTgcCoin->setIsPositiveDeltaR((*tgcCoin_it)->isPositiveDeltaR());
	      }
	      
	      rawHits->addMuonCalibRawHit( rawTgcCoin );
	    }
	  }
      }
    } // if ( m_doTGCCoin )
    
    if( m_doRPCCoin ){
      const Muon::RpcCoinDataContainer* rpcCoinContainer = 0;
      if( !evtStore()->contains<Muon::RpcCoinDataContainer>("RPC_triggerHits") ){
	ATH_MSG_DEBUG("RpcCoinDataContainer not contained in SG");
      }else if(evtStore()->retrieve(rpcCoinContainer,"RPC_triggerHits").isFailure()) {
	ATH_MSG_DEBUG("Cannot retrieve RpcCoinDataContainer");
      }else{
	
	ATH_MSG_DEBUG("Retrieved RpcCoinDataContainer " << rpcCoinContainer->numberOfCollections());
	
	Muon::RpcCoinDataContainer::const_iterator it = rpcCoinContainer->begin();
	Muon::RpcCoinDataContainer::const_iterator it_end = rpcCoinContainer->end();
	for( ; it!=it_end; ++it ) {
	  const Muon::RpcCoinDataCollection* col = *it;
	  if( !col ) continue;
	  Muon::RpcCoinDataCollection::const_iterator cit = col->begin();
	  Muon::RpcCoinDataCollection::const_iterator cit_end = col->end();
	  for( ;cit!=cit_end;++cit ){
	    if( !*cit ) continue;
	    const Muon::RpcCoinData& coinData = **cit;
	    MuonCalibRawRpcHit* rawRpc = createRawRpcHit(coinData);
	    if( !rawRpc ) continue;
	    bool lowPtCm  = coinData.isLowPtCoin();
	    MuonCalibRawRpcCoin* rawRpcCoin = new MuonCalibRawRpcCoin(*rawRpc,
								      coinData.ijk(),
								      coinData.threshold(),
								      coinData.overlap(),
								      coinData.parentCmId(),
								      coinData.parentPadId(),
								      coinData.parentSectorId(),
								      lowPtCm);
	    rawHits->addMuonCalibRawHit(rawRpcCoin);
	    
	    if(msgLvl(MSG::DEBUG)){
	      MuonCalibRawRpcCoin& hit = *rawRpcCoin;
	      ATH_MSG_DEBUG(" RpcCoinData: sector " << hit.parentSectorId() << " pad " << hit.parentPadId() 
		  << " CLHEP::cm " << hit.parentCmId() << " ijk " << hit.ijk() << " lowPT " << hit.lowPtCm() 
		  << " overlap " << hit.overlap() << " time " << hit.t() 
		  << " eta " << hit.globalPosition().eta()
		  << " phi " << hit.globalPosition().phi());
	    }
	    delete rawRpc;
	  }
	}
      }
    } // if( m_doRPCCoin )

    return rawHits;
  }
  
  const MuonCalibRawTriggerHitCollection* MuonCalibAlg::retrieveRawTriggerHits() const{
    // Rpc trigger hits
    const RpcPadContainer* rpcRDO   = nullptr;  
    const RpcPad*          rdoColl  = nullptr; 
    MuonCalibRawTriggerHitCollection* rawTriggerHits = new MuonCalibRawTriggerHitCollection();
    if( !evtStore()->contains<RpcPadContainer>("RPCPAD") ){
      ATH_MSG_DEBUG("RpcPadContainer not contained in SG");
    }else if(evtStore()->retrieve(rpcRDO,"RPCPAD").isFailure()) {
      ATH_MSG_DEBUG("Cannot retrieve RpcPadContainer");
    }else{
      ATH_MSG_VERBOSE("Retrieved RawRpcTriggerDataContainer " << rpcRDO->numberOfCollections());
      for (RpcPadContainer::const_iterator rdoColli = rpcRDO->begin(); rdoColli!=rpcRDO->end(); ++rdoColli){
	rdoColl = *rdoColli;
	// Now loop on the RDO
	if ( (rdoColl)->size() != 0 ) {      // number of cma/pad   
	  // for each pad, loop over cma
	  RpcPad::const_iterator it_cma = (rdoColl)->begin(); 
	  RpcPad::const_iterator it_cma_end = (rdoColl)->end();
	  for (; it_cma!=it_cma_end; ++it_cma) { 
	    const RpcCoinMatrix * cma = (*it_cma);
	    //for each cma loop over fired channels
	    ATH_MSG_DEBUG("Number of fired channels for this CM is " << (cma)->size());
	    RpcCoinMatrix::const_iterator it_chan = (*it_cma)->begin(); 
	    RpcCoinMatrix::const_iterator it_chan_end = (*it_cma)->end();
	    for (; it_chan!=it_chan_end; ++it_chan) { 
	      const RpcFiredChannel * fChannel = (*it_chan);
	      //Make a new trigger hit
	      MuonCalibRawRpcTriggerHit* rawRpcTriggerHit = new MuonCalibRawRpcTriggerHit( (rdoColl)->sector(),
											   (rdoColl)->onlineId(),
											   (rdoColl)->status(),
											   (rdoColl)->errorCode(),
											   cma->onlineId(),
											   cma->fel1Id(),
											   cma->febcId(),
											   cma->crc(),
											   fChannel->bcid(),
											   fChannel->time(),
											   fChannel->ijk(),
											   fChannel->channel(),
											   fChannel->ijk() == 7 ? fChannel->ovl() : -1,
											   fChannel->ijk() == 7 ? fChannel->thr() : -1 );
	      //and push it back!
	      rawTriggerHits->addMuonCalibRawTriggerHit(rawRpcTriggerHit);
	    }
	  }
	}
      } //for (RpcPadContainer::const_iterator rdoColli.... 
    } 
    ATH_MSG_VERBOSE("rawTriggerHits for this event has " << rawTriggerHits->numberOfMuonCalibRawRpcTriggerHits() << "  RPC trigger hits ");
    return rawTriggerHits;
  }
  
  const MuonCalibEvent* MuonCalibAlg::retrieveEvent() const {
    //Retrieve all ingredients needed to build an MuonCalibEvent
    MuonCalibEventInfo                 eventInfo = retrieveEventInfo();
    
    const MuonCalibEvent::MCPVec       patterns  = retrievePatterns();
    
    MuonCalibEvent* event = new MuonCalibEvent( patterns );

    const MuonCalibTriggerTimeInfo*    triggerTimeInfo = retrieveTriggerTimeInfo();
    if(triggerTimeInfo) event->setMuonCalibTriggerTimeInfo( triggerTimeInfo );
    
    if(m_doRawData){
      //After constructing the patternVector, use this to set occupancy values
      //of the raw hits in the MuonCalibRawHitCollection.
      const MuonCalibRawHitCollection*   rawHits = retrieveRawHits( patterns );
      event->setMuonCalibRawHitCollection( rawHits );
    } else event->setMuonCalibRawHitCollection( 0 );
    
    if( m_doRPCTriggerHits ){
      const MuonCalibRawTriggerHitCollection*   trighits = retrieveRawTriggerHits();
      event->setMuonCalibRawTriggerHitCollection( trighits );
    } else event->setMuonCalibRawTriggerHitCollection( 0 );
        
    if(m_doTruth){
      const MuonCalibTruthCollection*    truth     = retrieveTruthCollection();
      event->setMuonCalibTruthCollection( truth );
    } else event->setMuonCalibTruthCollection( 0 );
    
    if( m_doRpcSectorLogic ){
      const RpcSectorLogicContainer* slLogic = nullptr;
      if( !evtStore()->contains<RpcSectorLogicContainer>("") ){
	ATH_MSG_DEBUG(" RpcSectorLogicContainer not contained in SG");
      }else if(evtStore()->retrieve(slLogic).isFailure()) {
	ATH_MSG_DEBUG("Cannot retrieve RpcSectorLogicContainer");
      }
      event->setRpcSectorLogicContainer(slLogic);
      
      if( slLogic ){
	const RpcSectorLogicContainer& slContainer = *slLogic;
	// loop over container
	RpcSectorLogicContainer::const_iterator slit = slContainer.begin();
	RpcSectorLogicContainer::const_iterator slit_end = slContainer.end();
	for( ;slit!=slit_end;++slit ){
	  const RpcSectorLogic& slLogic = **slit;
	  if( slLogic.empty() ) continue;
	  ATH_MSG_DEBUG(" RpcSectorLogic: sector " << slLogic.sectorId()
	      << " felId " << slLogic.fel1Id()
	      << " bcId " << slLogic.bcid()
	      << " errCode " << slLogic.errorCode()
		<< " crc " << slLogic.crc());
	  RpcSectorLogic::const_iterator hitIt = slLogic.begin();
	  RpcSectorLogic::const_iterator hitIt_end = slLogic.end();
	  for( ;hitIt != hitIt_end; ++hitIt ){
	    const RpcSLTriggerHit& slHit = **hitIt;
	    ATH_MSG_DEBUG(" hit: rowinBcid " << slHit.rowinBcid() 
		<< " padId " << slHit.padId() 
		<< " ptid " << slHit.ptId()
		<< " roi " << slHit.roi()
		<< " outerPlane" << slHit.outerPlane()
		<< " overlapPhi " << slHit.overlapPhi() 
		<< " overlapEta " << slHit.overlapEta()
		<< " triggerBcid " << slHit.triggerBcid()
		  << " isInput" << slHit.isInput());
	  }
	}
      }
    }
    event->setMuonCalibEventInfo( eventInfo );
    return event;
  }
  
  MuonCalibRawRpcHit* MuonCalibAlg::createRawRpcHit( const Muon::RpcPrepData& prd ) const {
    MuonCalibRawRpcHit* rawRpcHit = new MuonCalibRawRpcHit();     	
    MuonFixedId fID = m_idToFixedIdTool->idToFixedId( prd.identify() ) ;
    rawRpcHit->setGlobalPosition( prd.globalPosition() );
    rawRpcHit->setId( fID );
    rawRpcHit->setT( prd.time() );
    // get detector element
    const MuonGM::RpcReadoutElement* detEl = prd.detectorElement();
    rawRpcHit->setWidth( detEl->StripWidth( m_idHelperSvc->rpcIdHelper().measuresPhi(prd.identify()) ) );
    rawRpcHit->setLength( detEl->StripLength(m_idHelperSvc->rpcIdHelper().measuresPhi(prd.identify())));
    return rawRpcHit;
  }
  
}//end namespace MuonCalib
  
