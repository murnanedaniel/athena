/*
  General-purpose view creation algorithm <bwynne@cern.ch>
  
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
*/

#include "EventViewCreatorAlgorithmWithJets.h"
#include "AthLinks/ElementLink.h"
#include "AthViews/ViewHelper.h"
#include "AthViews/View.h"
#include "DecisionHandling/TrigCompositeUtils.h"
#include "DecisionHandling/HLTIdentifier.h"


EventViewCreatorAlgorithmWithJets::EventViewCreatorAlgorithmWithJets( const std::string& name, ISvcLocator* pSvcLocator )
  : EventViewCreatorAlgorithm( name, pSvcLocator ) {}

EventViewCreatorAlgorithmWithJets::~EventViewCreatorAlgorithmWithJets() {}

StatusCode EventViewCreatorAlgorithmWithJets::initialize() {
  EventViewCreatorAlgorithm::initialize();

  ATH_CHECK( m_inViewJets.initialize() );

  return StatusCode::SUCCESS;
}

StatusCode EventViewCreatorAlgorithmWithJets::execute( const EventContext& context ) const { 
  auto outputHandles = decisionOutputs().makeHandles( context );     
  // make and store the views
  auto viewsHandle = SG::makeHandle( m_viewsKey ); 
  auto viewVector1 = std::make_unique< ViewContainer >();
  ATH_CHECK( viewsHandle.record(  std::move( viewVector1 ) ) );
  auto viewVector = viewsHandle.ptr();

  auto contexts = std::vector<EventContext>( );
  unsigned int viewCounter = 0;
  unsigned int conditionsRun = context.getExtension<Atlas::ExtendedEventContext>().conditionsRun();

  //  const TrigRoiDescriptor* previousRoI = 0;
  // mapping RoI with index of the View in the vector
  // This is used to link the same view to differnt decisions that come from the same RoI
  std::map <const TrigRoiDescriptor*, int> viewMap;
  std::map <const TrigRoiDescriptor*, int>::iterator itViewMap;
  int outputIndex = -1;
 // Loop over all input containers, which are of course TrigComposites, and request their features
  // this is the same as InputMaker, apart from the view creation. The loop can be splitted in two loops, to have one common part
  for ( auto inputKey: decisionInputs() ) {
    outputIndex++;
    auto inputHandle = SG::makeHandle( inputKey, context );
    if( not inputHandle.isValid() ) {
      ATH_MSG_DEBUG( "Got no decisions from input " << inputKey.key() );
      continue;
    }
    
    if( inputHandle->size() == 0 ) { // input filtered out
      ATH_MSG_ERROR( "Got 0 decisions from valid input "<< inputKey.key()<<". Is it expected?");
      return StatusCode::FAILURE;
    }
    ATH_MSG_DEBUG( "Got input " << inputKey.key() << " with " << inputHandle->size() << " elements" );
    
     // prepare output decisions
    TrigCompositeUtils::createAndStore(outputHandles[outputIndex]);
    TrigCompositeUtils::DecisionContainer* outputDecisions = outputHandles[outputIndex].ptr();
    
    const TrigRoiDescriptor* prevRoIDescriptor = nullptr;
    int inputCounter = -1;
    for ( auto inputDecision: *inputHandle ) {
      inputCounter++;      
      // Retrieve jets ...
      ATH_MSG_DEBUG( "Checking there are jets linked to decision object" );
      TrigCompositeUtils::LinkInfo< xAOD::JetContainer > jetELInfo = TrigCompositeUtils::findLink< xAOD::JetContainer >( inputDecision,m_jetsLink );
      ATH_CHECK( jetELInfo.isValid() );
      const xAOD::Jet *jet = *jetELInfo.link;
      ATH_MSG_DEBUG( "Placing xAOD::JetContainer " );
      ATH_MSG_DEBUG( "   -- pt="<< jet->p4().Et() <<" eta="<< jet->eta() << " phi="<< jet->phi() );

      // pull RoI descriptor
      TrigCompositeUtils::LinkInfo<TrigRoiDescriptorCollection> roiELInfo = TrigCompositeUtils::findLink<TrigRoiDescriptorCollection>(inputDecision, m_roisLink );
      ATH_CHECK( roiELInfo.isValid() );
      // associate this RoI to output decisions
      auto roiDescriptor = *roiELInfo.link;
      ATH_MSG_DEBUG( "Placing TrigRoiDescriptor " );
      ATH_MSG_DEBUG( "   " << *roiDescriptor );

      TrigCompositeUtils::Decision* newDecision = nullptr;      
      if ( prevRoIDescriptor != roiDescriptor ) {
	//make one TC decision output per input and connect to previous
	newDecision = TrigCompositeUtils::newDecisionIn( outputDecisions, name() );
	TrigCompositeUtils::linkToPrevious( newDecision, inputKey.key(), inputCounter );
	insertDecisions( inputDecision, newDecision );
	newDecision->setObjectLink( "initialRoI", roiELInfo.link );
	newDecision->setObjectLink( "jets", jetELInfo.link );
	prevRoIDescriptor = roiDescriptor;
      } else {
	newDecision = outputDecisions->back();
	newDecision->setObjectLink( "seedEnd", ElementLink<TrigCompositeUtils::DecisionContainer>( inputHandle.key(), inputCounter ) );
	insertDecisions( inputDecision, newDecision );
	ATH_MSG_DEBUG("No need to create another output decision object, just adding decision IDs");
      }
      
      // search for existing view
      itViewMap = viewMap.find(roiDescriptor);
      if ( itViewMap != viewMap.end() ) {
	int iview = itViewMap->second;
	newDecision->setObjectLink( "view", ElementLink< ViewContainer >(m_viewsKey.key(), iview ) ); //adding view to TC
	ATH_MSG_DEBUG( "Adding already mapped view " << iview << " in ViewVector , to new decision");
	//	need to check if this View has parent views? can we have more than one parent views?
      } else {
	
	ATH_MSG_DEBUG( "Positive decisions on RoI, preparing view" );
	
	// make the view
	ATH_MSG_DEBUG( "Making the View" );
	auto newView = ViewHelper::makeView( name()+"_view", viewCounter++, m_viewFallThrough ); //pointer to the view
	viewVector->push_back( newView );
	contexts.emplace_back( context );
	contexts.back().setExtension( Atlas::ExtendedEventContext( viewVector->back(), conditionsRun ) );
	
	// link decision to this view
	newDecision->setObjectLink( "view", ElementLink< ViewContainer >(m_viewsKey.key(), viewVector->size()-1 ));//adding view to TC
      	viewMap[roiDescriptor]=viewVector->size()-1;
	ATH_MSG_DEBUG( "Adding new view to new decision; storing view in viewVector component " << viewVector->size()-1 );

	ATH_CHECK( linkViewToParent( inputDecision, viewVector->back() ) );
	ATH_CHECK( placeRoIInView( roiDescriptor, viewVector->back(), contexts.back() ) );	
	ATH_CHECK( placeJetInView( jet, viewVector->back(), contexts.back() ) );
      }
    }
    
    ATH_MSG_DEBUG( "Recording output key " <<  decisionOutputs()[ outputIndex ].key() <<" of size "<< outputDecisions->size()  <<" at index "<< outputIndex);
  }

  ATH_MSG_DEBUG( "Launching execution in " << viewVector->size() << " views" );
  ATH_CHECK( ViewHelper::ScheduleViews( viewVector,           // Vector containing views
					m_viewNodeName,             // CF node to attach views to
					context,                    // Source context
					m_scheduler.get() ) );
  
  // report number of views, stored already when container was created
  ATH_MSG_DEBUG( "Store "<< viewsHandle->size() <<" Views");

  size_t validInputCount = countInputHandles( context );  
  size_t validOutputCount = 0;
  for ( auto outHandle: outputHandles ) {
    if( not outHandle.isValid() ) continue;
    validOutputCount++;
  }
  
  ATH_MSG_DEBUG("Produced " << validOutputCount << " decisions containers");
  if(validInputCount != validOutputCount ) {
    ATH_MSG_ERROR("Found " << validInputCount << " inputs and " << validOutputCount << " outputs");
    return StatusCode::FAILURE;
  }
  printDecisions( outputHandles );     
  return StatusCode::SUCCESS;
}

StatusCode EventViewCreatorAlgorithmWithJets::placeJetInView( const xAOD::Jet* theObject, SG::View* view, const EventContext& context ) const {
  // fill the Jet output collection  
  ATH_MSG_DEBUG( "Adding Jet To View : " << m_inViewJets.key() );
  auto oneObjectCollection = std::make_unique< ConstDataVector< xAOD::JetContainer > >();
  oneObjectCollection->clear( SG::VIEW_ELEMENTS ); 
  oneObjectCollection->push_back( theObject );

  //store in the view 
  auto handle = SG::makeHandle( m_inViewJets,context );
  ATH_CHECK( handle.setProxyDict( view ) );
  ATH_CHECK( handle.record( std::move( oneObjectCollection ) ) ); 
  return StatusCode::SUCCESS;
}

