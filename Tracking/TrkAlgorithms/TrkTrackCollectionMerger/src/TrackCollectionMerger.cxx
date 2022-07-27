/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
//   Implementation file for class Trk::TrackCollectionMerger
///////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////
// Version 1.0 11/26/2007 Thomas Koffas
///////////////////////////////////////////////////////////////////

#include "GaudiKernel/MsgStream.h"
#include "TrkPrepRawData/PrepRawData.h"
#include "TrkTrackCollectionMerger/TrackCollectionMerger.h"

///////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////

Trk::TrackCollectionMerger::TrackCollectionMerger
(const std::string& name, ISvcLocator* pSvcLocator  ) :
  AthAlgorithm(name, pSvcLocator ),
  m_createViewCollection(true),
  m_updateSharedHits(true),
  m_updateAdditionalInfo(false),
  m_doTrackOverlay(false)
{
  m_outtracklocation         = "CombinedInDetTracks"    ;

  declareProperty("TracksLocation",                 m_tracklocation          );
  declareProperty("OutputTracksLocation",           m_outtracklocation       ); 
  declareProperty("SummaryTool" ,                   m_trkSummaryTool         );
  declareProperty("CreateViewColllection" ,         m_createViewCollection   );
  declareProperty("UpdateSharedHits" ,              m_updateSharedHits);
  declareProperty("UpdateAdditionalInfo" ,          m_updateAdditionalInfo);
  declareProperty("DoTrackOverlay",                 m_doTrackOverlay);
}

///////////////////////////////////////////////////////////////////
// Initialisation
///////////////////////////////////////////////////////////////////

StatusCode Trk::TrackCollectionMerger::initialize()
{

  ATH_MSG_DEBUG("Initializing TrackCollectionMerger");
  ATH_CHECK(  m_tracklocation.initialize() );
  ATH_CHECK(  m_pileupTRT.initialize( m_doTrackOverlay ) );
  ATH_CHECK(  m_pileupPixel.initialize( m_doTrackOverlay ) );
  ATH_CHECK(  m_pileupSCT.initialize( m_doTrackOverlay ) );
  ATH_CHECK( m_outtracklocation.initialize() );
  if (not m_trkSummaryTool.name().empty()) ATH_CHECK( m_trkSummaryTool.retrieve() );
  if (not m_assoTool.name().empty()) ATH_CHECK( m_assoTool.retrieve() );
  ATH_CHECK( m_assoMapName.initialize( !m_assoMapName.key().empty() ));
  return StatusCode::SUCCESS;
}


///////////////////////////////////////////////////////////////////
// Execute
///////////////////////////////////////////////////////////////////
StatusCode 
Trk::TrackCollectionMerger::execute(){
  std::unique_ptr<TrackCollection> outputCol = m_createViewCollection ?
    std::make_unique<TrackCollection>(SG::VIEW_ELEMENTS) : std::make_unique<TrackCollection>();
  ATH_MSG_DEBUG("Number of Track collections " << m_tracklocation.size());
  
  // pre-loop to reserve enough space in the output collection
  std::vector<const TrackCollection*> trackCollections;
  trackCollections.reserve(m_tracklocation.size());
  size_t ttNumber = 0;
  for (auto& tcname : m_tracklocation){
    ///Retrieve tracks from StoreGate
    SG::ReadHandle<TrackCollection> trackCol (tcname);
    trackCollections.push_back(trackCol.cptr());
    ttNumber += trackCol->size();
  }
  std::unique_ptr<Trk::PRDtoTrackMap> pPrdToTrackMap(m_assoTool ? m_assoTool->createPRDtoTrackMap(): nullptr);
  // reserve the right number of entries for the output collection
  outputCol->reserve(ttNumber);
  // merging loop
  for(auto& tciter : trackCollections){
      // merge them in
    if(mergeTrack(tciter, pPrdToTrackMap.get(), outputCol.get()).isFailure()){
	     ATH_MSG_ERROR( "Failed to merge tracks! ");
      }
  }
  ATH_MSG_DEBUG("Size of combined tracks " << outputCol->size());
  if (m_trkSummaryTool){
    ATH_MSG_DEBUG("Update summaries");  
    // now loop over all tracks and update summaries with new shared hit counts
    // @TODO magic! tracks are now non-const !??
    const bool createTrackSummary = not (m_updateAdditionalInfo or m_updateSharedHits);
    for (Trk::Track* trk : *outputCol) {
      if (createTrackSummary){
         m_trkSummaryTool->computeAndReplaceTrackSummary(*trk, pPrdToTrackMap.get(), false /* DO NOT suppress hole search*/);
      } else {
        if (m_updateAdditionalInfo)  m_trkSummaryTool->updateAdditionalInfo(*trk);
        if (m_updateSharedHits) m_trkSummaryTool->updateSharedHitCount(*trk, pPrdToTrackMap.get());
      }
    }
  } else {
    ATH_MSG_WARNING("No track summary update performed because the TrackSummaryTool was not specified");  
  }
  SG::WriteHandle<TrackCollection> h_write(m_outtracklocation);
  ATH_CHECK(h_write.record(std::move(outputCol)));	     
  //
  if (!m_assoMapName.key().empty()) {
     SG::WriteHandle<Trk::PRDtoTrackMap> write_handle(m_assoMapName);
     if (write_handle.record( m_assoTool->reduceToStorableMap(std::move(pPrdToTrackMap))).isFailure()) {
        ATH_MSG_FATAL("Failed to add PRD to track association map.");
     }
  }
  //Print common event information
  ATH_MSG_DEBUG("Done !");  
  return StatusCode::SUCCESS;
}

///////////////////////////////////////////////////////////////////
// Finalize
///////////////////////////////////////////////////////////////////

StatusCode Trk::TrackCollectionMerger::finalize()
{
   return StatusCode::SUCCESS;
}



///////////////////////////////////////////////////////////////////
// Merge track collections and remove duplicates
///////////////////////////////////////////////////////////////////

StatusCode Trk::TrackCollectionMerger::mergeTrack(const TrackCollection* trackCol,
                                                  Trk::PRDtoTrackMap *pPrdToTrackMap,
                                                  TrackCollection* outputCol)
{
  // loop over tracks, accept them and add them into association tool
  if(trackCol && !trackCol->empty()) {
    ATH_MSG_DEBUG("Size of track collection " << trackCol->size());
    if (not pPrdToTrackMap) ATH_MSG_WARNING("No valid PRD to Track Map; was the association tool name missing?");
    // loop over tracks
    for(const auto *const rf: *trackCol){
      // add track into output
      // FIXME: const_cast
      // These objects are modified in the `Update summaries' section
      Trk::Track* newTrack = m_createViewCollection ? const_cast<Trk::Track*>(rf) : new Trk::Track(*rf);
      outputCol->push_back(newTrack);
      // add tracks into PRD tool
      if (m_assoTool and m_assoTool->addPRDs(*pPrdToTrackMap, *newTrack).isFailure())
	      ATH_MSG_WARNING( "Failed to add PRDs to map" );
    }
  }

  return StatusCode::SUCCESS;
}

