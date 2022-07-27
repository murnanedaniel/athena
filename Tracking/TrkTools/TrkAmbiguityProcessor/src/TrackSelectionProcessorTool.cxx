/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "TrackSelectionProcessorTool.h"
#include "TrackScoringTool.h"
#include "TrkToolInterfaces/IPRD_AssociationTool.h"
#include "TrkTrack/TrackCollection.h"
#include "AthContainers/ConstDataVector.h"
#include "GaudiKernel/MsgStream.h"
#include <map>

//==================================================================================================
Trk::TrackSelectionProcessorTool::TrackSelectionProcessorTool(const std::string& t, 
                const std::string& n,
                const IInterface*  p )
  :
  AthAlgTool(t,n,p),
  m_scoringTool("Trk::TrackScoringTool/TrackScoringTool"), 
  m_selectionTool("InDet::InDetAmbiTrackSelectionTool/InDetAmbiTrackSelectionTool")
{
  declareInterface<ITrackAmbiguityProcessorTool>(this);
  declareProperty("DropDouble"           , m_dropDouble         = true);
  declareProperty("ScoringTool"          , m_scoringTool);
  declareProperty("SelectionTool"        , m_selectionTool);
  declareProperty("DisableSorting"       , m_disableSorting     = false);
}
//==================================================================================================

Trk::TrackSelectionProcessorTool::~TrackSelectionProcessorTool()= default;
//==================================================================================================

StatusCode 
Trk::TrackSelectionProcessorTool::initialize(){
  StatusCode sc = AthAlgTool::initialize();
  if (sc.isFailure()) {
    ATH_MSG_FATAL( "AlgTool::initialise failed" );
    return StatusCode::FAILURE;
  }
  ATH_CHECK( m_assoMapName.initialize(!m_assoMapName.key().empty()));
  ATH_CHECK( m_assoTool.retrieve() );
  ATH_CHECK( m_scoringTool.retrieve());
  ATH_CHECK(m_selectionTool.retrieve());

  ATH_CHECK(m_clusterSplitProbContainerIn.initialize(!m_clusterSplitProbContainerIn.key().empty()));
  ATH_CHECK(m_clusterSplitProbContainerOut.initialize(!m_clusterSplitProbContainerOut.key().empty()));

  if (m_disableSorting) ATH_MSG_INFO( "Internal sorting disabled, using external ordering!" );    
  return sc;
}
//==================================================================================================

StatusCode 
Trk::TrackSelectionProcessorTool::finalize(){
  StatusCode sc = AlgTool::finalize(); 
  return sc;
}

//==================================================================================================

/** Do actual processing of event. Takes a track container, 
    and then returns the tracks which have been selected*/

const TrackCollection*  
Trk::TrackSelectionProcessorTool::process(const TrackCollection* tracksCol,
                                          Trk::PRDtoTrackMap *pPrdToTrackMap) const{
  //TODO: make sure the ownership; delete origin tracks from map?
  std::vector<const Track*> tracks;
  tracks.reserve(tracksCol->size());
  for(const auto *e: *tracksCol){
    tracks.push_back(e);
  }
  ATH_MSG_DEBUG ("Processing "<<tracks.size()<<" tracks");
  std::unique_ptr<Trk::PRDtoTrackMap> tmpPrdToTrackMap;
  if (!pPrdToTrackMap) {
     tmpPrdToTrackMap = m_assoTool->createPRDtoTrackMap();
     if (!m_assoMapName.key().empty()) {
        SG::ReadHandle<Trk::PRDtoTrackMap> inputPrdMap(m_assoMapName);
        if (!inputPrdMap.isValid()) {
           ATH_MSG_ERROR("Failed to retrieve prd to track map " << m_assoMapName.key() );
        } else {
           *tmpPrdToTrackMap = *inputPrdMap;
        }
     }
     pPrdToTrackMap = tmpPrdToTrackMap.get();
  }
  TrackScoreMap trackScoreTrackMap;
  //put tracks into maps etc
  addNewTracks(trackScoreTrackMap,*pPrdToTrackMap, tracks);
  // going to do simple algorithm for now:
  // - take track with highest score
  // - remove shared hits from all other tracks
  // - take next highest scoring tracks, and repeat 
  std::unique_ptr<ConstDataVector<TrackCollection> > result(std::make_unique<ConstDataVector<TrackCollection> >(SG::VIEW_ELEMENTS)); //TODO, old or new
  solveTracks(trackScoreTrackMap, *pPrdToTrackMap, *result);
  if (msgLvl(MSG::DEBUG)) dumpTracks(*result->asDataVector());
  return result.release()->asDataVector();
}


//==================================================================================================
void 
Trk::TrackSelectionProcessorTool::addNewTracks(TrackScoreMap &trackScoreTrackMap,
                                                    Trk::PRDtoTrackMap &prdToTrackMap,
                                                    const std::vector<const Track*> &tracks) const{
  ATH_MSG_DEBUG ("Number of tracks at Input: "<<tracks.size());
  PrdSignatureSet prdSigSet;
  TrackScore itrack=0;
  for (const Track*a_track : tracks )   {
    if(m_disableSorting) {
      // add track to map using ordering provided by the collection
      trackScoreTrackMap.insert( std::make_pair(itrack, TrackPtr(a_track)) );
      itrack++;
      continue;
    }
    bool reject = false;
    // only fitted tracks get hole search, input is not fitted
    TrackScore score = m_scoringTool->score( *a_track, true );
    // veto tracks with score 0
    if (score==0) { 
      ATH_MSG_DEBUG ("Track score is zero, reject it");
      reject = true;
    } else {
      if (m_dropDouble) {
        const std::vector<const Trk::PrepRawData*> & prds = m_assoTool->getPrdsOnTrack(prdToTrackMap, *a_track);
        // unfortunately PrepRawDataSet is not a set !
        PrdSignature prdSig;
        prdSig.insert( prds.begin(),prds.end() );
        // we try to insert it into the set, if we fail (pair.second), it then exits already
        if ( !(prdSigSet.insert(prdSig)).second ) {
          ATH_MSG_DEBUG ("Double track, reject it !");
          reject = true;
        } else {
          ATH_MSG_DEBUG ("Insert new track in PrdSignatureSet");
        }
      }
    }
    if (!reject) {
      // add track to map, map is sorted small to big ! set if fitted
      ATH_MSG_VERBOSE ("Track  ("<< a_track <<") has score "<<score);
      trackScoreTrackMap.insert( std::make_pair(-score, TrackPtr(a_track) ) );
    }
  }
  ATH_MSG_DEBUG ("Number of tracks in map:"<<trackScoreTrackMap.size());
}

void
Trk::TrackSelectionProcessorTool::solveTracks(TrackScoreMap &trackScoreTrackMap,
                                              Trk::PRDtoTrackMap &prdToTrackMap,
                                              ConstDataVector<TrackCollection> &result) const
{
  using namespace std;

  const EventContext& ctx = Gaudi::Hive::currentContext();
  SG::ReadHandle<Trk::ClusterSplitProbabilityContainer> splitProbContainerIn;
  if (!m_clusterSplitProbContainerIn.key().empty()) {
     splitProbContainerIn = SG::ReadHandle( m_clusterSplitProbContainerIn, ctx);
     if (!splitProbContainerIn.isValid()) {
        ATH_MSG_ERROR( "Failed to get input cluster split probability container "  << m_clusterSplitProbContainerIn.key());
     }
  }
  std::unique_ptr<Trk::ClusterSplitProbabilityContainer> splitProbContainerCleanup(!m_clusterSplitProbContainerIn.key().empty()
                                                                                      ? std::make_unique<ClusterSplitProbabilityContainer>(*splitProbContainerIn)
                                                                                      : std::make_unique<ClusterSplitProbabilityContainer>());
  SG::WriteHandle<Trk::ClusterSplitProbabilityContainer> splitProbContainerHandle;
  Trk::ClusterSplitProbabilityContainer *splitProbContainer;
  if (!m_clusterSplitProbContainerOut.key().empty()) {
     splitProbContainerHandle = SG::WriteHandle<Trk::ClusterSplitProbabilityContainer>( m_clusterSplitProbContainerOut, ctx);
     if (splitProbContainerHandle.record(std::move(splitProbContainerCleanup)).isFailure()) {
        ATH_MSG_FATAL( "Failed to record output cluster split probability container "  << m_clusterSplitProbContainerOut.key());
     }
     splitProbContainer=splitProbContainerHandle.ptr();
  }
  else {
     splitProbContainer=splitProbContainerCleanup.get();
  }

  ATH_MSG_VERBOSE ("Starting to solve tracks");
  // now loop as long as map is not empty
  while ( !trackScoreTrackMap.empty() ) {
    TrackScoreMap::iterator itnext = trackScoreTrackMap.begin();
    TrackPtr atrack( std::move(itnext->second) );
    TrackScore ascore( itnext->first);
    trackScoreTrackMap.erase(itnext);
    ATH_MSG_VERBOSE ("--- Trying next track "<<atrack.track()<<"\t with score "<<-ascore);
    std::unique_ptr<Trk::Track> cleanedTrack;
    const auto &[cleanedTrack_tmp, keepOriginal]  = m_selectionTool->getCleanedOutTrack( atrack.track() , -(ascore), *splitProbContainer, prdToTrackMap, -1, -1);
    cleanedTrack.reset(cleanedTrack_tmp);
    if (keepOriginal ){
      // track can be kept as identical to the input track
      ATH_MSG_DEBUG ("Accepted track "<<atrack.track()<<"\t has score "<<-(ascore));
      // add track to PRD_AssociationTool
      StatusCode sc = m_assoTool->addPRDs(prdToTrackMap,*atrack);
      if (sc.isFailure()) ATH_MSG_ERROR( "addPRDs() failed" );
      // add to output list
      result.push_back( atrack.track() );

    } else if ( !cleanedTrack ) {
      // track should be discarded
      ATH_MSG_DEBUG ("Track "<< atrack.track() << " doesn't meet the cuts of the AmbiTrack Selection tool");
    } else  {
      // delete cleaned track
      cleanedTrack.reset();
      // stripped down version cannot be handled discarding
      ATH_MSG_DEBUG("Selection tool returned a new track, cannot handle memory management of new track, deleting it. Check you configuration ");
    }
    // don't forget to drop track from map
  }
  ATH_MSG_DEBUG ("Finished, number of track on output: "<<result.size());
}

//==================================================================================================

void 
Trk::TrackSelectionProcessorTool::dumpTracks( const TrackCollection& tracks ) const{
  ATH_MSG_VERBOSE ("Dumping tracks in collection");
  int num=0;
  TrackScore totalScore = 0;
  TrackCollection::const_iterator it    = tracks.begin();
  TrackCollection::const_iterator itEnd = tracks.end();
  for (; it != itEnd ; ++it){
    // score track:
    const TrackScore score = m_scoringTool->score( **it, true );
    ATH_MSG_VERBOSE (num++<<"\tTrack :"<<*it<<"\tScore: "<<score);
    totalScore+=score;
  }
  ATH_MSG_DEBUG ("Total event score : "<<totalScore);
}
