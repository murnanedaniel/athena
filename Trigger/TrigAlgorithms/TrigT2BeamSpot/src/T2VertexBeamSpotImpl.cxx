/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

//============================================================
//
// $Id: T2VertexBeamSpotImpl.cxx 641309 2015-01-23 15:45:50Z smh $
//
// T2VertexBeamSpot.cxx, (c) ATLAS Detector software
// Trigger/TrigAlgorithms/TrigT2BeamSpot/T2VertexBeamSpot
//
// Online beam spot measurement and monitoring using
// L2 recontructed primary vertices.
//
// Authors : David W. Miller, Rainer Bartoldus,   
//           Su Dong
//
// Date : 11 May, 2008
//
//============================================================

// This algorithm
#include "T2VertexBeamSpotImpl.h"
#include "T2TrackClusterer.h"
#include "T2Timer.h"

// Specific to this algorithm
#include "TrigInDetEvent/TrigInDetTrack.h"
#include "TrigInDetEvent/TrigInDetTrackCollection.h"
#include "TrigInDetEvent/TrigVertex.h"
#include "TrigInDetEvent/TrigVertexCollection.h"
#include "TrigInDetToolInterfaces/ITrigPrimaryVertexFitter.h"

// Generic Trigger tools
#include "TrigNavigation/TriggerElement.h"
#include "TrigSteeringEvent/TrigRoiDescriptor.h"
#include "TrigTimeAlgs/TrigTimer.h"

// Generic tools
#include "CLHEP/Units/SystemOfUnits.h"
using CLHEP::GeV;

#include <string>
#include <sstream>
#include <cmath>

using std::string;
using std::ostringstream;
using std::vector;
using std::abs;

using namespace PESA;

//#define DEBUG if (msgLvl()<=MSG::DEBUG) msg() << MSG::DEBUG

namespace Beamspot
{ 
  class TrackPTSort
  {
  public:
    inline bool operator () (const TrigInDetTrack* trk1, const TrigInDetTrack* trk2) 
    { 
      return abs(trk1->param()->pT()) > abs(trk2->param()->pT());
    }
  };

  class TrackPhiSort
  {
  public:
    inline bool operator () (const TrigInDetTrack* trk1, const TrigInDetTrack* trk2) 
    { 
      return abs(trk1->param()->phi0()) > abs(trk2->param()->phi0());
    }
  };
}


namespace
{
  TrigVertex::AlgoId getPrimaryAlgoId( const TrigInDetTrack& track )
  {
    // Doing a linear search here is not ideal but still okay given that there is
    // only one strategy running in P1 at any given time, which we simply keep in the front.
    return
      ( track.algorithmId() == TrigInDetTrack::STRATEGY_B_ID ) ? TrigVertex::BSFULL_STRATEGY_B_ID :
      ( track.algorithmId() == TrigInDetTrack::STRATEGY_A_ID ) ? TrigVertex::BSFULL_STRATEGY_A_ID :
      ( track.algorithmId() == TrigInDetTrack::STRATEGY_F_ID ) ? TrigVertex::BSFULL_STRATEGY_F_ID :
      ( track.algorithmId() == TrigInDetTrack::SITRACKID     ) ? TrigVertex::BSFULLSITRACKID      :
      ( track.algorithmId() == TrigInDetTrack::IDSCANID      ) ? TrigVertex::BSFULLIDSCANID       :
      TrigVertex::NULLID;
  }

  TrigVertex::AlgoId getSplitAlgoId( const TrigInDetTrack& track )
  {
    return
      ( track.algorithmId() == TrigInDetTrack::STRATEGY_B_ID ) ? TrigVertex::BSSPLIT_STRATEGY_B_ID :
      ( track.algorithmId() == TrigInDetTrack::STRATEGY_A_ID ) ? TrigVertex::BSSPLIT_STRATEGY_A_ID :
      ( track.algorithmId() == TrigInDetTrack::STRATEGY_F_ID ) ? TrigVertex::BSSPLIT_STRATEGY_F_ID :
      ( track.algorithmId() == TrigInDetTrack::SITRACKID     ) ? TrigVertex::BSSPLITSITRACKID :
      ( track.algorithmId() == TrigInDetTrack::IDSCANID      ) ? TrigVertex::BSSPLITIDSCANID :
      TrigVertex::NULLID;
  }

  void stripTrackList( TrigVertexCollection& coll )
  {
    for ( TrigVertexCollection::iterator vertex = coll.begin(); vertex != coll.end(); ++vertex )
      {
        // Drop the track list from the vertex so we won't persist it.
        // FIXME: There is no interface for this yet, so for the moment we'll just clear the list.
        TrackInVertexList* vertexTracks = const_cast<TrackInVertexList*>( (*vertex)->tracks() );
        if ( vertexTracks )
          vertexTracks->clear();
      }
  }
}


void
T2VertexBeamSpotImpl::processTEs( const vector<vector<HLT::TriggerElement*> >& tes_in,
                                  TrigInDetTrackCollection& mySelectedTrackCollection )
{
  T2Timer t( m_timer[allTE] );
  //
  // Loop over all input TEs and process the ROIs they contain
  //
  m_TotalROIs       = 0;
  m_TotalTracks     = 0;
  m_TotalTracksPass = 0;
  m_TotalHiPTTracks = 0;

  for ( vector<vector<HLT::TriggerElement*> >::const_iterator iTE = tes_in.begin();
        iTE != tes_in.end(); ++iTE )
    {
      // Number of ROIs in this TE
      m_ROIperTE.push_back( iTE->size() );
      
      // If this TE has no ROIs, go to the next one
      if ( iTE->empty() )
        continue;
    
      // Event has ROI
      eventStage( hasROI );
    
      // Process the ROIs
      processROIs( *iTE, mySelectedTrackCollection );

      // Update counters
      m_TotalROIs         += m_ROIperTE.back(); // FIXME: spelling
      m_TotalTracks       += m_TracksPerTE.back();
      m_TotalTracksPass   += m_TracksPerTEPass.back();
      m_TotalHiPTTracks   += m_HiPTTracksPerTE.back();
    }

  // Reporting
  if (msgLvl()<=MSG::DEBUG) msg() << MSG::DEBUG
                                  << " Total Tracks              : " << m_TotalTracks << '\n'
                                  << " Tracks passing selection  : " << m_TotalTracksPass << '\n'
                                  << " Tracks in my container    : " << mySelectedTrackCollection.size() << '\n'
                                  << " Number passed seed Tracks : " << m_TotalHiPTTracks << endreq;
}


void
T2VertexBeamSpotImpl::processROIs( const HLT::TEVec& myTEVec,
                                   TrigInDetTrackCollection& mySelectedTrackCollection )

{
  T2Timer t( m_timer[allROI] );
  //
  // Loop over all ROIs in the given TEs and select tracks from
  // the collections they contain
  //
  unsigned int countTracksPerTE       = 0;
  unsigned int countTracksPassedPerTE = 0;
  unsigned int countHiPTTracksPerTE   = 0;
  
  const HLT::TEVec::const_iterator firstTE = myTEVec.begin();
  const HLT::TEVec::const_iterator lastTE  = myTEVec.end();

  for ( HLT::TEVec::const_iterator TEiterator = firstTE; TEiterator != lastTE; ++TEiterator )
    {
      if (msgLvl()<=MSG::DEBUG) msg() << MSG::DEBUG  << "TE iterator count = " << TEiterator - firstTE << endreq;

      // Get the ROI descriptor
      const TrigRoiDescriptor* roiDescriptor = 0;
      HLT::ErrorCode status = getFeature( *TEiterator, roiDescriptor );

      if ( status != HLT::OK )
        msg() << MSG::WARNING << "No RoIDescriptors for this Trigger Element! " << endreq;
      else
        if (msgLvl()<=MSG::DEBUG) msg() << MSG::DEBUG  << "Looking at ROI descriptor" << endreq;
        
      // Print out ROI descriptor info
      if ( roiDescriptor )
        {
          // Fill local variables for ROI reference position   
          m_RoIEta.push_back( roiDescriptor->eta()  ); 
          m_RoIPhi.push_back( roiDescriptor->phi()  ); 

          m_RoIID .push_back( roiDescriptor->roiId() );
      
          if (msgLvl()<=MSG::DEBUG)
            msg() << MSG::DEBUG 
                  << "Using inputTE(" << TEiterator - firstTE << ")->getId(): " << (*TEiterator)->getId()
                  << "; RoI ID = " << roiDescriptor->roiId() << ": Eta = " << roiDescriptor->eta()
                  << ", Phi = "    << roiDescriptor->phi()  << endreq;
        } 

      //  Get all track collections attached to this ROI
      vector<const TrigInDetTrackCollection*> vectorOfTrackCollections;
      status = getFeatures( *TEiterator, vectorOfTrackCollections );

      if ( status != HLT::OK )
        {
          if (msgLvl()<=MSG::DEBUG) msg() << MSG::DEBUG << " Failed to get InDetTrackCollections " << endreq;
        }
      else
        if (msgLvl()<=MSG::DEBUG) msg() << MSG::DEBUG << " Got " << vectorOfTrackCollections.size() << " InDetTrackCollections" << endreq;
        
      // Select tracks from this ROI's collections
      selectTracks( vectorOfTrackCollections, mySelectedTrackCollection );
        
      countTracksPerTE       += m_TracksPerROI.back();
      countTracksPassedPerTE += m_TracksPerROIPass.back();
      countHiPTTracksPerTE   += m_HiPTTracksPerROI.back();
    }
        
  // Update per-TE counters
  m_TracksPerTE    .push_back( countTracksPerTE       );
  m_TracksPerTEPass.push_back( countTracksPassedPerTE );
  m_HiPTTracksPerTE.push_back( countHiPTTracksPerTE   );
}


void
T2VertexBeamSpotImpl::selectTracks( const vector<const TrigInDetTrackCollection*>& vectorOfTrackCollections,
                                    TrigInDetTrackCollection& mySelectedTrackCollection )
{
  T2Timer t( m_timer[allTrack] );
  //
  // Loop over all tracks in the given track collections
  //
  unsigned int countTracksPerROI       = 0;
  unsigned int countTracksPassedPerROI = 0;
  unsigned int countHiPTTracksPerROI   = 0;

  // Loop over all track collections
  for ( std::vector<const TrigInDetTrackCollection*>::const_iterator
          pTrackColl = vectorOfTrackCollections.begin();
        pTrackColl != vectorOfTrackCollections.end(); ++pTrackColl )
    {
      // Loop over all tracks
      for ( TrigInDetTrackCollection::const_iterator trackIter = (*pTrackColl)->begin();
            trackIter != (*pTrackColl)->end(); ++trackIter )
        { 
          const TrigInDetTrack& track = **trackIter;

          // Make sure that the event has tracks
          eventStage( hasTracks );
            
          // Add to counters
          countTracksPerROI++;
            
          // Track Algorithm selection and track parameter selection (pT, Eta, Z0, D0, chi2/NDF)
          //    Algorithm ID's: SiTrack = 1, IdScan = 2, L2StarA = 5, L2StarB = 6
          //    (from /offline/Trigger/TrigEvent/TrigInDetEvent/TrigInDetEvent/TrigInDetTrack.h)
               
          // Only fill histograms for the individual algorithms separately if it's requested (saves time) 
          if ( m_saveTrackAlgs )
            {
              if      ( track.algorithmId() == TrigInDetTrack::STRATEGY_B_ID )
                {
                  if (msgLvl()<=MSG::DEBUG) msg() << MSG::DEBUG << "Looking at L2StarB track: " << track << endreq;
                  m_L2StarB.push_back( track );
                }
              else if ( track.algorithmId() == TrigInDetTrack::STRATEGY_A_ID )
                {
                  if (msgLvl()<=MSG::DEBUG) msg() << MSG::DEBUG << "Looking at L2StarA track: " << track << endreq;
                  m_L2StarA.push_back( track );
                }
              else if ( track.algorithmId() == TrigInDetTrack::STRATEGY_F_ID )
                {
                  if (msgLvl()<=MSG::DEBUG) msg() << MSG::DEBUG << "Looking at L2StarF track: " << track << endreq;
                  m_L2StarF.push_back( track );
                }
              else if ( track.algorithmId() == TrigInDetTrack::SITRACKID )
                {
                  if (msgLvl()<=MSG::DEBUG) msg() << MSG::DEBUG << "Looking at SiTrack track: " << track << endreq;
                  m_SiTrack.push_back( track );
                }
              else if ( track.algorithmId() == TrigInDetTrack::IDSCANID )
                {
                  if (msgLvl()<=MSG::DEBUG) msg() << MSG::DEBUG << "Looking at IDScan track: " << track << endreq;
                  m_IdScan.push_back( track );
                }
            }
          
          // Save Track parameters for the tracks we actually want to use
          if ( track.algorithmId() == m_TrackAlgoId )
            {
              const T2Track myTrack( track );

              // Check for passed track
              if ( isGoodTrack( myTrack ) )
                {
                  // Add this track to the set used to find a vertex
                  mySelectedTrackCollection.push_back( *trackIter );
    
                  // Fill accepted track histograms
                  m_trackPass.push_back( track );
                                    
                  countTracksPassedPerROI++;

                  // Check for high-pT track
                  if ( myTrack.Pt() > m_trackSeedPt )
                    {
                      countHiPTTracksPerROI++;
                    }
                        
                }
              else
                if (msgLvl()<=MSG::DEBUG) msg() << MSG::DEBUG << "Track failed selection: " << track << endreq;
            } 
        }
    }
        
  // Update per-ROI counters
  m_TracksPerROI    .push_back( countTracksPerROI       );
  m_TracksPerROIPass.push_back( countTracksPassedPerROI );
  m_HiPTTracksPerROI.push_back( countHiPTTracksPerROI   );
}


void
T2VertexBeamSpotImpl::reconstructVertices( TrigInDetTrackCollection& mySelectedTrackCollection,
                                           TrigVertexCollection& myVertexCollection,
                                           DataVector< TrigVertexCollection >& mySplitVertexCollections )
{
  T2Timer t( m_timer[allVertex] );
     
  // Make collections for vertex splitting (unsorted)
  TrigInDetTrackCollection mySplitTrackCollection( mySelectedTrackCollection );

  // Sort tracks by track pT
  {
    T2Timer t( m_timer[trackSort] );
    mySelectedTrackCollection.sort( Beamspot::TrackPTSort() );
  }

  // Prepare a track clustering algorithm with the given parameters
  //  T2TrackClusterer trackClusterer( m_trackClusDZ, m_trackSeedPt, m_weightSeed, m_vtxNTrkMax );

  // Create clusters from the track collection until all its tracks are used
  while ( ! mySelectedTrackCollection.empty() )
    {
      // Report
      if (msgLvl()<=MSG::DEBUG)
        {
          msg() << MSG::DEBUG << "Number of tracks remaining = " << mySelectedTrackCollection.size() << endreq;
          msg() << MSG::DEBUG << "pT of first (seed) track (GeV) = "
                << abs( (*mySelectedTrackCollection.begin())->param()->pT() )/GeV << endreq;
        }
      
      // Cluster tracks in z around the first (highest-pT track) in the collection, which is taken
      // as the seed.
      {
        T2Timer t( m_timer[zCluster] );
        m_trackClusterer->cluster( mySelectedTrackCollection );
      }

      // Sanity check:
      if ( m_trackClusterer->cluster().size() + m_trackClusterer->unusedTracks().size()
           != mySelectedTrackCollection.size() )
        {
          msg() << MSG::ERROR << "Track clustering check sum failed: "
                << "cluster().size()=" << m_trackClusterer->cluster().size()
                << " + unusedTracks().size()=" << m_trackClusterer->unusedTracks().size()
                << " != mySelectedTrackCollection.size()=" << mySelectedTrackCollection.size()
                << endreq;
        }

      // Continue with the remaining tracks - still pT-sorted
      mySelectedTrackCollection = m_trackClusterer->unusedTracks();
      // This always uses at least one track (the seed), so we are sure to reduce the track
      // selection and terminate the loop

      // Make sure we have enough tracks in the cluster
      if ( m_trackClusterer->cluster().size() < m_totalNTrkMin )
        {
          if (msgLvl()<=MSG::DEBUG) msg() << MSG::DEBUG << "Not enough tracks in cluster!" << endreq;
          continue;
        }
      
      // Event has a good cluster
      eventStage( hasCluster );
      
      m_NClusters++;

      // Monitor cluster properties
      m_clusterZ            .push_back( m_trackClusterer->seedZ0() );
      m_clusterNTracks      .push_back( m_trackClusterer->cluster().size() );
      m_clusterNTracksUnused.push_back( m_trackClusterer->unusedTracks().size() );

      for ( TrigInDetTrackCollection::const_iterator track = m_trackClusterer->cluster().begin();
            track !=  m_trackClusterer->cluster().end(); ++track )
        {
          // Monitor properties of tracks inside the cluster
          const double deltaZ0 = (*track)->param()->z0() - m_trackClusterer->seedZ0();
          const double z0Error = (*track)->param()->ez0();
          const double z0Pull  = ( z0Error > 0. ) ? deltaZ0 / z0Error : 0.;
          m_clusterDeltaZ0.push_back( deltaZ0 );
          m_clusterZ0Pull .push_back( z0Pull );
        }
 
      if (msgLvl()<=MSG::DEBUG) msg() << MSG::DEBUG << "Number of tracks remaining after cluster #(" << m_NClusters << ") = "
                                      << mySelectedTrackCollection.size() << endreq;
      
      if (msgLvl()<=MSG::DEBUG)
        {
          msg() << MSG::DEBUG << "Total number of tracks to fit    = " << m_trackClusterer->cluster().size() << endreq;
          msg() << MSG::DEBUG << "Average Z position (from trk Z0) = " << m_trackClusterer->seedZ0() << endreq;
          msg() << MSG::DEBUG << "Fitting tracks" << endreq;
        }
      
      // Fit a primary vertex to this cluster around its seed track
      TrigVertex* primaryVertex = 0;
      {
        T2Timer t( m_timer[vertexFit] );
        primaryVertex = m_primaryVertexFitter->fit( &m_trackClusterer->cluster(), m_trackClusterer->seedZ0() );
      }
      
      // Check to see if the fit succeeded / converged
      if ( ! primaryVertex )
        { 
          if (msgLvl()<=MSG::DEBUG) msg() << MSG::DEBUG << "Vertex fit failed" << endreq; 
          continue; 
        }
      
      if ( msgLvl()<=MSG::DEBUG ) msg() << MSG::DEBUG << "Successful vertex fit" << endreq;

      // Event has a vertex!
      eventStage( hasVertex );
      
      // Update vertex counter
      m_Nvtx++;
      
      // Extract beam spot parameters
      const T2BeamSpot beamSpot( m_beamCondSvc );

      if (msgLvl()<=MSG::DEBUG) msg() << MSG::DEBUG << "Beamspot from BeamCondSvc: " << beamSpot << endreq;

      // Compute vertex parameters with respect to beam spot position
      // Measure ctor time which includes calculating cumulative chi 2 probability
      //        T2Timer t( m_timer[cumChi2] );
      const T2Vertex myVertex( *primaryVertex, beamSpot, m_trackClusterer->seedZ0() );

      // Fill vertex histograms
      m_vertex.push_back( myVertex );

      if (msgLvl()<=MSG::DEBUG) msg() << MSG::VERBOSE << "Vertex fit: " << '\n' << myVertex << endreq;
      
      // Query for vertex selection
      const bool passVertex = isGoodVertex( myVertex );

      if ( ! passVertex )
        { 
          if (msgLvl()<=MSG::DEBUG) msg() << MSG::DEBUG << "Vertex failed selection" << endreq; 

          // Release the owned vertex
          delete primaryVertex;
          primaryVertex = 0;
          continue; 
        }

      // Tag vertex according to our first track's algorithm id
      primaryVertex->algorithmId( getPrimaryAlgoId( *mySelectedTrackCollection.front() ) );

      // Add primary vertex to collection
      myVertexCollection.push_back( primaryVertex ); // passes ownership to vertex collection

      // Event has good vertex
      eventStage( hasGoodVertex );
      
      // Fill accepted vertex histograms 
      m_vertexPass.push_back( myVertex );

      // Update good vertex counter
      m_NvtxPass++;    

      // Learn something more about the cluster that ended up in this vertex
      {
        const double deltaVertexZ =  m_trackClusterer->seedZ0() - myVertex.Z();
        m_clusterDeltaVertexZ.push_back( deltaVertexZ );
        //        m_clusterVertexZPull .push_back( m_trackClusterer->seedZ0() - myVertex.Z() );
      }
      //   clusterNTracksInVertex

      // If the vertex is good, and splits are requested, and we have enough tracks to split, split them
      if ( passVertex && m_nSplitVertices > 1 )
        {
          if ( m_splitWholeCluster )
            {
              // Alternative 1: Split the entire cluster of track into two
              mySplitTrackCollection = m_trackClusterer->cluster();
            }
          else
            {
              // Alternative 2: Split only the tracks that were successfully fit to a vertex
              // Ouch! tracks() returns a std:list while fit() expects a DataVector
              //          mySplitTrackCollection = primaryVertex->tracks();
              // No such luck: DataVector only defines default ctor
              //          mySplitTrackCollection = TrigInDetTrackCollection( primaryVertex->tracks()->begin(), primaryVertex->tracks()->end() );
              //          mySplitTrackCollection.clear();
              //          std::copy( primaryVertex->tracks()->begin(), primaryVertex->tracks()->end(), std::back_inserter( mySplitTrackCollection ) );
              mySplitTrackCollection.clear( SG::VIEW_ELEMENTS );
              mySplitTrackCollection.reserve( primaryVertex->tracks()->size() );
              for ( TrackInVertexList::const_iterator track = primaryVertex->tracks()->begin();
                    track != primaryVertex->tracks()->end(); ++track )
                mySplitTrackCollection.push_back( const_cast<TrigInDetTrack*>( *track ) );
            }

          if ( mySplitTrackCollection.size() >= m_nSplitVertices * m_totalNTrkMin )
            {
              // Split, optinally re-cluster, and fit separate vertices
              reconstructSplitVertices( mySplitTrackCollection, mySplitVertexCollections );
            }
          // Alternative 3: Split all the tracks and iterate with the remaining tracks
          //          mySplitTrackCollection = mySelectedTrackCollection;
        }

      if (msgLvl()<=MSG::DEBUG) msg() << MSG::DEBUG << "Vertex passed" << endreq;

      // Now check if this vertex passed the higher multiplicity cut to be used for per-BCID measurements
      const bool passVertexBCID = isGoodVertexBCID( myVertex );

      if ( passVertexBCID )
        {
          // Fill accepted per-BCID vertex histograms 
          m_vertexPassBCID.push_back( myVertex );

          // Update good per-BCID vertex counter
          m_NvtxPassBCID++;
        }
    }
}


void
T2VertexBeamSpotImpl::reconstructSplitVertices( TrigInDetTrackCollection& myFullTrackCollection,
                                                DataVector< TrigVertexCollection >& mySplitVertexCollections )
{
  T2Timer t( m_timer[allVertexSplit] );

  // Sort the tracks in phi for splitting in the most independent, uniform variable
  {
    T2Timer t( m_timer[trackSortSplit] );
    myFullTrackCollection.sort( Beamspot::TrackPhiSort() );
  }

  // Split the track collection (typically into halves)
  // This returns m_nSplitVertices (ideally) or fewer (if clustering fails) track collections
  m_trackManager.ResetKey( m_EventID % 2 - 1 );
  vector< TrigInDetTrackCollection > splitTrackCollections = m_trackManager.split( myFullTrackCollection );

  // Add a new track collection for the split vertices corresponding to this primary vertex
  // There can be anywhere between zero and m_nSplitVertices entries in the collection
  TrigVertexCollection* splitVertices = new TrigVertexCollection;
  mySplitVertexCollections.push_back( splitVertices ); // passes ownership

  // Loop over the split track collections to perform clustering and vertex fitting

  for ( vector< TrigInDetTrackCollection >::iterator tracks = splitTrackCollections.begin(); 
        tracks != splitTrackCollections.end(); ++tracks )
    {
      TrigVertex* splitVertex = 0;

      if ( m_reclusterSplit )
        {
          // Sort the tracks in pT for clustering around the highest-pT seed
          {
            T2Timer t( m_timer[trackSortSplit] );
            tracks->sort( Beamspot::TrackPTSort() );
          }

          // Cluster in z
          {
            T2Timer t( m_timer[zClusterSplit] );
            m_trackClusterer->cluster( *tracks );
          }

          // Check for a good cluster
          if ( m_trackClusterer->cluster().size() < m_totalNTrkMin )
            {
              continue;
            }

          if (msgLvl()<=MSG::DEBUG) msg() << MSG::DEBUG << "Checking splitCluster[" << tracks - splitTrackCollections.begin() << "]"
                                          << " which has " << m_trackClusterer->cluster().size() << " tracks" << endreq;

          // Fit a vertex to this split track cluster
          {
            T2Timer t( m_timer[vertexFitSplit] );
            splitVertex = m_primaryVertexFitter->fit( &m_trackClusterer->cluster(), m_trackClusterer->seedZ0() );
          }
        }
      else
        {
          // Alternatively: Fit a vertex to the unclustered split track
          // collection, using the seed Z0 from the primary vertex
          {
            T2Timer t( m_timer[vertexFitSplit] );
            splitVertex = m_primaryVertexFitter->fit( &(*tracks)                  , m_trackClusterer->seedZ0() );
          }
        }

      if ( splitVertex )
        {
          if (msgLvl()<=MSG::DEBUG) msg() << MSG::DEBUG << "Reconstructed a split vertex" << endreq;

          // Tag vertex according to our first track's algorithm id
          splitVertex->algorithmId( getSplitAlgoId( *tracks->front() ) );

          // Add split vertex to collection
          splitVertices->push_back( splitVertex );  // passes ownership to split vertex collection
        }
    }

  // Monitor split vertex distributions (if we found all of them)
  if ( m_nSplitVertices > 1 && splitVertices->size() == m_nSplitVertices )
    {
      if (msgLvl()<=MSG::DEBUG) msg() << MSG::DEBUG << "Split vertexing is ON."
                                      << "Attempting to split N = " << m_nSplitVertices << " vertices. " << endreq;
         
      // Store information on the first two vertices
      // There shouldn't be more, unless it's for systematic studies anyway
      const T2SplitVertex mySplitVertex( *(*splitVertices)[0], *(*splitVertices)[1] );

      // Fill split vertex histograms
      m_splitVertexPass.push_back( mySplitVertex );
    }
}


void
T2VertexBeamSpotImpl::createOutputTEs( TrigVertexCollection& myVertexCollection,
                                       DataVector< TrigVertexCollection >& mySplitVertexCollections,
                                       unsigned int type_out )
{
  T2Timer t( m_timer[allOutput] );

  //--------------------------------------------------
  //  create TEVec for creating output TE if necessary
  //--------------------------------------------------
  HLT::TEVec allTEs;

  // Save all events, or only those events which pass the Npv cuts (if activated)!
  if ( m_activateAllTE )
    {
      // Create an output TE seeded by the inputs
      HLT::TriggerElement* outputTE = config()->getNavigation()->addNode(allTEs, type_out);
      outputTE->setActiveState(true);
    }
  
  // Save event if 
  if ( ! m_activateAllTE && m_activateTE && m_NvtxPass > 0 )
    {
      // FIXME: Why do we need one TE per vertex?
      for ( unsigned i=0; i < m_NvtxPass; ++i )
        {
          // Create an output TE seeded by the inputs
          HLT::TriggerElement* outputTE = config()->getNavigation()->addNode(allTEs, type_out);
          outputTE->setActiveState(true);
        }
    }

  // Save all events, or only those events which pass the Npv cuts (if activated)
  if ( ! m_activateAllTE && m_activateSI &&
       m_minNpvTrigger <= m_NvtxPass && m_NvtxPass <= m_maxNpvTrigger )
    {
      // Create an output TE seeded by the inputs
      HLT::TriggerElement* outputTE = config()->getNavigation()->addNode(allTEs, type_out);
      outputTE->setActiveState(true);
    }
  
  // Save feature to output TE:
  //  if ( m_attachVertices && m_NvtxPass > 0 )
  if ( m_attachVertices )
    {
      if (msgLvl()<=MSG::DEBUG) msg() << MSG::DEBUG << "Saving primary vertex collection" << endreq; 

      // Save primary vertex collection
      HLT::TriggerElement* outputTE = config()->getNavigation()->addNode(allTEs, type_out);
      TrigVertexCollection* newColl = new TrigVertexCollection;
      newColl->clear( SG::VIEW_ELEMENTS );
      newColl->swap( myVertexCollection ); // swaps ownership
      // Strip list of tracks from the vertex so we don't persist them.
      stripTrackList( *newColl );
      const HLT::ErrorCode status = attachFeature( outputTE, newColl, m_vertexCollName );
      if ( status != HLT::OK )
        {
          msg() << MSG::ERROR << "Failed to attach primary vertex collection to output trigger element" << endreq;
        }

      if ( m_attachSplitVertices )
        {
          if (msgLvl()<=MSG::DEBUG) msg() << MSG::DEBUG << "Saving split vertex collections" << endreq; 

          // For each primary vertex we find, we save a group of
          // (typically two) split vertices in a separate collection
          // Thus we have:
          //
          //    VertexCollection    vector<VertexCollection>
          //    (primary vertex)    (split vertex collection)      (split vertex 1)  (split vertex 2) ...
          //  --------------------------------------------------------------------------------------------
          //    |      PV1           SplitVertexCollection_1  --->      SP1_1             SP1_2       ... 
          //    |      PV2           SplitVertexCollection_2  --->      SP2_1             SP2_2            
          //    v      PV3           SplitVertexCollection_3  --->      SP3_1             SP3_2            
          //           ...
          //

          HLT::TriggerElement* outputTE = config()->getNavigation()->addNode(allTEs, type_out);

          // FIXME: This is what we would like to be able to do:
          // Save vector of split vertex collections
          //          TrigVertexCollectionVector* newCollVec = new TrigVertexCollectionVector;
          //          newCollVec->clear( SG::VIEW_ELEMENTS );
          //          newCollVec->swap( mySplitVertexCollections ); // swaps ownership
          //          const HLT::ErrorCode status = attachFeature( outputTE, newCollVec, m_vertexCollName );
          //          if ( status != HLT::OK )
          //            {
          //              msg() << MSG::WARNING << "Write of Beamspot feature into outputTE failed" << endreq;
          //            }
          // FIXME: The above needs a new (container) object to be registered.
          // FIXME: For the time being we store the individual collections one by one.

          // Save split vertex collections
          for ( unsigned i=0; i < mySplitVertexCollections.size(); ++i )
            {
              ostringstream collName;
              collName << m_vertexCollName << "_" << i;
              TrigVertexCollection* newColl;
              // Wrestle back ownership of this collection from DataVector
              mySplitVertexCollections.swapElement( i, 0, newColl );
              // Strip list of tracks from the vertex so we don't persist them
              stripTrackList( *newColl );
              const HLT::ErrorCode status = attachFeature( outputTE, newColl, collName.str() );
              if ( status != HLT::OK )
                {
                  msg() << MSG::ERROR << "Failed to attach split vertex collection " << i
                        << " of " << mySplitVertexCollections.size() << " to output trigger element" << endreq;
                }
            }
        }
    }
}


//
// Selection criteria
//

bool
T2VertexBeamSpotImpl::isGoodTrack( const T2Track& track ) const
{
  // Vertex quality cuts
  return
    (
     abs( track.Pt() )  >= m_minTrackPt       &&
     track.SiHits()     >= m_minSiHits        &&
     track.PIXHits()    >= m_minPIXHits       &&
     track.SCTHits()    >= m_minSCTHits       &&
     track.TRTHits()    >= m_minTRTHits       &&
     track.NDF()        >= m_minTrackNDF      &&
     abs( track.D0() )  <= m_maxTrackD0       &&
     abs( track.Z0() )  <= m_maxTrackZ0       &&
     track.D0err()      <= m_maxTrackD0err    &&
     track.Z0err()      <= m_maxTrackZ0err    &&
     abs( track.Eta() ) <= m_maxTrackEta      &&
     track.Qual()       >= m_minTrackQual     && 
     track.Qual()       <= m_maxTrackQual     &&
     track.Chi2Prob()   >= m_minTrackChi2Prob &&
     true
     );
}


bool
T2VertexBeamSpotImpl::isGoodVertex( const T2Vertex& vertex ) const
{
  // Vertex quality cuts
  return
    (
     vertex.NTrks()       >= m_vtxNTrkMin     &&  // FIXME: harmonize spelling
     abs( vertex.X() )    <= m_vtxXposMax     &&
     abs( vertex.Y() )    <= m_vtxYposMax     &&
     abs( vertex.Z() )    <= m_vtxZposMax     &&
     abs( vertex.Xerr() ) <= m_vtxXerrMax     &&  // FIXME: do we have negative errors?
     abs( vertex.Yerr() ) <= m_vtxYerrMax     &&
     abs( vertex.Zerr() ) <= m_vtxZerrMax     &&
     vertex.Mass()        >= m_vtxMassMin     &&
     vertex.SumPt()       >= m_vtxSumPtMin    &&
     vertex.Qual()        >= m_vtxQualMin     &&
     vertex.Qual()        <= m_vtxQualMax     &&
     vertex.Chi2Prob()    >= m_vtxChi2ProbMin &&
     true
     );
}


bool
T2VertexBeamSpotImpl::isGoodVertexBCID( const T2Vertex& vertex ) const
{
  // Track quality cuts
  return
    (
     vertex.NTrks()       >= m_vtxBCIDNTrkMin &&  // FIXME: harmonize spelling
     true
     );
}

//
// Helpers
//

void
T2VertexBeamSpotImpl::resetMonitoredVariables()
{
  // Event
  m_RunID            = 0;
  m_EventID          = 0;
  m_LumiBlock        = 0;
  m_bcid             = 0;

  // Tracks
  m_TotalTracks      = 0;
  m_TotalTracksPass  = 0;
  m_TotalHiPTTracks  = 0;
  m_TotalROIs        = 0;
  m_TotalTEs         = 0;

  // Vertices
  m_NClusters        = 0;
  m_Nvtx             = 0;
  m_NvtxPass         = 0;
  m_NvtxPassBCID     = 0;
}
  

bool
T2VertexBeamSpotImpl::eventStage( Statistics stage )
{
  if ( ! m_eventStageFlag[ stage ] )
    {
      m_eventStageFlag[ stage ] = true;
      m_eventStage.push_back( stage );
      return true;
    }
  return false;
}
