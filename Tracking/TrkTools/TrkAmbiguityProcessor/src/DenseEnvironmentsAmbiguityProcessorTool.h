/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef DenseEnvironmentsAmbiguityProcessorTool_H
#define DenseEnvironmentsAmbiguityProcessorTool_H

// turn on debugging ? uncomment this
// #define SIMPLEAMBIGPROCDEBUGCODE

#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/IIncidentSvc.h"

#include "TrkToolInterfaces/ITrackAmbiguityProcessorTool.h"
#include "TrkEventPrimitives/TrackScore.h"
#include "TrkFitterInterfaces/ITrackFitter.h"
#include "TrkToolInterfaces/IAmbiTrackSelectionTool.h"
#include "InDetPrepRawData/PixelGangedClusterAmbiguities.h"


//need to include the following, since its a typedef and can't be forward declared.
#include "TrkTrack/TrackCollection.h"
#include "TrkTrack/TrackSeedMap.h"
#include "TrkParameters/TrackParameters.h"

#include <map>
#include <set>
#include <vector>
#include <functional>


#ifdef SIMPLEAMBIGPROCDEBUGCODE
  // --------------- DEBUG CODE
  #include "HepMC/GenEvent.h"
  #include "TrkTruthData/PRD_MultiTruthCollection.h"
  typedef std::map<const Trk::Track*, const Trk::Track*> TrackCollectionConnection;
  CLASS_DEF( TrackCollectionConnection , 148639440 , 1 )
#endif



class AtlasDetectorID;
class PixelID;

namespace InDet{
  class IPixelClusterSplitProbTool;
  class PixelCluster;
  class SCT_Cluster;
}

namespace Trk {

  class ITrackScoringTool;
  class IPRD_AssociationTool;
  class ITruthToTrack;

  class DenseEnvironmentsAmbiguityProcessorTool : public AthAlgTool, 
                                                  virtual public IIncidentListener, 
                                                  virtual public ITrackAmbiguityProcessorTool
    {
    public:

      // public types
      typedef std::pair< const Track*, const bool > TrackFittedPair;
      typedef std::multimap< TrackScore, TrackFittedPair > TrackScoreMap;
    
      typedef std::set<const PrepRawData*> PrdSignature;
      typedef std::set<PrdSignature> PrdSignatureSet;

      // default methods
      DenseEnvironmentsAmbiguityProcessorTool(const std::string&,const std::string&,const IInterface*);
      virtual ~DenseEnvironmentsAmbiguityProcessorTool ();
      virtual StatusCode initialize();
      virtual StatusCode finalize  ();

      /**Returns a processed TrackCollection from the passed 'tracks'
     @param tracks collection of tracks which will have ambiguities resolved. Will not be 
     modified.
     The tracks will be refitted if no fitQuality is given at input.
     @return new collections of tracks, with ambiguities resolved. Ownership is passed on 
     (i.e. client handles deletion)*/
      virtual TrackCollection*  process(const TrackCollection* tracks);
    
      /** statistics output */
      virtual void statistics();

      /** handle for incident service */
      void handle(const Incident& inc) ;

    private:
      
      void reset();
    
      /**Add passed TrackCollection, and Trk::PrepRawData from tracks to caches
         @param tracks the TrackCollection is looped over, 
         and each Trk::Track is added to the various caches. 
         The Trk::PrepRawData from each Trk::Track are added to the IPRD_AssociationTool*/
      void addNewTracks(const TrackCollection* tracks);

      void addTrack(const Track* track, const bool fitted);

      void solveTracks();
    
      /** add subtrack to map */
      void addSubTrack( const std::vector<const TrackStateOnSurface*>& tsos);

      /** refit track */
      void refitTrack( const Trk::Track* track);

      /** refit PRDs */
      const Track* refitPrds( const Track* track);

      /** refit ROTs corresponding to PRDs */
      const Track* refitRots( const Track* track);

      /** print out tracks and their scores for debugging*/
      void dumpTracks(const TrackCollection& tracks);
      
      
      /**  Find SiS Tracks that share hits in the track score map*/
      void overlapppingTracks();
     
      /**  Update pixel split information based using the fitted track*/    
      void updatePixelSplitInformation(std::map< const InDet::PixelCluster*, const Trk::TrackParameters* >& setOfClustersOnTrack);
      void updatePixelSplitInformationForCluster(const std::pair<const InDet::PixelCluster* const,
                                                                 const Trk::TrackParameters*> & clusterTrkPara );
      
 
      void updateSCT_SplitInformation(std::map< const InDet::SCT_Cluster*, const Trk::TrackParameters* >& setOfClustersOnTrack);

      
    
      // private data members

      /** brem recovery mode with brem fit ? */
      bool  m_tryBremFit;
      bool  m_caloSeededBrem;
      float m_pTminBrem;

      /** by default drop double tracks before refit*/
      bool m_dropDouble;

      /** by default tracks at input get refitted */
      bool m_forceRefit;

      /** by default refit tracks from PRD */
      bool m_refitPrds;

      /** suppress Hole Search */ 
      bool m_suppressHoleSearch;

      /** suppress Track Fit */ 
      bool m_suppressTrackFit;

      /** control material effects (0=non-interacting, 1=pion, 2=electron, 3=muon, 4=pion) read in as an integer 
      read in as an integer and convert to particle hypothesis */
      int m_matEffects;
      Trk::ParticleHypothesis m_particleHypothesis;
   
      /** IncidentSvc to catch begining of event and end of event */   
      ServiceHandle<IIncidentSvc>                           m_incidentSvc;    
      
      /**Scoring tool
         This tool is used to 'score' the tracks, i.e. to quantify what a good track is.
         @todo The actual tool that is used should be configured through job options*/
      ToolHandle<ITrackScoringTool> m_scoringTool;


      /** refitting tool - used to refit tracks once shared hits are removed. 
          Refitting tool used is configured via jobOptions.*/
      ToolHandle<ITrackFitter> m_fitterTool;
        
      /** selection tool - here the decision which hits remain on a track and
          which are removed are made */
      ToolHandle<IAmbiTrackSelectionTool> m_selectionTool;

      /** recalculate split prob tool **/
      ToolHandle<InDet::IPixelClusterSplitProbTool>       m_splitProbTool; 

      /**Association tool - used to work out which (if any) PRDs are shared between 
       tracks */ 
      ToolHandle<Trk::IPRD_AssociationTool> m_assoTool;
      
      /** unsorted container of track and track scores.*/
      TrackScoreMap m_trackScoreTrackMap;
    
      /** signature map to drop double track. */
      PrdSignatureSet m_prdSigSet;

      /**Tracks that will be passed out of AmbiProcessor. 
         Recreated anew each time process() is called*/ 
      TrackCollection* m_finalTracks;

      /**NN split sprob cut for 2 particle clusters */      
      float m_sharedProbCut;

      /**NN split sprob cut for 3 particle clusters */      
      float m_sharedProbCut2;

      /** monitoring statistics */
      int m_Nevents;
      std::vector<int> m_Ncandidates, m_NcandScoreZero, m_NcandDouble,
      m_NscoreOk,m_NscoreZeroBremRefit,m_NscoreZeroBremRefitFailed,m_NscoreZeroBremRefitScoreZero,m_NscoreZero,
      m_Naccepted,m_NsubTrack,m_NnoSubTrack,m_NacceptedBrem,
      m_NbremFits,m_Nfits,m_NrecoveryBremFits,m_NgoodFits,m_NfailedFits;
      /** internal monitoring: categories for counting different types of extension results*/
      enum StatIndex {iAll = 0, iBarrel = 1, iTransi = 2, iEndcap = 3};
      std::vector<float>  m_etabounds;           //!< eta intervals for internal monitoring

      /** helper for monitoring and validation: does success/failure counting */
      void increment_by_eta(std::vector<int>&,const Track*,bool=true);


      mutable InDet::PixelGangedClusterAmbiguities*       m_splitClusterMap;      //!< the actual split map         
      std::string                                         m_splitClusterMapName; //!< split cluster ambiguity map


      const PixelID *m_PixelHelper;

//==================================================================================================
//
//   FROM HERE EVERYTHING IS DEBUGGING CODE !!!
//
//==================================================================================================

#ifdef SIMPLEAMBIGPROCDEBUGCODE
    const PRD_MultiTruthCollection   * m_truthPIX;
    const PRD_MultiTruthCollection   * m_truthSCT;  
    std::string                        m_truth_locationPixel    ;
    std::string                        m_truth_locationSCT      ;  
#endif


#ifdef SIMPLEAMBIGPROCDEBUGCODE
//==================================================================================================
// PART 2 : Output statistics
//==================================================================================================

    private:

      std::set<const Trk::Track*> m_trueTracks;
      std::map<const Trk::Track*, const Trk::Track*> m_trackHistory;

      void findTrueTracks(const TrackCollection* recTracks);
      void keepTrackOfTracks(const Trk::Track* oldTrack, const Trk::Track* newTrack);
      void produceInputOutputConnection();
 
      std::string m_resolvedTrackConnection;
      std::string m_truthCollection;
      int n_trueFitFails;
      int n_fitFails;
      int numOutliersDiff;
      int numOutliersBefore;
      int numOutliersAfter;
      int numFirstFitLost;
      int numSecondFitLost;
      int numSharedTruth;
      int truthBefore;
      int truthAfter;

      bool isSharedTrack( const Track* track);
      bool isTrueTrack( const Track* track);

      const PixelID* m_pixelId;

      void addTrackToMap(Trk::Track* Tr);
      void findSharedTrueTracks(const TrackCollection* recTracks);    
      void prdTruth(const Track* track);
      void tsosTruth(const Track* track);
      const Track* origTrack( const Track* track);

      std::map<const Trk::PrepRawData*, const Trk::Track*> m_tracksShared;

      StatusCode getBremTruth();
      double originalMomentum( const HepMC::GenEvent* genEvent );
      double momentumLostByBrem( const HepMC::GenEvent* genEvent ) const;
      const std::vector<double> fractionOfIncidentMomentumLostPerVertex( const HepMC::GenEvent* genEvent ) const;
      const std::vector<Amg::Vector3D> positionsOfBremVertices( const HepMC::GenEvent* genEvent ) const;
      bool vertexAssociatedWithOriginalTrack( HepMC::GenVertex* genVertex) const;

      std::string                        m_generatedEventCollectionName; 
      Trk::ITruthToTrack*                m_truthToTrack         ;
      const PRD_MultiTruthCollection   * m_truthTRT               ;    
      std::string                        m_truth_locationTRT      ;

#endif // DebugCode

  };


} //end ns


#endif // TrackAmbiguityProcessorTool_H

