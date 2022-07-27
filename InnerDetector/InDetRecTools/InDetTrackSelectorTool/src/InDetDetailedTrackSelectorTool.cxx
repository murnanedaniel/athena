/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "InDetTrackSelectorTool/InDetDetailedTrackSelectorTool.h"
// forward declares
#include "TrkToolInterfaces/ITrackSummaryTool.h"
#include "TrkExInterfaces/IExtrapolator.h"
#include "VxVertex/Vertex.h"
#include "VxVertex/RecVertex.h"
#include "TrkEventPrimitives/FitQuality.h"
#include "TrkTrack/Track.h"
#include "TrkTrackSummary/TrackSummary.h"
#include "TrkParticleBase/TrackParticleBase.h"
#include "InDetRecToolInterfaces/ITrtDriftCircleCutTool.h"
#include "InDetRecToolInterfaces/IInDetTestPixelLayerTool.h"

#include "xAODTracking/TrackParticle.h"
#include "xAODTracking/Vertex.h"


// normal includes
#include "CLHEP/GenericFunctions/CumulativeChiSquare.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "TrkSurfaces/PerigeeSurface.h"
#include "GeoPrimitives/GeoPrimitives.h"
#include "EventPrimitives/EventPrimitives.h"
#include <cmath>

using CLHEP::GeV;
using CLHEP::mm;

namespace InDet
{

  // ---------------------------------------------------------------------
  InDetDetailedTrackSelectorTool::InDetDetailedTrackSelectorTool(const std::string& t, const std::string& n, const IInterface*  p)
    : AthAlgTool(t,n,p)
    , m_trackSumTool("Trk::TrackSummaryTool", this)
    , m_extrapolator("Trk::Extrapolator", this)
    , m_trtDCTool("InDet::InDetTrtDriftCircleCutTool", this)
    , m_inDetTestPixelLayerTool("", this)
    , m_trackSumToolAvailable(true)
    , m_usePtDependentCuts(false)
    , m_useEventInfoBs(false)


  {
    declareInterface<ITrackSelectorTool>(this);
    declareProperty("pTMin"            , m_pTMin             = 1.*GeV);
    declareProperty("pMin"             , m_pMin              = 0.);
    declareProperty("IPd0Max"          , m_IPd0Max           = 2.*mm);
    declareProperty("IPz0Max"          , m_IPz0Max           = 1.5*mm);
    declareProperty("z0Max"            , m_z0Max             = 9999.*mm);
    declareProperty("sigIPd0Max"       , m_sigIPd0Max        = 999.*mm);
    declareProperty("sigIPz0Max"       , m_sigIPz0Max        = 999.*mm);
    declareProperty("d0significanceMax", m_d0significanceMax = -1.);
    declareProperty("z0significanceMax", m_z0significanceMax = -1.);
    declareProperty("etaMax"           , m_etaMax            = 9999.);
    
    declareProperty("useTrackSummaryInfo", m_useTrackSummaryInfo = true);
    
    declareProperty("nHitBLayer"       , m_nHitBLayer        = 1);
    declareProperty("nHitPix"          , m_nHitPix           = 2);
    declareProperty("nHitBLayerPlusPix", m_nHitBLayerPlusPix = 0);
    declareProperty("nHitSct"          , m_nHitSct           = 0);
    declareProperty("nHitSi"           , m_nHitSi            = 7);
    
    declareProperty("nHitPixPhysical", m_nHitPixPhysical = 0);
    declareProperty("nHitSiPhysical",  m_nHitSiPhysical  = 3);
    
    declareProperty("nHitTrt"                            , m_nHitTrt                          = 0);
    declareProperty("nHitTrtPlusOutliers"                , m_nHitTrtPlusOutliers              = 0);
    declareProperty("nHitTrtHighE"                       , m_nHitTrtHighE                     = 0);
    declareProperty("nHitTrtPlusOutliersHighE"           , m_nHitTrtPlusOutliersHighE         = 0);
    declareProperty("nHitTrtHighEFractionMax"            , m_nHitTrtHighEFraction             = 999);
    declareProperty("nHitTrtHighEFractionWithOutliersMax", m_nHitTrtHighEFractionWithOutliers = 999);
    
    declareProperty("useSharedHitInfo", m_useSharedHitInfo = false);
    declareProperty("nSharedBLayer"   , m_nSharedBLayer    = 0);
    declareProperty("nSharedPix"      , m_nSharedPix       = 0);
    declareProperty("nSharedSct"      , m_nSharedSct       = 1);
    declareProperty("nSharedSi"       , m_nSharedSi        = 999);
    
    declareProperty("nHoles"      , m_nHoles       = 999);
    declareProperty("nDoubleHoles", m_nDoubleHoles = 999);
    declareProperty("nHolesPixel" , m_nHolesPix    = 999);
    declareProperty("nHolesSct"   , m_nHolesSct    = 999);
    
    declareProperty("useTrackQualityInfo", m_useTrackQualityInfo = true);
    declareProperty("fitChi2"            , m_fitChi2             = 99999.);
    declareProperty("fitProb"            , m_fitProb             = -1.);
    declareProperty("fitChi2OnNdfMax"    , m_fitChi2OnNdfMax     = 999.);
    declareProperty("TrtMaxEtaAcceptance", m_TrtMaxEtaAcceptance = 999.);
    
    declareProperty("usePreselectionCuts", m_usePreselectionCuts = false);
    declareProperty("d0MaxPreselection"  , m_d0MaxPreselection   = 10.);
    
    declareProperty("useEtaDepententMinHitTrt", m_useEtaDepententMinHitTrt = false);
    declareProperty("scaleMinHitTrt"          , m_scaleMinHitTrt           = 1.);
    declareProperty("addToMinHitTrt"          , m_addToMinHitTrt           = 0);
    
    declareProperty("useEtaDepententMinHitTrtWithOutliers", m_useEtaDepententMinHitTrtWithOutliers = false);
    declareProperty("scaleMinHitTrtWithOutliers"          , m_scaleMinHitTrtWithOutliers           = 1.);
    declareProperty("addToMinHitTrtWithOutliers"          , m_addToMinHitTrtWithOutliers           = 0);
    
    declareProperty("TrackSummaryTool"   , m_trackSumTool);
    declareProperty("Extrapolator"       , m_extrapolator);
    declareProperty("TrtDCCutTool"       , m_trtDCTool);
    declareProperty("InDetTestPixelLayerTool", m_inDetTestPixelLayerTool);
    
   
    declareProperty("UsePtDependentCuts", m_usePtDependentCuts = false);  
    declareProperty("PtBenchmarks"      , m_ptBenchmarks);  
    declareProperty("SCTCutValues"      , m_nSCTValues);
    declareProperty("UseEventInfoBS"    , m_useEventInfoBs);
  }

  Trk::Vertex* InDetDetailedTrackSelectorTool::getBeamSpot(const EventContext& ctx) const
  {
    if(m_useEventInfoBs){
      SG::ReadHandle<xAOD::EventInfo> evt(m_eventInfo_key, ctx);
      if (evt.isValid()) {
        InDet::BeamSpotData temp(evt->beamStatus(), evt->beamPosX(), evt->beamPosY(), evt->beamPosZ(),
                                 evt->beamPosSigmaX(), evt->beamPosSigmaY(), evt->beamPosSigmaZ(),
                                 evt->beamTiltXZ(), evt->beamTiltYZ(), evt->beamPosSigmaXY());
        return new Trk::RecVertex(temp.beamVtx());
      } else {
        ATH_MSG_WARNING( " Cannot get beamSpot center from xAOD::EventInfo. Using (0,0,0)... " );
        return new Trk::Vertex(Amg::Vector3D(0,0,0));
      }
    }else{
      SG::ReadCondHandle<InDet::BeamSpotData> beamSpotHandle { m_beamSpotKey, ctx };
      if (beamSpotHandle.isValid()) {
        return new Trk::RecVertex(beamSpotHandle->beamVtx());
      } else {
        ATH_MSG_WARNING( " Cannot get beamSpot center from BeamSpotData. Using (0,0,0)... " );
        return new Trk::Vertex(Amg::Vector3D(0,0,0));
      }
    }
  }

  // ---------------------------------------------------------------------
  InDetDetailedTrackSelectorTool::~InDetDetailedTrackSelectorTool()
  = default;
  
  // ---------------------------------------------------------------------
  StatusCode  
  InDetDetailedTrackSelectorTool::initialize(){
    m_trackSumToolAvailable=false;
    if (!m_trackSumTool.empty()) {  
	    if(m_trackSumTool.retrieve().isFailure()){
	      ATH_MSG_INFO(" Unable to retrieve. OK if running on AOD. "<<m_trackSumTool);
	    }else{
	      ATH_MSG_INFO("Track summary tool retrieved");
	      m_trackSumToolAvailable=true;
	    }
    } 
    ATH_CHECK( m_extrapolator.retrieve() );
    ATH_CHECK(m_beamSpotKey.initialize(!m_useEventInfoBs));
    ATH_CHECK(m_eventInfo_key.initialize(m_useEventInfoBs));
    if (m_useEtaDepententMinHitTrt || m_useEtaDepententMinHitTrtWithOutliers){
	    if(m_trtDCTool.empty()) {
	      ATH_MSG_ERROR(" Eta delendent cut on number of TRT hits requested but TrtDCCutTool not specified. ");
	      return StatusCode::FAILURE;
	    } else if(m_trtDCTool.retrieve().isFailure()) {
	      ATH_MSG_ERROR(" Unable to retrieve tool "<<m_trtDCTool);
	      return StatusCode::FAILURE;
	    }
	    ATH_MSG_INFO("Retrieved tool "<<m_trtDCTool);
	    if(m_useEtaDepententMinHitTrt){
	      ATH_MSG_INFO("Using eta dependent cut on number of TRT hits.");
	    }
	    if(m_useEtaDepententMinHitTrtWithOutliers){
	      ATH_MSG_INFO("Using eta dependent cut on number of TRT hits + outliers.");
	    }
    }else{
      m_trtDCTool.disable();
    }

    // Read handle for AtlasFieldCacheCondObj
    ATH_CHECK( m_fieldCacheCondObjInputKey.initialize() );

    ATH_MSG_INFO("Using cuts on the number of Silicon hits");
    if(m_usePtDependentCuts) {
	    //checking whether sizes of cuts and pt interval expressed in vectors match
	    if( m_ptBenchmarks.size() != m_nSCTValues.size()){
	      ATH_MSG_ERROR( "Number of cuts DOES NOT match the number of intervals to apply. Please check jobOptions. ");
	      return StatusCode::FAILURE;
	    } else if (m_ptBenchmarks.empty()){
	      ATH_MSG_ERROR( "Zero vectors for number of cuts and pt intervals. Please check jobOptions. ");
	      return StatusCode::FAILURE;
	    }//end of vector size protection block   
    }//end of memory protection 
    return StatusCode::SUCCESS;
  }
    
  // ---------------------------------------------------------------------
  StatusCode InDetDetailedTrackSelectorTool::finalize()
  {
    ATH_MSG_DEBUG( "Finalize successful" );
    return StatusCode::SUCCESS;
  }

  // ---------------------------------------------------------------------
  bool 
  InDetDetailedTrackSelectorTool::decision(const Trk::Track& track,const Trk::Vertex* vertex) const{
    int nHitTrt = m_nHitTrt;
    int nHitTrtPlusOutliers = m_nHitTrtPlusOutliers;
    const Trk::Perigee* perigeeBeforeExtrapolation=dynamic_cast<const Trk::Perigee*>(track.perigeeParameters());
    if (perigeeBeforeExtrapolation && m_usePreselectionCuts){
	    bool preselectionDecision=preselectionBeforeExtrapolation(*perigeeBeforeExtrapolation);
	    if (!preselectionDecision) {
	      ATH_MSG_DEBUG("Track rejected because of preselection decision!");
	      return false;
	    }
    } else if (m_usePreselectionCuts){
	    ATH_MSG_INFO( " Preselection was requested but cannot be made since no Perigee in Track is available. This is not an error." );
    }
    const Trk::Vertex* myVertex=vertex;
    //in case no Vertex is provided by the user, beam position will be used if available
    if (myVertex==nullptr) {
      myVertex = getBeamSpot(Gaudi::Hive::currentContext());
    }
    Trk::PerigeeSurface perigeeSurface(myVertex->position());
    const Trk::TrackParameters *firstmeaspar=nullptr;
    for (unsigned int i=0;i<track.trackParameters()->size();i++){
      if ( (*track.trackParameters())[i]->covariance() && !dynamic_cast<const Trk::Perigee*>((*track.trackParameters())[i])) {
        firstmeaspar=(*track.trackParameters())[i];
        break;
      }
    }
    if (!firstmeaspar) {
      //assumes perigeeParameters exist...
      //no track selection if firstmeas + perigee does not exist !
      firstmeaspar=track.perigeeParameters();
      if (!firstmeaspar){
	      ATH_MSG_WARNING( " First measurment on track is missing. Using perigee Parameters, but they are missing: 0 pointer! Track selection failed " );
	      //clean up vertex
	      if (myVertex!=vertex) {
	        delete myVertex;
	        myVertex=nullptr;
	      }
	      return false;
	    }
    }
    const Trk::TrackParameters* extrapolatedParameters= m_extrapolator->extrapolate(Gaudi::Hive::currentContext(),
                                                                                    *firstmeaspar,
                                                                                    perigeeSurface,
                                                                                    Trk::anyDirection,
                                                                                    true,
                                                                                    track.info().particleHypothesis() ).release(); 
    const Trk::Perigee* extrapolatedPerigee = extrapolatedParameters ? dynamic_cast<const Trk::Perigee*>(extrapolatedParameters) : nullptr; 
    if (!extrapolatedPerigee || !extrapolatedPerigee->covariance() ) {
      ATH_MSG_WARNING( "Track Selector failed to extrapolate track to the vertex: " << myVertex->position() );
      if (extrapolatedParameters) {
        ATH_MSG_WARNING( "The return object of the extrapolator was not a perigee even if a perigeeSurface was used!" );
        delete extrapolatedParameters;
        extrapolatedParameters=nullptr;
      }
    }
    
    //decision based on the track parameters 
    const Trk::RecVertex* recVertex = dynamic_cast<const Trk::RecVertex*>(myVertex);
    bool dec = decision(extrapolatedPerigee, recVertex ? &recVertex->covariancePosition() : nullptr );
    if (myVertex!=vertex) {
      delete myVertex;
      myVertex=nullptr;
    }
    bool isInTrtAcceptance=true;
    if (!extrapolatedPerigee || std::fabs(extrapolatedPerigee->momentum().eta())>m_TrtMaxEtaAcceptance) {
      isInTrtAcceptance=false;
    }
    if (extrapolatedPerigee!=track.perigeeParameters()) {
      delete extrapolatedPerigee;
      extrapolatedPerigee=nullptr;
    }
    if(!dec) { 
      ATH_MSG_DEBUG("Track rejected because of perigee parameters!");
      return false;
    }
    if (m_useTrackQualityInfo) {
      const Trk::FitQuality*  TrkQuality=track.fitQuality();
      if (TrkQuality==nullptr) {
        ATH_MSG_WARNING( "Requested cut on track quality was not possible. Track has no FitQuality object attached. Selection failed." );
        return false;
      }
      if (!decision(TrkQuality)) {
        return false;
      }
    }
    if (m_useTrackSummaryInfo) {
      //number of hits, silicon hits, b-layer
      // first ask track for summary
      std::unique_ptr<Trk::TrackSummary> summaryUniquePtr;
      const Trk::TrackSummary* summary = track.trackSummary();
      if (m_trackSumToolAvailable && summary == nullptr) {
        summaryUniquePtr = m_trackSumTool->summary(track);
        summary = summaryUniquePtr.get();
      }
      if (nullptr==summary ) {
        ATH_MSG_FATAL( "Track preselection: cannot create a track summary (but useTrackSummary is true). Selection failed." );
        return false;
      }
      // get the minimum nimber of TRT hits based on eta of the track
      if(m_useEtaDepententMinHitTrt) {
        nHitTrt = m_trtDCTool->minNumberDCs( (*track.trackParameters())[0] );
        if(m_addToMinHitTrt!=0){
          nHitTrt += m_addToMinHitTrt;
        }else{
          nHitTrt = (int)((double)nHitTrt*m_scaleMinHitTrt);
        }
      }
      // get the minimum nimber of TRT hits + outliers based on eta of the track
      if(m_useEtaDepententMinHitTrtWithOutliers) {
        nHitTrtPlusOutliers = m_trtDCTool->minNumberDCs( (*track.trackParameters())[0] );
        if(m_addToMinHitTrtWithOutliers!=0){
          nHitTrtPlusOutliers += m_addToMinHitTrtWithOutliers;
        }else{
          nHitTrtPlusOutliers = (int)((double)nHitTrtPlusOutliers*m_scaleMinHitTrtWithOutliers);
        }
      }
      if (!decision(summary,m_useSharedHitInfo,isInTrtAcceptance, perigeeBeforeExtrapolation,
                    nHitTrt, nHitTrtPlusOutliers)) {
        return false;
      }
    }
    return true;  
  }

  // ---------------------------------------------------------------------
  bool 
  InDetDetailedTrackSelectorTool::decision(const Trk::TrackParticleBase& track,const Trk::Vertex* vertex) const{
    int nHitTrt = m_nHitTrt;
    int nHitTrtPlusOutliers = m_nHitTrtPlusOutliers;
    const Trk::TrackParameters* definintParameters=&(track.definingParameters());
    const Trk::Perigee* perigeeBeforeExtrapolation=dynamic_cast<const Trk::Perigee*>(definintParameters);
    if (perigeeBeforeExtrapolation && m_usePreselectionCuts) {
      bool preselectionDecision=preselectionBeforeExtrapolation(*perigeeBeforeExtrapolation);  
      if (!preselectionDecision) {
        ATH_MSG_DEBUG("Track rejected because of preselection decision!");
        return false;
      }
    } else if (m_usePreselectionCuts) {
      ATH_MSG_WARNING( " Preselection was requested but cannot be made since the Perigee is not the defining Parameter of the TrackParticle. This is not an error." );
    }
    bool isInTrtAcceptance=true;
    if (!perigeeBeforeExtrapolation || std::fabs(perigeeBeforeExtrapolation->momentum().eta())>m_TrtMaxEtaAcceptance) {
      isInTrtAcceptance=false;
    }
    if (m_useTrackQualityInfo) {
      const Trk::FitQuality*  TrkQuality=track.fitQuality();
      if (TrkQuality==nullptr) {
        ATH_MSG_WARNING( "Requested cut on track quality was not possible. TrackParticleBase has no FitQuality object attached. Selection failed." );
        return false;
      }
      if (!decision(TrkQuality)) {
        return false;
      }
    }
    if (m_useTrackSummaryInfo) {
      //number of hits, silicon hits, b-layer
      const Trk::TrackSummary* summary = track.trackSummary();
      if (nullptr==summary ) { 
        ATH_MSG_WARNING( "Track preselection: cannot create a track summary (but useTrackSummary is true). Selection failed." );
        return false;
      }
      if(m_useEtaDepententMinHitTrt) {
        nHitTrt = m_trtDCTool->minNumberDCs( (track.trackParameters())[0] );
        if(m_addToMinHitTrt!=0){
          nHitTrt += m_addToMinHitTrt;
        }else{
          nHitTrt = (int)((double)nHitTrt*m_scaleMinHitTrt);
        }
      }
      if(m_useEtaDepententMinHitTrtWithOutliers) {
        nHitTrtPlusOutliers = m_trtDCTool->minNumberDCs( (track.trackParameters())[0] );
        if(m_addToMinHitTrtWithOutliers!=0){
          nHitTrtPlusOutliers += m_addToMinHitTrtWithOutliers;
        }else{
          nHitTrtPlusOutliers = (int)((double)nHitTrtPlusOutliers*m_scaleMinHitTrtWithOutliers);
        }
      }
      if ((!perigeeBeforeExtrapolation) or
          (!decision(summary, m_useSharedHitInfo, isInTrtAcceptance, perigeeBeforeExtrapolation,
                     nHitTrt, nHitTrtPlusOutliers))) {
	      return false;
      }      
    }
    const Trk::Perigee* extrapolatedPerigee=dynamic_cast<const Trk::Perigee*>(definintParameters);
    const Trk::Vertex* myVertex=vertex;
    if (vertex==nullptr) {
      myVertex = getBeamSpot(Gaudi::Hive::currentContext());
    }
    Trk::PerigeeSurface perigeeSurface(myVertex->position());
    const Trk::TrackParameters *firstmeaspar=nullptr;
    for (unsigned int i=0;i<track.trackParameters().size();i++) {
      if ((track.trackParameters())[i]->covariance() &&
	     !dynamic_cast<const Trk::Perigee*>((track.trackParameters())[i])) {
        firstmeaspar=(track.trackParameters())[i];
        break;
      }
    }
    if (!firstmeaspar) {
      if (!extrapolatedPerigee || !extrapolatedPerigee->covariance() ) {
	      ATH_MSG_DEBUG( " Track Paraemters at first measurement not found. Perigee not found. Cannot do TrackSelection..." );
	      if (myVertex!=vertex) {
	        delete myVertex;
	        myVertex=nullptr;
	      }
	      return false;
      }
      //using perigee instead of firstmeasurement, since first measurement was not found...
      firstmeaspar=&(track.definingParameters());
    }
    
    ATH_MSG_VERBOSE ("Input to extrapolation: "    << *firstmeaspar);
    ATH_MSG_VERBOSE ("Extrapolating to position: " << myVertex->position()[0] << " , " <<
		     myVertex->position()[1] << " , " << myVertex->position()[2]);
    const Trk::TrackParameters* extrapolatedParameters= firstmeaspar ?
      m_extrapolator->extrapolate(Gaudi::Hive::currentContext(),
                                  *firstmeaspar,
                                  perigeeSurface,
                                  Trk::anyDirection,
                                  true,Trk::pion ).release() : nullptr;
    extrapolatedPerigee = extrapolatedParameters ? dynamic_cast<const Trk::Perigee*>(extrapolatedParameters) : nullptr; 
    if (extrapolatedPerigee==nullptr || !extrapolatedPerigee->covariance()) {
      ATH_MSG_WARNING( "Track Selector failed to extrapolate track to the vertex: " << myVertex->position() );
      if (extrapolatedParameters) {
	      ATH_MSG_WARNING( "The return object of the extrapolator was not a perigee even if a perigeeSurface was used!" );
	      delete extrapolatedParameters;
	      extrapolatedParameters = nullptr;
      }
    }
    if (extrapolatedParameters) ATH_MSG_VERBOSE ("Result: " << *extrapolatedParameters);
    const Trk::RecVertex* recVertex = dynamic_cast<const Trk::RecVertex*>(myVertex);
    bool dec = decision(extrapolatedPerigee, recVertex ? &recVertex->covariancePosition() : nullptr );
    if (myVertex!=vertex) {
      delete myVertex;
      myVertex=nullptr;
    }
    if (extrapolatedPerigee!=&(track.definingParameters())) {
      delete extrapolatedPerigee;
      extrapolatedPerigee=nullptr;
    }
    if(!dec) {
      ATH_MSG_DEBUG("Track rejected because of perigee parameters!");
      return false;
    }
    return true;
  }

  Amg::Vector3D InDetDetailedTrackSelectorTool::getPosOrBeamSpot(const xAOD::Vertex* vertex) const
  {
    if(vertex) return vertex->position();
    if(m_useEventInfoBs){
      SG::ReadHandle<xAOD::EventInfo> evt(m_eventInfo_key);
      if (evt.isValid()) {
        InDet::BeamSpotData temp(evt->beamStatus(), evt->beamPosX(), evt->beamPosY(), evt->beamPosZ(),
                                 evt->beamPosSigmaX(), evt->beamPosSigmaY(), evt->beamPosSigmaZ(),
                                 evt->beamTiltXZ(), evt->beamTiltYZ(), evt->beamPosSigmaXY());
        return temp.beamVtx().position();
      } else {
        ATH_MSG_WARNING( " Cannot get beamSpot center from xAOD::EventInfo. Using (0,0,0)... " );
        return Amg::Vector3D(0,0,0);
      }
    }else{
      SG::ReadCondHandle<InDet::BeamSpotData> beamSpotHandle { m_beamSpotKey };
      if (beamSpotHandle.isValid()) {
        return beamSpotHandle->beamVtx().position();
      } else {
        ATH_MSG_WARNING( " Cannot get beamSpot center from BeamSpotData. Using (0,0,0)... " );
        return Amg::Vector3D(0,0,0);
      }
    }
  }

  // ---------------------------------------------------------------------
  bool 
  InDetDetailedTrackSelectorTool::decision(const xAOD::TrackParticle& tp,const xAOD::Vertex* vertex) const
  {
    int nHitTrt = m_nHitTrt;
    int nHitTrtPlusOutliers = m_nHitTrtPlusOutliers;
    
    const Trk::Perigee& perigee=tp.perigeeParameters();
    if (m_usePreselectionCuts && !preselectionBeforeExtrapolation(perigee)) {
      ATH_MSG_DEBUG("Track rejected because of preselection decision!");
      return false;
    }

    if (m_useTrackQualityInfo && !decision(tp.chiSquared(),tp.numberDoF())) {
      ATH_MSG_DEBUG("Track rejected because of bad fit quality!");
      return false;
    }

    if (m_useTrackSummaryInfo) {
      //number of hits, silicon hits, b-layer

      if(m_useEtaDepententMinHitTrt) {
        nHitTrt = m_trtDCTool->minNumberDCs( &perigee );
        if(m_addToMinHitTrt!=0){
          nHitTrt += m_addToMinHitTrt;
        }else{
          nHitTrt = (int)((double)nHitTrt*m_scaleMinHitTrt);
        }
      }
      if(m_useEtaDepententMinHitTrtWithOutliers) {
        nHitTrtPlusOutliers = m_trtDCTool->minNumberDCs( &perigee );
        if(m_addToMinHitTrtWithOutliers!=0){
          nHitTrtPlusOutliers += m_addToMinHitTrtWithOutliers;
        }else{
          nHitTrtPlusOutliers = (int)((double)nHitTrtPlusOutliers*m_scaleMinHitTrtWithOutliers);
        }
      }
      int nb = getCount(tp,xAOD::numberOfInnermostPixelLayerHits );
      int np = getCount(tp,xAOD::numberOfPixelHits );
      int npd = getCount(tp,xAOD::numberOfPixelDeadSensors );
      int ns  = getCount(tp,xAOD::numberOfSCTHits );
      int nhp = getCount(tp,xAOD::numberOfPixelHoles );
      int nhs = getCount(tp,xAOD::numberOfSCTHoles );
      int ndhs = getCount(tp,xAOD::numberOfSCTDoubleHoles);
      bool eiph = (getCount(tp,xAOD::expectInnermostPixelLayerHit)==1);

      //**-----------------------------------------------------------------------
      if(m_usePtDependentCuts) {
	      double pt = tp.pt();
	      unsigned int it = 0;
        for(; it< m_ptBenchmarks.size()-1; ++it ) {
          if(pt>m_ptBenchmarks[it] && pt <=m_ptBenchmarks[it+1] && ns < m_nSCTValues[it]) {
            ATH_MSG_DEBUG("Track rejected because of Pt-Dependent SCT Hit cut (CAREFUL! Excludes dead modules)") ;
            return false;
          }
        }//end of  pt intervals loop
	
        //now cutting all the rest by the last value in the vector   
        if(pt>m_ptBenchmarks[it+1] && ns < m_nSCTValues[it+1]) {
          ATH_MSG_DEBUG("Track rejected because of Pt-Dependent SCT Hit cut (CAREFUL! Excludes dead modules)") ;
          return false;
        }
      }
      
      //*--------------------------------------------------------------------------------

      //normal cuts in all their variety    
    
      if(nb == 0 && nb < m_nHitBLayer) {
        ATH_MSG_DEBUG("Track rejected because of nHitBLayer "<<nb<<" < "<<m_nHitBLayer);
	if (eiph) {
          ATH_MSG_DEBUG("and track rejected because at least one hit is expected in the innermost pixel layer") ;
          return false;
        }else  ATH_MSG_DEBUG("recovered track as no b-layer expected") ;
      }//end of checking the b-layer
     
      if(np+npd < m_nHitPix) {
        ATH_MSG_DEBUG("Track rejected because of nHitPix "<<np+npd<<" < "<<m_nHitPix);
        return false;
      }
    
      if(np < m_nHitPixPhysical) {
        ATH_MSG_DEBUG("Track rejected because of nHitPixPhysical "<<np<<" < "<<m_nHitPixPhysical);
        return false;
      }
    
      int nsd = getCount(tp,xAOD::numberOfSCTDeadSensors);
      if(ns+nsd < m_nHitSct) {
        ATH_MSG_DEBUG("Track rejected because of nHitSct "<<ns+nsd<<" < "<<m_nHitSct);
        return false;
      }
    
      if(np+ns+npd+nsd < m_nHitSi) {
        ATH_MSG_DEBUG("Track rejected because of nHitSi "<<np+npd+ns+nsd<<" < "<<m_nHitSi);
        return false;
      }
    
      if(np+ns < m_nHitSiPhysical) {
        ATH_MSG_DEBUG("Track rejected because of nHitSiPhysical "<<np+ns<<" < "<<m_nHitSiPhysical);
        return false;
      }
    
      if (nb+np+npd< m_nHitBLayerPlusPix){
	ATH_MSG_DEBUG("Track rejected because of nHitBLayerPlusPix "<<nb+np+npd<<" < "<<m_nHitBLayerPlusPix);
	return false;
      }

      // Cuts on number of Holes
      if (nhp+nhs > m_nHoles){
	ATH_MSG_DEBUG("Track rejected because of nHolesPixPlusSCT "<<nhp+nhs<<" > "<<m_nHoles);
	return false;
      }
   
      if(ndhs > m_nDoubleHoles){
	ATH_MSG_DEBUG("Track rejected because of nDoubleHolesSCT "<<ndhs<<" > "<<m_nDoubleHoles);
	return false;
      }
   
      if(nhp > m_nHolesPix){
        ATH_MSG_DEBUG("Track rejected because of nHolesPix "<<nhp<<" > "<<m_nHolesPix);
        return false;
      }
   
      if (nhs > m_nHolesSct){
        ATH_MSG_DEBUG("Track rejected because of nHolesSct "<<nhs<<" > "<<m_nHolesSct);
        return false;
      }
   
      if (std::fabs(tp.eta())>m_TrtMaxEtaAcceptance) {
        int nh = getCount(tp,xAOD::numberOfTRTHits);
        if(nh < nHitTrt) {
          ATH_MSG_DEBUG("Track rejected because of nHitTrt "<<nh<<" < "<<nHitTrt);
          return false;
        }

      int nhh = getCount(tp, xAOD::numberOfTRTHits ) + getCount(tp, xAOD::numberOfTRTOutliers );
      if (nhh<nHitTrtPlusOutliers) {
        ATH_MSG_DEBUG("Track rejected because of nHitTrtPlusOutliers "<<nhh<<" < "<<nHitTrtPlusOutliers);
        return false;
      }

      int nhthits= getCount(tp,xAOD::numberOfTRTHighThresholdHits);
      if (nhthits<m_nHitTrtHighE) {
        ATH_MSG_DEBUG("Track rejected because of nHitTrtHighE "<<nhthits<<" < "<<m_nHitTrtHighE);
        return false;
      }
      
      int nhthitsWithOutliers= getCount(tp,xAOD::numberOfTRTHighThresholdHits) + getCount(tp,xAOD::numberOfTRTHighThresholdOutliers);
      if (nhthitsWithOutliers<m_nHitTrtPlusOutliersHighE) {
        ATH_MSG_DEBUG("Track rejected because of nHitTrtPlusOutliersHighE "<<nhthitsWithOutliers<<" < "<<m_nHitTrtPlusOutliersHighE);
        return false;
      }

      if ( getCount(tp, xAOD::numberOfTRTHits )>0) {
        double nhe = getCount(tp,xAOD::numberOfTRTHighThresholdHits);
        nhe /= getCount(tp, xAOD::numberOfTRTHits);
        if(nhe > m_nHitTrtHighEFraction ) {
          ATH_MSG_DEBUG("Track rejected because of nHitTrtHighEFraction "<<nhe<<" < "<<m_nHitTrtHighEFraction);
          return false;
        }
      }
      
      if ( getCount(tp, xAOD::numberOfTRTHits ) + getCount(tp, xAOD::numberOfTRTOutliers ) > 0 ) {
        double nheh = (double)(getCount(tp,xAOD::numberOfTRTHighThresholdHits) + getCount(tp,xAOD::numberOfTRTHighThresholdOutliers))/
                (double)(getCount(tp, xAOD::numberOfTRTHits) + getCount(tp, xAOD::numberOfTRTOutliers ) );
        if(nheh<0.) nheh=0.;
        if (nheh>1.) nheh=1.;
        if(nheh > m_nHitTrtHighEFractionWithOutliers ) {
          ATH_MSG_DEBUG("Track rejected because of nHitTrtHighEFractionWithOutliers "<<nheh<<" < "<<m_nHitTrtHighEFractionWithOutliers);
          return false;
        }
      }
      }
      if (m_useSharedHitInfo) {
        int nbs = getCount(tp,xAOD::numberOfInnermostPixelLayerSharedHits);
        if (nbs>1) nbs=1;
        if(nbs>m_nSharedBLayer) {
          ATH_MSG_DEBUG("Track rejected because of nSharedBLayer "<<nbs<<" < "<<m_nSharedBLayer);
          return false;
        }

	int nps = getCount(tp,xAOD::numberOfPixelSharedHits);
	if(nps>m_nSharedPix) {
	  ATH_MSG_DEBUG("Track rejected because of nSharedPix "<<nps<<" < "<<m_nSharedPix);
	  return false;
	}

	int nss = getCount(tp,xAOD::numberOfSCTSharedHits);
	if(nss > m_nSharedSct) {
	  ATH_MSG_DEBUG("Track rejected because of nSharedSct "<<nss<<" < "<<m_nSharedSct);
	  return false;
	}

	int nst = nps + nss;
	if(nst>m_nSharedSi) {
	  ATH_MSG_DEBUG("Track rejected because of nSharedSi "<<nst<<" < "<<m_nSharedSi);
	  return false;
	}
      }
    }
    Trk::PerigeeSurface perigeeSurface( getPosOrBeamSpot(vertex) );
    
    const Trk::TrackParameters* extrapolatedParameters= m_extrapolator->extrapolate(
      Gaudi::Hive::currentContext(),
      perigee,perigeeSurface,
      Trk::anyDirection,true,Trk::pion).release();
    const Trk::Perigee* extrapolatedPerigee = extrapolatedParameters ? dynamic_cast<const Trk::Perigee*>(extrapolatedParameters) : nullptr; 
    if (extrapolatedPerigee==nullptr) {
      ATH_MSG_WARNING( "Extrapolation to the vertex failed: " << perigeeSurface << std::endl << perigee );
      if (extrapolatedParameters!=nullptr) {
        ATH_MSG_WARNING( "The return object of the extrapolator was not a perigee even if a perigeeSurface was used!" );
        delete extrapolatedParameters;
        extrapolatedParameters=nullptr;
      }
      return false;
    }
    bool dec = false;
    if( vertex ){
      // for now copy the position error
      AmgSymMatrix(3) vertexError = vertex->covariancePosition();
      dec = decision(extrapolatedPerigee,&vertexError);
    }else{
      dec = decision(extrapolatedPerigee,nullptr);
    }

    delete extrapolatedPerigee;

    if(!dec) {
      ATH_MSG_DEBUG("Track rejected because of perigee parameters!");
      return false;
    }
    
    return true;
  }

  
  // ---------------------------------------------------------------------
  bool InDetDetailedTrackSelectorTool::decision(const Trk::Perigee* track,const AmgSymMatrix(3)* covariancePosition) const {

    // checking pointer first
    if(nullptr==track || !track->covariance()) { 
      ATH_MSG_WARNING( "Decision on measured perigee: Zero pointer to measured perigee passed. Selection failed." );
      return false;
    }

    const AmgVector(5)& perigeeParms = track->parameters();
    
    // only check pt if mag. field is on
    const EventContext& ctx = Gaudi::Hive::currentContext();
    SG::ReadCondHandle<AtlasFieldCacheCondObj> readHandle{m_fieldCacheCondObjInputKey, ctx};
    const AtlasFieldCacheCondObj* fieldCondObj{*readHandle};
    if (fieldCondObj == nullptr) {
      ATH_MSG_ERROR("execute: Failed to retrieve AtlasFieldCacheCondObj with key " << m_fieldCacheCondObjInputKey.key());
      return false;
    }
    MagField::AtlasFieldCache fieldCache;
    fieldCondObj->getInitializedCache (fieldCache);

    if (fieldCache.solenoidOn()){//B field
      if (perigeeParms[Trk::qOverP] == 0.) {
	ATH_MSG_DEBUG("Track rejected because of qOverP == 0.");
	return false;
      }
      double p = std::fabs(1./perigeeParms[Trk::qOverP]);
      if (p<m_pMin) {
	ATH_MSG_DEBUG("Track rejected because of p " << p << " < " << m_pMin);
	return false;
      }
      double pt = p*std::sin(perigeeParms[Trk::theta]);
      if (pt<m_pTMin) {
	ATH_MSG_DEBUG("Track rejected because of pt " << pt << " < " << m_pTMin);
	return false;
      }
    }
    
    if (std::fabs(perigeeParms[Trk::d0]) > m_IPd0Max) {
      ATH_MSG_DEBUG("Track rejected because of fabs(d0) " << std::fabs(perigeeParms[Trk::d0]) << " > " << m_IPd0Max);
      return false;
    }
    
    if (std::fabs(perigeeParms[Trk::z0]*std::sin(perigeeParms[Trk::theta])) > m_IPz0Max) {
      ATH_MSG_DEBUG("Track rejected because of fabs(z0*sin(theta)) " << std::fabs(perigeeParms[Trk::z0]*std::sin(perigeeParms[Trk::theta])) << " > " << m_IPz0Max);
      return false;
    }
    
    if (std::fabs(perigeeParms[Trk::z0]) > m_z0Max) {
      ATH_MSG_DEBUG("Track rejected because of fabs(z0) " << std::fabs(perigeeParms[Trk::z0]) << " > " << m_z0Max);
      return false;
    }

    if (sqrt( (*track->covariance())(Trk::z0,Trk::z0) )*std::sin(perigeeParms[Trk::theta])>m_sigIPz0Max) {
      ATH_MSG_DEBUG("Track rejected because of err(z0)*sin(theta) " << sqrt( (*track->covariance())(Trk::z0,Trk::z0) )*std::sin(perigeeParms[Trk::theta]) << " > " << m_sigIPz0Max);
      return false;
    }
    
    if (sqrt( (*track->covariance())(Trk::d0,Trk::d0) )>m_sigIPd0Max) {
      ATH_MSG_DEBUG("Track rejected because of err(d0) " << sqrt( (*track->covariance())(Trk::d0,Trk::d0) ) << " > " << m_sigIPd0Max);
      return false;
    }
    
    if (m_d0significanceMax>0 || m_z0significanceMax>0) {
      
      double sinTheta = std::sin(perigeeParms[Trk::theta]);
      double cosTheta = std::cos(perigeeParms[Trk::theta]);
      double d0wrtPriVtx = perigeeParms[Trk::d0];
      double deltaZ = perigeeParms[Trk::z0];
      double z0wrtPriVtx = deltaZ*sinTheta;
      double testtrackSigD0 = sqrt( (*track->covariance())(Trk::d0,Trk::d0) );
      double testtrackSigZ0 = sqrt( (*track->covariance())(Trk::z0,Trk::z0) );
      double testtrackSigTh = sqrt( (*track->covariance())(Trk::theta,Trk::theta) );
      // error on IP:
      double trackPhi = perigeeParms[Trk::phi];
      double dIPdx = std::sin(trackPhi);
      double dIPdy = -std::cos(trackPhi);
      double DD0 = testtrackSigD0*testtrackSigD0;
      double newD0Err=0;
      if (covariancePosition) {
        double DXX = dIPdx*dIPdx* (*covariancePosition)(0,0);
        double DYY = dIPdy*dIPdy* (*covariancePosition)(1,1);
        double DXY = 2.*dIPdx*dIPdy* (*covariancePosition)(0,1);
        newD0Err = DD0 + DXX + DYY + DXY;
      } else {
	newD0Err = DD0;
      }
      
      double d0ErrwrtPriVtx = (newD0Err>0 ? sqrt(newD0Err) : -10e-9);
      
      if (d0ErrwrtPriVtx<0) {
	ATH_MSG_WARNING( " error on d0 is negative: numeric error... (not expected. please report!)" );
      }
      
      if (m_d0significanceMax>0) {
	if (std::fabs(d0wrtPriVtx/d0ErrwrtPriVtx)>m_d0significanceMax) {
	  ATH_MSG_DEBUG("Track rejected because of fabs(d0wrtPriVtx/d0ErrwrtPriVtx) " << std::fabs(d0wrtPriVtx/d0ErrwrtPriVtx) << " > " << m_d0significanceMax);
	  return false;
	}
      }
      
      if (m_z0significanceMax>0) {

	// error on zIP:
	double dZIPdTheta = deltaZ*cosTheta;
	double dZIPdz0 = sinTheta;
	double dZIPdzV = -sinTheta;
	double DTheta2 = dZIPdTheta*dZIPdTheta*testtrackSigTh*testtrackSigTh;
	double DZ02 = dZIPdz0*dZIPdz0*testtrackSigZ0*testtrackSigZ0;
	double DThetaZ0 = 2.*dZIPdTheta*dZIPdz0*(*track->covariance())(Trk::theta,Trk::z0);
	double newZ0Err(0);
	if (covariancePosition) {
	  double DZV2 = dZIPdzV*dZIPdzV* (*covariancePosition)(2,2);
	  newZ0Err = DTheta2 + DZ02 + DZV2 + DThetaZ0;
	} else {
	  newZ0Err = DTheta2 + DZ02 + DThetaZ0;
	}
	
	double z0ErrwrtPriVtx = (newZ0Err>0 ? sqrt(newZ0Err) : -10e-9);
        
	if (z0ErrwrtPriVtx<0) {
	  ATH_MSG_WARNING( " error on z0 is negative: numeric error... (not expected. please report!)" );
	}

	if (std::fabs(z0wrtPriVtx/z0ErrwrtPriVtx)>m_z0significanceMax) {
	  ATH_MSG_DEBUG("Track rejected because of fabs(z0wrtPriVtx/z0ErrwrtPriVtx) " << std::fabs(z0wrtPriVtx/z0ErrwrtPriVtx) << " > " << m_z0significanceMax);
	  return false;
	}
      }
      
    }
    
    if (std::fabs(track->momentum().eta())>m_etaMax) {
      ATH_MSG_DEBUG("Track rejected because of fabs(eta) " << std::fabs(track->momentum().eta()) << " > " << m_etaMax);
      return false;
    }
    
    return true;
  }      
  
  // ---------------------------------------------------------------------
  bool InDetDetailedTrackSelectorTool::decision(const Trk::FitQuality*  trkQuality) const
  {
    if(nullptr == trkQuality) {
      ATH_MSG_WARNING( "Null FitQuality pointer passed. No track Quality cut possible. Selection failed." );
      return false;
    }
    return decision(trkQuality->chiSquared(),trkQuality->numberDoF());
  }

  bool InDetDetailedTrackSelectorTool::decision(double chi2, int ndf ) const{

    double proba = 1.;

    if(ndf>0 && chi2>=0.) {
      Genfun::CumulativeChiSquare myCumulativeChiSquare(ndf);
      proba = 1.-myCumulativeChiSquare(chi2);
    }
    
    if(chi2>m_fitChi2) {
      ATH_MSG_DEBUG("Track rejected because of chi2 "<<chi2<<" > "<<m_fitChi2);
      return false;
    }

    if(proba<m_fitProb) {
      ATH_MSG_DEBUG("Track rejected because of fit probability "<<proba<<" > "<<m_fitProb);
      return false;
    }
    if(ndf>0 && chi2/double(ndf)>m_fitChi2OnNdfMax) {
      ATH_MSG_DEBUG("Track rejected because of chi2/ndof "<<chi2/double(ndf)<<" > "<<m_fitChi2OnNdfMax);
      return false;
    }
    
    return true;
  }


  // ---------------------------------------------------------------------
  bool InDetDetailedTrackSelectorTool::decision(const Trk::TrackSummary* summary,bool useSharedHitInfo,bool useTrtHitInfo,
						const Trk::Perigee * track,
                                                const int nHitTrt, const int nHitTrtPlusOutliers) const
  {
    if (summary==nullptr) {
      ATH_MSG_WARNING( "Null TrackSummary pointer passed. Selection failed." );
      return false;
    }
    
    int nb = summary->get(Trk::numberOfInnermostPixelLayerHits);
    
    if(nb<0) nb=0; 
    
    int np  = summary->get(Trk::numberOfPixelHits);
    if(np<0) np=0;
    
    int npd = summary->get(Trk::numberOfPixelDeadSensors);
    if(npd<0) npd=0;
    
    int ns = summary->get(Trk::numberOfSCTHits);
    if(ns<0) ns=0;
    
    int nhp = summary->get(Trk::numberOfPixelHoles);
    if (nhp < 0) nhp = 0;
    
    int nhs = summary->get(Trk::numberOfSCTHoles);
    if (nhs < 0) nhs = 0;
    
    int ndhs = summary->get(Trk::numberOfSCTDoubleHoles);
    if (ndhs < 0) ndhs = 0;
    
    //**-----------------------------------------------------------------------

    if(m_usePtDependentCuts) {
      if (!track) {
         return false;
      }
      const AmgVector(5)& perigeeParms = track->parameters();
      double p = std::fabs(1./perigeeParms[Trk::qOverP]);  
      double pt = p*std::sin(perigeeParms[Trk::theta]);
     
      unsigned int it = 0;
      for(; it< m_ptBenchmarks.size()-1; ++it ) {
	if(pt>m_ptBenchmarks[it] && pt <=m_ptBenchmarks[it+1] && ns < m_nSCTValues[it]) {
	  ATH_MSG_DEBUG("Track rejected because of Pt-Dependent SCT Hit cut (CAREFUL! Excludes dead modules)") ;
	  return false;
	}
      }//end of  pt intervals loop
    
      //now cutting all the rest by the last value in the vector   
      if(pt>m_ptBenchmarks[it+1] && ns < m_nSCTValues[it+1]) {
	ATH_MSG_DEBUG("Track rejected because of Pt-Dependent SCT Hit cut (CAREFUL! Excludes dead modules)") ;
	return false;
      }
     
    }
   
    //*--------------------------------------------------------------------------------

    //normal cuts in all their variety    
    
    if(nb == 0 && nb < m_nHitBLayer) {
      ATH_MSG_DEBUG("Track rejected because of nHitBLayer "<<nb<<" < "<<m_nHitBLayer);
      if(m_inDetTestPixelLayerTool.empty()) {
	ATH_MSG_DEBUG("and no blayer tool configured, so will not try to recover track");
	return false;
      } else if (m_inDetTestPixelLayerTool->expectHitInInnermostPixelLayer(track)) {
	ATH_MSG_DEBUG("and track rejected because at least one hit is expected in the innermost pixel layer") ;
	return false;
      }else  ATH_MSG_DEBUG("recovered track as no b-layer expected") ;
    }//end of checking the b-layer
     
    if(np+npd < m_nHitPix) {
      ATH_MSG_DEBUG("Track rejected because of nHitPix "<<np+npd<<" < "<<m_nHitPix);
      return false;
    }
    
    if(np < m_nHitPixPhysical) {
      ATH_MSG_DEBUG("Track rejected because of nHitPixPhysical "<<np<<" < "<<m_nHitPixPhysical);
      return false;
    }
    
    int nsd = summary->get(Trk::numberOfSCTDeadSensors);
    if(nsd<0)
      nsd=0;
    
    if(ns+nsd < m_nHitSct) 
      {
	ATH_MSG_DEBUG("Track rejected because of nHitSct "<<ns+nsd<<" < "<<m_nHitSct);
	return false;
      }
    
    if((np+ns+npd+nsd) < m_nHitSi) 
      {
	ATH_MSG_DEBUG("Track rejected because of nHitSi "<<np+npd+ns+nsd<<" < "<<m_nHitSi);
	return false;
      }
    
    if((np+ns) < m_nHitSiPhysical) 
      {
	ATH_MSG_DEBUG("Track rejected because of nHitSiPhysical "<<np+ns<<" < "<<m_nHitSiPhysical);
	return false;
      }
    
    if ((nb+np+npd)< m_nHitBLayerPlusPix) 
      {
	ATH_MSG_DEBUG("Track rejected because of nHitBLayerPlusPix "<<nb+np+npd<<" < "<<m_nHitBLayerPlusPix);
	return false;
      }

    // Cuts on number of Holes

    if ((nhp+nhs) > m_nHoles)
      {
	ATH_MSG_DEBUG("Track rejected because of nHolesPixPlusSCT "<<nhp+nhs<<" > "<<m_nHoles);
	return false;
      }
   
    if (ndhs > m_nDoubleHoles)
      {
	ATH_MSG_DEBUG("Track rejected because of nDoubleHolesSCT "<<ndhs<<" > "<<m_nDoubleHoles);
	return false;
      }
   
    if (nhp > m_nHolesPix)
      {
	ATH_MSG_DEBUG("Track rejected because of nHolesPix "<<nhp<<" > "<<m_nHolesPix);
	return false;
      }
   
    if (nhs > m_nHolesSct)
      {
	ATH_MSG_DEBUG("Track rejected because of nHolesSct "<<nhs<<" > "<<m_nHolesSct);
	return false;
      }
   
    if (useTrtHitInfo) {
 
      int nh = summary->get(Trk::numberOfTRTHits);
      if(nh<0) nh=0;
      if(nh < nHitTrt) {
	ATH_MSG_DEBUG("Track rejected because of nHitTrt "<<nh<<" < "<<nHitTrt);
	return false;
      }

      int nhh = summary->get( Trk::numberOfTRTHits ) + summary->get( Trk::numberOfTRTOutliers );
      if (nhh<0) nhh=0;
      if (nhh<nHitTrtPlusOutliers) {
	ATH_MSG_DEBUG("Track rejected because of nHitTrtPlusOutliers "<<nhh<<" < "<<nHitTrtPlusOutliers);
	return false;
      }

      int nhthits=summary->get(Trk::numberOfTRTHighThresholdHits);
      if (nhthits<0) nhthits=0;
      if (nhthits<m_nHitTrtHighE) {
	ATH_MSG_DEBUG("Track rejected because of nHitTrtHighE "<<nhthits<<" < "<<m_nHitTrtHighE);
	return false;
      }
      
      int nhthitsWithOutliers=summary->get(Trk::numberOfTRTHighThresholdHits) + summary->get(Trk::numberOfTRTHighThresholdOutliers);
      if (nhthitsWithOutliers<0) nhthitsWithOutliers=0;
      if (nhthitsWithOutliers<m_nHitTrtPlusOutliersHighE) {
	ATH_MSG_DEBUG("Track rejected because of nHitTrtPlusOutliersHighE "<<nhthitsWithOutliers<<" < "<<m_nHitTrtPlusOutliersHighE);
	return false;
      }

      if (summary->get( Trk :: numberOfTRTHits )>0) {
	double nhe = (double)summary->get(Trk::numberOfTRTHighThresholdHits) / (double)summary->get( Trk::numberOfTRTHits );
	if(nhe<0.) nhe=0.;
	if(nhe > m_nHitTrtHighEFraction ) {
	  ATH_MSG_DEBUG("Track rejected because of nHitTrtHighEFraction "<<nhe<<" < "<<m_nHitTrtHighEFraction);
	  return false;
	}
      }
      
      if ( summary->get( Trk :: numberOfTRTHits ) + summary->get( Trk :: numberOfTRTOutliers ) > 0 ) {
	double nheh = (double)(summary->get(Trk::numberOfTRTHighThresholdHits) + summary->get(Trk::numberOfTRTHighThresholdOutliers))/
	  (double)(summary->get( Trk::numberOfTRTHits) + summary->get( Trk :: numberOfTRTOutliers ) );
	if(nheh<0.) nheh=0.;
	if (nheh>1.) nheh=1.;
	if(nheh > m_nHitTrtHighEFractionWithOutliers ) {
	  ATH_MSG_DEBUG("Track rejected because of nHitTrtHighEFractionWithOutliers "<<nheh<<" < "<<m_nHitTrtHighEFractionWithOutliers);
	  return false;
	}
      }
    }

    if (useSharedHitInfo) {

      int nbs = summary->get(Trk::numberOfInnermostPixelLayerSharedHits);
      if(nbs < 0) nbs = 0;
      if (nbs>1) nbs=1;
      if(nbs>m_nSharedBLayer) {
	ATH_MSG_DEBUG("Track rejected because of nSharedBLayer "<<nbs<<" < "<<m_nSharedBLayer);
	return false;
      }

      int nps = summary->get(Trk::numberOfPixelSharedHits);
      if(nps < 0) nps = 0;
      if(nps>m_nSharedPix) {
	ATH_MSG_DEBUG("Track rejected because of nSharedPix "<<nps<<" < "<<m_nSharedPix);
	return false;
      }

      int nss = summary->get(Trk::numberOfSCTSharedHits);
      if(nss < 0) nss = 0;
      if(nss > m_nSharedSct) {
	ATH_MSG_DEBUG("Track rejected because of nSharedSct "<<nss<<" < "<<m_nSharedSct);
	return false;
      }

      int nst = nps + nss;
      if(nst>m_nSharedSi) {
	ATH_MSG_DEBUG("Track rejected because of nSharedSi "<<nst<<" < "<<m_nSharedSi);
	return false;
      }
    }
    
    return true;
    
  }

  // ---------------------------------------------------------------------
  bool InDetDetailedTrackSelectorTool::preselectionBeforeExtrapolation(const Trk::Perigee & myPerigee) const
  {
    const AmgVector(5)& perigeeParms = myPerigee.parameters();
    
    // only check pt if mag. field is on
    const EventContext& ctx = Gaudi::Hive::currentContext();
    SG::ReadCondHandle<AtlasFieldCacheCondObj> readHandle{m_fieldCacheCondObjInputKey, ctx};
    const AtlasFieldCacheCondObj* fieldCondObj{*readHandle};
    if (fieldCondObj == nullptr) {
      ATH_MSG_ERROR("execute: Failed to retrieve AtlasFieldCacheCondObj with key " << m_fieldCacheCondObjInputKey.key());
      return false;
    }
    MagField::AtlasFieldCache fieldCache;
    fieldCondObj->getInitializedCache (fieldCache);

    if (fieldCache.solenoidOn()){//B field
      if (perigeeParms[Trk::qOverP] == 0.) {
	ATH_MSG_DEBUG("Track rejected because of perigee qOverP == 0.");
	return false;
      }
      double p = std::fabs(1./perigeeParms[Trk::qOverP]);
      if (p<m_pMin) {
	ATH_MSG_DEBUG("Track rejected because of p " << p << " < " << m_pMin);
	return false;
      }
      double pt = p*std::sin(perigeeParms[Trk::theta]);
      if (pt<m_pTMin) {
	ATH_MSG_DEBUG("Track rejected because of pt " << pt << " < " << m_pTMin);
	return false;
      }
    }

    if (std::fabs(perigeeParms[Trk::d0]) > m_d0MaxPreselection) {
      ATH_MSG_DEBUG("Track rejected because of fabs(d0) "<<std::fabs(perigeeParms[Trk::d0])<<" < "<<m_d0MaxPreselection);
      return false;
    }

    return true;
  }

} //end of namespace definitions
