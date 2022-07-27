/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "InDetTrackSelectorTool/InDetConversionTrackSelectorTool.h"
#include "VxVertex/Vertex.h"
#include "TrkTrack/Track.h"
#include "TrkParticleBase/TrackParticleBase.h"
// normal includes
#include "TrkTrackSummary/TrackSummary.h"
#include "TrkSurfaces/PerigeeSurface.h"
#include "VxVertex/RecVertex.h"
#include "TrkParameters/TrackParameters.h"
#include "xAODTracking/Vertex.h"
#include <cmath>


namespace InDet
{

 InDetConversionTrackSelectorTool::InDetConversionTrackSelectorTool(const std::string& t, const std::string& n, const IInterface*  p)
 :AthAlgTool(t,n,p),
  m_maxSiD0   (35.),
  m_maxTrtD0 (100.),
  m_maxSiZ0  (200.),
  m_maxTrtZ0 (1200),
  m_minPt    (500.),
  m_trRatio1  (0.5), /** eProbHT cut for Si trks, ntrt<=15 **/
  m_trRatio2  (0.1), /** eProbHT cut for Si trks, ntrt >15 **/
  m_trRatio3 (0.05), /** eProbHT cut for Si trks, ntrt >25 **/
  m_trRatioTRT(0.1), /** eProbHT cut for all TRT tracks **/
  m_trRatioV0  (1.),
  m_sD0_Si     (2.),
  m_sD0_Trt   (0.5),
  m_sZ0_Trt    (3.),
  m_isConv(true)
 {
   // There is room for 10 bins, for future development.
   // The default below represents one eta bin between 0 and 999
   // and no cuts applied.
   m_TRTTrksEtaBins.clear();
   m_TRTTrksBinnedRatioTRT.clear();
   for (unsigned int i=0;i<10;++i) {
     m_TRTTrksEtaBins.push_back(999);
     m_TRTTrksBinnedRatioTRT.push_back(0);
   }

   declareInterface<ITrackSelectorTool>(this);
   declareProperty("maxSiD0",              m_maxSiD0);
   declareProperty("maxTrtD0",             m_maxTrtD0);
   declareProperty("maxSiZ0",              m_maxSiZ0);
   declareProperty("maxTrtZ0",             m_maxTrtZ0);
   declareProperty("minPt",                m_minPt);
   declareProperty("RatioCut1",            m_trRatio1); /** eProbHT cut for Si trks, ntrt<=15 **/
   declareProperty("RatioCut2",            m_trRatio2); /** eProbHT cut for Si trks, ntrt >15 **/
   declareProperty("RatioCut3",            m_trRatio3); /** eProbHT cut for Si trks, ntrt >25 **/
   declareProperty("RatioTRT",             m_trRatioTRT); /** eProbHT cut for all TRT tracks **/
   // See InDet Trt Track Scoring Tool for nTRT cut on TRT-only tracks
   declareProperty("TRTTrksEtaBins"       , m_TRTTrksEtaBins); /* expects 10 eta bins (set unused bins to e.g. 999) */
   declareProperty("TRTTrksBinnedRatioTRT"          , m_TRTTrksBinnedRatioTRT); /* expects 10 values */
   declareProperty("RatioV0",              m_trRatioV0);
   declareProperty("significanceD0_Si",    m_sD0_Si);
   declareProperty("significanceD0_Trt",   m_sD0_Trt);
   declareProperty("significanceZ0_Trt",   m_sZ0_Trt);
   declareProperty("IsConversion",         m_isConv);
   declareProperty("PIDonlyForXe",         m_PIDonlyForXe = false,
     "Only check TRT PID if all hits are Xe hits");
 }


 StatusCode  InDetConversionTrackSelectorTool::initialize()
 {
   ATH_CHECK( m_extrapolator.retrieve() );
   ATH_CHECK(m_beamSpotKey.initialize());

   return StatusCode::SUCCESS;
 }


 unsigned int InDetConversionTrackSelectorTool::getEtaBin(const Trk::Perigee& perigee) const
 {
   // Find the correct bin for applying eta-dependent cuts

   double tanThetaOver2 = std::tan( perigee.parameters()[Trk::theta] / 2.);
   double abs_eta = (tanThetaOver2 == 0) ? 999.0 : std::fabs( std::log(tanThetaOver2) );

   for (unsigned int i=0;i<m_TRTTrksEtaBins.size();++i) {
     if (abs_eta < m_TRTTrksEtaBins[i]) {
       return i;
     }
   }
   return m_TRTTrksEtaBins.size()-1;
 }


 bool InDetConversionTrackSelectorTool::decision(const Trk::Track& track,const Trk::Vertex* vx) const
 {
   const EventContext& ctx = Gaudi::Hive::currentContext();
   bool pass = false;
   const Trk::Perigee* perigee=dynamic_cast<const Trk::Perigee*>(track.perigeeParameters());
   const bool vertexSuppliedByUser{vx!=nullptr};
   const Trk::Vertex* myVertex=vx;
   //in case no Vertex is provided by the user, beam position will be used if available
   if (not vertexSuppliedByUser) {
     SG::ReadCondHandle<InDet::BeamSpotData> beamSpotHandle{ m_beamSpotKey,
                                                             ctx };
     if (beamSpotHandle.isValid()) {
       myVertex=new Trk::RecVertex(beamSpotHandle->beamVtx());
     } else {
       ATH_MSG_WARNING(" Cannot get beamSpot center from BeamSpotData. Using (0,0,0)... " );
       myVertex=new Trk::Vertex(Amg::Vector3D(0,0,0));
     }
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
       if (not vertexSuppliedByUser) delete myVertex;
       return false;
     }
   }

   const Trk::TrackParameters* extrapolatedParameters =
     m_extrapolator->extrapolate(ctx,
                                 *firstmeaspar,
                                 perigeeSurface,
                                 Trk::anyDirection,
                                 true,
                                 track.info().particleHypothesis()).release();
   perigee = extrapolatedParameters
               ? dynamic_cast<const Trk::Perigee*>(extrapolatedParameters)
               : nullptr;
   if (perigee==nullptr || !perigee->covariance() ) {
     ATH_MSG_WARNING( "Track Selector failed to extrapolate track to the vertex: " << myVertex->position() );
     if (extrapolatedParameters!=nullptr) {
       ATH_MSG_WARNING( "The return object of the extrapolator was not a perigee even if a perigeeSurface was used!");
       delete extrapolatedParameters;
       if (not vertexSuppliedByUser) delete myVertex;
       return false;
     }
     if (not vertexSuppliedByUser) delete myVertex;
     return false;
   }

   double qOverP = perigee->parameters()[Trk::qOverP];
   double pt = std::fabs(1./qOverP)*std::sin(perigee->parameters()[Trk::theta]);
   double d0 = perigee->parameters()[Trk::d0];
   double z0 = perigee->parameters()[Trk::z0];
   const Trk::TrackSummary* tSum = track.trackSummary();
   if(tSum){
     double ratioTrk = 1.0;
     int nclus  = tSum->get(Trk::numberOfPixelHits) + tSum->get(Trk::numberOfSCTHits);
     bool isSilicon = (nclus > 0);
     int nTrtHits       = tSum->get(Trk::numberOfTRTHits);
     int nTrtOutliers   = tSum->get(Trk::numberOfTRTOutliers);
     int ntrt           = nTrtHits + nTrtOutliers;
     int nTrtXenonHits  = tSum->get(Trk::numberOfTRTXenonHits);
     if(m_isConv) {
       if(ntrt > 0 && (!m_PIDonlyForXe || nTrtXenonHits==ntrt) ){
         // only check TRT PID if m_PIDonlyForXe is false or all TRT hits are Xenon hits
          ATH_MSG_FATAL( "eProbabilityHT not available for Trk::Tracks only xAOD::TrackParticle objects" );
       }

       // Start of track cuts
       if ( pt >= m_minPt ) {

         // Silicon track cuts
         if ( isSilicon && (std::fabs(d0)<=m_maxSiD0) && (fabs(z0)<=m_maxSiZ0) ) {
           if((ntrt<=15 && ratioTrk>=m_trRatio1) ||
              (ntrt>15 && ntrt<=25 && ratioTrk>=m_trRatio2) ||
              (ntrt>25 && ratioTrk>=m_trRatio3)) pass = true;
         }

         // TRT-only track cuts
         if ( (not isSilicon) && (std::fabs(d0)<=m_maxTrtD0) && (std::fabs(z0)<=m_maxTrtZ0) ) {

           unsigned int eta_bin = getEtaBin(*perigee);
           // TRT-only Tracks: eProbabilityHT cut below
           // See InDet Trt Track Scoring Tool for nTRT cut on TRT-only tracks
           double trRatioTRT = std::max( m_trRatioTRT, m_TRTTrksBinnedRatioTRT[eta_bin] );

           if ( ratioTrk >= trRatioTRT ) pass = true; // TRT Track cuts
         }
       } // end of track cuts for isConv

     } else {
       //The cuts below are necessary for the V0 track selection
       const AmgSymMatrix(5)& err = *perigee->covariance();
       double sd0sq = err(0,0);
       double sd0 = (sd0sq>0.)?std::sqrt(sd0sq):0.;
       double sz0sq = err(1,1);
       double sz0 = (sz0sq>0.)?std::sqrt(sz0sq):0.;
       if(nclus == 0){
	       if(std::fabs(d0)>=m_sD0_Trt*sd0 && std::fabs(d0)<=m_maxTrtD0 && std::fabs(z0)<=m_sZ0_Trt*sz0 && pt>=m_minPt) pass = true;
       }else{
	       if(std::fabs(d0)>=m_sD0_Si*sd0 && std::fabs(z0)<=m_maxSiZ0 && pt>=m_minPt) pass = true;
       }
       ratioTrk = 1.0;
       if(ntrt > 0 && (!m_PIDonlyForXe || nTrtXenonHits==ntrt) ) { // only check TRT PID if m_PIDonlyForXe is false or all TRT hits are Xenon hits
          ATH_MSG_FATAL( "eProbabilityHT not available for Trk::Tracks only xAOD::TrackParticle objects" );
       }
       if(ratioTrk>m_trRatioV0) pass = false;
     }

   } else pass = false;
   if (not vertexSuppliedByUser) delete myVertex;
   if (perigee!=track.perigeeParameters()) delete perigee;
   return pass;
 }

 bool InDetConversionTrackSelectorTool::decision(const Trk::TrackParticleBase& track,const Trk::Vertex* vx) const
 {
   const EventContext& ctx = Gaudi::Hive::currentContext();
   bool pass = false;
   const Trk::TrackParameters* definintParameters=&(track.definingParameters());
   const Trk::Perigee* perigee=dynamic_cast<const Trk::Perigee*>(definintParameters);
   const Trk::Vertex* myVertex=vx;
   const bool vertexSuppliedByUser{vx!=nullptr};
   //in case no Vertex is provided by the user, beam position will be used if available
   if (not vertexSuppliedByUser) {
     SG::ReadCondHandle<InDet::BeamSpotData> beamSpotHandle{ m_beamSpotKey,
                                                             ctx };
     if (beamSpotHandle.isValid()) {
       myVertex=new Trk::RecVertex(beamSpotHandle->beamVtx());
     } else {
       ATH_MSG_WARNING( " Cannot get beamSpot center from BeamSpotData. Using (0,0,0)... " );
       myVertex=new Trk::Vertex(Amg::Vector3D(0,0,0));
     }
   }

   Trk::PerigeeSurface perigeeSurface(myVertex->position());
   const Trk::TrackParameters *firstmeaspar=nullptr;
   for (unsigned int i=0;i<track.trackParameters().size();i++){
     if ( (track.trackParameters())[i]->covariance() && !dynamic_cast<const Trk::Perigee*>((track.trackParameters())[i])) {
       firstmeaspar=(track.trackParameters())[i];
       break;
     }
   }
   if (!firstmeaspar) {
     //using perigee instead of firstmeasurement, since first measurement was not found...
     firstmeaspar=&(track.definingParameters());
     if (!firstmeaspar){
       ATH_MSG_WARNING( " Track Paraemters at first measurement not found. Perigee not found. Cannot do TrackSelection..." );
       if (not vertexSuppliedByUser) delete myVertex;
       return false;
     }
   }
   const Trk::TrackParameters* extrapolatedParameters =
     m_extrapolator->extrapolate(
       ctx, *firstmeaspar, perigeeSurface, Trk::anyDirection, true, Trk::pion).release();
   perigee = extrapolatedParameters
               ? dynamic_cast<const Trk::Perigee*>(extrapolatedParameters)
               : nullptr;
   if (perigee == nullptr || !perigee->covariance()) {
     ATH_MSG_WARNING( "Track Selector failed to extrapolate track to the vertex: " << myVertex->position() );
     if (extrapolatedParameters!=nullptr) {
       ATH_MSG_WARNING( "The return object of the extrapolator was not a perigee even if a perigeeSurface was used!" );
       delete extrapolatedParameters;
       if (not vertexSuppliedByUser) delete myVertex;
       return false;
     }
     if (not vertexSuppliedByUser) delete myVertex;
     return false;
   }

   double qOverP = perigee->parameters()[Trk::qOverP];
   double pt = std::fabs(1./qOverP)*std::sin(perigee->parameters()[Trk::theta]);
   double d0 = perigee->parameters()[Trk::d0];
   double z0 = perigee->parameters()[Trk::z0];
   const Trk::TrackSummary* tSum = track.trackSummary();
   if(tSum){
     double ratioTrk = 1.0;
     int nclus  = tSum->get(Trk::numberOfPixelHits) + tSum->get(Trk::numberOfSCTHits);
     bool isSilicon = (nclus > 0);
     int nTrtHits       = tSum->get(Trk::numberOfTRTHits);
     int nTrtOutliers   = tSum->get(Trk::numberOfTRTOutliers);
     int ntrt           = nTrtHits + nTrtOutliers;
     int nTrtXenonHits  = tSum->get(Trk::numberOfTRTXenonHits);

     if(m_isConv){
       if(ntrt > 0 && (!m_PIDonlyForXe || nTrtXenonHits==ntrt) ){
         // only check TRT PID if m_PIDonlyForXe is false or all TRT hits are Xenon hits
         ATH_MSG_FATAL( "eProbabilityHT not available for Trk::TrackParticleBase only xAOD::TrackParticle objects" );
       }

       // Start of track cuts
       if ( pt >= m_minPt ) {

         // Silicon track cuts
         if ( isSilicon && (std::fabs(d0)<=m_maxSiD0) && (fabs(z0)<=m_maxSiZ0) ) {
           if((ntrt<=15 && ratioTrk>=m_trRatio1) ||
              (ntrt>15 && ntrt<=25 && ratioTrk>=m_trRatio2) ||
              (ntrt>25 && ratioTrk>=m_trRatio3)) pass = true;
         }

         // TRT-only track cuts
         if ( (not isSilicon) && (std::fabs(d0)<=m_maxTrtD0) && (std::fabs(z0)<=m_maxTrtZ0) ) {

           unsigned int eta_bin = getEtaBin(*perigee);
           // TRT-only Tracks: eProbabilityHT cuts below
           // See InDet Trt Track Scoring Tool for nTRT cut on TRT-only tracks
           double trRatioTRT = std::max( m_trRatioTRT, m_TRTTrksBinnedRatioTRT[eta_bin] );

           if ( ratioTrk >= trRatioTRT ) pass = true; // TRT Track cuts
         }
       } // end of track cuts for isConv

     } else {
       //The cuts below are necessary for the V0 track selection
       const AmgSymMatrix(5)& err = *perigee->covariance();
       double sd0sq = err(0,0);
       double sd0 = (sd0sq>0.)?std::sqrt(sd0sq):0.;
       double sz0sq = err(1,1);
       double sz0 = (sz0sq>0.)?std::sqrt(sz0sq):0.;
       if(nclus == 0){
	       if(std::fabs(d0)>=m_sD0_Trt*sd0 && std::fabs(d0)<= m_maxTrtD0 && std::fabs(z0)<=m_sZ0_Trt*sz0 && pt>=m_minPt) pass = true;
       }else{
	       if(std::fabs(d0)>=m_sD0_Si*sd0 && std::fabs(z0)<=m_maxSiZ0 && pt>=m_minPt) pass = true;
       }

       ratioTrk = 1.0;
       if(ntrt > 0 && (!m_PIDonlyForXe || nTrtXenonHits==ntrt) ) {// only check TRT PID if m_PIDonlyForXe is false or all TRT hits are Xenon hits
          ATH_MSG_FATAL( "eProbabilityHT not available for Trk::TrackParticleBase only xAOD::TrackParticle objects" );
       }
       if(ratioTrk>m_trRatioV0) pass = false;
     }
   } else pass = false;
   if (not vertexSuppliedByUser) delete myVertex;
   if (perigee!=&(track.definingParameters())) delete perigee;

   return pass;
 }

 Amg::Vector3D
 InDetConversionTrackSelectorTool::getPosOrBeamSpot(
   const EventContext& ctx,
   const xAOD::Vertex* vertex) const
 {
   if (vertex) {
     return vertex->position();
   }
   SG::ReadCondHandle<InDet::BeamSpotData> beamSpotHandle{ m_beamSpotKey, ctx };
   if (beamSpotHandle.isValid()) {
     return beamSpotHandle->beamVtx().position();
   } else {
     return Amg::Vector3D(0, 0, 0);
   }
 }

   // ---------------------------------------------------------------------
  bool InDetConversionTrackSelectorTool::decision(const xAOD::TrackParticle& tp,const xAOD::Vertex* vertex) const
  {
    const EventContext& ctx = Gaudi::Hive::currentContext();
    bool pass = false;
    const Trk::Perigee& perigee=tp.perigeeParameters();
    // in case no Vertex is provided by the user, beam position will be used if
    // available
    Trk::PerigeeSurface perigeeSurface(getPosOrBeamSpot(ctx, vertex));
    const Trk::TrackParameters* extrapolatedParameters =
      m_extrapolator->extrapolate(
        ctx, perigee, perigeeSurface, Trk::anyDirection, false, Trk::pion).release();
    if (extrapolatedParameters == nullptr) {
      ATH_MSG_WARNING("Extrapolation to the vertex failed: " << perigeeSurface
                                                             << "\n"
                                                             << perigee);
      return false;
    }
    double qOverP = perigee.parameters()[Trk::qOverP];
    double pt = std::fabs(1./qOverP)*std::sin(perigee.parameters()[Trk::theta]);
    double d0 = extrapolatedParameters->parameters()[Trk::d0];
    double z0 = extrapolatedParameters->parameters()[Trk::z0];

    double ratioTrk = 1.0;
    int nclus  = getCount(tp,xAOD::numberOfPixelHits) + getCount(tp,xAOD::numberOfSCTHits);
    bool isSilicon = (nclus > 0);
    int nTrtHits       = getCount(tp,xAOD::numberOfTRTHits);
    int nTrtOutliers   = getCount(tp,xAOD::numberOfTRTOutliers);
    int ntrt           = nTrtHits + nTrtOutliers;
    int nTrtXenonHits  = getCount(tp,xAOD::numberOfTRTXenonHits);

    if(m_isConv){
     float temp(0);
     if(ntrt > 0 && (!m_PIDonlyForXe || nTrtXenonHits==ntrt) ){
       // only check TRT PID if m_PIDonlyForXe is false or all TRT hits are Xenon hits
       ratioTrk = tp.summaryValue(temp,xAOD::eProbabilityHT)  ?  temp: 0 ;
     }

     // Start of track cuts
     if ( pt >= m_minPt ) {

       // Silicon track cuts
       if ( isSilicon && (std::fabs(d0)<=m_maxSiD0) && (fabs(z0)<=m_maxSiZ0) ) {
         if((ntrt<=15 && ratioTrk>=m_trRatio1) ||
            (ntrt>15 && ntrt<=25 && ratioTrk>=m_trRatio2) ||
            (ntrt>25 && ratioTrk>=m_trRatio3)) pass = true;
       }

       // TRT-only track cuts
       if ( (not isSilicon) && (std::fabs(d0)<=m_maxTrtD0) && (std::fabs(z0)<=m_maxTrtZ0) ) {

         unsigned int eta_bin = getEtaBin(perigee);
         // TRT-only Tracks: eProbabilityHT cuts below
         // See InDet Trt Track Scoring Tool for nTRT cut on TRT-only tracks
         double trRatioTRT = std::max( m_trRatioTRT, m_TRTTrksBinnedRatioTRT[eta_bin] );

         if ( ratioTrk >= trRatioTRT ) pass = true; // TRT Track cuts
       }
     } // end of track cuts for isConv

    } else {
     //The cuts below are necessary for the V0 track selection
     const AmgSymMatrix(5)& err = *perigee.covariance();
     double sd0sq = err(0,0);
     double sd0 = (sd0sq>0.)?std::sqrt(sd0sq):0.;
     double sz0sq = err(1,1);
     double sz0 = (sz0sq>0.)?std::sqrt(sz0sq):0.;
     if(nclus == 0){
       if(std::fabs(d0)>=m_sD0_Trt*sd0 && std::fabs(d0)<= m_maxTrtD0 && std::fabs(z0)<=m_sZ0_Trt*sz0 && pt>=m_minPt) pass = true;
     }else{
       if(std::fabs(d0)>=m_sD0_Si*sd0 && std::fabs(z0)<=m_maxSiZ0 && pt>=m_minPt) pass = true;
     }
     ratioTrk = 1.0;
     float temp(0);
     if(ntrt > 0 && (!m_PIDonlyForXe || nTrtXenonHits==ntrt) ) // only check TRT PID if m_PIDonlyForXe is false or all TRT hits are Xenon hits
      ratioTrk = tp.summaryValue(temp,xAOD::eProbabilityHT)  ?  temp: 0 ;
     if(ratioTrk>m_trRatioV0) pass = false;
    }


    delete extrapolatedParameters;
    return pass;
  }
}//end of namespace definitions
