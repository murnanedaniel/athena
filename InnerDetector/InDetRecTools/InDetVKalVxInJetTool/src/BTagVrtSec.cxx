/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

// Header include
#include "InDetVKalVxInJetTool/InDetVKalVxInJetTool.h"
#include "GeoPrimitives/GeoPrimitivesHelpers.h"
//-------------------------------------------------
// Other stuff
#include  "AnalysisUtils/AnalysisMisc.h"
#include  "TrkVKalVrtFitter/TrkVKalVrtFitter.h"

#include "TProfile.h"
#include "TH2D.h"
#include  "TMath.h"

///----------------------------------------------------------------------------------------
///  getVrtSec returns the vector results with the following
///   1) Vertex mass
///   2) Vertex/jet energy fraction
///   3) Number of initially selected 2tr vertices
///   4) Number of selected for vertexing tracks in jet 
///   5) Number of track in secondary vertex
///   6) 3D SV-PV significance with sign
///   7) Jet energy used in (2) calculation 
///   8) Minimal distance between vertex and any material layer
///   9) Transverse vertex/jet energy fraction. Jet pT independent.
///   10) "Product" variable
///   11) "Boost" variable
///
///  Author:  Vadim Kostyukhin (vadim.kostyukhin@cern.ch)
///---------------------------------------------------------------------------------------- 
namespace InDet{

  float median(std::vector<float> &Vec){
    int N=Vec.size();
    if(N==1) return Vec[0];
    if(N==2) return (Vec[0]+Vec[1])/2.;
    if(N==3) return Vec[1];
    if(N>3){
      std::vector<float> tmp(Vec);
      int N05m=(N-1)/2, N05p=N/2;
      //std::sort(tmp.begin(),tmp.end());  //can use nth_element instead of completely sorting, it's quicker
      if(N05m==N05p){ 
         std::nth_element(tmp.begin(),tmp.begin()+N05m,tmp.end());
         return tmp[N05m];
      } else { 
         std::nth_element(tmp.begin(),tmp.begin()+N05m,tmp.end());
         std::nth_element(tmp.begin(),tmp.begin()+N05p,tmp.end());
         return (tmp[N05m]+tmp[N05p])/2.;
      }
    }
    return 0.;
  }
  //
  // Probability-like 2-track B-vertex score function based on track scores.
  // It is distributed approximetely flat in [0,1]
  inline float getVrtScore(int i, int j, std::vector<std::vector<float>> & trkScore){
     if(i==j)return 0.;
     float res = trkScore.at(i)[0] *  trkScore.at(j)[0];
     res  *=    (1.-trkScore[i][1])*(1.-trkScore[j][1]);
     res  *=    (1.-trkScore[i][2])*(1.-trkScore[j][2]);
     return res;
  }

  const double c_vrtBCMassLimit=5500.;  // Mass limit to consider a vertex not coming from B,C-decays


  xAOD::Vertex* InDetVKalVxInJetTool::getVrtSec(const std::vector<const xAOD::TrackParticle*>& inpTrk,
                                                const xAOD::Vertex                           & primVrt,
                                                const TLorentzVector                         & jetDir,
                                                std::vector<double>                          & results,
                                                std::vector<const xAOD::TrackParticle*>      & listSecondTracks,
                                                std::vector<const xAOD::TrackParticle*>      & trkFromV0,
                                                int & nRefPVTrk,
                                                compatibilityGraph_t                         & compatibilityGraph)
  const
  {

      ATH_MSG_DEBUG("getVrtSec() called with xAOD::TrackParticle=" <<inpTrk.size());
      if( inpTrk.size() < 2 ) { return nullptr;} // 0,1 track => nothing to do!

      std::vector<double> errorMatrix,Impact,ImpactError;
      Amg::Vector3D          fitVertex;
      std::vector< std::vector<double> > TrkAtVrt; 
      TLorentzVector    Momentum;
      double Signif3D=0., Chi2=0.;

      std::vector<const xAOD::TrackParticle*> selectedTracks(0);
      results.clear();        
      listSecondTracks.clear();
      nRefPVTrk=0;

      float evtWgt=1.;
      const EventContext & ctx = Gaudi::Hive::currentContext();
      SG::ReadHandle<xAOD::EventInfo> eventInfo{m_eventInfoKey, ctx};
      if (eventInfo.isValid()) {
        if(eventInfo->hasBeamSpotWeight()) evtWgt *= eventInfo->beamSpotWeight();
      } else ATH_MSG_DEBUG("No event info object found!");

      nRefPVTrk = selGoodTrkParticle( inpTrk, primVrt, jetDir, selectedTracks);
      if(m_fillHist){
        Hists& h = getHists();
        h.m_hb_ntrkjet->Fill( (double)selectedTracks.size(), evtWgt);
        h.m_pr_NSelTrkMean->Fill(jetDir.Pt(),(double)selectedTracks.size());
      }
      long int NTracks = (int) (selectedTracks.size());
      if( NTracks < 2 ) { return nullptr;} // 0,1 selected track => nothing to do!

      if(msgLvl(MSG::DEBUG))msg(MSG::DEBUG) << "Number of selected tracks inside jet= " <<NTracks << endmsg;
      
      TLorentzVector MomentumJet = totalMom(selectedTracks);
      if(m_fillHist){
        Hists& h = getHists();
        h.m_hb_jmom->Fill( MomentumJet.E(), evtWgt);
      }


//--------------------------------------------------------------------------------------------
//                    Initial xAOD::TrackParticle list ready

      std::vector<const xAOD::TrackParticle*>  tracksForFit;
      std::vector<double> inpMass(NTracks,m_massPi);
      int Vrt2TrackNumber = select2TrVrt(selectedTracks, tracksForFit, primVrt, jetDir, inpMass, nRefPVTrk,
                   trkFromV0,listSecondTracks,
                   compatibilityGraph,evtWgt);
//
//--- Cleaning
// 
      if( trkFromV0.size() > 1) {
        removeDoubleEntries(trkFromV0);
        AnalysisUtils::Sort::pT (&trkFromV0); 
      }
//---
      removeDoubleEntries(listSecondTracks);
      AnalysisUtils::Sort::pT (&listSecondTracks);
      for(const auto *iv0 : trkFromV0){ auto itf=std::find(selectedTracks.begin(),selectedTracks.end(),iv0);
                                 if(itf!=selectedTracks.end())  selectedTracks.erase(itf);}
//---
      if(msgLvl(MSG::DEBUG))msg(MSG::DEBUG)<<" Found different xAOD tracks in pairs="<< listSecondTracks.size()<<endmsg;
      if(listSecondTracks.size() < 2 ) return nullptr;

//--Ranking of selected tracks
      std::vector<float> trkRank(0);
      for(const auto *tk : listSecondTracks) trkRank.push_back( m_trackClassificator->trkTypeWgts(tk, primVrt, jetDir)[0] );
      while( median(trkRank)<0.3 && trkRank.size()>3 ) {
        int Smallest= std::min_element(trkRank.begin(),trkRank.end()) - trkRank.begin();
        removeEntryInList(listSecondTracks,trkRank,Smallest);
      }
//
//-----------------------------------------------------------------------------------------------------
//            Secondary track list is ready
//            Now common vertex fit
//
      double Signif3DP=0, Signif3DS=0;

      Chi2 =  fitCommonVrt( listSecondTracks, trkRank, primVrt, jetDir, inpMass, fitVertex, errorMatrix, Momentum, TrkAtVrt);
      if( Chi2 < 0 && listSecondTracks.size()>2 ) { // Vertex not reconstructed. Try to remove one track with biggest pt.
        double tpmax=0; int ipmax=-1;
        for(int it=0; it<(int)listSecondTracks.size(); it++) if(tpmax<listSecondTracks[it]->pt()){tpmax=listSecondTracks[it]->pt(); ipmax=it;}
        if(ipmax>=0)removeEntryInList(listSecondTracks,trkRank,ipmax);
        Chi2 =  fitCommonVrt( listSecondTracks, trkRank, primVrt, jetDir, inpMass, fitVertex, errorMatrix, Momentum, TrkAtVrt);
        if(msgLvl(MSG::DEBUG))msg(MSG::DEBUG)<<" Second fitCommonVrt try="<< Chi2<<" Ntrk="<<listSecondTracks.size()<<endmsg;
      }
      if(msgLvl(MSG::DEBUG))msg(MSG::DEBUG)<<" fitCommonVrt result="<< Chi2<<" Ntrk="<<listSecondTracks.size()<<endmsg;
//
      if( Chi2 < 0) return nullptr;
//
// Check jet tracks not in secondary vertex
      std::map<double,const xAOD::TrackParticle*> AdditionalTracks;
      vrtVrtDist(primVrt, fitVertex, errorMatrix, Signif3D);
      if(Signif3D>8.){
        int hitL1=0, nLays=0, hitIBL=0, hitBL=0; 
        for (const auto *i_ntrk : selectedTracks) {
          if(  find( listSecondTracks.begin(), listSecondTracks.end(), i_ntrk) != listSecondTracks.end() ) continue; // Track is used already
          std::vector<float> trkScore=m_trackClassificator->trkTypeWgts(i_ntrk, primVrt, jetDir);
          if( trkScore[0] < 0.1) continue; //Remove very low track HF score
          Signif3DS = m_fitSvc->VKalGetImpact(i_ntrk, fitVertex         , 1, Impact, ImpactError);
          if( Signif3DS > 10.) continue;
          getPixelLayers(i_ntrk , hitIBL , hitBL, hitL1, nLays );
          if( hitIBL<=0 && hitBL<=0 ) continue;                  // No IBL and BL pixel hits => non-precise track
          Signif3DP = m_fitSvc->VKalGetImpact(i_ntrk, primVrt.position(), 1, Impact, ImpactError);
          if(Signif3DP<1.)continue;
          if(m_fillHist){
            Hists& h = getHists();
            h.m_hb_diffPS->Fill( Signif3DP-Signif3DS, evtWgt);
          }
          if(Signif3DP-Signif3DS>4.0) AdditionalTracks[Signif3DP-Signif3DS]=i_ntrk;
        }
      }
//
// Add found tracks and refit
//
      if( !AdditionalTracks.empty()){
        while (AdditionalTracks.size()>3) AdditionalTracks.erase(AdditionalTracks.begin());//Tracks are in increasing DIFF order.
        for (auto atrk : AdditionalTracks) listSecondTracks.push_back(atrk.second);        //3tracks with max DIFF are selected
        trkRank.clear();
        for(const auto *tk : listSecondTracks) trkRank.push_back( m_trackClassificator->trkTypeWgts(tk, primVrt, jetDir)[0] );
        Chi2 =  fitCommonVrt( listSecondTracks, trkRank, primVrt, jetDir, inpMass, fitVertex, errorMatrix, Momentum, TrkAtVrt);
        if(msgLvl(MSG::DEBUG))msg(MSG::DEBUG)<<" Added track fitCommonVrt output="<< Chi2<<endmsg;
        if( Chi2 < 0) return nullptr;
      }
//
//  Saving of results
//
//
      if( listSecondTracks.size()==2 ){         // If there are 2 only tracks
        if(msgLvl(MSG::DEBUG))msg(MSG::DEBUG)<<" Start Ntr=2 vertex check"<<endmsg;
        int Charge=0;
        for (const auto *i_ntrk : listSecondTracks) { Charge +=  (int) i_ntrk->charge();}
        vrtVrtDist(primVrt, fitVertex, errorMatrix, Signif3D);
// Check track pixel hit patterns vs vertex position.
        if(m_useVertexCleaningPix){
          if(!check2TrVertexInPixel(listSecondTracks[0],listSecondTracks[1],fitVertex,errorMatrix)) return nullptr;
        }
// Check track first measured points vs vertex position.
        if(m_useVertexCleaningFMP){
          float hitR1  = listSecondTracks[0]->radiusOfFirstHit();
          float hitR2  = listSecondTracks[1]->radiusOfFirstHit();
          float vrErr  = vrtRadiusError(fitVertex, errorMatrix);
          if(std::abs(hitR1-hitR2)>25.) return nullptr;                                 // Hits in different pixel layers
          if( fitVertex.perp()-std::min(hitR1,hitR2) > 2.*vrErr) return nullptr; // Vertex is behind hit in pixel 
        }
//--------
//
        if(m_fillHist){
          Hists& h = getHists();
          if(Charge){
            h.m_hb_totmass2T1->Fill(Momentum.M(),evtWgt);
          }
          else{
            h.m_hb_totmass2T0->Fill(Momentum.M(),evtWgt);
          }
        }
        if( !Charge && std::abs(Momentum.M()-m_massK0)<15. ) {       // Final rejection of K0
          trkFromV0.push_back(listSecondTracks[0]);
          trkFromV0.push_back(listSecondTracks[1]);
          if( trkFromV0.size() > 1) {
             removeDoubleEntries(trkFromV0);
             AnalysisUtils::Sort::pT (&trkFromV0);
          }
          return nullptr;
        }
        if(msgLvl(MSG::DEBUG))msg(MSG::DEBUG)<<" Ntr=2 vertex check passed"<<endmsg;
      }


      double jetVrtDir = projSV_PV(fitVertex,primVrt,jetDir);
      if(msgLvl(MSG::DEBUG))msg(MSG::DEBUG)<<"Combined SV neg.dir="<<jetVrtDir<<endmsg;
      if(  m_getNegativeTag )
         { if( jetVrtDir>0. )   return nullptr; }
      else if( m_getNegativeTail )
         { ; }
      else 
         { if( jetVrtDir<0. ) return nullptr; } 

      double xvt=fitVertex.x(); double yvt=fitVertex.y();
      double Dist2DBP=std::hypot( xvt-m_beampipeX, yvt-m_beampipeY); 
      double Dist2DBL=std::hypot( xvt-m_xLayerB, yvt-m_yLayerB); 
      double Dist2DL1=std::hypot( xvt-m_xLayer1, yvt-m_yLayer1);
      double Dist2DL2=std::hypot( xvt-m_xLayer2, yvt-m_yLayer2);
      double minDstMat=39.9;
      minDstMat=TMath::Min(minDstMat,std::abs(Dist2DBP-m_beampipeR));
      minDstMat=TMath::Min(minDstMat,std::abs(Dist2DBL-m_rLayerB));
      minDstMat=TMath::Min(minDstMat,std::abs(Dist2DL1-m_rLayer1));
      minDstMat=TMath::Min(minDstMat,std::abs(Dist2DL2-m_rLayer2));
      if(m_existIBL) minDstMat=std::min(minDstMat,std::abs(Dist2DL2-m_rLayer3));  // 4-layer pixel detector
 
 
      vrtVrtDist(primVrt, fitVertex, errorMatrix, Signif3D);
      if(jetVrtDir < 0) Signif3D = -Signif3D;


      results.push_back(Momentum.M());                             //1st
      double eRatio = Momentum.E()/MomentumJet.E();
      results.push_back(  eRatio<0.99999 ? eRatio : 0.99999);      //2nd
      results.push_back((double)Vrt2TrackNumber);                  //3rd
      results.push_back((double)NTracks);                          //4th
      results.push_back((double)listSecondTracks.size());          //5th
      results.push_back(Signif3D);                                 //6th
      results.push_back(MomentumJet.E());                          //7th
      results.push_back(minDstMat);                                //8th
      double nRatio = Momentum.Et(jetDir.Vect())/std::sqrt(MomentumJet.Perp());   nRatio /= (nRatio+4.); 
      results.push_back( nRatio );                                 //9th   New transverse energy ration
      results.push_back((Momentum.M()-2.*m_massPi)*eRatio/m_massB);           //10th   "Product" variable
      results.push_back((Momentum.Pt()/Momentum.M())*(m_massB/jetDir.Pt()) ); //11th   "Boost" variable

      if(m_fillHist){
          // Find highest track Pt with respect to jet direction
          double trackPt, trackPtMax=0.;
          for (int tr=0; tr<(int)listSecondTracks.size(); tr++) {
            if(listSecondTracks[tr]->pt()/jetDir.Pt()>0.5)continue;
            trackPt=pTvsDir(Amg::Vector3D(jetDir.X(),jetDir.Y(),jetDir.Z()) , TrkAtVrt[tr]);
            if(trackPt>trackPtMax)trackPtMax=trackPt;
          }
          Hists& h = getHists();
          h.m_hb_rNdc->Fill( fitVertex.perp(), evtWgt);
          h.m_hb_trkPtMax->Fill( trackPtMax, evtWgt);
          h.m_pr_effVrt->Fill((float)nRefPVTrk,1.);              
          h.m_pr_effVrtEta->Fill( jetDir.Eta(),1.);
          h.m_hb_mom->Fill( MomentumJet.E(), evtWgt);
          h.m_hb_ratio->Fill( results[1], evtWgt); 
          h.m_hb_totmass->Fill( results[0], evtWgt); 
          h.m_hb_nvrt2->Fill( results[2], evtWgt);
          h.m_hb_sig3DTot->Fill( Signif3D, evtWgt);
          h.m_hb_dstToMat->Fill( minDstMat, evtWgt);
          float R=jetDir.DeltaR(TLorentzVector(fitVertex.x()-primVrt.x(),fitVertex.y()-primVrt.y(),
                                               fitVertex.z()-primVrt.z(), 1.e4));
          h.m_hb_deltaRSVPV->Fill( R, evtWgt);
          if(h.m_curTup){
            h.m_curTup->TotM=Momentum.M();
            if(!m_multiVertex){
               h.m_curTup->nNVrt=1;
               h.m_curTup->NVrtNT[0]     =listSecondTracks.size();
               h.m_curTup->NVrtDist2D[0] =fitVertex.perp();
               h.m_curTup->NVrtSig3D[0]  =Signif3D;
               h.m_curTup->NVrtM[0]      =Momentum.M();
               h.m_curTup->NVrtChi2[0]   =Chi2;
               h.m_curTup->NVrtMaxW[0]   =eRatio;
               h.m_curTup->NVrtDR[0]     =R;
            }
          }
       }

//-------------------------------------------------------------------------------------
//Return xAOD::Vertex
       xAOD::Vertex * tmpVertex=new (std::nothrow) xAOD::Vertex();
       if(!tmpVertex)return nullptr;
       tmpVertex->makePrivateStore();
       tmpVertex->setPosition(fitVertex);
       std::vector<float> floatErrMtx; floatErrMtx.resize(errorMatrix.size());
       for(int i=0; i<(int)errorMatrix.size(); i++) floatErrMtx[i]=errorMatrix[i];
       tmpVertex->setCovariance(floatErrMtx);
       tmpVertex->setFitQuality(Chi2, (float)(listSecondTracks.size()*2.-3.));

       std::vector<Trk::VxTrackAtVertex> & tmpVTAV=tmpVertex->vxTrackAtVertex();    tmpVTAV.clear();
       for(int ii=0; ii<(int)listSecondTracks.size(); ii++) {
         AmgSymMatrix(5) CovMtxP;
         CovMtxP.setIdentity(); 
         Trk::Perigee * tmpMeasPer  =  new (std::nothrow) Trk::Perigee( 0.,0., TrkAtVrt[ii][0], TrkAtVrt[ii][1], TrkAtVrt[ii][2],
                                                                Trk::PerigeeSurface(fitVertex), std::move(CovMtxP) );
         tmpVTAV.emplace_back( 1., tmpMeasPer );
         ElementLink<xAOD::TrackParticleContainer> TEL;  TEL.setElement( listSecondTracks[ii] );
         const xAOD::TrackParticleContainer* cont = (const xAOD::TrackParticleContainer* ) (listSecondTracks[ii]->container() );
         TEL.setStorableObject(*cont);
         tmpVertex->addTrackAtVertex(TEL,1.);
       }
       return tmpVertex;


  }




//
//--------------------------------------------------------
//   Template routine for global secondary vertex fitting
//

  template <class Track>
  double InDetVKalVxInJetTool::fitCommonVrt(std::vector<const Track*>& listSecondTracks,
 				  std::vector<float>        & trkRank,
                                  const xAOD::Vertex        & primVrt,
 	                          const TLorentzVector      & jetDir,
                                  std::vector<double>       & inpMass, 
	                          Amg::Vector3D             & fitVertex,
                                  std::vector<double>       & errorMatrix,
	                          TLorentzVector            & Momentum,
		        std::vector< std::vector<double> >  & TrkAtVrt)
  const
 {
      if(msgLvl(MSG::DEBUG))msg(MSG::DEBUG) << "fitCommonVrt() called " <<listSecondTracks.size()<< endmsg;
//preparation
      StatusCode sc;
      std::vector<double> Chi2PerTrk;
      long int           Charge;
      double             Chi2 = 0.;
      Amg::Vector3D      tmpVertex;

      int Outlier=1, i=0;
//
// Start of fit
//
      std::unique_ptr<Trk::IVKalState> state = m_fitSvc->makeState();
      m_fitSvc->setMassInputParticles( inpMass, *state );            // Use pions masses
      sc=VKalVrtFitFastBase(listSecondTracks,fitVertex,*state);          /* Fast crude estimation */
      if(sc.isFailure() || fitVertex.perp() > m_rLayer2*2. ) {    /* No initial estimation */ 
         m_fitSvc->setApproximateVertex(primVrt.x(), primVrt.y(), primVrt.z(),*state); /* Use as starting point */
      } else {
         m_fitSvc->setApproximateVertex(fitVertex.x(),fitVertex.y(),fitVertex.z(),*state); /*Use as starting point*/
      }
//fit itself
      int NTracksVrt = listSecondTracks.size(); double FitProb=0.;
      std::vector<double> trkFitWgt(0);
      for (i=0; i < NTracksVrt-1; i++) {
         if(m_RobustFit)m_fitSvc->setRobustness(m_RobustFit, *state);
         else m_fitSvc->setRobustness(0, *state);
         sc=VKalVrtFitBase(listSecondTracks,fitVertex,Momentum,Charge,
                           errorMatrix,Chi2PerTrk,TrkAtVrt,Chi2,
                           *state, true);
         if(sc.isFailure() ||  Chi2 > 1000000. ) { return -10000.;}    // No fit
         if(m_RobustFit){
           sc=GetTrkFitWeights(trkFitWgt, *state);
           if(sc.isFailure()){ return -10000.;}    // No weights
	   Outlier=std::min_element(trkFitWgt.begin(),trkFitWgt.end())-trkFitWgt.begin();
         }else{
           Outlier = findMax( Chi2PerTrk, trkRank ); 
         }
	 FitProb=TMath::Prob( Chi2, 2*listSecondTracks.size()-3);
	 if(listSecondTracks.size() == 2 )              break;         // Only 2 tracks left
//////////////////////////////
         double signif3Dproj=vrtVrtDist( primVrt, fitVertex, errorMatrix, jetDir);
         if(signif3Dproj<0 && (!m_getNegativeTail) && (!m_getNegativeTag)){
	   double maxDst=-1.e12; int maxT=-1; double minChi2=1.e12;
	   for(int it=0; it<(int)listSecondTracks.size(); it++){
              std::vector<const Track*> tmpList(listSecondTracks);
              tmpList.erase(tmpList.begin()+it);
              sc=VKalVrtFitBase(tmpList,tmpVertex,Momentum,Charge,errorMatrix,Chi2PerTrk,TrkAtVrt,Chi2,*state,true);
              if(sc.isFailure())continue;
              signif3Dproj=vrtVrtDist( primVrt, tmpVertex, errorMatrix, jetDir);
	      if(signif3Dproj>maxDst  && maxDst<10. ){maxDst=signif3Dproj; maxT=it; minChi2=Chi2;}
	      else if(signif3Dproj>0. && maxDst>10. && Chi2<minChi2) {minChi2=Chi2; maxT=it;}
	   }
	   if(maxT>=0){ Outlier=maxT;   removeEntryInList(listSecondTracks,trkRank,Outlier);
                        m_fitSvc->setApproximateVertex(fitVertex.x(),fitVertex.y(),fitVertex.z(),*state);
                        if(msgLvl(MSG::DEBUG))msg(MSG::DEBUG)<<" Remove negative outlier="<< maxT<<" from "
                            <<listSecondTracks.size()+1<<" tracks"<<endmsg;
			continue;}
	 }
////////////////////////////////////////
//	 if( Momentum.m() < c_vrtBCMassLimit) {
//	   if( Chi2PerTrk[Outlier] < m_secTrkChi2Cut && FitProb > 0.001)  break;  // Solution found
//	 }
	 if( FitProb > 0.001) {
	   if( Momentum.M() <c_vrtBCMassLimit) {
	     if( Chi2PerTrk[Outlier] < m_secTrkChi2Cut*m_chiScale[listSecondTracks.size()<10?listSecondTracks.size():10])  break;  // Solution found
           } else {
	     double minM=1.e12; int minT=-1; double minChi2=1.e12;
	     for(int it=0; it<(int)listSecondTracks.size(); it++){
                std::vector<const Track*> tmpList(listSecondTracks);
                tmpList.erase(tmpList.begin()+it);
                sc=VKalVrtFitBase(tmpList,tmpVertex,Momentum,Charge,errorMatrix,Chi2PerTrk,TrkAtVrt,Chi2,*state,true);
                if(sc.isFailure())continue;
		if(projSV_PV(tmpVertex,primVrt,jetDir)<0.)continue; // Drop negative direction 
	        Chi2 += trkRank[it];                                // Remove preferably non-HF-tracks
		if(Momentum.M()<minM  && minM>c_vrtBCMassLimit){minM=Momentum.M(); minT=it; minChi2=Chi2;}
		else if(Momentum.M()<c_vrtBCMassLimit && minM<c_vrtBCMassLimit && Chi2<minChi2){minChi2=Chi2; minT=it;}
	     }
             if(msgLvl(MSG::DEBUG))msg(MSG::DEBUG)<<" Big mass. Remove trk="<<minT<<" New mass="<<minM<<" New Chi2="<<minChi2<<endmsg;
	     if(minT>=0)Outlier=minT;
	   }
	 }
         if(msgLvl(MSG::DEBUG))msg(MSG::DEBUG)<<"SecVrt remove trk="<<Outlier<<" from "<< listSecondTracks.size()<<" tracks"<<endmsg;
         removeEntryInList(listSecondTracks,trkRank,Outlier);
         m_fitSvc->setApproximateVertex(fitVertex.x(),fitVertex.y(),fitVertex.z(),*state); /*Use as starting point*/
      }
//--
      if(msgLvl(MSG::DEBUG))msg(MSG::DEBUG)<<" SecVrt fit converged. Ntr="<< listSecondTracks.size()<<" Chi2="<<Chi2
         <<" Chi2_trk="<<Chi2PerTrk[Outlier]<<" Prob="<<FitProb<<" M="<<Momentum.M()<<" Dir="<<projSV_PV(fitVertex,primVrt,jetDir)<<endmsg;
//--
      if( listSecondTracks.size()==2 ){
	 if( Momentum.M() > c_vrtBCMassLimit || FitProb < 0.001 || Chi2PerTrk[Outlier] > m_secTrkChi2Cut) return -10000.;  
      } 
//
//-- To kill remnants of conversion
      double Dist2D=std::sqrt(fitVertex.x()*fitVertex.x()+fitVertex.y()*fitVertex.y());
      if( listSecondTracks.size()==2  && (Dist2D > 20.) && Charge==0 ) {
        double mass_EE   =  massV0( TrkAtVrt,m_massE,m_massE);
        if( mass_EE < 40. ) return -40.;
      }
//-- Test creation of Trk::Track
//      if( listSecondTracks.size()==2 && Charge==0 ) {      
//        std::vector<double> VKPerigee,CovPerigee;
//        sc=m_fitSvc->VKalVrtCvtTool(fitVertex,Momentum,errorMatrix,0,VKPerigee,CovPerigee,*state);
//        if(sc.isSuccess())const Trk::Track* TT = m_fitSvc->CreateTrkTrack(VKPerigee,CovPerigee,*state);
//      }
//      if( listSecondTracks.size()==2  && (Dist2D > 20.) ) {    // Protection against fake vertices    
//         std::vector<double> Impact,ImpactError;double ImpactSignif;
//         ImpactSignif = m_fitSvc->VKalGetImpact(listSecondTracks[0], primVrt.position(), 1, Impact, ImpactError);
//         if( ImpactSignif<4. && ImpactSignif>-3. ) return -41.;
//         ImpactSignif = m_fitSvc->VKalGetImpact(listSecondTracks[1], primVrt.position(), 1, Impact, ImpactError);
//         if( ImpactSignif<4. && ImpactSignif>-3. ) return -41.;
//      }
      return Chi2;
}







//
//
//--------------------------------------------------------
//   Template routine for 2track secondary vertices selection
//

    template <class Track>
    int InDetVKalVxInJetTool::select2TrVrt(std::vector<const Track*>  & selectedTracks,
                                  std::vector<const Track*>            & tracksForFit,
                                  const xAOD::Vertex                   & primVrt,
                                  const TLorentzVector                 & jetDir,
                                  std::vector<double>                  & inpMass, 
                                  int                                  & nRefPVTrk,
                                  std::vector<const Track*>   & trkFromV0,
                                  std::vector<const Track*>   & listSecondTracks,
                                  compatibilityGraph_t        & compatibilityGraph,
                                  float evtWgt)
    const
    {
      std::vector<double> Chi2PerTrk,VKPerigee,CovPerigee;
      std::vector<double> Impact,ImpactError;
      double             Signif3D, Dist2D, jetVrtDir;
      long int           Charge;
      int i,j;
      StatusCode sc;
      Vrt2Tr tmpVrt;
      std::vector<Vrt2Tr> all2TrVrt(0);
      int NTracks = (int) (selectedTracks.size());

//
//  Impact parameters with sign calculations
//
      double SignifR=0.,SignifZ=0.;
      std::vector<int> hitIBL(NTracks,0), hitBL(NTracks,0);
      std::vector<double> TrkSig3D(NTracks);
      std::vector< std::vector<float> > trkScore(NTracks);
      AmgVector(5) tmpPerigee; tmpPerigee.setZero();
      for (i=0; i<NTracks; i++) {
         TrkSig3D[i] = m_fitSvc->VKalGetImpact(selectedTracks[i], primVrt.position(), 1, Impact, ImpactError);
         tmpPerigee = getPerigee(selectedTracks[i])->parameters(); 
         if( sin(tmpPerigee[2]-jetDir.Phi())*Impact[0] < 0 ){ Impact[0] = -std::abs(Impact[0]);}
                                                        else{ Impact[0] =  std::abs(Impact[0]);}
         if(  (tmpPerigee[3]-jetDir.Theta())*Impact[1] < 0 ){ Impact[1] = -std::abs(Impact[1]);}
                                                        else{ Impact[1] =  std::abs(Impact[1]);}
         SignifR = Impact[0]/ std::sqrt(ImpactError[0]);
         SignifZ = Impact[1]/ std::sqrt(ImpactError[2]);
         int hL1=0, nLays=0; getPixelLayers(selectedTracks[i] , hitIBL[i] , hitBL[i], hL1, nLays );
         //----
         trkScore[i]=m_trackClassificator->trkTypeWgts(selectedTracks[i], primVrt, jetDir);
         if(m_fillHist){
            Hists& h = getHists();
            h.m_hb_impactR->Fill( SignifR, evtWgt); 
            h.m_hb_impactZ->Fill( SignifZ, evtWgt); 
            h.m_hb_impactRZ->Fill(SignifR, SignifZ, evtWgt); 
            h.m_hb_impact->Fill( TrkSig3D[i], evtWgt);
            if(i<DevTuple::maxNTrk && h.m_curTup){
                 h.m_curTup->etatrk[i]=selectedTracks[i]->eta();
                 h.m_curTup->p_prob[i]=rankBTrk(selectedTracks[i]->pt(),jetDir.Pt(),0.);
                 h.m_curTup->s_prob[i]=rankBTrk(0.,0.,TrkSig3D[i]); 
                 h.m_curTup->SigR[i]=SignifR; h.m_curTup->SigZ[i]=SignifZ; 
                 h.m_curTup->d0[i]=Impact[0]; h.m_curTup->Z0[i]=Impact[1];
                 h.m_curTup->idMC[i]=getG4Inter(selectedTracks[i]); 
                 if(getIdHF(selectedTracks[i]))h.m_curTup->idMC[i]=2;
                 if(getMCPileup(selectedTracks[i]))h.m_curTup->idMC[i]=3;
                 h.m_curTup->wgtB[i]=trkScore[i][0]; h.m_curTup->wgtL[i]=trkScore[i][1]; h.m_curTup->wgtG[i]=trkScore[i][2]; 
                 h.m_curTup->sig3D[i]=TrkSig3D[i];
                 h.m_curTup->chg[i]=tmpPerigee[4]<0. ? 1: -1;
                 h.m_curTup->ibl[i]=hitIBL[i];
                 h.m_curTup->bl[i]=hitBL[i];
                 h.m_curTup->fhitR[i]=selectedTracks[i]->radiusOfFirstHit();
                 TLorentzVector TLV=selectedTracks[i]->p4();
                 h.m_curTup->pTvsJet[i]=TLV.Perp(jetDir.Vect());
                 TLorentzVector normJ;  normJ.SetPtEtaPhiM(1.,jetDir.Eta(),jetDir.Phi(),0.);
                 h.m_curTup->prodTJ[i]=std::sqrt(TLV.Dot(normJ));
                 h.m_curTup->nVrtT[i]=0;
            }
         }
      }
      if(m_fillHist){
        Hists& h = getHists();
        h.m_curTup->ptjet=jetDir.Perp();
        h.m_curTup->etajet=jetDir.Eta();
        h.m_curTup->phijet=jetDir.Phi();
        h.m_curTup->nTrkInJet=std::min(NTracks,DevTuple::maxNTrk);
      };

      listSecondTracks.reserve(2*NTracks);                 // Reserve memory for single vertex


      Amg::Vector3D iniVrt(0.,0.,0.);
      for (i=0; i<NTracks-1; i++) {
         for (j=i+1; j<NTracks; j++) {
             if(m_multiWithPrimary) {  // For multi-vertex with primary one search
                if(trkScore[i][2] > 0.75)continue;  //---- Use classificator to remove Pileup+Interactions only
                if(trkScore[j][2] > 0.75)continue;
             }else{
                if(trkScore[i][0]==0.)continue;  //Explicitly remove non-classified tracks
                if(trkScore[j][0]==0.)continue;
                if(getVrtScore(i,j,trkScore) < m_cutBVrtScore) continue;
             }
             int badTracks = 0;                                       //Bad tracks identification 
             tracksForFit.resize(2);
             std::unique_ptr<Trk::IVKalState> state = m_fitSvc->makeState();
             m_fitSvc->setMassInputParticles( inpMass, *state );     // Use pion masses for fit
             tracksForFit[0]=selectedTracks[i];
             tracksForFit[1]=selectedTracks[j];
             sc=VKalVrtFitFastBase(tracksForFit,tmpVrt.fitVertex,*state);              /* Fast crude estimation*/
             if( sc.isFailure() || tmpVrt.fitVertex.perp() > m_rLayer2*2. ) {   /* No initial estimation */ 
                iniVrt=primVrt.position();
                if( m_multiWithPrimary ) iniVrt.setZero(); 
             } else {
                jetVrtDir = projSV_PV(tmpVrt.fitVertex,primVrt,jetDir);
                if( m_multiWithPrimary ) jetVrtDir=std::abs(jetVrtDir); /* Always positive when primary vertex is seeked for*/ 
                if( jetVrtDir>0. ) iniVrt=tmpVrt.fitVertex;                /* Good initial estimation */ 
                else               iniVrt=primVrt.position();
             }
             m_fitSvc->setApproximateVertex(iniVrt.x(), iniVrt.y(), iniVrt.z(),*state);
             tmpVrt.i=i; tmpVrt.j=j;
             m_fitSvc->setRobustness(6, *state);
             sc=VKalVrtFitBase(tracksForFit,tmpVrt.fitVertex, tmpVrt.momentum, Charge,
                               tmpVrt.errorMatrix, tmpVrt.chi2PerTrk, tmpVrt.trkAtVrt, tmpVrt.chi2,
                               *state, true);
             if(sc.isFailure())                       continue;          /* No fit */ 
             if(tmpVrt.chi2 > m_sel2VrtChi2Cut)       continue;          /* Bad Chi2 */
             if(std::abs(tmpVrt.fitVertex.z())> 650.)     continue;  // definitely outside of Pixel detector
             Dist2D=tmpVrt.fitVertex.perp(); 
             if(Dist2D    > 180. )             continue;  // can't be from B decay
             double vrErr  = vrtRadiusError(tmpVrt.fitVertex, tmpVrt.errorMatrix);
             if(m_useFrozenVersion){
               if(vrErr>1.5&&getVrtScore(i,j,trkScore) < 4.*m_cutBVrtScore) continue;
             }
//---------------
             double mass_PiPi =  tmpVrt.momentum.M();
             if(mass_PiPi > m_vrt2TrMassLimit)      continue;  // can't be from B decay
             vrtVrtDist(primVrt, tmpVrt.fitVertex, tmpVrt.errorMatrix, Signif3D);
             tmpVrt.signif3D=Signif3D;
             vrtVrtDist2D(primVrt, tmpVrt.fitVertex, tmpVrt.errorMatrix, tmpVrt.signif2D);
//---
             TVector3 SVmPV(tmpVrt.fitVertex.x()-primVrt.x(),tmpVrt.fitVertex.y()-primVrt.y(),tmpVrt.fitVertex.z()-primVrt.z());
             tmpVrt.dRSVPV=jetDir.DeltaR(TLorentzVector(SVmPV, 1.)); //DeltaR SV-PV vs jet
             if(tmpVrt.dRSVPV > m_coneForTag ) continue;  // SV is outside of the jet cone
//---
             jetVrtDir = SVmPV.Dot(jetDir.Vect());
             double vPos=SVmPV.Dot(tmpVrt.momentum.Vect())/tmpVrt.momentum.Rho();
             if((!m_multiWithPrimary) &&(!m_getNegativeTail) && (!m_getNegativeTag) &&  jetVrtDir<0. )  continue; /* secondary vertex behind primary*/
             if(vPos<-100.) continue;                                              /* Secondary vertex is too far behind primary*/
//
// Check track pixel hit patterns vs vertex position.
             if(m_useVertexCleaningPix && !check2TrVertexInPixel(selectedTracks[i],selectedTracks[j],tmpVrt.fitVertex,tmpVrt.errorMatrix)) continue;
// Check track first measured points vs vertex position.
             if(m_useVertexCleaningFMP){
               float ihitR  = selectedTracks[i]->radiusOfFirstHit();
               float jhitR  = selectedTracks[j]->radiusOfFirstHit();
               if(std::abs(ihitR-jhitR)>25.) continue;                            // Hits in different pixel layers
               if( tmpVrt.fitVertex.perp()-std::min(ihitR,jhitR) > 2.*vrErr) continue; // Vertex is behind hit in pixel 
             }
//--------
//
             double signif3Dproj=vrtVrtDist( primVrt, tmpVrt.fitVertex, tmpVrt.errorMatrix, jetDir);
             tmpVrt.signif3DProj=signif3Dproj;
             if(m_fillHist){ 
                Hists& h = getHists();
                double Signif3DSign=Signif3D; if(jetVrtDir<0) Signif3DSign=-Signif3D;
                h.m_hb_signif3D->Fill( Signif3DSign, evtWgt);
                h.m_hb_sig3DNtr->Fill(signif3Dproj,  evtWgt);
             }

             //if( m_multiWithPrimary || m_multiVertex) { // For multivertex
             //  add_edge(i,j,compatibilityGraph);
             //} 
             if( m_multiWithPrimary )   continue;   /* Multivertex with primary one. All below is not needed */
//
//  Check if V0 or material interaction on Pixel layer is present
//
             if( Charge == 0 && Signif3D>8. && mass_PiPi<900.) {
               double mass_PPi  =  massV0( tmpVrt.trkAtVrt,m_massP,m_massPi);
               double mass_EE   =  massV0( tmpVrt.trkAtVrt,m_massE,m_massE);
               if(m_fillHist && !m_multiVertex){
                 Hists& h = getHists();
                 h.m_hb_massEE->Fill( mass_EE, evtWgt);
               } 
               if(       mass_EE <  40.)  { 
                 badTracks = 3;
               }else{
                 if(m_fillHist && !m_multiVertex){
                   Hists& h = getHists();
                   h.m_hb_massPiPi->Fill( mass_PiPi, evtWgt);
                 } /* Total mass with input particles masses*/
                 if(m_fillHist && !m_multiVertex){
                   Hists& h = getHists();
                   h.m_hb_massPPi->Fill( mass_PPi, evtWgt);
                 } 
                 if( std::abs(mass_PiPi-m_massK0) < 22. )  badTracks = 1;
                 if( std::abs(mass_PPi-m_massLam) <  8. )  badTracks = 2;
               }

//
//  Creation of V0 tracks
//
               if(badTracks){
                  std::vector<double> inpMassV0;
                  //Reset VKalVrt settings
                  state = m_fitSvc->makeState();
                                                              //matrix are calculated
                  if( badTracks == 1 ) {  // K0 case
                    inpMassV0.push_back(m_massPi);inpMassV0.push_back(m_massPi);
                    m_fitSvc->setMassInputParticles( inpMassV0, *state );
                    m_fitSvc->setMassForConstraint(m_massK0, *state);
                    m_fitSvc->setCnstType(1, *state);       // Set mass  constraint
                  }
                  if( badTracks == 2 ) {  // Lambda case
                    if( std::abs(1./tmpVrt.trkAtVrt[0][2]) > std::abs(1./tmpVrt.trkAtVrt[1][2]) ) {
                            inpMassV0.push_back(m_massP);inpMassV0.push_back(m_massPi);
                    }else{  inpMassV0.push_back(m_massPi);inpMassV0.push_back(m_massP); }
                    m_fitSvc->setMassInputParticles( inpMassV0, *state );
                    m_fitSvc->setMassForConstraint(m_massLam, *state);
                    m_fitSvc->setCnstType(1, *state);       // Set mass  constraint
                  }
                  if( badTracks == 3 ) {  // Gamma case
                    inpMassV0.push_back(m_massE);inpMassV0.push_back(m_massE);
                    m_fitSvc->setMassInputParticles( inpMassV0, *state );
                    m_fitSvc->setCnstType(12, *state);       // Set 3d angular constraint
                  }
                  m_fitSvc->setApproximateVertex(tmpVrt.fitVertex.x(),tmpVrt.fitVertex.y(),tmpVrt.fitVertex.z(),*state);
                  TLorentzVector MomentumV0;
                  Amg::Vector3D  fitVertexV0;
                  std::vector< std::vector<double> > TrkAtVrtV0; 
                  std::vector<double> errorMatrixV0;
                  double Chi2V0;
                  sc=VKalVrtFitBase(tracksForFit, fitVertexV0, MomentumV0, Charge,
                                    errorMatrixV0,Chi2PerTrk,TrkAtVrtV0,Chi2V0,
                                    *state, true);
                  if(sc.isSuccess()) {
                    sc=m_fitSvc->VKalVrtCvtTool(fitVertexV0,MomentumV0,errorMatrixV0,0,VKPerigee,CovPerigee,*state);
                    if(sc.isSuccess()) {
                      const Trk::Track* TT = m_fitSvc->CreateTrkTrack(VKPerigee,CovPerigee,*state); 
                      double ImpactSignifV0=m_fitSvc->VKalGetImpact(TT, primVrt.position(), 0, Impact, ImpactError, *state);
                      if(m_fillHist){
                        Hists& h = getHists();
                        h.m_hb_impV0->Fill( ImpactSignifV0, evtWgt);
                      }
                      if(ImpactSignifV0>3.0 ) badTracks=0;
                      delete TT;
                    } else { badTracks=0;}
                 }  // else { badTracks=0;}
               }
             }
//
//  Check interactions on material layers
//
	     if(m_useITkMaterialRejection){
	       double xvt = tmpVrt.fitVertex.x(); double yvt = tmpVrt.fitVertex.y();
	       float zvt = tmpVrt.fitVertex.z(); float Rvt = std::hypot(xvt,yvt);
	       int bin = m_ITkPixMaterialMap->FindBin(zvt,Rvt);
	       if(m_ITkPixMaterialMap->GetBinContent(bin)>0) badTracks = 4;
	     }
	     else{
	       float minWgtI = std::min(trkScore[i][2],trkScore[j][2]);
	       if(minWgtI >0.50 && Dist2D > m_beampipeR-vrtRadiusError(tmpVrt.fitVertex, tmpVrt.errorMatrix) ) badTracks = 4;
	     }

//
//-----------------------------------------------
             tmpVrt.badVrt=badTracks;          //
             all2TrVrt.push_back(tmpVrt);      //
//-----------------------------------------------
             if(m_fillHist){
               Hists& h = getHists();
               h.m_hb_r2d->Fill( tmpVrt.fitVertex.perp(), evtWgt);
             }
//
//  Creation of tracks from V0 list
//
             if( badTracks ){
                trkFromV0.push_back(selectedTracks[i]);
                trkFromV0.push_back(selectedTracks[j]);
             }
         }
      } 
//============================================================================
//-- Save results
      listSecondTracks.clear();
      std::map<float,int> trkHF;
      for( auto &vv : all2TrVrt){ 
        if( m_multiWithPrimary || m_multiVertex) add_edge(vv.i,vv.j,compatibilityGraph);
        trkHF[trkScore[vv.i][0]]=vv.i; trkHF[trkScore[vv.j][0]]=vv.j;
      }
      for( auto it : trkHF) { listSecondTracks.push_back(selectedTracks[it.second]); }
//-Debug
      if( m_fillHist ) {
        Hists& h = getHists();
        if ( h.m_curTup ){ 
         h.m_curTup->ewgt=evtWgt;
         for( auto &it : trkHF) { h.m_curTup->itHF[h.m_curTup->NTHF++]=it.second; }
         for( auto &vv : all2TrVrt){ h.m_curTup->nVrtT[vv.i]++; h.m_curTup->nVrtT[vv.j]++; }
         fillVrtNTup(all2TrVrt);
        }
      }
//
//--------------------------------------------------------------------
//-- Post-selection checks 
//--------------------------------------------------------------------
      if(!listSecondTracks.empty() ){ 
        if(m_fillHist){
          Hists& h = getHists();
          h.m_pr_effVrt2tr->Fill((float)nRefPVTrk,1.);
          h.m_pr_effVrt2trEta->Fill( jetDir.Eta(),1.);
        }
      } else if(listSecondTracks.empty()) { if(m_fillHist){
          Hists& h = getHists();
          h.m_pr_effVrt2tr->Fill((float)nRefPVTrk,0.);
          h.m_pr_effVrt2trEta->Fill(jetDir.Eta(),0.);
        }}
      return all2TrVrt.size();
   }




   template <class Track>
   bool InDetVKalVxInJetTool::check2TrVertexInPixel( const Track* p1, const Track* p2,
                                              Amg::Vector3D &fitVertex, std::vector<double> & vrtErr)
   const
   {
	int blTrk[2]={0,0};
	int blP[2]={0,0};
	int l1Trk[2]={0,0};
	int l1P[2]={0,0};
	int l2Trk[2]={0,0};
	int nLays[2]={0,0};
        getPixelLayers( p1, blTrk[0] , l1Trk[0], l2Trk[0], nLays[0] );
        getPixelLayers( p2, blTrk[1] , l1Trk[1], l2Trk[1], nLays[1] );    // Very close to PV. Both b-layer hits are mandatory.
        getPixelProblems(p1, blP[0], l1P[0] );
        getPixelProblems(p2, blP[1], l1P[1] );
        double xvt=fitVertex.x(),  yvt=fitVertex.y(); 
        double radiusError=vrtRadiusError(fitVertex, vrtErr);
        double Dist2DBL=std::hypot(xvt-m_xLayerB, yvt-m_yLayerB);
        if      (Dist2DBL < m_rLayerB-radiusError) {       //----------------------------------------- Inside B-layer
          if( blTrk[0]==0 && blTrk[1]==0) return false;  // No b-layer hits at all, but all expected
	  if( blTrk[0]<1  && l1Trk[0]<1 ) return false;
	  if( blTrk[1]<1  && l1Trk[1]<1 ) return false;
          if(  nLays[0]           <2 )    return false;  // Less than 2 layers on track 0
          if(  nLays[1]           <2 )    return false;  // Less than 2 layers on track 1
	  return true;
        }else if(Dist2DBL > m_rLayerB+radiusError){      //----------------------------------------- Outside b-layer
          if( blTrk[0]>0 && blP[0]==0 && blTrk[1]>0 && blP[1]==0 ) return false;  // Good hit in b-layer is present
        }
// 
// L1 and L2 are considered only if vertex is in acceptance
//
	if(std::abs(fitVertex.z())<400.){
          double Dist2DL1=std::hypot(xvt-m_xLayer1, yvt-m_yLayer1);
          double Dist2DL2=std::hypot(xvt-m_xLayer2, yvt-m_yLayer2);
          if      (Dist2DL1 < m_rLayer1-radiusError) {   //------------------------------------------ Inside 1st-layer
	     if( l1Trk[0]==0 && l1Trk[1]==0 )     return false;  // No L1 hits at all
             if( l1Trk[0]<1  && l2Trk[0]<1  )     return false;  // Less than 1 hits on track 0
             if( l1Trk[1]<1  && l2Trk[1]<1  )     return false;  // Less than 1 hits on track 1
             return true;
          }else if(Dist2DL1 > m_rLayer1+radiusError) {  //------------------------------------------- Outside 1st-layer
	     if( l1Trk[0]>0 && l1P[0]==0 && l1Trk[1]>0 && l1P[1]==0 )       return false;  //  Good L1 hit is present
          }
          
          if      (Dist2DL2 < m_rLayer2-radiusError) {  //------------------------------------------- Inside 2nd-layer
	     if( (l2Trk[0]+l2Trk[1])==0 )  return false;           // At least one L2 hit must be present
          }else if(Dist2DL2 > m_rLayer2+radiusError) {  
	  //   if( (l2Trk[0]+l2Trk[1])>0  )  return false;           // L2 hits are present
	  }     
        } else {
	  int d0Trk[2]={0,0}; 
	  int d1Trk[2]={0,0}; 
	  int d2Trk[2]={0,0}; 
          getPixelDiscs( p1, d0Trk[0] , d1Trk[0], d2Trk[0] );
          getPixelDiscs( p2, d0Trk[1] , d1Trk[1], d2Trk[1] );
          if( d0Trk[0]+d1Trk[0]+d2Trk[0] ==0 )return false;
          if( d0Trk[1]+d1Trk[1]+d2Trk[1] ==0 )return false;
        }
        return true;
   }

}  //end of namespace
