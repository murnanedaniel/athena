/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

////////////////////////////////////////////////////////////////////////////////
// TrigL2LowPtTrackFitter tool
// -------------------------------
// ATLAS Collaboration
//
// 17.06.2010 Package created
//
// Author: Dmitry Emeliyanov, RAL
// e-mail: D.Emeliyanov@rl.ac.uk
//
////////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <iostream>
#include <string.h>


#include "TrkSurfaces/Surface.h"
#include "TrkSurfaces/PlaneSurface.h"
#include "TrkParameters/TrackParameters.h"
#include "TrkPrepRawData/PrepRawData.h"
#include "TrkDistributedKalmanFilter/TrkBaseNode.h"
#include "TrkDistributedKalmanFilter/TrkFilteringNodes.h"
#include "TrkDistributedKalmanFilter/TrkTrackState.h"
#include "TrkDistributedKalmanFilter/TrkPlanarSurface.h"

#include "TrkRIO_OnTrack/RIO_OnTrack.h"

#include "TrigInDetTrackFitter/TrigL2LowPtTrackFitter.h"

TrigL2LowPtTrackFitter::TrigL2LowPtTrackFitter(const std::string& t, 
						 const std::string& n,
						 const IInterface*  p ): 
  AthAlgTool(t,n,p),
  m_recalibrate(false),
  m_fastExtrapolator("TrigL2FastExtrapolationTool"),
  m_ROTcreator("Trk::RIO_OnTrackCreator/InDetTrigBroadInDetRotCreator"),
  m_extrapolator("Trk::Extrapolator/InDetTrigExtrapolator")
{
  declareInterface< ITrigL2TrackFittingTool >( this );  
  declareProperty( "fastExtrapolator", m_fastExtrapolator, "TrigL2FastExtrapolationTool");
  declareProperty( "useROTs",m_recalibrate = true);
  declareProperty( "ROTcreator", m_ROTcreator, "Trk::RIO_OnTrackCreator/InDetTrigBroadInDetRotCreator");
  declareProperty( "TrackExtrapolatorTool",m_extrapolator, "Trk::Extrapolator/InDetTrigExtrapolator");
}

StatusCode TrigL2LowPtTrackFitter::initialize()
{
  StatusCode sc = AthAlgTool::initialize();

  if(m_recalibrate)
    {
      sc=m_ROTcreator.retrieve();
      if(sc.isFailure())
	{
	  ATH_MSG_ERROR("Could not retrieve RIO_OnTrack creator tool "<<m_ROTcreator);
	  return StatusCode::FAILURE;
	}
    }
  sc=m_fastExtrapolator.retrieve();
  if(sc.isFailure())
    {
      ATH_MSG_ERROR("Could not retrieve extrapolation tool"<<m_fastExtrapolator);
      return sc;
    }
  sc=m_extrapolator.retrieve();
  if(sc.isFailure())
    {
      ATH_MSG_ERROR("Could not retrieve extrapolation tool"<<m_extrapolator);
      return sc;
    }

  ATH_MSG_INFO("TrigL2LowPtTrackFitter initialized ");
  return sc;
}

StatusCode TrigL2LowPtTrackFitter::finalize()
{
  StatusCode sc = AthAlgTool::finalize(); 
  return sc;
}

void TrigL2LowPtTrackFitter::recalibrateFilteringNode(Trk::TrkBaseNode* pN, Trk::TrkTrackState* pTS)
{
  if(pTS->getSurface()==NULL) return;

  AmgSymMatrix(5)* pM = new AmgSymMatrix(5);

  for(int i=0;i<5;i++) 
    for(int j=0;j<5;j++)
      (*pM)(i,j)=pTS->getTrackCovariance(i,j);
  
  const Trk::PlaneSurface& pTrkSB = dynamic_cast<const Trk::PlaneSurface&>(pN->getPrepRawData()->detectorElement()->surface());
  Trk::TrackParameters* pTP=new Trk::AtaPlane(pTS->getTrackState(0),pTS->getTrackState(1),
					      pTS->getTrackState(2),pTS->getTrackState(3),
					      pTS->getTrackState(4),pTrkSB,
					      pM);

  const Trk::RIO_OnTrack* pRIO = m_ROTcreator->correct(*(pN->getPrepRawData()),*pTP);
  pN->updateWithRIO(pRIO);
  
  delete pTP;
  delete pRIO;
}

Trk::TrkTrackState* TrigL2LowPtTrackFitter::extrapolateOffline(Trk::TrkTrackState* pTS, 
                                                               Trk::TrkPlanarSurface* pSB,
                                                               Trk::TrkPlanarSurface* pSE,
                                                               int dir)
{
  //1. create starting parameters
  Trk::TrackParameters* pTP=NULL;
  Trk::TrkTrackState* pTE=NULL;

  if(pSB==NULL)
    {
      // 1a. MeasuredPerigee

      AmgSymMatrix(5)* pM = new AmgSymMatrix(5);
      pM->setZero();

      for(int i=0;i<5;i++) 
	for(int j=0;j<5;j++)
	  (*pM)(i,j)=pTS->getTrackCovariance(i,j);
      const Trk::PerigeeSurface perSurf;
      pTP=new Trk::Perigee(pTS->getTrackState(0),pTS->getTrackState(1),
			   pTS->getTrackState(2),pTS->getTrackState(3),
			   pTS->getTrackState(4),perSurf, pM);
    }
  else
    {
      const Trk::PlaneSurface* pTrkSB = dynamic_cast<const Trk::PlaneSurface*>(pSB->getTrkSurface());
      AmgSymMatrix(5)* pM = new AmgSymMatrix(5);
      pM->setZero();

      for(int i=0;i<5;i++) 
	for(int j=0;j<5;j++)
	  (*pM)(i,j)=pTS->getTrackCovariance(i,j);

      pTP=new Trk::AtaPlane(pTS->getTrackState(0),pTS->getTrackState(1),
			    pTS->getTrackState(2),pTS->getTrackState(3),
			    pTS->getTrackState(4),*pTrkSB,pM);
    }

  // 2. Extrapolation

  const Trk::TrackParameters* predPar = NULL;
  if(dir>0)
    {
      predPar = m_extrapolator->extrapolate(*pTP,*pSE->getTrkSurface(),
					    Trk::alongMomentum,false,
					    Trk::pion);
    }
  else
    {
      if(pSE!=NULL)
	predPar = m_extrapolator->extrapolate(*pTP,*pSE->getTrkSurface(),
					      Trk::oppositeMomentum,false,
					      Trk::pion);
      else
	{
	  Trk::PerigeeSurface perSurf;
	  predPar = m_extrapolator->extrapolate(*pTP,perSurf,
						Trk::oppositeMomentum,false,
						Trk::pion);
	}
    }

  if(predPar!=NULL)
    {
      if(pSE!=NULL)
	{
	  const Trk::AtaPlane* pTPE = dynamic_cast<const Trk::AtaPlane*>(predPar);
	  if(pTPE!=NULL) {
	    // 4. Create new TrackState
	    double Re[5],Ge[5][5];
	    Re[0]=pTPE->parameters()[Trk::locX];Re[1]=pTPE->parameters()[Trk::locY];
	    Re[2]=pTPE->parameters()[Trk::phi];Re[3]=pTPE->parameters()[Trk::theta];
	    Re[4]=pTPE->parameters()[Trk::qOverP];
	    
	    for(int i=0;i<5;i++) 
	      for(int j=0;j<5;j++)
		Ge[i][j]=(*pTPE->covariance())(i,j);
	    
	    pTE=new Trk::TrkTrackState(pTS);
	    pTE->setTrackState(Re);
	    pTE->setTrackCovariance(Ge);
	    pTE->attachToSurface(pSE);
	  }
	  else pTE=NULL;
	}
      else
	{
	  const Trk::Perigee* pTPE = dynamic_cast<const Trk::Perigee*>(predPar);
	  if(pTPE!=NULL) {
	    // 4. Create new TrackState
	    double Re[5],Ge[5][5];
	    Re[0]=pTPE->parameters()[Trk::d0];Re[1]=pTPE->parameters()[Trk::z0];
	    Re[2]=pTPE->parameters()[Trk::phi];Re[3]=pTPE->parameters()[Trk::theta];
	    Re[4]=pTPE->parameters()[Trk::qOverP];
	    
	    for(int i=0;i<5;i++) 
	      for(int j=0;j<5;j++)
		Ge[i][j]=(*pTPE->covariance())(i,j);
	    
	    pTE=new Trk::TrkTrackState(pTS);
	    pTE->setTrackState(Re);
	    pTE->setTrackCovariance(Ge);
	  }
	  else pTE=NULL;
	}
      delete predPar;
    }
  delete pTP;
  return pTE;
}

Trk::TrkTrackState* TrigL2LowPtTrackFitter::fit(Trk::TrkTrackState* pTS, std::vector<Trk::TrkBaseNode*>& vpTrkNodes,
						bool runSmoother)
{

  bool OK=true;
  Trk::TrkPlanarSurface *pSB=NULL;
  
  // 1. Forward filter

  for(std::vector<Trk::TrkBaseNode*>::iterator pnIt=vpTrkNodes.begin();pnIt!=vpTrkNodes.end();++pnIt)
    {
      Trk::TrkPlanarSurface* pSE=(*pnIt)->getSurface();
      Trk::TrkTrackState* pNS=m_fastExtrapolator->extrapolate(pTS,pSB,pSE,false);
      pSB=pSE;
      if(pNS!=NULL)
	{
	  // recalibrateFilteringNode((*pnIt),pNS);

	  (*pnIt)->validateMeasurement(pNS);
		ATH_MSG_DEBUG("dChi2="<<(*pnIt)->getChi2());
	  (*pnIt)->updateTrackState(pNS);
	  delete pTS;
	  pTS=pNS;
	  double Pt=sin(pTS->getTrackState(3))/pTS->getTrackState(4);
	  if(fabs(Pt)<200.0)
	    {
				ATH_MSG_DEBUG("Estimated Pt is too low "<<Pt<<" - skipping fit");
	      delete pTS;
	      OK=false;break;
	    }
	}
      else
	{
	  delete pTS;OK=false;break;
	}
    }
  if(!OK) return NULL;

  if(!runSmoother) return pTS;

  // 2. Backward filter
  // reset covariance
  double Gk[5][5];
  memset(&Gk[0][0],0,sizeof(Gk));
  Gk[0][0]=100.0;Gk[1][1]=100.0;Gk[2][2]=0.01;Gk[3][3]=0.01;Gk[4][4]=1e-6;
  pTS->setTrackCovariance(Gk);

  for(std::vector<Trk::TrkBaseNode*>::reverse_iterator pnrIt=vpTrkNodes.rbegin();pnrIt!=vpTrkNodes.rend();++pnrIt)
    {
      Trk::TrkPlanarSurface* pSE=(*pnrIt)->getSurface();
      Trk::TrkTrackState* pNS=NULL;
      if(pSE!=pSB)
	{		  
	  double C1,C2,dist=0.0;
	  for(int i=0;i<3;i++)
	    {
	      C1=pSB->getCenter()[i];C2=pSE->getCenter()[i];
	      dist+=(C2-C1)*(C2-C1);
	    }
	  dist=sqrt(dist);
	  if(dist>100.0)
	    pNS=extrapolateOffline(pTS,pSB,pSE,-1);
	  else
	    pNS=m_fastExtrapolator->extrapolate(pTS,pSB,pSE,false);
	}
      else 
	pNS=new Trk::TrkTrackState(pTS);

      pSB=pSE;
      if(pNS!=NULL)
	{
	  if(m_recalibrate)
	    {
	      recalibrateFilteringNode((*pnrIt),pNS);
	    }
	  (*pnrIt)->validateMeasurement(pNS);
    ATH_MSG_DEBUG("dChi2="<<(*pnrIt)->getChi2());
	  (*pnrIt)->updateTrackState(pNS);
	  delete pTS;
	  pTS=pNS;
	  double Pt=sin(pTS->getTrackState(3))/pTS->getTrackState(4);
	  if(fabs(Pt)<200.0)
	    {
				ATH_MSG_DEBUG("Estimated Pt is too low "<<Pt<<" - skipping fit");
	      delete pTS;
	      OK=false;break;
	    }
	}
      else
	{
	  delete pTS;OK=false;break;
	}
    }
  if(!OK) return NULL;

  // 3. Extrapolating back to perigee
  Trk::TrkTrackState* pNS=extrapolateOffline(pTS,pSB,NULL,-1);
  delete pTS;pTS=pNS;

  return pTS;
}
