/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

////////////////////////////////////////////////////////////////////////////////
// TrigL2VertexFitter tool
// -------------------------------
// ATLAS Collaboration
//
// 04.07.2006 Package created
//
// a new version of LVL2 vertex fitter based on a Kalman filter with
// decorrelating measurement noise transformation
// 
// this tool will supersede TrigVertexFitter
//
// Author: Dmitry Emeliyanov, RAL
// e-mail: D.Emeliyanov@rl.ac.uk
//
////////////////////////////////////////////////////////////////////////////////


#include <cmath>
#include <iostream>

#include "TrigInDetEvent/TrigL2Vertex.h"
#include "TrigInDetEvent/TrigInDetTrack.h"
#include "TrigInDetEvent/TrigInDetTrackCollection.h"
#include "TrigTimeAlgs/TrigTimerSvc.h"

#include "TrigInDetToolInterfaces/ITrigVertexingTool.h"
#include "TrigVertexFitter/TrigL2VertexFitter.h"

TrigL2VertexFitter::TrigL2VertexFitter(const std::string& t, 
				       const std::string& n,
				       const IInterface*  p ): AthAlgTool(t,n,p), 
							       m_vertexingTool("TrigVertexingTool")
{
  declareInterface< ITrigL2VertexFitter >( this );
  declareProperty( "numberOfIterations", m_numIter=5);
  declareProperty( "maxChi2Increase", m_maxChi2Increase=1000.0);
  declareProperty( "VertexingTool", m_vertexingTool, "TrigVertexingTool");
}

StatusCode TrigL2VertexFitter::initialize()
{
  StatusCode sc = AlgTool::initialize();

  MsgStream athenaLog(msgSvc(), name());

  sc = m_vertexingTool.retrieve();
  if ( sc.isFailure() ) {
    athenaLog << MSG::FATAL << "Unable to locate " <<m_vertexingTool<<endreq;
    return sc;
  } 
  else athenaLog << MSG::INFO << "TrigVertexingTool retrieved"<< endreq;

  athenaLog << MSG::INFO << "Number of iterations is set to " << m_numIter << endreq; 
  ITrigTimerSvc* timerSvc;
  StatusCode scTime = service( "TrigTimerSvc", timerSvc);
  if(scTime.isFailure()) 
    {
      athenaLog << MSG::INFO<< "Unable to locate Service TrigTimerSvc " << endreq;
      m_timers = false;
    } 
  else m_timers = true;  
  if(m_timers) 
    {
      m_timer[0] = timerSvc->addItem("TrigL2Vertex");
      m_timer[1] = timerSvc->addItem("TL2VFindDCA");
      m_timer[2] = timerSvc->addItem("TL2VPrep");
      m_timer[3] = timerSvc->addItem("TL2VChi2");
      m_timer[4] = timerSvc->addItem("TL2VUpd");
      m_timer[5] = timerSvc->addItem("TL2VChi2C");
      m_timer[6] = timerSvc->addItem("TL2VUpdC");
    }
  athenaLog << MSG::INFO << "TrigL2VertexFitter constructed "<< endreq;
  return sc;
}

StatusCode TrigL2VertexFitter::finalize()
{
  StatusCode sc = AlgTool::finalize(); 
  return sc;
}

TrigL2VertexFitter::~TrigL2VertexFitter()
{

}

StatusCode TrigL2VertexFitter::fit(TrigL2Vertex* pV)
{
  int i,iter;
  bool fitFailed=false;

  MsgStream athenaLog(msgSvc(), name());
  int outputLevel = msgSvc()->outputLevel( name() );

  if(m_timers) m_timer[0]->start();

  if (outputLevel <= MSG::DEBUG) 
    {
      athenaLog<<MSG::DEBUG<<"vertex has "<<pV->m_getTracks()->size()<<" tracks"<<endreq;
    }

  if(pV->m_getTracks()->size()<2)
    {
      athenaLog<<MSG::WARNING<<"vertex has less than 2 tracks - fit failed"<<endreq;
      if(m_timers) m_timer[0]->stop();
      return StatusCode::FAILURE;
    }

  std::list<TrigVertexFitInputTrack*>::iterator it,it1,it2;

  if (outputLevel <= MSG::VERBOSE) 
    {
      athenaLog<<MSG::VERBOSE<<"Track list:"<<endreq;
      for(it=pV->m_getTracks()->begin();it!=pV->m_getTracks()->end();++it)
	{
	  athenaLog<<MSG::VERBOSE<<(*it);
	}
    }

  double V[3],dist=-999.0;
  it1=pV->m_getTracks()->begin();it2=it1;++it2;

  if(m_timers) m_timer[1]->start();
  dist=m_vertexingTool->FindClosestApproach((*it1),(*it2),V);
  if(m_timers) m_timer[1]->stop();

  if(dist<-1.0)
    {
      athenaLog<<MSG::DEBUG<<"FindClosestApproach failed - fit failed"<<endreq;
      if(m_timers) m_timer[0]->stop();
      return StatusCode::FAILURE;
    }

  if (outputLevel <= MSG::DEBUG) 
    {
      if ((*it1)->m_getTrigTrack() && (*it2)->m_getTrigTrack()) {
        athenaLog<<MSG::DEBUG<<"Track 1 AlgId="<<(*it1)->m_getTrigTrack()->algorithmId()<<endreq;
        athenaLog<<MSG::DEBUG<<"Track 2 AlgId="<<(*it2)->m_getTrigTrack()->algorithmId()<<endreq;
      }
      athenaLog<<MSG::DEBUG<<"Min dist "<<dist<<" x="<<V[0]<<" y="<<V[1]<<" z="<<V[2]<<endreq;
    }

  if(m_timers) m_timer[2]->start();
  pV->m_prepareForFit();
  for(it=pV->m_getTracks()->begin();it!=pV->m_getTracks()->end();++it)
    {
      (*it)->m_initializeVertex(pV);
    }
  for(i=0;i<3;i++) pV->m_getParametersVector()[i]=V[i];
  if(m_timers) m_timer[2]->stop();

  if(m_timers) 
    {
      m_timer[3]->start();
      m_timer[3]->pause();
      m_timer[4]->start();
      m_timer[4]->pause();
    }

  for(iter=0;iter<m_numIter;iter++)
  {
    pV->m_reset();
    pV->m_Gk[0][0]=100.0;
    pV->m_Gk[1][1]=100.0;
    pV->m_Gk[2][2]=400.0;
    for(it=pV->m_getTracks()->begin();it!=pV->m_getTracks()->end();++it)
      {
	if(m_timers) m_timer[3]->resume();
	double dchi2=(*it)->m_getChi2Distance(pV);
	if(m_timers) m_timer[3]->pause();

	if (outputLevel <= MSG::VERBOSE) 
	  {
	    athenaLog<<MSG::VERBOSE<<"Track "<<(*it)->m_getIndex()<<" dchi2="<<dchi2<<endreq;
	  }
	if(std::isnan(dchi2)||(dchi2<0.0)||(dchi2>m_maxChi2Increase))
	  {
	    fitFailed=true;
	    break;
	  }
	pV->m_addChi2(dchi2);
	pV->m_addNdof(2);
	if(m_timers) m_timer[4]->resume();
	(*it)->m_updateVertex(pV);
	if(m_timers) m_timer[4]->pause();
	if (outputLevel <= MSG::VERBOSE) 
	  {
	    athenaLog<<MSG::VERBOSE<<"Updated vertex"<<endreq;
	    athenaLog<<pV;
	  }
      }
    if(fitFailed) break;
  }
  if(m_timers) 
    {
      m_timer[0]->stop();
      m_timer[1]->stop();
      m_timer[2]->stop();
      m_timer[3]->stop();
      m_timer[4]->stop();
    }
  if(!fitFailed) pV->m_setStatus(1);
  else return StatusCode::FAILURE;
  return StatusCode::SUCCESS;
}

StatusCode TrigL2VertexFitter::fitWithConstraints(TrigL2Vertex* pV)
{
  int i,iter;
  bool fitFailed=false;

  MsgStream athenaLog(msgSvc(), name());
  int outputLevel = msgSvc()->outputLevel( name() );
  if(m_timers) m_timer[0]->start();
  if (outputLevel <= MSG::DEBUG) 
    {
      athenaLog<<MSG::DEBUG<<"vertex has "<<pV->m_getTracks()->size()<<" tracks"<<endreq;
    }

  if(pV->m_getTracks()->size()<2)
    {
      athenaLog<<MSG::WARNING<<"vertex has less than 2 tracks - fit failed"<<endreq;
      if(m_timers) m_timer[0]->stop();
      return StatusCode::FAILURE;
    }

  std::list<TrigVertexFitInputTrack*>::iterator it,it1,it2;

  if (outputLevel <= MSG::VERBOSE) 
    {
      athenaLog<<MSG::VERBOSE<<"Track list:"<<endreq;
      for(it=pV->m_getTracks()->begin();it!=pV->m_getTracks()->end();++it)
	{
	  athenaLog<<MSG::VERBOSE<<(*it);
	}
    }

  double V[3],dist=-999.0;
  it1=pV->m_getTracks()->begin();it2=it1;++it2;

  if(m_timers) m_timer[1]->start();
  dist=m_vertexingTool->FindClosestApproach((*it1),(*it2),V);
  if(m_timers) m_timer[1]->stop();

  if(dist<-1.0)
    {
      athenaLog<<MSG::DEBUG<<"vertex has 2 identical tracks - fit failed"<<endreq;
      if(m_timers) m_timer[0]->stop();
      return StatusCode::FAILURE;
    }

  if (outputLevel <= MSG::DEBUG) 
    {
      athenaLog<<MSG::DEBUG<<"Min dist "<<dist<<" x="<<V[0]<<" y="<<V[1]<<" z="<<V[2]<<endreq;
    }
  if(m_timers) m_timer[2]->start();
  pV->m_prepareForFit();
  for(it=pV->m_getTracks()->begin();it!=pV->m_getTracks()->end();++it)
    {
      (*it)->m_initializeVertex(pV);
    }
  for(i=0;i<3;i++) pV->m_getParametersVector()[i]=V[i];
  if(m_timers) m_timer[2]->stop();
  if(m_timers) 
    {
      m_timer[3]->start();
      m_timer[3]->pause();
      m_timer[4]->start();
      m_timer[4]->pause();
    }
  for(iter=0;iter<m_numIter;iter++)
    {
      pV->m_reset();
      pV->m_Gk[0][0]=100.0;
      pV->m_Gk[1][1]=100.0;
      pV->m_Gk[2][2]=400.0;
      
      for(it=pV->m_getTracks()->begin();it!=pV->m_getTracks()->end();++it)
	{
	  m_timer[3]->resume();
	  double dchi2=(*it)->m_getChi2Distance(pV);
	  m_timer[3]->pause();
	  if (outputLevel <= MSG::VERBOSE) 
	    {
	      athenaLog<<MSG::VERBOSE<<"Track "<<(*it)->m_getIndex()<<" dchi2="<<dchi2<<endreq;
	    }
	  if(std::isnan(dchi2)||(dchi2<0.0)||(dchi2>m_maxChi2Increase))
	    {
	      fitFailed=true;
	      break;
	    }
	  pV->m_addChi2(dchi2);
	  pV->m_addNdof(2);
	  if(m_timers) m_timer[4]->resume();
	  (*it)->m_updateVertex(pV);
	  if(m_timers) m_timer[4]->pause();
	  if (outputLevel <= MSG::VERBOSE) 
	    {
	      athenaLog<<MSG::VERBOSE<<"Updated vertex"<<endreq;
	      athenaLog<<pV;
	    }
	}
      if(fitFailed) break;
      for(std::list<TrigVertexFitConstraint*>::iterator itC=pV->m_getConstraints()->begin();
	  itC!=pV->m_getConstraints()->end();++itC)
	{
	  if(m_timers) m_timer[5]->resume();
	  double dchi2=(*itC)->m_getChi2Distance(pV);
	  if(m_timers) m_timer[5]->pause();
	  if (outputLevel <= MSG::VERBOSE) 
	    {
	      athenaLog<<MSG::VERBOSE<<"Constraint "<<(*itC)->m_getValue()<<" dchi2="<<dchi2<<endreq;
	    }
	  if((dchi2<0.0)||std::isnan(dchi2))
	    {
	      fitFailed=true;
	      break;
	    }
	  pV->m_addChi2(dchi2);
	  pV->m_addNdof(1);
	  if(m_timers) m_timer[6]->resume();
	  (*itC)->m_updateVertex(pV);
	  if(m_timers) m_timer[6]->pause();
	  if (outputLevel <= MSG::VERBOSE) 
	    {
	      athenaLog<<MSG::VERBOSE<<"Updated vertex"<<endreq;
	      athenaLog<<pV;
	    }
	}
      if(fitFailed) break;
    }
  if(m_timers) 
    {
      m_timer[0]->stop();
      m_timer[1]->stop();
      m_timer[2]->stop();
      m_timer[3]->stop();
      m_timer[4]->stop();
      m_timer[5]->stop();
      m_timer[6]->stop();
    }
  if(!fitFailed) pV->m_setStatus(1);
  else return StatusCode::FAILURE;
  return StatusCode::SUCCESS;
}


