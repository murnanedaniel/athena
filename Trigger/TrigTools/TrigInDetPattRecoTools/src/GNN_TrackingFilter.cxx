/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include<iostream>
#include<cmath>
#include<cstring>

#include "GNN_TrackingFilter.h"

#include<algorithm>
#include<list>


void TrigFTF_GNN_EdgeState::initialize(TrigFTF_GNN_EDGE* pS) {

  m_initialized = true;

  m_J = 0.0;
  m_vs.clear();


  //n2->n1

  float dx = pS->m_n1->m_sp.x() - pS->m_n2->m_sp.x();
  float dy = pS->m_n1->m_sp.y() - pS->m_n2->m_sp.y();
  float L  = sqrt(dx*dx + dy*dy);

  m_s = dy/L;
  m_c = dx/L;

  //transform for extrapolation and update
  // x' =  x*m_c + y*m_s
  // y' = -x*m_s + y*m_c

  m_refY = pS->m_n2->m_sp.r();
  m_refX = pS->m_n2->m_sp.x()*m_c + pS->m_n2->m_sp.y()*m_s;

  //X-state: y, dy/dx, d2y/dx2

  m_X[0] = -pS->m_n2->m_sp.x()*m_s + pS->m_n2->m_sp.y()*m_c;
  m_X[1] = 0.0;
  m_X[2] = 0.0;

  //Y-state: z, dz/dr

  m_Y[0] = pS->m_n2->m_sp.z();
  m_Y[1] = (pS->m_n1->m_sp.z() - pS->m_n2->m_sp.z())/(pS->m_n1->m_sp.r() - pS->m_n2->m_sp.r());

  memset(&m_Cx[0][0], 0, sizeof(m_Cx));
  memset(&m_Cy[0][0], 0, sizeof(m_Cy));

  m_Cx[0][0] = 0.25;
  m_Cx[1][1] = 0.001;
  m_Cx[2][2] = 0.001;

  m_Cy[0][0] = 1.5;
  m_Cy[1][1] = 0.001;

}

void TrigFTF_GNN_EdgeState::clone(const TrigFTF_GNN_EDGE_STATE& st) {

  memcpy(&m_X[0], &st.m_X[0], sizeof(m_X));
  memcpy(&m_Y[0], &st.m_Y[0], sizeof(m_Y));
  memcpy(&m_Cx[0][0], &st.m_Cx[0][0], sizeof(m_Cx));
  memcpy(&m_Cy[0][0], &st.m_Cy[0][0], sizeof(m_Cy));
  m_refX = st.m_refX;
  m_refY = st.m_refY;
  m_c = st.m_c;
  m_s = st.m_s;
  m_J = st.m_J;
  m_vs.clear();
  m_vs.reserve(st.m_vs.size());
  std::copy(st.m_vs.begin(), st.m_vs.end(), std::back_inserter(m_vs));

  m_initialized = true;
}

TrigFTF_GNN_TrackingFilter::TrigFTF_GNN_TrackingFilter(const std::vector<TRIG_INDET_SI_LAYER>& g, std::vector<TrigFTF_GNN_EDGE>& sb) : m_geo(g), m_segStore(sb) {

}

void TrigFTF_GNN_TrackingFilter::followTrack(TrigFTF_GNN_EDGE* pS, TrigFTF_GNN_EDGE_STATE& output) {

  
  if(pS->m_level == -1) return;//already collected

  m_globalStateCounter = 0;

  //create track state

  TrigFTF_GNN_EDGE_STATE* pInitState = &m_stateStore[m_globalStateCounter++];
  
  pInitState->initialize(pS);
  
  m_stateVec.clear();
  
  //recursive branching and propagation

  propagate(pS, *pInitState);


  if(m_stateVec.empty()) return;

  std::sort(m_stateVec.begin(), m_stateVec.end(), TrigFTF_GNN_EDGE_STATE::Compare());

  TrigFTF_GNN_EDGE_STATE* best = (*m_stateVec.begin());


  output.clone(*best);
  
  m_globalStateCounter = 0;
 
}

void TrigFTF_GNN_TrackingFilter::propagate(TrigFTF_GNN_EDGE* pS, TrigFTF_GNN_EDGE_STATE& ts) {

  if(m_globalStateCounter >= MAX_EDGE_STATE) return;
  
  TrigFTF_GNN_EDGE_STATE* p_new_ts = &m_stateStore[m_globalStateCounter++];
  
  TrigFTF_GNN_EDGE_STATE& new_ts = *p_new_ts;
  new_ts.clone(ts);

  new_ts.m_vs.push_back(pS);
  
  bool accepted = update(pS, new_ts); //update using n1 of the segment

  if(!accepted) return;//stop further propagation
  
  int level = pS->m_level;

  std::list<TrigFTF_GNN_EDGE*> lCont;

  for(int nIdx=0;nIdx<pS->m_nNei;nIdx++) {//loop over the neighbours of this segment
    unsigned int nextSegmentIdx = pS->m_vNei[nIdx];
    
    TrigFTF_GNN_EDGE* pN = &(m_segStore.at(nextSegmentIdx));
    
    if(pN->m_level == -1) continue;//already collected
    
    if(pN->m_level == level - 1) {

      lCont.push_back(pN);
    }
  }

  if(lCont.empty()) {//the end of chain

    //store in the vector
    if(m_globalStateCounter < MAX_EDGE_STATE) {

      if(m_stateVec.empty()) {//add the first segment state
        TrigFTF_GNN_EDGE_STATE* p = &m_stateStore[m_globalStateCounter++];
        p->clone(new_ts);
        m_stateVec.push_back(p);
      }
      else {//compare with the best and add
        float best_so_far = (*m_stateVec.begin())->m_J;
        if(new_ts.m_J > best_so_far) {
          TrigFTF_GNN_EDGE_STATE* p = &m_stateStore[m_globalStateCounter++];
          p->clone(new_ts);
          m_stateVec.push_back(p);
        }
      }
    }
  } 
  else {//branching
    int nBranches = 0;
    for(std::list<TrigFTF_GNN_EDGE*>::iterator sIt = lCont.begin();sIt!=lCont.end();++sIt, nBranches++) {
      propagate((*sIt), new_ts);//recursive call
    }
  }
}

bool TrigFTF_GNN_TrackingFilter::update(TrigFTF_GNN_EDGE* pS, TrigFTF_GNN_EDGE_STATE& ts) {

  const float sigma_t = 0.0003;
  const float sigma_w = 0.00009;

  const float sigmaMS = 0.016;

  const float sigma_x = 0.25;//was 0.22
  const float sigma_y = 2.5;//was 1.7

  const float weight_x = 0.5;
  const float weight_y = 0.5;

  const float maxDChi2_x = 60.0;//35.0;
  const float maxDChi2_y = 60.0;//31.0;

  const float add_hit = 14.0;

  if(ts.m_Cx[2][2] < 0.0 || ts.m_Cx[1][1] < 0.0 || ts.m_Cx[0][0] < 0.0) {
    std::cout<<"Negative cov_x"<<std::endl;
  }

  if(ts.m_Cy[1][1] < 0.0 || ts.m_Cy[0][0] < 0.0) {
    std::cout<<"Negative cov_y"<<std::endl;
  }

  //add ms.

  ts.m_Cx[2][2] += sigma_w*sigma_w;
  ts.m_Cx[1][1] += sigma_t*sigma_t;

  int type1 = getLayerType(pS->m_n2->m_sp.layer());

  float t2 = type1 == 0 ? 1.0 + ts.m_Y[1]*ts.m_Y[1] : 1.0 + 1.0/(ts.m_Y[1]*ts.m_Y[1]); 
  float s1 = sigmaMS*t2;
  float s2 = s1*s1;

  s2 *= sqrt(t2);

  ts.m_Cy[1][1] += s2;  

  //extrapolation

  float X[3], Y[2];
  float Cx[3][3], Cy[2][2];

  float refX, refY, mx, my;

  float x, y, z, r;

  x = pS->m_n1->m_sp.x();
  y = pS->m_n1->m_sp.y();
  z = pS->m_n1->m_sp.z();
  r = pS->m_n1->m_sp.r();

  refX =  x*ts.m_c + y*ts.m_s;
  mx   = -x*ts.m_s + y*ts.m_c;//measured X[0]
  refY = r;
  my   = z;//measured Y[0]

  float A = refX - ts.m_refX;
  float B = 0.5*A*A;
  float dr = refY - ts.m_refY;

  X[0] = ts.m_X[0] + ts.m_X[1]*A + ts.m_X[2]*B;
  X[1] = ts.m_X[1] + ts.m_X[2]*A;
  X[2] = ts.m_X[2];

  Cx[0][0] = ts.m_Cx[0][0] + 2*ts.m_Cx[0][1]*A + 2*ts.m_Cx[0][2]*B + A*A*ts.m_Cx[1][1] + 2*A*B*ts.m_Cx[1][2] + B*B*ts.m_Cx[2][2];
  Cx[0][1] = Cx[1][0] = ts.m_Cx[0][1] + ts.m_Cx[1][1]*A + ts.m_Cx[1][2]*B + ts.m_Cx[0][2]*A + A*A*ts.m_Cx[1][2]  + A*B*ts.m_Cx[2][2];
  Cx[0][2] = Cx[2][0] = ts.m_Cx[0][2] + ts.m_Cx[1][2]*A + ts.m_Cx[2][2]*B;   

  Cx[1][1] = ts.m_Cx[1][1] + 2*A*ts.m_Cx[1][2] + A*A*ts.m_Cx[2][2];
  Cx[1][2] = Cx[2][1] =  ts.m_Cx[1][2] + ts.m_Cx[2][2]*A;

  Cx[2][2] = ts.m_Cx[2][2];

  Y[0] = ts.m_Y[0] + ts.m_Y[1]*dr;
  Y[1] = ts.m_Y[1];

  Cy[0][0] = ts.m_Cy[0][0] + 2*ts.m_Cy[0][1]*dr + dr*dr*ts.m_Cy[1][1];
  Cy[0][1] = Cy[1][0] = ts.m_Cy[0][1] + dr*ts.m_Cy[1][1];
  Cy[1][1] = ts.m_Cy[1][1];

  //chi2 test

  float resid_x = mx - X[0];
  float resid_y = my - Y[0];

  float CHx[3] = {Cx[0][0], Cx[0][1], Cx[0][2]};
  float CHy[2] = {Cy[0][0], Cy[0][1]};


  float sigma_rz = 0.0;

  int type = getLayerType(pS->m_n1->m_sp.layer());

  if(type == 0) {//barrel TO-DO: split into barrel Pixel and barrel SCT
    sigma_rz = sigma_y*sigma_y;
  }
  else {
    sigma_rz = sigma_y*ts.m_Y[1];
    sigma_rz = sigma_rz*sigma_rz;
  }

  float Dx = 1.0/(Cx[0][0] + sigma_x*sigma_x);

  float Dy = 1.0/(Cy[0][0] + sigma_rz);

  float dchi2_x = resid_x*resid_x*Dx;
  float dchi2_y = resid_y*resid_y*Dy;


  if(dchi2_x > maxDChi2_x || dchi2_y > maxDChi2_y) {
    return false;
  }

  ts.m_J += add_hit - dchi2_x*weight_x - dchi2_y*weight_y;

  //state update

  float Kx[3] = {Dx*Cx[0][0], Dx*Cx[0][1], Dx*Cx[0][2]};
  float Ky[2] = {Dy*Cy[0][0], Dy*Cy[0][1]};
  
  for(int i=0;i<3;i++) ts.m_X[i] = X[i] + Kx[i]*resid_x;
  for(int i=0;i<2;i++) ts.m_Y[i] = Y[i] + Ky[i]*resid_y;

  for(int i=0;i<3;i++) {
    for(int j=0;j<3;j++) {
      ts.m_Cx[i][j] = Cx[i][j] - Kx[i]*CHx[j];
    }
  }

  for(int i=0;i<2;i++) {
    for(int j=0;j<2;j++) {
      ts.m_Cy[i][j] = Cy[i][j] - Ky[i]*CHy[j];
    }
  }
  ts.m_refX = refX;
  ts.m_refY = refY;
  return true;
}

int TrigFTF_GNN_TrackingFilter::getLayerType(int l) {
  return m_geo.at(l).m_type;
}
