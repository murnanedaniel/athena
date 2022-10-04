/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include <cmath>
#include <iostream>
#include <list>
#include <vector>
#include <algorithm>
#include <iterator>
#include "TrigInDetEvent/TrigSiSpacePointBase.h"
#include "TrigInDetPattRecoEvent/TrigInDetTriplet.h"
#include "TrigInDetPattRecoEvent/TrigInDetSiLayer.h"
#include "GNN_DataStorage.h"
#include "GNN_TrackingFilter.h"

#include "TrigInDetPattRecoTools/TrigTrackSeedGenerator_ITk.h"
#include "TrigInDetPattRecoTools/FasTrackConnector.h"
#include "IRegionSelector/IRoiDescriptor.h"
#include "IRegionSelector/RoiUtil.h"
#include "InDetPrepRawData/PixelCluster.h"


TrigTrackSeedGeneratorITk::TrigTrackSeedGeneratorITk(const TrigCombinatorialSettings& tcs) 
  : m_settings(tcs), 
    m_minDeltaRadius(2.0),
    m_zTol(3.0) {
  
  m_maxDeltaRadius = m_settings.m_doublet_dR_Max;
  m_phiSliceWidth = 2*M_PI/m_settings.m_nMaxPhiSlice;

  m_storage = new TrigFTF_GNN_DataStorage(*m_settings.m_geo, m_phiSliceWidth);

  //mult scatt. variance for doublet matching
  const double radLen = 0.036;
  m_CovMS = std::pow((13.6/m_settings.m_tripletPtMin),2)*radLen;
  const double ptCoeff = 0.29997*1.9972/2.0;// ~0.3*B/2 - assumes nominal field of 2*T
  m_minR_squ = m_settings.m_tripletPtMin*m_settings.m_tripletPtMin/std::pow(ptCoeff,2);
  m_dtPreCut = std::tan(2.0*m_settings.m_tripletDtCut*sqrt(m_CovMS));

}

TrigTrackSeedGeneratorITk::~TrigTrackSeedGeneratorITk() {
  delete m_storage;
  m_storage = nullptr;
}

void TrigTrackSeedGeneratorITk::loadSpacePoints(const std::vector<TrigSiSpacePointBase>& vSP) {

  int nPixels = 0;
  int nStored = 0;

  for(std::vector<TrigSiSpacePointBase>::const_iterator it = vSP.begin();it != vSP.end();++it) {
 
    bool isPixel = (*it).isPixel();

    if(!isPixel) continue;

    nPixels++;

    int result = m_storage->addSpacePoint((*it));
    
    if(result == 0) nStored++;
    
  }
  m_storage->sortByPhi();
  m_storage->generatePhiIndexing(1.5*m_phiSliceWidth);

}

void TrigTrackSeedGeneratorITk::createSeeds(const IRoiDescriptor* roiDescriptor) {

  const int MaxEdges = 1500000;

  const float cut_dphi_min      =-0.012;
  const float cut_dphi_max      = 0.012;
  const float cut_dcurv_min     =-0.001;
  const float cut_dcurv_max     = 0.001;
  const float cut_tau_ratio_min =-0.007;
  const float cut_tau_ratio_max = 0.007;

  //1. loop over stages

  int currentStage = 0;

  std::map<unsigned int, int> n2StageMap;

  const FASTRACK_CONNECTOR& conn = *(m_settings.m_conn);

  std::vector<TrigFTF_GNN_EDGE> edgeStorage;
  
  edgeStorage.resize(MaxEdges);

  int nEdges = 0;


  for(std::map<int, std::vector<FASTRACK_CONNECTION*> >::const_iterator it = conn.m_connMap.begin();it!=conn.m_connMap.end();++it, currentStage++) {
    
    const std::vector<FASTRACK_CONNECTION*> & vConn = (*it).second;

    //2. loop over links

    for(std::vector<FASTRACK_CONNECTION*>::const_iterator cIt=vConn.begin();cIt!=vConn.end();++cIt) {
      
      unsigned int src = (*cIt)->m_src;//n2 : the new connectors
      unsigned int dst = (*cIt)->m_dst;//n1

      const TrigFTF_GNN_LAYER* pL1 = m_settings.m_geo->getTrigFTF_GNN_LayerByKey(dst);
      const TrigFTF_GNN_LAYER* pL2 = m_settings.m_geo->getTrigFTF_GNN_LayerByKey(src);

      if (pL1==nullptr) {
	continue; 
      }
      if (pL2==nullptr) {
	continue; 
      }


      //3. loops over eta-bins

      int nSrcBins = pL2->m_bins.size();
      int nDstBins = pL1->m_bins.size();

      for(int b1=0;b1<nDstBins;b1++) {//loop over bins in Layer 1

	const TrigFTF_GNN_ETA_BIN& B1 = m_storage->getEtaBin(pL1->m_bins.at(b1));

	if(B1.empty()) continue;

	for(int b2=0;b2<nSrcBins;b2++) {//loop over bins in Layer 2

          const TrigFTF_GNN_ETA_BIN& B2 = m_storage->getEtaBin(pL2->m_bins.at(b2));

          if(B2.empty()) continue;

	  //calculated delta Phi for rb1 ---> rb2 extrapolation

	  float deltaPhi = 0.5*m_phiSliceWidth;

	  //loop over nodes in Layer 1

	  unsigned int first_it = 0;

	  for(std::vector<TrigFTF_GNN_NODE*>::const_iterator n1It = B1.m_vn.begin();n1It!=B1.m_vn.end();++n1It) {

	    TrigFTF_GNN_NODE& n1 = *(*n1It);

	    if(n1.m_in.size() == MAX_SEG_PER_NODE) continue;


	    float phi1 = n1.m_sp.phi();
	    float r1 = n1.m_sp.r();
	    float x1 = n1.m_sp.x(); 
	    float y1 = n1.m_sp.y(); 
	    float z1 = n1.m_sp.z();
	    
	    //sliding window phi1 +/- deltaPhi

	    float minPhi = phi1 - deltaPhi;
	    float maxPhi = phi1 + deltaPhi;

	    for(unsigned int n2PhiIdx = first_it; n2PhiIdx<B2.m_vPhiNodes.size();n2PhiIdx++) {
	      
	      float phi2 = B2.m_vPhiNodes.at(n2PhiIdx).first;
	      
	      if(phi2 < minPhi) {
		first_it = n2PhiIdx;
		continue;
	      }
	      if(phi2 > maxPhi) break;

	      TrigFTF_GNN_NODE& n2 = *B2.m_vn.at(B2.m_vPhiNodes.at(n2PhiIdx).second);
	      
	      if(n2.m_out.size() == MAX_SEG_PER_NODE) continue;
	      if(n2.isFull()) continue;
	      

	      float r2 = n2.m_sp.r();
	      float dr = r2 - r1;
	      
	      if(dr < m_minDeltaRadius) {

		continue;
	      }
	      
	      float z2 = n2.m_sp.z();

	      float dz = z2 - z1;
	      float tau = dz/dr;
	      float ftau = fabs(tau);
	      if (ftau > 36.0) {
		
		continue;
	      }
	      
	      if(ftau < n1.m_minCutOnTau) continue;
	      if(ftau < n2.m_minCutOnTau) continue;
	      if(ftau > n1.m_maxCutOnTau) continue;
	      if(ftau > n2.m_maxCutOnTau) continue;
	      
	      
	      
	      float z0 = z1 - r1*tau;
	      
	      if (m_settings.m_doubletFilterRZ) {
				
		if (!RoiUtil::contains( *roiDescriptor, z0, tau)) {
		  continue;
		}
	      }

	      float dx = n2.m_sp.x() - x1;
	      float dy = n2.m_sp.y() - y1;
	      
	      float L2 = 1/(dx*dx+dy*dy);
                  
	      float D = (n2.m_sp.y()*x1 - y1*n2.m_sp.x())/(r1*r2);

	      float kappa = D*D*L2; 
	      
	      if(kappa > 0.4/m_minR_squ) {// was 0.4
		continue;
	      }

	      float curv = D*sqrt(L2);//curvature

	      float df = curv*r2;
	      float df2 = df*df;
	      float dPhi2 = df*(1 + df2*(0.1667 + 0.075*df2));//asinf
                  
	      float ef = curv*r1;
	      float ef2 = ef*ef;
	      float dPhi1 = ef*(1 + ef2*(0.1667 + 0.075*ef2));//asinf
              
	      if(nEdges < MaxEdges) {

		TrigFTF_GNN_EDGE* pE = &(edgeStorage.at(nEdges));

		float* params = pE->m_p;//exp(-eta), curvature, phi1, phi2
	      
		params[0] = std::sqrt(1+tau*tau)-tau;
		params[1] = curv;
		params[2] = phi1 + dPhi1;
		params[3] = phi2 + dPhi2;
	      
		pE->initialize(&n1, &n2);
		pE->m_n1->addIn(nEdges);
		pE->m_n2->addOut(nEdges);
				
		nEdges++;
		
	      }
	    }
	  }
	}
      }
    }
  }


  std::vector<const TrigFTF_GNN_NODE*> vNodes;

  m_storage->getConnectingNodes(vNodes);

  if(vNodes.empty()) return;

  int nNodes = vNodes.size();

  int nLinks = 0;
    
  for(int nodeIdx=0;nodeIdx<nNodes;nodeIdx++) {
      
    const TrigFTF_GNN_NODE* pN = vNodes.at(nodeIdx);

    std::vector<std::pair<float, int > > in_sort, out_sort;
    in_sort.resize(pN->m_in.size());
    out_sort.resize(pN->m_out.size());

    for(int inIdx = 0;inIdx<static_cast<int>(pN->m_in.size());inIdx++) {
      int inEdgeIdx = pN->m_in.at(inIdx);
      TrigFTF_GNN_EDGE* pS = &(edgeStorage.at(inEdgeIdx));
      in_sort[inIdx].second  = inEdgeIdx;
      in_sort[inIdx].first = pS->m_p[0];
    }
    for(int outIdx = 0;outIdx<static_cast<int>(pN->m_out.size());outIdx++) {
      int outEdgeIdx = pN->m_out.at(outIdx);
      TrigFTF_GNN_EDGE* pS = &(edgeStorage.at(outEdgeIdx));
      out_sort[outIdx].second  = outEdgeIdx;
      out_sort[outIdx].first = pS->m_p[0];
    }

    std::sort(in_sort.begin(), in_sort.end());
    std::sort(out_sort.begin(), out_sort.end());

    unsigned int last_out = 0;

    for(unsigned int in_idx=0;in_idx<in_sort.size();in_idx++) {//loop over incoming edges

      int inEdgeIdx = in_sort[in_idx].second;

      TrigFTF_GNN_EDGE* pS = &(edgeStorage.at(inEdgeIdx));

      pS->m_nNei  = 0;
      float tau1  = pS->m_p[0];
      float uat_1 = 1.0/tau1;
      float curv1 = pS->m_p[1];
      float Phi1  = pS->m_p[2];

      for(unsigned int out_idx = last_out;out_idx<out_sort.size();out_idx++) {

	int outEdgeIdx = out_sort[out_idx].second;

	TrigFTF_GNN_EDGE* pNS = &(edgeStorage.at(outEdgeIdx));
	
	
	float tau2 = pNS->m_p[0];
	float tau_ratio = tau2*uat_1 - 1.0;
	
	if(tau_ratio < cut_tau_ratio_min) {
	  last_out = out_idx;
	  continue;
	}
	if(tau_ratio > cut_tau_ratio_max) break;


	float dPhi =  pNS->m_p[3] - Phi1;
        
	if(dPhi<-M_PI) dPhi += 2*M_PI;
	else if(dPhi>M_PI) dPhi -= 2*M_PI;
	
	if(dPhi < cut_dphi_min || dPhi > cut_dphi_max) {
	  continue;
	}
	
	float curv2 = pNS->m_p[1];
	float dcurv = curv2-curv1;
        
	
	if(dcurv < cut_dcurv_min || dcurv > cut_dcurv_max) {
	  continue;
	}
		
	pS->m_vNei[pS->m_nNei++] = outEdgeIdx;
	nLinks++;
	if(pS->m_nNei >= N_SEG_CONNS) break;
      }
    }
  }
  

  const int maxIter = 15;

  int maxLevel = 0;

  int iter = 0;

  
  std::vector<TrigFTF_GNN_EDGE*> v_old;
  
  for(int edgeIndex=0;edgeIndex<nEdges;edgeIndex++) {

    TrigFTF_GNN_EDGE* pS = &(edgeStorage.at(edgeIndex));
    if(pS->m_nNei == 0) continue;
    v_old.push_back(pS);//TO-DO: increment level for segments as they already have at least one neighbour
  }

  for(;iter<maxIter;iter++) {


    //generate proposals
    std::vector<TrigFTF_GNN_EDGE*> v_new;
    v_new.clear();

    for(std::vector<TrigFTF_GNN_EDGE*>::iterator sIt = v_old.begin();sIt != v_old.end();++sIt) {
      
      TrigFTF_GNN_EDGE* pS = (*sIt);
      
      int next_level = pS->m_level;
          
      for(int nIdx=0;nIdx<pS->m_nNei;nIdx++) {
            
        unsigned int nextEdgeIdx = pS->m_vNei[nIdx];
            
        TrigFTF_GNN_EDGE* pN = &(edgeStorage.at(nextEdgeIdx));
            
        if(pS->m_level == pN->m_level) {
          next_level = pS->m_level + 1;
          v_new.push_back(pS);
          break;
        }
      }
      
      pS->m_next = next_level;//proposal
    }
  
    //update

    int nChanges = 0;
      
    for(std::vector<TrigFTF_GNN_EDGE*>::iterator it=v_new.begin();it!=v_new.end();++it) {
        
      TrigFTF_GNN_EDGE* pS = (*it); 
      if(pS->m_next != pS->m_level) {
        nChanges++;
        pS->m_level = pS->m_next;
        if(maxLevel < pS->m_level) maxLevel = pS->m_level;
      }
    }

    if(nChanges == 0) break;


    v_old = std::move(v_new);
    v_new.clear();
  }


  int minLevel = 3;//a triplet + 2 confirmation
  

  std::vector<TrigFTF_GNN_EDGE*> vSeeds;

  for(int edgeIndex=0;edgeIndex<nEdges;edgeIndex++) {
    TrigFTF_GNN_EDGE* pS = &(edgeStorage.at(edgeIndex));

    if(pS->m_level < minLevel) continue;

    vSeeds.push_back(pS);
  }

  m_triplets.clear();

  std::sort(vSeeds.begin(), vSeeds.end(), TrigFTF_GNN_EDGE::CompareLevel());

  if(vSeeds.empty()) return;


  //backtracking

  TrigFTF_GNN_TRACKING_FILTER tFilter(m_settings.m_layerGeometry, edgeStorage);

  INTERNAL_TRIPLET_BUFFER output;

  for(std::vector<TrigFTF_GNN_EDGE*>::iterator sIt=vSeeds.begin();sIt!=vSeeds.end();++sIt) {

    TrigFTF_GNN_EDGE* pS = (*sIt);
    

    if(pS->m_level == -1) continue;

    TrigFTF_GNN_EDGE_STATE rs(false);
    
    tFilter.followTrack((*sIt), rs);

    if(!rs.m_initialized) {
      continue;
    }

    
    
    if(static_cast<int>(rs.m_vs.size()) < minLevel) continue;

    std::vector<const TrigSiSpacePointBase*> vSP;
    
    for(std::vector<TrigFTF_GNN_EDGE*>::reverse_iterator sIt=rs.m_vs.rbegin();sIt!=rs.m_vs.rend();++sIt) {
      
      if((*sIt)->m_level < minLevel-1) 
	(*sIt)->m_level = -1;//mark as collected
      
      if(sIt == rs.m_vs.rbegin()) {
	vSP.push_back(&(*sIt)->m_n1->m_sp);
	vSP.push_back(&(*sIt)->m_n2->m_sp);
      }
      else {
	vSP.push_back(&(*sIt)->m_n2->m_sp);
      }
    }

    if(vSP.size()<3) continue;

    //making triplets

    unsigned int nTriplets = 0;

    for(unsigned int idx_m = 1;idx_m < vSP.size()-1;idx_m++) {
      const TrigSiSpacePointBase& spM = *vSP.at(idx_m);
      const double pS_r = spM.r();
      const double pS_x = spM.x();
      const double pS_y = spM.y();
      const double cosA = pS_x/pS_r;
      const double sinA = pS_y/pS_r;
      
      for(unsigned int idx_o = idx_m+1; idx_o < vSP.size(); idx_o++) {
	const TrigSiSpacePointBase& spO = *vSP.at(idx_o);
	
	double dx = spO.x() - pS_x;
	double dy = spO.y() - pS_y;
	double R2inv = 1.0/(dx*dx+dy*dy);
	double xn = dx*cosA + dy*sinA;
	double yn =-dx*sinA + dy*cosA;

	const double uo = xn*R2inv;
	const double vo = yn*R2inv;

	for(unsigned int idx_i = 0; idx_i < idx_m; idx_i++) {
	  const TrigSiSpacePointBase& spI = *vSP.at(idx_i);
	

	  dx = spI.x() - pS_x;
	  dy = spI.y() - pS_y;
	  R2inv = 1.0/(dx*dx+dy*dy);
	  
	  xn = dx*cosA + dy*sinA;
	  yn =-dx*sinA + dy*cosA;
    
	  const double ui = xn*R2inv;
	  const double vi = yn*R2inv;

	  //1. pT estimate
    
	  const double du = uo - ui;
	  if(du==0.0) continue;
	  const double A = (vo - vi)/du;
	  const double B = vi - A*ui;
	  const double R_squ = (1 + A*A)/(B*B);
	  
	  if(R_squ < 0.5*m_minR_squ) {
	    continue;
	  }

	  //2. d0 cut
    
	  const double fabs_d0 = std::fabs(pS_r*(B*pS_r - A));
	  
	  if(fabs_d0 > m_settings.m_tripletD0Max) {
	    continue;
	  }

	  //3. phi0 cut
	  
	  if (!roiDescriptor->isFullscan()) {
	    const double uc = 2*B*pS_r - A;
	    const double phi0 = atan2(sinA - uc*cosA, cosA + uc*sinA);
	    if ( !RoiUtil::containsPhi( *roiDescriptor, phi0 ) ) {
	      continue;
	    }
	  }

	  //4. add new triplet

	  const double Q = fabs_d0*fabs_d0;

	  TrigInDetTriplet* t = new TrigInDetTriplet(spI, spM, spO, Q);

	  output.push_back(std::pair<double, TrigInDetTriplet*>(Q,t));

	  nTriplets++;
	  if(nTriplets >= m_settings.m_maxTripletBufferLength) break;
	}
	if(nTriplets >= m_settings.m_maxTripletBufferLength) break;
      }
      if(nTriplets >= m_settings.m_maxTripletBufferLength) break;
    }
  }

  if(!output.empty()) storeTriplets(output);

  
}


void TrigTrackSeedGeneratorITk::createSeedsZv() {


}

void TrigTrackSeedGeneratorITk::storeTriplets(INTERNAL_TRIPLET_BUFFER& tripletVec) {
  for(INTERNAL_TRIPLET_BUFFER::iterator it=tripletVec.begin();it!=tripletVec.end();++it) {
    double Q = (*it).first;
    if (m_settings.m_LRTmode) {
      // In LRT mode penalize pixels in Triplets
      if((*it).second->s1().isPixel()) Q+=1000;
      if((*it).second->s2().isPixel()) Q+=1000;
      if((*it).second->s3().isPixel()) Q+=1000;
    } else {
      // In normal (non LRT) mode penalise SSS by 1000, PSS (if enabled) and PPS by 10000
      if((*it).second->s3().isSCT()) {
	Q += (*it).second->s1().isSCT() ? 1000.0 : 10000.0;
      } 
    }
    m_triplets.push_back(std::pair<double, TrigInDetTriplet*>(Q, (*it).second));
  }
}

void TrigTrackSeedGeneratorITk::getSeeds(std::vector<TrigInDetTriplet>& vs) {


  vs.clear();
  vs.reserve(m_triplets.size());
  std::sort(m_triplets.begin(), m_triplets.end(), 
    [](const std::pair<float, TrigInDetTriplet>& A, const std::pair<float, TrigInDetTriplet>& B) {
      return A.first > B.first;
    }
  );
  for(INTERNAL_TRIPLET_BUFFER::reverse_iterator it=m_triplets.rbegin();it!=m_triplets.rend();++it) {
    vs.push_back((*it).second);
  }
  m_triplets.clear();


}
