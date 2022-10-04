/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "TrigInDetEvent/TrigSiSpacePointBase.h"
#include "TrigInDetPattRecoEvent/TrigInDetSiLayer.h"
#include "TrigInDetPattRecoTools/GNN_Geometry.h"

#include<cmath>
#include<cstring>
#include<algorithm>

TrigFTF_GNN_Layer::TrigFTF_GNN_Layer(const TRIG_INDET_SI_LAYER& ls, float ew, int bin0) : m_layer(ls), m_etaBinWidth(ew) {

  if(m_layer.m_type == 0) {//barrel
    m_r1 = m_layer.m_refCoord;
    m_r2 = m_layer.m_refCoord;
    m_z1 = m_layer.m_minBound;
    m_z2 = m_layer.m_maxBound;
  }
  else {//endcap
    m_r1 = m_layer.m_minBound;
    m_r2 = m_layer.m_maxBound;
    m_z1 = m_layer.m_refCoord;
    m_z2 = m_layer.m_refCoord;
  }

  float t1   = m_z1/m_r1;
  float eta1 = -log(sqrt(1+t1*t1)-t1);


  float t2   = m_z2/m_r2;
  float eta2 = -log(sqrt(1+t2*t2)-t2);

  m_minEta = eta1;
  m_maxEta = eta2;
  
  if(m_maxEta < m_minEta) {
    m_minEta = eta2;
    m_maxEta = eta1;
  }

  m_maxEta += 1e-6;//increasing them slightly to avoid range_check exceptions
  m_minEta -= 1e-6;

  float deltaEta = m_maxEta - m_minEta;

  int binCounter = bin0;

  if(deltaEta < m_etaBinWidth) {
    m_nBins = 1;
    m_bins.push_back(binCounter++);
    m_etaBin = deltaEta;
    if(m_layer.m_type == 0) {//barrel
      m_radii.push_back(m_layer.m_refCoord);
      m_minBinCoord.push_back(m_layer.m_minBound);
      m_maxBinCoord.push_back(m_layer.m_maxBound);
    }
    else {//endcap
      float eta = 0.5*(m_maxEta + m_minEta);
      float r = m_layer.m_refCoord*std::exp(-eta)/(1-std::exp(-eta)*std::exp(-eta));
      m_radii.push_back(r);
      m_minBinCoord.push_back(m_layer.m_minBound);
      m_maxBinCoord.push_back(m_layer.m_maxBound);
    }
  }
  else {
    float nB = static_cast<int>(deltaEta/m_etaBinWidth);
    m_nBins = nB;
    if(deltaEta - m_etaBinWidth*nB > 0.5*m_etaBinWidth) m_nBins++;

    m_etaBin = deltaEta/m_nBins;

    if(m_nBins == 1) {
      m_bins.push_back(binCounter++);
      if(m_layer.m_type == 0) {//barrel
	m_radii.push_back(m_layer.m_refCoord);
	m_minBinCoord.push_back(m_layer.m_minBound);
	m_maxBinCoord.push_back(m_layer.m_maxBound);
      }
      else {//endcap
	float eta = 0.5*(m_maxEta + m_minEta);
	float r = 2*m_layer.m_refCoord*std::exp(-eta)/(1-std::exp(-eta)*std::exp(-eta));
	m_radii.push_back(r);
	m_minBinCoord.push_back(m_layer.m_minBound);
	m_maxBinCoord.push_back(m_layer.m_maxBound);
      }
    }
    else {

      float eta = m_minEta+0.5*m_etaBin;

      for(int i=1;i<=m_nBins;i++) {

	m_bins.push_back(binCounter++);
	
	float e1 = eta - 0.5*m_etaBin;
	float e2 = eta + 0.5*m_etaBin;
	
	if(m_layer.m_type == 0) {//barrel
	  m_radii.push_back(m_layer.m_refCoord);
	  float z1 = 0.5*m_layer.m_refCoord*(std::exp(e1)-std::exp(-e1));
	  m_minBinCoord.push_back(z1);
	  float z2 = 0.5*m_layer.m_refCoord*(std::exp(e2)-std::exp(-e2));
	  m_maxBinCoord.push_back(z2);
	}
	else {//endcap
	  float r = 2*m_layer.m_refCoord*std::exp(-eta)/(1-std::exp(-eta)*std::exp(-eta));
	  m_radii.push_back(r);
	  r = 2*m_layer.m_refCoord*std::exp(-e1)/(1-std::exp(-e1)*std::exp(-e1));
	  m_minBinCoord.push_back(r);
	  r = 2*m_layer.m_refCoord*std::exp(-e2)/(1-std::exp(-e2)*std::exp(-e2));
	  m_maxBinCoord.push_back(r);
	}
	
	eta += m_etaBin;
      }
    }
  }
}

bool TrigFTF_GNN_Layer::verifyBin(const struct TrigFTF_GNN_Layer* pL, int b1, int b2, float min_z0, float max_z0) const {

  float z1min = m_minBinCoord.at(b1);
  float z1max = m_maxBinCoord.at(b1);
  float r1 = m_layer.m_refCoord;

  if(m_layer.m_type == 0 && pL->m_layer.m_type == 0) {//barrel -> barrel

    float tol = 5.0;

    float min_b2 = pL->m_minBinCoord.at(b2);
    float max_b2 = pL->m_maxBinCoord.at(b2);

    float r2 = pL->m_layer.m_refCoord;

    float A = r2/(r2-r1);
    float B = r1/(r2-r1);

    float z0_min = z1min*A - max_b2*B;
    float z0_max = z1max*A - min_b2*B;
    
    if(z0_max < min_z0-tol || z0_min > max_z0+tol) return false;

    return true;
  }

  if(m_layer.m_type == 0 && pL->m_layer.m_type != 0) {//barrel -> endcap

    float tol = 10.0;

    float z2 = pL->m_layer.m_refCoord;
    float r2max = pL->m_maxBinCoord.at(b2);
    float r2min = pL->m_minBinCoord.at(b2);

    if(r2max <= r1) return false;

    if(r2min <= r1) {
      r2min = r1 + 1e-3;
    }

    float z0_max = 0.0;
    float z0_min = 0.0;
    
    if(z2 > 0) { 
      z0_max = (z1max*r2max - z2*r1)/(r2max-r1);//
      z0_min = (z1min*r2min - z2*r1)/(r2min-r1);
    }
    else {
      z0_max = (z1max*r2min - z2*r1)/(r2min-r1);
      z0_min = (z1min*r2max - z2*r1)/(r2max-r1);//
    }

    if(z0_max < min_z0-tol || z0_min > max_z0+tol) return false;
    return true;
  }

  
  if(m_layer.m_type != 0 && pL->m_layer.m_type != 0) {//endcap -> endcap
    
     
    float r1min = m_minBinCoord.at(b1);
    float r2max = pL->m_maxBinCoord.at(b2);

    if(r2max <= r1min) return false;
    
     return true;
   }
  

   return true;
}


int TrigFTF_GNN_Layer::getEtaBin(float zh, float rh) const {
  
  if(m_bins.size() == 1) return m_bins.at(0);

  float t1   = zh/rh;
  float eta = -std::log(std::sqrt(1+t1*t1)-t1);
  

  int idx = static_cast<int>((eta-m_minEta)/m_etaBin);

  
  if(idx < 0) idx = 0;
  if(idx >= static_cast<int>(m_bins.size())) idx = static_cast<int>(m_bins.size())-1;

  return m_bins.at(idx);//index in the global storage
}

float TrigFTF_GNN_Layer::getBinRadius(int idx) const {
  if(idx >= static_cast<int>(m_radii.size())) idx = idx-1;
  if(idx >= static_cast<int>(m_radii.size())){
    std::cout << "TrigFTF_GNN_Layer::getBinRadius Elliot idx = " << idx << ">= static_cast<int>(m_radii.size()) = " << static_cast<int>(m_radii.size()) << " returning 0" << std::endl;
    return 0;
  }
  return m_radii.at(idx);
}

TrigFTF_GNN_Layer::~TrigFTF_GNN_Layer() {
  m_bins.clear();
}

TrigFTF_GNN_Geometry::TrigFTF_GNN_Geometry(const std::vector<TRIG_INDET_SI_LAYER>& layers, const FASTRACK_CONNECTOR* conn) : m_nEtaBins(0) {

  const float min_z0 = -168.0;
  const float max_z0 =  168.0;

  m_etaBinWidth = conn->m_etaBin;

  std::cout<<"TrigFTF_GNN Geometry::initialize, delta_eta="<<m_etaBinWidth<<std::endl;

  for(std::vector<TRIG_INDET_SI_LAYER>::const_iterator it = layers.begin();it!=layers.end();++it) {
    const TrigFTF_GNN_Layer* pL = addNewLayer((*it), m_nEtaBins);
    m_nEtaBins += pL->num_bins();
  }
 
  std::cout<<"total number of TrigFTF_GNN layers "<<m_layMap.size()<<std::endl;

  int nBins = 0;
  
  for(std::vector<TrigFTF_GNN_LAYER*>::iterator it =  m_layArray.begin();it!=m_layArray.end();++it) {
    nBins += (*it)->num_bins();
  }

  std::cout<<"Total "<<nBins<<" eta bins created out of "<<m_nEtaBins<<std::endl;

  //calculating bin tables in the connector...

  for(std::map<int, std::vector<FASTRACK_CONNECTION*> >::const_iterator it = conn->m_connMap.begin();it!=conn->m_connMap.end();++it) {

    const std::vector<FASTRACK_CONNECTION*>& vConn = (*it).second;

    for(std::vector<FASTRACK_CONNECTION*>::const_iterator cIt=vConn.begin();cIt!=vConn.end();++cIt) {
      
      unsigned int src = (*cIt)->m_src;//n2 : the new connectors
      unsigned int dst = (*cIt)->m_dst;//n1
      
      const TrigFTF_GNN_LAYER* pL1 = getTrigFTF_GNN_LayerByKey(dst);
      const TrigFTF_GNN_LAYER* pL2 = getTrigFTF_GNN_LayerByKey(src);
      
      if (pL1==nullptr) {
	std::cout << " skipping invalid dst layer " << dst << std::endl; 
	  continue; 
      }
      if (pL2==nullptr) {
	std::cout << " skipping invalid src layer " << src << std::endl; 
	  continue; 
      }
      int nSrcBins = pL2->m_bins.size();
      int nDstBins = pL1->m_bins.size();
      
      (*cIt)->m_binTable.resize(nSrcBins*nDstBins, 0);

      for(int b1=0;b1<nDstBins;b1++) {//loop over bins in Layer 1
	for(int b2=0;b2<nSrcBins;b2++) {//loop over bins in Layer 2
	  if(!pL1->verifyBin(pL2, b1, b2, min_z0, max_z0)) continue;
	  int address = b1 + b2*nDstBins;
	  (*cIt)->m_binTable.at(address) = 1;
	}
      }
    }
  }
}

TrigFTF_GNN_Geometry::~TrigFTF_GNN_Geometry() {
  for(std::vector<TrigFTF_GNN_LAYER*>::iterator it =  m_layArray.begin();it!=m_layArray.end();++it) {
    delete (*it);
  }
  m_layMap.clear();m_layArray.clear();
}

const TrigFTF_GNN_LAYER* TrigFTF_GNN_Geometry::getTrigFTF_GNN_LayerByKey(unsigned int key) const {
  std::map<unsigned int, TrigFTF_GNN_LAYER*>::const_iterator it = m_layMap.find(key);
  if(it == m_layMap.end()) {
    std::cout << " TrigFTF_GNN_LAYER* TrigFTF_GNN_Geometry::getTrigFTF_GNN_LayerByKey("<<key<<") returning nullptr" << std::endl; 
    return nullptr;
  }
  return (*it).second;
}

const TrigFTF_GNN_LAYER* TrigFTF_GNN_Geometry::getTrigFTF_GNN_LayerByIndex(int idx) const {
  return m_layArray.at(idx);
}



const TrigFTF_GNN_LAYER* TrigFTF_GNN_Geometry::addNewLayer(const TRIG_INDET_SI_LAYER& l, int bin0) {

  unsigned int layerKey = l.m_subdet;

  float ew = m_etaBinWidth;
  
  TrigFTF_GNN_LAYER* pHL = new TrigFTF_GNN_LAYER(l, ew, bin0);
  
  m_layMap.insert(std::pair<unsigned int, TrigFTF_GNN_LAYER*>(layerKey, pHL));
  m_layArray.push_back(pHL);
  return pHL;
}
