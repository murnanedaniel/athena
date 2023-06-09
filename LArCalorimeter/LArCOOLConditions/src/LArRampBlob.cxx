/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#include "LArCOOLConditions/LArRampBlob.h" 
#include "AthenaPoolUtilities/CondAttrListCollection.h"
#include "CoralBase/Blob.h"

LArRampBlob::LArRampBlob()
  :m_nChannels(0),
   m_nPoints(0)
{}

LArRampBlob::~LArRampBlob() 
{}


void LArRampBlob::readBlob(const CondAttrListCollection* attrList, MsgStream& msg) {
  m_nChannels=0;
  m_nPoints=0;
  m_pRamp.clear();

  if (!attrList) return;

  CondAttrListCollection::const_iterator gainIt=attrList->begin();
  CondAttrListCollection::const_iterator gainIt_e=attrList->end();
  
  m_pRamp.resize(attrList->size());
  msg << MSG::DEBUG << "Found data for " << attrList->size() << " gains." << endreq;
  
  int blobSize=0;  

  for(;gainIt!=gainIt_e;++gainIt) {
    const unsigned gain=gainIt->first;
    if (gain>=attrList->size() || gain>2) {
      msg << MSG::ERROR << "Found unexpected COOL-channel (=gain) number:" << gain << endreq;
      return; //ERROR
    }
    const coral::AttributeList& attr=gainIt->second;
    const coral::Blob& rampBlob = attr["RampVec"].data<coral::Blob>();
    if (blobSize==0) blobSize=rampBlob.size();
    if (m_nPoints==0) m_nPoints=attr["nPoints"].data<unsigned>();
    
    //Sanity checks:
    if (blobSize!=rampBlob.size()) {
      msg << MSG::ERROR << "Unequal blob size (" << blobSize << "/" 
	       << rampBlob.size() << ")" <<endreq;
      return;
    }
    if (m_nPoints!=attr["nPoints"].data<unsigned>()) {
      msg << MSG::ERROR << "Unequal polynom degree (" << m_nPoints << "/" 
	  << attr["nPoints"].data<unsigned>() << ")" << endreq;
      return;
    }
    
    m_pRamp[gain]=static_cast<const float*>(rampBlob.startingAddress());
  }// end loop over COOL channels

  
  if (m_nPoints==0) {
    msg << MSG::ERROR << "Number of points is zero!" << endreq;
    return;
  }
  m_nChannels=blobSize/(sizeof(float)*m_nPoints);
  msg << MSG::DEBUG << "Found data for " << m_nChannels << endreq;
  return;
}

