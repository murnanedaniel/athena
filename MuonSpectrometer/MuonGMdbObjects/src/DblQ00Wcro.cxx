/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

/***************************************************************************
 DB data - Muon Station components
 -----------------------------------------
 ***************************************************************************/

#include "MuonGMdbObjects/DblQ00Wcro.h"
#include "RDBAccessSvc/IRDBRecordset.h"
#include "RDBAccessSvc/IRDBAccessSvc.h"
#include "AmdcDb/AmdcDb.h"
#include "AmdcDb/AmdcDbRecord.h"

#include <iostream>
#include <sstream>

namespace MuonGM
{

  DblQ00Wcro::DblQ00Wcro(IRDBAccessSvc *pAccessSvc, const std::string & GeoTag, const std::string & GeoNode):
    m_nObj(0) {

    IRDBRecordset_ptr wcro = pAccessSvc->getRecordsetPtr(getName(),GeoTag, GeoNode);

    if(wcro->size()>0) {
    m_nObj = wcro->size();
    m_d = new WCRO[m_nObj];
    if (m_nObj == 0) std::cerr<<"NO Wcro banks in the MuonDD Database"<<std::endl;

    size_t i=0;
    while(i<wcro->size()) {
        m_d[i].version     = (*wcro)[i]->getInt("VERS");    
        m_d[i].jsta        = (*wcro)[i]->getInt("JSTA");
        m_d[i].num         = (*wcro)[i]->getInt("NUM");
        m_d[i].heightness     = (*wcro)[i]->getFloat("HEIGHTNESS");
        m_d[i].largeness      = (*wcro)[i]->getFloat("LARGENESS");
        m_d[i].thickness      = (*wcro)[i]->getFloat("THICKNESS");
	i++;
    }
  }
  else {
    m_d = new WCRO[0];
    std::cerr<<"NO Wcro banks in the MuonDD Database"<<std::endl;
  }
}

DblQ00Wcro::DblQ00Wcro(AmdcDb* wcro) :
    m_nObj(0) {
  IRDBRecordset_ptr pIRDBRecordset = wcro->getRecordsetPtr(std::string(getObjName()),"Amdc");
  std::vector<IRDBRecord*>::const_iterator it = pIRDBRecordset->begin();

  m_nObj = pIRDBRecordset->size();
  m_d = new WCRO[m_nObj];
  if (m_nObj == 0) std::cerr<<"NO Wcro banks in the AmdcDbRecord"<<std::endl;

  const AmdcDbRecord* pAmdcDbRecord = dynamic_cast<const AmdcDbRecord*>((*it));
  if (pAmdcDbRecord == nullptr){
    std::cerr << "No way to cast in AmdcDbRecord for " << getObjName() << std::endl;
    return;
  }
  
  std::vector< std::string> VariableList = pAmdcDbRecord->getVariableList();
  int ItemTot = VariableList.size() ;
  for(int Item=0 ; Item<ItemTot ; Item++){
    std::string DbVar = VariableList[Item];
  }

  int i = -1;
  it = pIRDBRecordset->begin();
  for( ; it<pIRDBRecordset->end(); ++it){
     pAmdcDbRecord = dynamic_cast<const AmdcDbRecord*>((*it));
     if(pAmdcDbRecord == nullptr){
       std::cerr << "No way to cast in AmdcDbRecord for " << getObjName() << std::endl;
       return;
     }

     i = i + 1;

     m_d[i].version = (*it)->getInt("VERS");    
     m_d[i].jsta = (*it)->getInt("JSTA");
     m_d[i].num = (*it)->getInt("NUM");
     m_d[i].heightness = (*it)->getFloat("HEIGHTNESS");
     m_d[i].largeness = (*it)->getFloat("LARGENESS");
     m_d[i].thickness = (*it)->getFloat("THICKNESS");
  }
}

DblQ00Wcro::~DblQ00Wcro()
{
    delete [] m_d;
}

} // end of namespace MuonGM
