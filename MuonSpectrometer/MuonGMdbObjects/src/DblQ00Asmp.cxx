/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

/***************************************************************************
 DB data - Muon Station components
 -----------------------------------------
 ***************************************************************************/

#include "MuonGMdbObjects/DblQ00Asmp.h"
#include "RDBAccessSvc/IRDBRecordset.h"
#include "RDBAccessSvc/IRDBAccessSvc.h"
#include "AmdcDb/AmdcDb.h"
#include "AmdcDb/AmdcDbRecord.h"

#include <iostream>

namespace MuonGM
{

  DblQ00Asmp::DblQ00Asmp(IRDBAccessSvc *pAccessSvc, const std::string & GeoTag, const std::string & GeoNode):
    m_nObj(0) {

    IRDBRecordset_ptr asmp = pAccessSvc->getRecordsetPtr(getName(),GeoTag, GeoNode);

    if(asmp->size()>0) {
    m_nObj = asmp->size();
    m_d = new ASMP[m_nObj];
    if (m_nObj == 0) std::cerr<<"NO Asmp banks in the MuonDD Database"<<std::endl;

    size_t i=0;
    while(i<asmp->size()) {
        m_d[i].version     = (*asmp)[i]->getInt("VERS");    
        m_d[i].indx        = (*asmp)[i]->getInt("INDX");
        m_d[i].n           = (*asmp)[i]->getInt("N");          
        m_d[i].jtyp        = (*asmp)[i]->getInt("JTYP");       
        i++;
    }
  }
  else {
    m_d = new ASMP[0];
    std::cerr<<"NO Asmp banks in the MuonDD Database"<<std::endl;
  }
}

DblQ00Asmp::DblQ00Asmp(AmdcDb* asmp) :
    m_nObj(0) {
  IRDBRecordset_ptr pIRDBRecordset = asmp->getRecordsetPtr(std::string(getObjName()),"Amdc");
  std::vector<IRDBRecord*>::const_iterator it = pIRDBRecordset->begin();

  m_nObj = pIRDBRecordset->size();
  m_d = new ASMP[m_nObj];
  if (m_nObj == 0) std::cerr<<"NO Asmp banks in the AmdcDbRecord"<<std::endl;

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
     m_d[i].indx =  (*it)->getInt("INDX");
     m_d[i].n = (*it)->getInt("N");
     m_d[i].jtyp = (*it)->getInt("JTYP");
  }
}

DblQ00Asmp::~DblQ00Asmp()
{
    delete [] m_d;
}

} // end of namespace MuonGM
