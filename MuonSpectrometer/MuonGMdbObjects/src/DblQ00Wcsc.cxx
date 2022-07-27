/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

/***************************************************************************
 DB data - Muon Station components
 -----------------------------------------
 ***************************************************************************/

#include "MuonGMdbObjects/DblQ00Wcsc.h"
#include "RDBAccessSvc/IRDBRecordset.h"
#include "RDBAccessSvc/IRDBAccessSvc.h"
#include "AmdcDb/AmdcDb.h"
#include "AmdcDb/AmdcDbRecord.h"

#include <iostream>
#include <sstream>
#include <stdexcept>

namespace MuonGM
{

  DblQ00Wcsc::DblQ00Wcsc(IRDBAccessSvc *pAccessSvc, const std::string & GeoTag, const std::string & GeoNode):
    m_nObj(0) {

    IRDBRecordset_ptr wcsc = pAccessSvc->getRecordsetPtr(getName(),GeoTag, GeoNode);

    if(wcsc->size()>0) {
    m_nObj = wcsc->size();
    m_d = new WCSC[m_nObj];
    if (m_nObj == 0) std::cerr<<"NO Wcsc banks in the MuonDD Database"<<std::endl;

    size_t i=0;
    while(i<wcsc->size()) {
        m_d[i].version     = (*wcsc)[i]->getInt("VERS");    
        m_d[i].jsta        = (*wcsc)[i]->getInt("JSTA");    
        m_d[i].laycsc      = (*wcsc)[i]->getInt("LAYCSC");
        m_d[i].ttotal      = (*wcsc)[i]->getFloat("TTOTAL");
        m_d[i].tnomex      = (*wcsc)[i]->getFloat("TNOMEX"); 
        m_d[i].tlag10      = (*wcsc)[i]->getFloat("TLAG10"); 
        m_d[i].wispa       = (*wcsc)[i]->getFloat("WISPA"); 
        m_d[i].dancat      = (*wcsc)[i]->getFloat("DANCAT"); 
        m_d[i].pcatre      = (*wcsc)[i]->getFloat("PCATRE"); 
        m_d[i].gstrip      = (*wcsc)[i]->getFloat("GSTRIP"); 
        m_d[i].wrestr      = (*wcsc)[i]->getFloat("WRESTR"); 
        m_d[i].wflstr      = (*wcsc)[i]->getFloat("WFLSTR"); 
        m_d[i].trrwas      = (*wcsc)[i]->getFloat("TRRWAS"); 
        m_d[i].wroxa       = (*wcsc)[i]->getFloat("WROXA"); 
        m_d[i].groxwi      = (*wcsc)[i]->getFloat("GROXWI"); 
        m_d[i].wgasba      = (*wcsc)[i]->getFloat("WGASBA"); 
        m_d[i].tgasba      = (*wcsc)[i]->getFloat("TGASBA"); 
        m_d[i].wgascu      = (*wcsc)[i]->getFloat("WGASCU"); 
        m_d[i].tgascu      = (*wcsc)[i]->getFloat("TGASCU"); 
        m_d[i].wfixwi      = (*wcsc)[i]->getFloat("WFIXWI"); 
        m_d[i].tfixwi      = (*wcsc)[i]->getFloat("TFIXWI"); 
        m_d[i].pba1wi      = (*wcsc)[i]->getFloat("PBA1WI"); 
        m_d[i].pba2wi      = (*wcsc)[i]->getFloat("PBA2WI"); 
        m_d[i].pba3wi      = (*wcsc)[i]->getFloat("PBA3WI"); 
        m_d[i].psndco      = (*wcsc)[i]->getFloat("PSNDCO");
        m_d[i].azcat = 0.;
        float  azcat = 0.;
        try {
            azcat = (*wcsc)[i]->getFloat("AZCAT");
            m_d[i].azcat =   azcat;
        }
        catch (const std::runtime_error&)
        {
            std::cerr<<" azcat field does not exists !"<<std::endl;
            m_d[i].azcat =   0.;
        }
       i++;
    }
  }
  else {
    m_d = new WCSC[0];
    std::cerr<<"NO Wcsc banks in the MuonDD Database"<<std::endl;
  }
}

DblQ00Wcsc::DblQ00Wcsc(AmdcDb* wcsc) :
    m_nObj(0) {
  IRDBRecordset_ptr pIRDBRecordset = wcsc->getRecordsetPtr(std::string(getObjName()),"Amdc");
  std::vector<IRDBRecord*>::const_iterator it = pIRDBRecordset->begin();

  m_nObj = pIRDBRecordset->size();
  m_d = new WCSC[m_nObj];
  if (m_nObj == 0) std::cerr<<"NO Wcsc banks in the AmdcDbRecord"<<std::endl;

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
     m_d[i].laycsc = (*it)->getInt("LAYCSC");
     m_d[i].ttotal = (*it)->getFloat("TTOTAL");
     m_d[i].tnomex = (*it)->getFloat("TNOMEX"); 
     m_d[i].tlag10 = (*it)->getFloat("TLAG10"); 
     m_d[i].wispa = (*it)->getFloat("WISPA"); 
     m_d[i].dancat = (*it)->getFloat("DANCAT"); 
     m_d[i].pcatre = (*it)->getFloat("PCATRE"); 
     m_d[i].gstrip = (*it)->getFloat("GSTRIP"); 
     m_d[i].wrestr = (*it)->getFloat("WRESTR"); 
     m_d[i].wflstr = (*it)->getFloat("WFLSTR"); 
     m_d[i].trrwas = (*it)->getFloat("TRRWAS"); 
     m_d[i].wroxa = (*it)->getFloat("WROXA"); 
     m_d[i].groxwi = (*it)->getFloat("GROXWI"); 
     m_d[i].wgasba = (*it)->getFloat("WGASBA"); 
     m_d[i].tgasba = (*it)->getFloat("TGASBA"); 
     m_d[i].wgascu = (*it)->getFloat("WGASCU"); 
     m_d[i].tgascu = (*it)->getFloat("TGASCU"); 
     m_d[i].wfixwi = (*it)->getFloat("WFIXWI"); 
     m_d[i].tfixwi = (*it)->getFloat("TFIXWI"); 
     m_d[i].pba1wi = (*it)->getFloat("PBA1WI"); 
     m_d[i].pba2wi = (*it)->getFloat("PBA2WI"); 
     m_d[i].pba3wi = (*it)->getFloat("PBA3WI"); 
     m_d[i].psndco = (*it)->getFloat("PSNDCO");
     m_d[i].azcat = 0.;
     float azcat = 0.;
     if((*it)->getFloat("AZCAT") != 999999999999.) 
     {
       azcat = (*it)->getFloat("AZCAT");
       m_d[i].azcat =   azcat;
     }
     else
     {
       std::cerr<<" azcat field does not exists !"<<std::endl;
       m_d[i].azcat = 0.;
     }
  }
}

DblQ00Wcsc::~DblQ00Wcsc()
{
    delete [] m_d;
}

} // end of namespace MuonGM
