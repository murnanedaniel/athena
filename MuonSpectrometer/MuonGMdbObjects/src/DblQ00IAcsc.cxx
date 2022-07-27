/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

/***************************************************************************
 DB data - Muon Station components
 -----------------------------------------
 ***************************************************************************/

#include "MuonGMdbObjects/DblQ00IAcsc.h"
#include "RDBAccessSvc/IRDBRecordset.h"
#include "RDBAccessSvc/IRDBAccessSvc.h"
#include "AmdcDb/AmdcDb.h"
#include "AmdcDb/AmdcDbRecord.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdlib>

namespace MuonGM
{
DblQ00IAcsc::DblQ00IAcsc()
{
    m_nObj = 0;
    m_d = nullptr;    
}

  DblQ00IAcsc::DblQ00IAcsc(IRDBAccessSvc *pAccessSvc, const std::string & GeoTag, const std::string & GeoNode):
    m_nObj(0) {

    IRDBRecordset_ptr iacsc = pAccessSvc->getRecordsetPtr(getName(),GeoTag, GeoNode);

    if(iacsc->size()>0) {
    
    m_nObj = iacsc->size();
    m_d = new IACSC[m_nObj];
    if (m_nObj == 0) std::cerr<<"NO IAcsc banks in the MuonDD Database"<<std::endl;

    size_t i=0;
    while(i<iacsc->size()) {
	
        m_d[i].version        = (*iacsc)[i]->getInt("VERS");    
        m_d[i].line           = i; 
        m_d[i].jff            = (*iacsc)[i]->getInt("JFF");
        m_d[i].jzz            = (*iacsc)[i]->getInt("JZZ");
        m_d[i].job            = (*iacsc)[i]->getInt("JOB");
        m_d[i].wireLayer      = (*iacsc)[i]->getInt("JLAY");
        m_d[i].tras           = 10.*(*iacsc)[i]->getFloat("TRAS"); // I lines in mm, but ISZT in cm
        m_d[i].traz           = 10.*(*iacsc)[i]->getFloat("TRAZ"); // I lines in mm, but ISZT in cm
        m_d[i].trat           = 10.*(*iacsc)[i]->getFloat("TRAT"); // I lines in mm, but ISZT in cm
        m_d[i].rots           = (*iacsc)[i]->getFloat("ROTS");
        m_d[i].rotz           = (*iacsc)[i]->getFloat("ROTZ");
        m_d[i].rott           = (*iacsc)[i]->getFloat("ROTT");
        sprintf(m_d[i].type,"%s",(*iacsc)[i]->getString("TYP").c_str());
        i++;
    }
  }
  else {
    m_d = new IACSC[0];
    std::cerr<<"NO IAcsc banks in the MuonDD Database"<<std::endl;
  }
}

DblQ00IAcsc::DblQ00IAcsc(AmdcDb* iacsc) :
    m_nObj(0) {
  IRDBRecordset_ptr pIRDBRecordset = iacsc->getRecordsetPtr("ISZT","Amdc");
  std::vector<IRDBRecord*>::const_iterator it = pIRDBRecordset->begin();

  m_nObj = pIRDBRecordset->size();
  m_d = new IACSC[m_nObj];
  if (m_nObj == 0) std::cerr<<"NO IAcsc banks in the AmdcDbRecord"<<std::endl;

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
     m_d[i].line = i; 
     m_d[i].jff = (*it)->getInt("JFF");
     m_d[i].jzz = (*it)->getInt("JZZ");
     m_d[i].job = (*it)->getInt("JOB");
     m_d[i].wireLayer = (*it)->getInt("JLAY");
     m_d[i].tras = 10.*(*it)->getFloat("TRAS");
     m_d[i].traz = 10.*(*it)->getFloat("TRAZ");
     m_d[i].trat = 10.*(*it)->getFloat("TRAT");
     m_d[i].rots = (*it)->getFloat("ROTS");
     m_d[i].rotz = (*it)->getFloat("ROTZ");
     m_d[i].rott = (*it)->getFloat("ROTT");
     sprintf(m_d[i].type,"%s",(*it)->getString("TYP").c_str());
  }
}

DblQ00IAcsc::DblQ00IAcsc(const std::string& asciiFileName) {
  std::cerr<<"IAcsc with asciiFileName = : <"<<asciiFileName<<"> "<<std::endl;
  // open file and count number of lines
  m_nObj=0;
  std::ifstream iacscFile(asciiFileName.c_str());
  if (!iacscFile.is_open()) 
    std::cerr<<" bad ascii file: "<<asciiFileName<<std::endl;
  
  
  m_nObj = std::count(std::istreambuf_iterator<char>(iacscFile),
		      std::istreambuf_iterator<char>(),'\n');
  std::cout<<"Number of lines in the CSc Internal A-line file <"<<asciiFileName<<"> is "<< m_nObj <<std::endl;
  
  
  m_d = new IACSC[m_nObj];
  if (m_nObj == 0) std::cerr<<"NO IAcsc banks in "<<asciiFileName<<std::endl;
  
  int j=0;

  // close and reopen file for input
  iacscFile.close();
  iacscFile.open(asciiFileName.c_str());

  char AlineMarker;
  while ( iacscFile 
          >> AlineMarker 
	  >> m_d[j].type
	  >> m_d[j].jff
	  >> m_d[j].jzz
	  >> m_d[j].job
	  >> m_d[j].wireLayer
	  >> m_d[j].tras
	  >> m_d[j].traz
	  >> m_d[j].trat
	  >> m_d[j].rots
	  >> m_d[j].rotz
	  >> m_d[j].rott
	  )
  {  
      std::cout<<" IAcsc:: line "<<j+1<<" --- jtyp, jff, jzz, job, w-layer "<<m_d[j].type<<" "
	     <<m_d[j].jff<<" "<<m_d[j].jzz  <<" "
	     <<m_d[j].job<<" "<<m_d[j].wireLayer  <<std::endl;
      m_d[j].line = j+1;
      j++;
  }
  

  if (j!=(int)m_nObj) { 
    std::cerr<<"problem with DblQ00IAcsc: j="<<j<<" m_nObj="<<(int)m_nObj<<std::endl; 
  }  

}

DblQ00IAcsc::~DblQ00IAcsc()
{
    if  (m_nObj > 0) delete [] m_d;
}

void DblQ00IAcsc::WriteIAcscToAsciiFile(const std::string& filename)
{
  std::ofstream iacscFile;
  iacscFile.open(filename.c_str());
  for (int j=0;j<(int)m_nObj;j++) {
    iacscFile
        <<"A "
        << m_d[j].type        <<" " 
        << m_d[j].jff         <<" " 
        << m_d[j].jzz         <<" " 
        << m_d[j].job         <<"  "
        << m_d[j].wireLayer   <<"  "
        << m_d[j].tras  <<" "  // here mm !
        << m_d[j].traz  <<" "  // here mm !
        << m_d[j].trat  <<" "  // here mm !
        << m_d[j].rots  <<" " 
        << m_d[j].rotz  <<" " 
        << m_d[j].rott  <<" " 
        << "\n";
  }
  iacscFile.close();  
}

} // end of namespace MuonGM
