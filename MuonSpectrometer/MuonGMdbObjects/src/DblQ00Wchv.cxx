/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

/***************************************************************************
 DB data - Muon Station components
 -----------------------------------------
 Copyright (C) 2004 by ATLAS Collaboration
 ***************************************************************************/

//<doc><file>	$Id: DblQ00Wchv.cxx,v 1.4 2007-02-12 17:33:50 stefspa Exp $
//<version>	$Name: not supported by cvs2svn $

//<<<<<< INCLUDES                                                       >>>>>>

#include "MuonGMdbObjects/DblQ00Wchv.h"
#include "RDBAccessSvc/IRDBQuery.h"
#include <iostream>
#include <sstream>
//#include <stdio>

//<<<<<< PRIVATE DEFINES                                                >>>>>>
//<<<<<< PRIVATE CONSTANTS                                              >>>>>>
//<<<<<< PRIVATE TYPES                                                  >>>>>>
//<<<<<< PRIVATE VARIABLE DEFINITIONS                                   >>>>>>
//<<<<<< PUBLIC VARIABLE DEFINITIONS                                    >>>>>>
//<<<<<< CLASS STRUCTURE INITIALIZATION                                 >>>>>>
//<<<<<< PRIVATE FUNCTION DEFINITIONS                                   >>>>>>
//<<<<<< PUBLIC FUNCTION DEFINITIONS                                    >>>>>>
//<<<<<< MEMBER FUNCTION DEFINITIONS                                    >>>>>>

namespace MuonGM
{

DblQ00Wchv::DblQ00Wchv(IRDBQuery* m_wchv)
 : m_nObj(0)
{
  if(m_wchv) {
    m_wchv->execute();
    m_nObj = m_wchv->size();
    m_d = new WCHV[m_nObj];
    if (m_nObj == 0) std::cerr<<"NO Wchv banks in the MuonDD Database"<<std::endl;

    int i=0;
    while(m_wchv->next()) {
        m_d[i].version     = m_wchv->data<int>("WCHV_DATA.VERS");    
        m_d[i].jsta        = m_wchv->data<int>("WCHV_DATA.JSTA");
        m_d[i].num         = m_wchv->data<int>("WCHV_DATA.NUM");
        m_d[i].heightness     = m_wchv->data<float>("WCHV_DATA.HEIGHTNESS");
        m_d[i].largeness      = m_wchv->data<float>("WCHV_DATA.LARGENESS");
        m_d[i].thickness      = m_wchv->data<float>("WCHV_DATA.THICKNESS");
//         std::cerr<<i<<" type, iz, iphi, z, r, s "<<m_d[i].type<<" "<<m_d[i].iz <<" ";
//         for(unsigned int j=0; j<8; j++)std::cerr<<m_d[i].iphi[j];
//         std::cerr<<" "<<m_d[i].z<<" "<<m_d[i].r<<" "<<m_d[i].s    <<std::endl;
        i++;
    }
    m_wchv->finalize();
  }
  else {
    m_d = new WCHV[0];
    std::cerr<<"NO Wchv banks in the MuonDD Database"<<std::endl;
  }
}
    
DblQ00Wchv::~DblQ00Wchv()
{
    delete [] m_d;
}

} // end of namespace MuonGM
