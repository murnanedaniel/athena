/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

/***************************************************************************
 DB data - Muon Station components
 -----------------------------------------
 Copyright (C) 2004 by ATLAS Collaboration
 ***************************************************************************/

//<doc><file>	$Id: DblQ00Atln.cxx,v 1.4 2007-02-12 17:33:50 stefspa Exp $
//<version>	$Name: not supported by cvs2svn $

//<<<<<< INCLUDES                                                       >>>>>>

#include "MuonGMdbObjects/DblQ00Atln.h"
#include "RDBAccessSvc/IRDBQuery.h"
#include <iostream>
#include <stdio.h>

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

DblQ00Atln::DblQ00Atln(IRDBQuery* m_atln)
 : m_nObj(0)
{
  if(m_atln) {
    m_atln->execute();
    m_nObj = m_atln->size();
    m_d = new ATLN[m_nObj];
    if (m_nObj == 0) std::cerr<<"NO Atln banks in the MuonDD Database"<<std::endl;

     unsigned fieldVers(1);
     unsigned fieldI(2);
     unsigned fieldIcovol(3);
     unsigned fieldZpovol(4);
     unsigned fieldWidvol(5);
     unsigned fieldNamvol(6);
     unsigned fieldJsta(7);

    int i=0;
    while(m_atln->next()) {
        m_d[i].version     = m_atln->data<int>(fieldVers);    
        m_d[i].i           = m_atln->data<int>(fieldI);          
        m_d[i].icovol      = m_atln->data<int>(fieldIcovol);
        m_d[i].zpovol      = m_atln->data<float>(fieldZpovol);
        m_d[i].widvol      = m_atln->data<float>(fieldWidvol);
        sprintf(m_d[i].namvol,"%s",m_atln->data<std::string>(fieldNamvol).c_str());
        m_d[i].jsta        = m_atln->data<int>(fieldJsta);          
        i++;
    }
    m_atln->finalize();
  }
  else {
    m_d = new ATLN[0];
    std::cerr<<"NO Atln banks in the MuonDD Database"<<std::endl;
  }
}
    
DblQ00Atln::~DblQ00Atln()
{
    delete [] m_d;
}

} // end of namespace MuonGM
