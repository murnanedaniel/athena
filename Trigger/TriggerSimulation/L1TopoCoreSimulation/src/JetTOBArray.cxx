//
//  ClusterTOBArray.cxx
//  TopoCore
//
//  Created by Joerg Stelzer on 11/17/12.
//  Copyright (c) 2012 Joerg Stelzer. All rights reserved.
//

#include "L1TopoCoreSimulation/JetTOBArray.h"

void
TCS::JetTOBArray::push_back(const TCS::JetTOB& tob) {
   m_data.push_back(JetTOB::createOnHeap(tob));
}


void
TCS::JetTOBArray::print(std::ostream &o) const {
   for(const_iterator tob = begin(); tob != end(); ++tob) {
      if( tob!=begin() ) o << std::endl;
      o << **tob;
   }
}
