//
//  Decision.cxx
//  TopoCore
//
//  Created by Joerg Stelzer on 11/16/12.
//  Copyright (c) 2012 Joerg Stelzer. All rights reserved.
//

#include "L1TopoCoreSimulation/Decision.h"
#include "L1TopoCoreSimulation/Exception.h"

void
TCS::Decision::setBit(unsigned int index, bool set) {
   if(index>=m_nBits) {
      TCS_EXCEPTION("Trying to set bit " << index << ", while max index allowed is " << m_nBits-1)
   }
	 if(set) { 
      m_decision |= (0x1<<index);
	 } else {
	    m_decision &= ~(0x1<<index);
	 }
}


std::ostream &
operator<<(std::ostream & o, const TCS::Decision & dec) {
   unsigned int bit = dec.nBits();
   if(dec.nBits()==1) {
      o << "decision (1 bit, position " << dec.firstBit() << "): ";
   } else {
      o << "decision (" << dec.nBits() << " bits, position " << dec.firstBit() << " - " << dec.lastBit() << "): ";
   }
   while(bit!=0) {
      o << (dec.bit(--bit)?1:0);
   }
   return o;
}
