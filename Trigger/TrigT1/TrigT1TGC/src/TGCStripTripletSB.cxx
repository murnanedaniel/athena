/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "TrigT1TGC/TGCStripTripletSB.h"
#include "TrigT1TGC/TGCPatchPanelOut.h"
#include <iostream>
#include <cstdlib>

namespace LVL1TGCTrigger {

TGCStripTripletSB::TGCStripTripletSB()
 : TGCSlaveBoard()
{}

void TGCStripTripletSB::createSlaveBoardOut()
{
  if ( m_slaveBoardOut!=0 ) {
    delete m_slaveBoardOut;
    m_slaveBoardOut =0;
  } 

  if(m_coincidenceOut!=0){
    m_slaveBoardOut = new  TGCSlaveBoardOut(this, m_bid);
    m_slaveBoardOut->clear();
    m_slaveBoardOut->setNumberOfData(s_NumberOfStripTripletSBData);

    // fill SlaveBoardOut.
    // select largest R hit in each sections.
    int lengthOfSection = m_lengthOfCoincidenceOut/s_NumberOfStripTripletSBData;
    int i,j;
#ifdef TGCDEBUG
    std::cout <<" lengthOfCoincidenceOut= "<< m_lengthOfCoincidenceOut<<std::endl
         <<" lengthOfSection= "<<lengthOfSection<<std::endl;
#endif
    for( i=0; i<s_NumberOfStripTripletSBData; i+=1){// i=3:d 2:c 1:b 0:a, 7:d 6:c 5:b 4:a
      m_slaveBoardOut->setHit(i,false);
      for( j=0; j<lengthOfSection; j+=1) {
        if(m_coincidenceOut->getChannel(j+i*lengthOfSection)){
          m_slaveBoardOut->setPos(i,j);
          m_slaveBoardOut->setHit(i,true);
          break;
        }
      }
      if(m_slaveBoardOut->getHit(i)){
        m_slaveBoardOut->setbPos(i, m_slaveBoardOut->getPos(i));
      }
    }
  }
}

void TGCStripTripletSB::doCoincidence()
{
  const TGCHitPattern* pattern[2];
  pattern[0] = m_patchPanelOut->getHitPattern(0);
  pattern[1] = m_patchPanelOut->getHitPattern(1);

  int length;
  if(pattern[0]!=0){
    length = pattern[0]->getLength();
  }else if(pattern[1]!=0){
    length = pattern[1]->getLength();
  }else{
    length = -1;
  }
  if(m_coincidenceOut!=0) delete m_coincidenceOut;
  m_coincidenceOut =0;

  if(length>0){
    m_lengthOfCoincidenceOut = 2*length;
    m_coincidenceOut = new TGCHitPattern(m_lengthOfCoincidenceOut);

    // rearrange bit pattern for coincidence.
    bool* b = new bool [m_lengthOfCoincidenceOut];

    int j;
    for( j=0; j<m_lengthOfCoincidenceOut; j+=1){
      b[j]=false;
    }

    for( int i=0; i<length/2; i+=1){
      if(pattern[0]!=0){
	b[2*i+1]   = pattern[0]->getChannel(i+length/2);// smaller in phi
	b[2*i] = pattern[0]->getChannel(i);
      }
      if(pattern[1]!=0){
	b[length+2*i+1]   = pattern[1]->getChannel(i+length/2);// smaller in phi
	b[length+2*i] = pattern[1]->getChannel(i); 
      }
    }

    // perform 1/2 coincidence
    int block;
    for( block=0; block<2; block+=1){
      int base=block*length;
/*
        b &  c
       +b & !a & !c
       +a &  c & !b
       +c & !b & !d
       +b &  d & !c
*/

      int i=base;
      m_coincidenceOut->setChannel(i,( b[i] & !b[i+1] ));

      i=base+1;
      m_coincidenceOut->setChannel(i,(( b[i-1] &  b[i] )|
				      ( b[i-1] & !b[i] )|
				      ( b[i]   & !b[i-1] & !b[i+1] )|
				      ( b[i-1] &  b[i+1] & !b[i] )));
	
      for( i=base+2; i<base+length-1; i+=1){
        m_coincidenceOut->setChannel(i,(( b[i-1] &  b[i] )|
					( b[i-1] & !b[i-2] & !b[i] )|
					( b[i-2] &  b[i]   & !b[i-1] )|
					( b[i]   & !b[i-1] & !b[i+1] )|
					( b[i-1] &  b[i+1] & !b[i] )));
      }
      i=base+length-1;
      m_coincidenceOut->setChannel(i,(( b[i-1] &  b[i] )|
				      ( b[i-1] & !b[i-2] & !b[i] )|
				      ( b[i-2] &  b[i]   & !b[i-1] )|
				      ( b[i]   & !b[i-1] )));

#ifdef TGCCOUT
      std::cout << "StripTripletCoincidence OUT ";
      m_coincidenceOut->printb();
      std::cout << std::endl;
#endif

    }

    delete [] b;
  }
}

} //end of namespace bracket
