/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef MUONDIGIT_CONTAINER_p1_H
#define MUONDIGIT_CONTAINER_p1_H

#include <vector>
#include <string>
#include "MuonEventTPCnv/MuonDigitContainer/MM_Digit_p1.h"
#include "MuonEventTPCnv/MuonDigitContainer/STGC_Digit_p1.h"

namespace Muon {
/** Class to contain the Muon Prep Raw Data. */
template <class DIGIT>
class MuonDigitContainer_p1   
{
 public:
/// Default constructor
  MuonDigitContainer_p1 () {};
  
  //Contrary to PRDs, store collection info directly here. Don't bother with XXXDigitCollection_p1 etc
  
  /** Identifier of the collections.
  @todo Ideally we wouldn't write this out, and would get the Id from the hash, but this doesn't seem to be easy to do with the current IdHelper.*/ 
  std::vector<unsigned int> m_collectionId;

  /** Hash Identifier of this collection */ 
  std::vector<unsigned int> m_collectionHashId;    

  /** Collection size into master collection
      Note I use a uint16_t. If a collection has more than 2^16 elements (unlikely) it 
      has to be split into two separate collections*/
  std::vector<uint16_t>  m_size;
  
  std::vector<DIGIT>                           m_digits;
   /** The delta identifiers of the PRD i.e. we only store the difference from the last one. This saves space, compared to always writing the full 32-bit word*/
  std::vector<unsigned short>                  m_digitDeltaId;
};

typedef MuonDigitContainer_p1<MM_Digit_p1> MM_DigitContainer_p1;
typedef MuonDigitContainer_p1<STGC_Digit_p1> STGC_DigitContainer_p1;

} 

#endif
