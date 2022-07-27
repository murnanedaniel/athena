/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef RAWDATACOLLECTIONREADOUT_H
#define RAWDATACOLLECTIONREADOUT_H

#include "ALFA_ReadOut.h"
#include "ALFA_RawDataReadOut.h"
#include "ALFA_RawDataReadOut_charge.h"
#include <stdint.h>

#include "AthenaBaseComps/AthAlgTool.h"

#include "eformat/SourceIdentifier.h"
#include "eformat/Version.h"
#include "ByteStreamData/RawEvent.h"

#include "ALFA_RawEv/ALFA_RawDataContainer.h"
#include "ALFA_RawEv/ALFA_RawDataCollection.h"
#include "ALFA_RawEv/ALFA_RawData.h"

#include "ALFA_RawEv/ALFA_RawDataContainer_charge.h"
#include "ALFA_RawEv/ALFA_RawDataCollection_charge.h"
#include "ALFA_RawEv/ALFA_RawData_charge.h"

// Decoding methods for a block
// S. Diglio


class ALFA_RawDataCollectionReadOut : public ALFA_ReadOut
{

 private:
  
  // Raw data word and word header
  uint32_t m_dataWord;
  uint16_t m_wordHeader;
  uint16_t m_wordHeader2;
  
  
  // Decoded data
  uint32_t m_lvl1Id;   // Lvl1 Id  bit 23-0 BOB
  uint32_t m_clwc;   // MB word count  bit 15-0 LWC: 2*18+3*22+LWC+BOL+TLP+TWC=0x6a
  uint16_t m_ftyp;   // flag type
  
  uint16_t m_ecnt_LWC;   // Event counter bit 19-16 LWC
  uint16_t m_ecnt_TWC;   // Event counter bit 23-12 TWC
  uint16_t m_ecnt_EOB;   // Event counter bit 15-0  EOB
  
  uint16_t m_ctwc;   // Trailer word count bit 11-0 EOT
  
  uint16_t m_BOLflag; // Bol flags bit 23-16 BOL
  
  uint16_t m_MrodNum; // MROD number  bit 15-4 BOL
  uint16_t m_RODinput; // ROD input bit 3-0 BOL
  uint32_t m_KaptId; // Kapton number  bit 19-0 TLP
  
  // Word header Position and values
  // The word header occupies the 8 highest bits of each   word 
  static constexpr uint16_t s_headerPos  = 24;
  static constexpr uint16_t s_headerBits = 0xff;
  
  // Word header Position and values
  // The word header occupies the 12 highest bits of each   word 
  static constexpr uint16_t s_headerPos2  = 20;
  static constexpr uint16_t s_headerBits2 = 0xfff;
   
  // Beginning of buffer
  static constexpr uint16_t s_BOBvalue = 0x80;
  // Link Word Count: first word of a   link block
  static constexpr uint16_t s_LWCvalue   = 0x810;
  // Beginning of Link: signals which   link is giving its data
  static constexpr uint16_t s_BOLvalue   = 0x18;
  // TDC Link Present
  static constexpr uint16_t s_TLPvalue   = 0x890;
  // Trailer Word Count: last word of a   link block
  static constexpr uint16_t s_TWCvalue   = 0x8a;
  // End of Block (end of data from the  s)
  static constexpr uint16_t s_EOBvalue   = 0xf0;
  
 public:
  ALFA_RawDataCollectionReadOut();

  void decodeWord(uint32_t dataWord);

  bool is_BOB() const {return m_wordHeader == s_BOBvalue;};  
  bool is_LWC() const {return m_wordHeader2 == s_LWCvalue;};  
  bool is_BOL() const {return m_wordHeader == s_BOLvalue;};
  bool is_TLP() const {return m_wordHeader2 == s_TLPvalue;};
  bool is_TWC() const {return m_wordHeader == s_TWCvalue;}; 
  bool is_EOB() const {return m_wordHeader == s_EOBvalue;}; 

  // Methods to access the decoded information
  uint32_t lvl1Id() const {return m_lvl1Id;}
  uint32_t clwc()   const {return m_clwc;}
  uint16_t BOLflag() const {return m_BOLflag;}
  uint16_t ftyp()  const {return m_ftyp;}
 
  uint16_t ecnt_LWC()  const {return m_ecnt_LWC;}
  uint16_t ecnt_TWC()  const {return m_ecnt_TWC;}
  uint16_t ecnt_EOB()  const {return m_ecnt_EOB;}
  uint16_t ctwc()  const  {return m_ctwc;}
  uint16_t MrodNum() const {return m_MrodNum;}
  uint16_t RODinput() const {return m_RODinput;}
  uint16_t KaptId() const {return m_KaptId;}
  
 private:

  // Private methods
  void setZero();

  // Methods for the decoding
  //ALFA_RawDataCollectionReadOut * m_ALFA_RawDataCollectionReadOut;
  //ALFA_RawDataReadOut * m_ALFA_RawDataReadOut;

};

#endif  // RAWDATACOLLECTIONREADOUT_H




