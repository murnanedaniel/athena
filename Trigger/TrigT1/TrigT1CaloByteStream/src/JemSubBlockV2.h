/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TRIGT1CALOBYTESTREAM_JEMSUBBLOCKV2_H
#define TRIGT1CALOBYTESTREAM_JEMSUBBLOCKV2_H

#include <stdint.h>
#include <vector>

#include "L1CaloSubBlock.h"

namespace LVL1BS {

class JemJetElement;

/** Sub-Block class for JEM data post LS1.
 *
 *  Based on "ATLAS Level-1 Calorimeter Trigger Read-out Driver"
 *           Version X.xxx                                          //<<== CHECK
 *
 *  @author Peter Faulkner
 */

class JemSubBlockV2 : public L1CaloSubBlock {

 public:
   JemSubBlockV2();
   ~JemSubBlockV2();

   /// Clear all data
   void clear();

   /// Store JEM header
   void setJemHeader(int version, int format, int slice, int crate,
                     int module, int timeslices);
   /// Store jet element data
   void fillJetElement(int slice, const JemJetElement& jetEle);
   /// Store energy subsum data
   void setEnergySubsums(int slice, unsigned int ex,
	                 unsigned int ey, unsigned int et);

   /// Return jet element for given channel
   JemJetElement jetElement(int slice, int channel) const;
   /// Return energy subsum Ex
   unsigned int ex(int slice)      const;
   /// Return energy subsum Ey
   unsigned int ey(int slice)      const;
   /// Return energy subsum Et
   unsigned int et(int slice)      const;
   /// Return number of timeslices
   int  timeslices()               const;

   /// Pack data
   bool pack();
   /// Unpack data
   bool unpack();

 private:
   /// JEM header word ID
   static const int      s_wordIdVal       = 0xc;
   /// Data word length
   static const int      s_wordLength      = 32;
   //  Jet Element data word bit positions and masks
   static const int      s_dataIdBit       = 30;
   static const int      s_jeWordId        = 0x1;
   static const uint32_t s_dataIdMask      = 0x3;
   //  Energy subsum data bit positions and masks
   static const int      s_energyWordId    = 0x2;
   static const int      s_sourceIdBit     = 28;
   static const uint32_t s_sourceIdMask    = 0x3;
   static const int      s_etId            = 2;
   static const int      s_exEyId          = 3;
   static const int      s_exBit           = 0;
   static const int      s_eyBit           = 14;
   static const int      s_etBit           = 14;
   static const uint32_t s_exMask          = 0x3fff;
   static const uint32_t s_eyMask          = 0x3fff;
   static const uint32_t s_etMask          = 0x3fff;
   //  Neutral format data lengths
   static const int      s_pairsPerPin       = 3;
   static const int      s_jetElementBits    = 9;
   static const int      s_jePaddingBits     = 22;
   static const int      s_energyBits        = 14;
   static const int      s_bunchCrossingBits = 12;
   static const int      s_energyPaddingBits = 12;
   static const int      s_glinkBitsPerSlice = 67;

   int sourceId(uint32_t word)        const;
   int dataId(uint32_t word)          const;
   int index(int slice, int channels) const;
   void resize(std::vector<uint32_t>& vec, int channels);

   /// Pack neutral data
   bool packNeutral();
   /// Pack uncompressed data
   bool packUncompressed();
   /// Unpack neutral data
   bool unpackNeutral();
   /// Unpack uncompressed data
   bool unpackUncompressed();

   /// Jet element data
   std::vector<uint32_t> m_jeData;
   /// Energy subsum data
   std::vector<uint32_t> m_energySubsums;
   /// Number of jet element channels
   int m_channels;
   /// Number of energy data words
   int m_energyWords;

};

inline int JemSubBlockV2::sourceId(const uint32_t word) const
{
  return (word >> s_sourceIdBit) & s_sourceIdMask;
}

inline int JemSubBlockV2::dataId(const uint32_t word) const
{
  return (word >> s_dataIdBit) & s_dataIdMask;
}

} // end namespace

#endif
