/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/


// System include(s):
#include <stdexcept>

// xAOD include(s):
#include "xAODCore/AuxStoreAccessorMacros.h"

// Local include(s):
#include "xAODTrigger/jFexLRJetRoI.h"

namespace xAOD {

  const float jFexLRJetRoI_v1::s_tobEtScale = 200.;
  const float jFexLRJetRoI_v1::s_towerEtaWidth = 0.1;
  const float jFexLRJetRoI_v1::s_towerPhiWidth = 0.1;
  const float jFexLRJetRoI_v1::s_minEta = -4.9;

   jFexLRJetRoI_v1::jFexLRJetRoI_v1()
     : SG::AuxElement() {
   }
   void jFexLRJetRoI_v1::initialize( uint8_t jFexNumber, uint8_t fpgaNumber, uint32_t word0) {
 
     setWord0( word0 );
     setjFexNumber( jFexNumber );
     setfpgaNumber(fpgaNumber);
     setTobEt(unpackEtTOB());
     setEta( unpackEtaIndex() );
     setPhi( unpackPhiIndex() ); 
     setSatFlag(unpackSaturationIndex());
     setGlobalEta(getGlobalEta());
     setGlobalPhi(getGlobalPhi());
    
   //include in future when xTOB in jFEX has been implemented.

   // If the object is a TOB then the isTOB should be true.
   // For xTOB default is false, but should be set if a matching TOB is found 
   // if (type() == TOB) setIsTOB(1);
   // else               setIsTOB(0);

      return;
   }

   //----------------
   /// Raw data words
   //----------------

   AUXSTORE_PRIMITIVE_SETTER_AND_GETTER( jFexLRJetRoI_v1, uint32_t, word0, setWord0 )
   AUXSTORE_PRIMITIVE_SETTER_AND_GETTER( jFexLRJetRoI_v1, uint8_t, jFexNumber, setjFexNumber )
   AUXSTORE_PRIMITIVE_SETTER_AND_GETTER( jFexLRJetRoI_v1, uint8_t, fpgaNumber, setfpgaNumber)   
   /// Extracted from data words
   AUXSTORE_PRIMITIVE_SETTER_AND_GETTER( jFexLRJetRoI_v1, uint16_t, tobEt, setTobEt )
   AUXSTORE_PRIMITIVE_SETTER_AND_GETTER( jFexLRJetRoI_v1, uint8_t, iEta, setEta )
   AUXSTORE_PRIMITIVE_SETTER_AND_GETTER( jFexLRJetRoI_v1, uint8_t, iPhi, setPhi )
   AUXSTORE_PRIMITIVE_SETTER_AND_GETTER( jFexLRJetRoI_v1, uint8_t, satFlag, setSatFlag)

   //global coordinates, stored for future use but not sent to L1Topo
   AUXSTORE_PRIMITIVE_SETTER_AND_GETTER(jFexLRJetRoI_v1, int8_t, globalEta, setGlobalEta)
   AUXSTORE_PRIMITIVE_SETTER_AND_GETTER(jFexLRJetRoI_v1, uint8_t, globalPhi, setGlobalPhi)

   //-----------------
   /// Methods to decode data from the TOB/RoI and return to the user
   //-----------------

  //include in future when xTOB in jFEX has been implemented.
   
   /// TOB or xTOB?
   //jFexLRJetRoI_v1::ObjectType jFexLRJetRoI_v1::type() const {
   //if (Word1() == 0) return TOB;
   //else              return xTOB;
   //}

   //Hardware coordinate elements  

   //Raw ET on TOB scale (200 MeV/count)
    unsigned int jFexLRJetRoI_v1::unpackEtTOB() const{
     //Data content = TOB
     return (word0() >> s_etBit) & s_etMask;

    } 

   //Return an eta index
   unsigned int jFexLRJetRoI::unpackEtaIndex() const {
     return (word0() >> s_etaBit) & s_etaMask;
   }
   //Return a phi index
   unsigned int jFexLRJetRoI::unpackPhiIndex() const {
     return (word0() >> s_phiBit) & s_phiMask;
   }

   //Return sat flag
   unsigned int jFexLRJetRoI::unpackSaturationIndex() const{
     return (word0() >> s_satBit) & s_satMask;
   }

   /// Methods that require combining results or applying scales

   /// ET on TOB scale
   unsigned int jFexLRJetRoI_v1::et() const {
    //Return the TOB Et in a 200 MeV scale
     return tobEt();
   }

   /// Returns the local coordinated within the FPGA core area
   unsigned int jFexLRJetRoI_v1::eta() const{
      return iEta();
   }

  unsigned int jFexLRJetRoI_v1::phi() const {
     return iPhi();
   }
  //Global coords
  int8_t jFexLRJetRoI_v1::getGlobalEta() const{
     int8_t globalEta = iEta() + (8*(jFexNumber() -3) -1); 
     return globalEta; 
  }


  uint8_t jFexLRJetRoI_v1::getGlobalPhi() const{
     uint8_t globalPhi = iPhi() + (fpgaNumber() * 16);
     return globalPhi;
  }


} // namespace xAOD


