// Dear emacs, this is -*- c++ -*-

/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#ifndef XAODTRIGGER_VERSIONS_GFEXJETROI_V1_H
#define XAODTRIGGER_VERSIONS_GFEXJETROI_V1_H

// System include(s):
extern "C" {
#   include <stdint.h>
}
#include <vector>
#include <string>

// xAOD include(s):
#include "AthContainers/AuxElement.h"

namespace xAOD {

   /// Class describing properties of a LVL1 gFEX jet Trigger Object (TOB) 
   /// in the xAOD format.

   class gFexJetRoI_v1 : public SG::AuxElement {

   public:
      /// Default constructor
      gFexJetRoI_v1();

      /// Initialise the object with its most important properties: only the word for gFEX
      void initialize( uint32_t word );

      /// Object types
      enum ObjectType {
         gBlock   = 0, ///< This object is a TOB (32 bit word)
         gJet     = 1, ///< This object is a TOB (32 bit word)
         gRho     = 2  ///< This object is a TOB (32 bit word)
      };

      /// The "raw" 32-bit word describing the object candidate
      uint32_t word() const;
      /// Set the "raw" 32-bit words describing the object candidate
      void setWord( uint32_t value );


      /// TOB ET (decoded from TOB, stored for convenience)
      uint16_t tobEt() const;    /// getter for integer ET on TOB scale (3.2 GeV/count)
      void     setTobEt( uint16_t value); /// setter for the above
      unsigned int unpackEtIndex( ) const; /// retrieves the Et index from the 32-bit word
      float    etMax() const; /// floating point value (GeV, TOB scale)
      float    etMin() const; /// floating point value (GeV, TOB scale)

      /// Eta Coordinates (decoded from TOB, stored for convenience)
      uint8_t iEta() const;  /// getter for integer eta index (0-63) 
      void    setEta( uint8_t value); /// setter for the above 
      unsigned int unpackEtaIndex( ) const; /// retrieves the Eta index from the 32-bit word
      float etaMin() const; /// Floating point 
      float etaMax() const; /// Floating point 

      /// Phi coordinates
      uint8_t iPhi() const; /// Getter for integer phi index (0-32) --> check numbers for gFEX
      void  setPhi( uint8_t value); /// Setter for the above
      unsigned int unpackPhiIndex( ) const; /// retrieves the phi index from the 32-bit word
      float phiMin() const; /// Minimum value of phi corresponding to phi index. Range is always 0.2
      float phiMax() const; /// Minimum value of phi corresponding to phi index. Range is always 0.2

      /// TOB status: set to 1 if TOB Et exceeds TOB threshold (gBlocks & gJets). Status is set to 1 if Rho calculation is valid
      char status() const;
      void setStatus( char value ) ; 
      unsigned int unpackStatus( ) const; /// retrieves the Status info from the 32-bit word
      /// Energy saturation: if any gTower is saturated within gBlock and gJet, this bit is set. Always 0 for Rho.
      char saturated() const;
      void setSaturated( char value ) ; 
      unsigned int unpackSaturated( ) const; /// retrieves the Saturated info from the 32-bit word

      int gFexType () const;
      void setgFexType ( int type) ;
      int unpackType( ) const;
      ///Identification of object type with flags
      bool isgBlock() const;
      //void setIsgBlock( char value);
      bool isgJet() const;
      //void setIsgJet( char value);
      bool isgRho() const;
      //void setIsgRho( char value);
      bool isLeadingJet() const;
  




   private:


      /// Constants used in converting to ATLAS units
      static const float s_tobEtScale;
      static const float s_centralPhiWidth;
      static const float s_forwardPhiWidth;
      static const std::vector<float> s_EtaPosition; 


      /** Constants used in decoding TOB words
          For TOB word format changes these can be replaced
          by arrays in the _v2 object so that different 
          versions can be decoded by one class */

      //  Data locations within word
      static const int s_saturBit        = 31;
      static const int s_phiBit          = 26;
      static const int s_etaBit          = 20;
      static const int s_etBit           =  8;
      static const int s_statusBit       =  7;
      static const int s_resBit          =  5;
      static const int s_tobIDBit        =  0;

      //  Data masks
      static const int s_saturMask       = 0x1;
      static const int s_phiMask         = 0x1f;
      static const int s_etaMask         = 0x3f;
      static const int s_etMask          = 0xfff;
      static const int s_statusMask      = 0x1;
      static const int s_resMask         = 0x3;
      static const int s_tobIDMask       = 0x1f;

   }; // class gFexJetRoI_v1

} // namespace xAOD

// Declare the inheritance of the type:
#include "xAODCore/BaseInfo.h"
SG_BASE( xAOD::gFexJetRoI_v1, SG::AuxElement );

#endif // XAODTRIGGER_VERSIONS_GFEXJETROI_V1_H

