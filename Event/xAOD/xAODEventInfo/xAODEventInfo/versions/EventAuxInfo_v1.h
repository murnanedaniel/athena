// Dear emacs, this is -*- c++ -*-

/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// $Id: EventAuxInfo_v1.h 594273 2014-04-25 14:12:33Z krasznaa $
#ifndef XAODEVENTINFO_VERSIONS_EVENTAUXINFO_V1_H
#define XAODEVENTINFO_VERSIONS_EVENTAUXINFO_V1_H

// System include(s):
extern "C" {
#   include <stdint.h>
}
#include <vector>
#include <string>

// EDM include(s):
#include "xAODCore/AuxInfoBase.h"

namespace xAOD {

   /// Auxiliary information about the event
   ///
   /// This object describes the "static" auxiliary information about
   /// an event.
   ///
   /// @author Attila Krasznahorkay <Attila.Krasznahorkay@cern.ch>
   ///
   /// $Revision: 594273 $
   /// $Date: 2014-04-25 16:12:33 +0200 (Fri, 25 Apr 2014) $
   ///
   class EventAuxInfo_v1 : public AuxInfoBase {

   public:
      /// Default constructor
      EventAuxInfo_v1();

   private:
      /// @name Basic event information
      /// @{
      uint32_t runNumber;
      unsigned long long eventNumber;
      uint32_t lumiBlock;
      uint32_t timeStamp;
      uint32_t timeStampNSOffset;
      uint32_t bcid;
      uint32_t detectorMask0;
      uint32_t detectorMask1;
      /// @}

      /// @name Event type information
      /// @{
      std::vector< std::pair< std::string, std::string > > detDescrTags;
      uint32_t eventTypeBitmask;
      /// @}

      /// @name Trigger related information
      /// @{
      uint32_t statusElement;
      uint32_t extendedLevel1ID;
      uint16_t level1TriggerType;
      std::vector< std::string > streamTagNames;
      std::vector< std::string > streamTagTypes;
      std::vector< char > streamTagObeysLumiblock;
      /// @}

      /// @name Pileup information
      /// @{
      float actualInteractionsPerCrossing;
      float averageInteractionsPerCrossing;
      /// @}

      /// @name Detector flags
      /// @{
      uint32_t pixelFlags;
      uint32_t sctFlags;
      uint32_t trtFlags;
      uint32_t larFlags;
      uint32_t tileFlags;
      uint32_t muonFlags;
      uint32_t forwardDetFlags;
      uint32_t coreFlags;
      uint32_t backgroundFlags;
      uint32_t lumiFlags;
      /// @}

      /// @name Beam spot information
      /// @{
      float beamPosX;
      float beamPosY;
      float beamPosZ;
      float beamPosSigmaX;
      float beamPosSigmaY;
      float beamPosSigmaZ;
      float beamPosSigmaXY;
      float beamTiltXZ;
      float beamTiltYZ;
      uint32_t beamStatus;
      /// @}

   }; // class EventAuxInfo_v1

} // namespace xAOD

// Set up a CLID for the object:
#ifndef XAOD_STANDALONE
#include "SGTools/CLASS_DEF.h"
CLASS_DEF( xAOD::EventAuxInfo_v1, 38853911, 1 )
#include "SGTools/BaseInfo.h"
SG_BASE( xAOD::EventAuxInfo_v1, xAOD::AuxInfoBase );
#endif // not XAOD_STANDALONE

#endif // XAODEVENTINFO_VERSIONS_EVENTAUXINFO_V1_H
