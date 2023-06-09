/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
// ITrkExtrapolator.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef ISF_TRACKINGINTERFACES_ITRKEXTRAPOLATOR_H
#define ISF_TRACKINGINTERFACES_ITRKEXTRAPOLATOR_H 1

// Gaudi
#include "GaudiKernel/IAlgTool.h"
 

namespace ISF {

  // forward declarations
  class ISFParticle;


  static const InterfaceID IID_ITrkExtrapolator("ITrkExtrapolator", 1, 0);
   
  /**
   @class ITrkExtrapolator
      
   An Athena AlgTool wrapper for the Tracking Extrapolator engine.
   
   The aim of this interface is to minimize compile-time dependencies
   of ISF code on Tracking code.

       
   @author Elmar.Ritsch -at- cern.ch
   */
     
  class ITrkExtrapolator : virtual public IAlgTool {
    public:
      /** virtual destructor */
      ~ITrkExtrapolator() { ; }

      /** AlgTool interface methods */
      static const InterfaceID& interfaceID() { return IID_ITrkExtrapolator; }
      
      /** Extrapolate the given ISFParticle */
      virtual ISF::ISFParticle* extrapolate( const ISF::ISFParticle &particle ) const = 0;
  };

} // end of namespace

#endif // ISF_TRACKINGINTERFACES_ITRKEXTRAPOLATOR_H
