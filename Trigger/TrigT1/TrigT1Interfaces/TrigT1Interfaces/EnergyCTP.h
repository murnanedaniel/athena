// Dear emacs, this is -*- c++ -*-
/***************************************************************************
                         EnergyCTP.h  -  description
                            -------------------
   begin                : Friday May 05 2002
   copyright            : (C) 2002 by moyse
   email                : e.moyse@qmul.ac.uk
***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef TRIGT1INTERFACES_ENERGYCTP_H
#define TRIGT1INTERFACES_ENERGYCTP_H

// Gaudi kernel stuff.
#include "GaudiKernel/DataObject.h"

namespace LVL1 {

   /**
    *  @short "Energy" input class to the CTP simulation
    *
    *          This class defines the Energy CTP
    *          which is generated by the Energy Trigger.
    *
    * @author E. Moyse
    *
    * $Revision: 187728 $
    * $Date: 2009-05-27 18:18:06 +0200 (Wed, 27 May 2009) $
    */
   class EnergyCTP : public DataObject {

   public:
      EnergyCTP( unsigned int cableword0 = 0 );
      ~EnergyCTP();

      /** return the data
       *
       * <code>|P|4b Et Sum Map|8b EtMiss Hit Map|0|</code>
       */
      unsigned int cableWord0() const;

   private:
      const unsigned int m_cableWord0;

   }; // class EnergyCTP

} // namespace LVL1

#ifndef CLIDSVC_CLASSDEF_H
#include "CLIDSvc/CLASS_DEF.h"
#endif
CLASS_DEF( LVL1::EnergyCTP, 6254, 1 )

#endif // TRIGT1INTERFACES_ENERGYCTP_H
