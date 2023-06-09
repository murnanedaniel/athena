// Dear emacs, this is -*- c++ -*-
#ifndef TRIGT1INTERFACES_LUCID_CTP_H
#define TRIGT1INTERFACES_LUCID_CTP_H

#include <stdint.h>
#include <string>

namespace LVL1 {

   /** @class LucidCTP
    *  @short Lucid input class to the CTP simulation
    *
    *         A StoreGate class to contain the output status of the
    *         level 1 LUCID trigger simulation for input into the CTP
    *         simulation.  This class contains two trigger bits in one
    *         32 bit unsigned int.
    *
    * @author Jacob Groth-Jensen <jacob.groth-jensen@hep.lu.se>
    *
    * $Revision: 187728 $
    * $Date: 2009-05-27 18:18:06 +0200 (Wed, 27 May 2009) $
    */
   class LucidCTP {
   public:
      LucidCTP( uint32_t word0 = 0 );

      /**
       * Returns an unsigned integer trigger word containing two trigger
       * bits.
       */
      uint32_t cableWord0(void) const {
         return m_cableWord0;
      }

      //! dump raw object content to string
      const std::string dump() const;

      //! print object content in a human readable form to string
      const std::string print() const;

   private:
      //! A data member to contain two trigger bits
      const uint32_t m_cableWord0;

   }; // class LucidCTP

} // namespace LVL1

#ifndef CLIDSVC_CLASSDEF_H
#include "CLIDSvc/CLASS_DEF.h"
#endif

CLASS_DEF( LVL1::LucidCTP , 48467911 , 1 )

#endif // TRIGT1INTERFACES_LUCID_CTP_H

