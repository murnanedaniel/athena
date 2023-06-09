  
 //***************************************************************************
 //  *                                                                         *
 //  *   This program is free software; you can redistribute it and/or modify  *
  // *   it under the terms of the GNU General Public License as published by  *
//   *   the Free Software Foundation; either version 2 of the License, or     *
 //  *   (at your option) any later version.                                   *
 //  *                                                                         *
 //  ***************************************************************************/
  
#ifndef CMXJetTobCnv_p1_H
#define CMXJetTobCnv_p1_H

// Gaudi/Athena include(s):
#include "AthenaPoolCnvSvc/T_AthenaPoolTPConverter.h"

// Trigger include(s):
#include "TrigT1CaloEvent/CMXJetTob.h"

// Local include(s):
#include "TrigT1EventTPCnv/CMXJetTob_p1.h"

//using namespace LVL1;

/** T/P convertor for TrigT1Calo CMXJetTob class */

class CMXJetTobCnv_p1 : public T_AthenaPoolTPCnvBase< LVL1::CMXJetTob, CMXJetTob_p1 > {

public:
  CMXJetTobCnv_p1() : m_CMXJetTobCnv(0) {}

  virtual void persToTrans( const CMXJetTob_p1* persObj, LVL1::CMXJetTob* transObj, MsgStream &log );
  virtual void transToPers( const LVL1::CMXJetTob* transObj, CMXJetTob_p1* persObj, MsgStream &log );

protected:
  CMXJetTobCnv_p1 *m_CMXJetTobCnv ;

}; // class CMXJetTobCnv_p1

#endif
