
 //***************************************************************************
 //                           CPMHitsCnv_p1.h  -  description
 //                              -------------------
 //     begin                : 21 04 2008
 //     copyright            : (C) 2008 by Alan Watson
 //     email                : Alan.Watson@cern.ch
 //  ***************************************************************************/
  
 //***************************************************************************
 //  *                                                                         *
 //  *   This program is free software; you can redistribute it and/or modify  *
  // *   it under the terms of the GNU General Public License as published by  *
//   *   the Free Software Foundation; either version 2 of the License, or     *
 //  *   (at your option) any later version.                                   *
 //  *                                                                         *
 //  ***************************************************************************/
  
#ifndef CPMHitsCnv_p1_H
#define CPMHitsCnv_p1_H 

// Gaudi/Athena include(s):
#include "AthenaPoolCnvSvc/T_AthenaPoolTPConverter.h"

// Trigger include(s):
#include "TrigT1CaloEvent/CPMHits.h"

// Local include(s):
#include "TrigT1EventTPCnv/CPMHits_p1.h"

//using namespace LVL1;

/** T/P convertor for TrigT1Calo CPMHits class */

class CPMHitsCnv_p1 : public T_AthenaPoolTPCnvBase< LVL1::CPMHits, CPMHits_p1 > {

public:
  CPMHitsCnv_p1() : m_CPMHitsCnv(0) {}

  virtual void persToTrans( const CPMHits_p1* persObj, LVL1::CPMHits* transObj, MsgStream &log );
  virtual void transToPers( const LVL1::CPMHits* transObj, CPMHits_p1* persObj, MsgStream &log );

protected:
  CPMHitsCnv_p1 *m_CPMHitsCnv ;

}; // class TriggerTowerCnv_p1

#endif
